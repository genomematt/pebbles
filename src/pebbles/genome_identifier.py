import pysam
import pandas as pd
import io
import requests
from pathlib import Path
from typing import Dict, Any, Optional, List, Set, Tuple
import importlib.resources


class GenomeIdentifier:
    """
    Identifies genomic assemblies from BAM headers of pre-mapped reads and
    reconciles sequence naming conventions (UCSC/Ensembl/GenBank) to
    versioned NCBI RefSeq IDs.
    """

    def __init__(self, cache_dir: Optional[str] = None):
        self.cache_dir = Path(cache_dir) if cache_dir else Path.home() / ".genome_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.registry = self._generate_manifest()

    def _generate_manifest(self) -> List[Tuple[str, str]]:
        """Internal registry of supported GCF accessions and assembly aliases."""
        manifest = []
        manifest.extend([(f"GCF_000001405.{26 + i}", f"GRCh38.p{i}") for i in range(15)])
        manifest.extend([(f"GCF_000001405.{13 + i}", f"GRCh37.p{i}") for i in range(14)])
        manifest.append(("GCF_000001635.27", "GRCm39"))
        manifest.extend([(f"GCF_000001635.{20 + i}", f"GRCm38.p{i}") for i in range(7)])
        return manifest

    def fetch_assembly_report(self, acc: str, aid: str, force_update: bool = False) -> pd.DataFrame:
        """Retrieves assembly report following priority: Cache -> Bundled -> Live."""
        cache_path = self.cache_dir / f"{acc}.report.txt"
        if cache_path.exists() and not force_update:
            return self._parse_report(cache_path)

        try:
            pkg_data = importlib.resources.files("genome_identifier.data.assemblies") / f"{acc}.report.txt"
            if pkg_data.exists() and not force_update:
                return self._parse_report(Path(str(pkg_data)))
        except (ImportError, AttributeError):
            pass

        return self.refresh_cache(acc, aid)

    def refresh_cache(self, acc: str, aid: str) -> pd.DataFrame:
        """Downloads the report and updates the local cache."""
        prefix = acc.replace("_", "/")[:15]
        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{acc}_{aid}/{acc}_{aid}_assembly_report.txt"

        response = requests.get(url, timeout=30)
        response.raise_for_status()

        cache_path = self.cache_dir / f"{acc}.report.txt"
        cache_path.write_text(response.text, encoding='utf-8')
        return self._parse_report(cache_path)

    def _parse_report(self, path: Path) -> pd.DataFrame:
        df = pd.read_csv(path, sep='\t', comment='#', header=None, low_memory=False)
        df.columns = ['name', 'role', 'molecule', 'genbank', 'rel', 'refseq', 'unit', 'length', 'ucsc', 'seq_acc']
        return df

    def identify_assembly(self, bam_path: str) -> Tuple[str, str]:
        """
        Detects the assembly version by calculating the intersection of
        sequence lengths. Robust against added contaminants or decoys.
        """
        with pysam.AlignmentFile(bam_path, "rb") as sam:
            bam_lengths = set(sam.lengths)

        best_acc, best_aid, max_match = None, None, 0

        for acc, aid in self.registry:
            try:
                df = self.fetch_assembly_report(acc, aid)
                report_lengths = set(df['length'].astype(int))
                # Jaccard-style intersection logic
                match_count = len(bam_lengths.intersection(report_lengths))
                if match_count > max_match:
                    max_match = match_count
                    best_acc, best_aid = acc, aid
            except Exception:
                continue

        if not best_acc:
            raise ValueError("No matching assembly found for BAM sequences.")
        return best_acc, best_aid

    def get_id_translation_table(self, bam_path: str) -> Dict[str, str]:
        """Creates an alias-to-RefSeq translation dictionary."""
        acc, aid = self.identify_assembly(bam_path)
        df = self.fetch_assembly_report(acc, aid)

        translation_table = {}
        for _, row in df.iterrows():
            refseq_id = str(row['refseq'])
            if refseq_id == 'na': continue
            for alias in [row['molecule'], row['genbank'], row['ucsc'], row['name']]:
                if pd.notna(alias) and str(alias).lower() != 'na':
                    translation_table[str(alias)] = refseq_id
        return translation_table

    def create_refseq_header(self, input_bam: str, output_sam: str):
        """Generates a header-only SAM with versioned RefSeq identifiers."""
        acc, aid = self.identify_assembly(input_bam)
        df = self.fetch_assembly_report(acc, aid)

        valid_lengths = {str(row['refseq']): int(row['length']) for _, row in df.iterrows() if row['refseq'] != 'na'}
        id_lookup = self.get_id_translation_table(input_bam)

        with pysam.AlignmentFile(input_bam, "rb") as sam:
            header = sam.header.to_dict()
            new_sq = []
            for sq in header.get('SQ', []):
                name, length = sq['SN'], sq['LN']
                translated_id = id_lookup.get(name, name)

                # Check for length-validated sequence matches
                if translated_id in valid_lengths and valid_lengths[translated_id] == length:
                    new_sq.append({'SN': translated_id, 'LN': length})
                else:
                    # Treat as unidentifiable or contaminant (e.g., EBV)
                    new_sq.append({'SN': name, 'LN': length})
            header['SQ'] = new_sq

        with pysam.AlignmentFile(output_sam, "wh", header=header) as out:
            pass