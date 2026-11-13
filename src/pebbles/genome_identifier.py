import pysam
import pandas as pd
import io
import requests
from pathlib import Path
from typing import Dict, Any, Optional, List, Set, Tuple
import importlib.resources

pysam.set_verbosity(0) # to silence HTSlib .bai errors


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

        # Human GRCh38: p0 - p14
        manifest.extend([(f"GCF_000001405.{26 + i}", f"GRCh38" if i == 0 else f"GRCh38.p{i}") for i in range(15)])

        # Human GRCh37: Filtered to verified available versions
        grch37_versions = [
            (13, "GRCh37"), (14, "GRCh37.p2"), (17, "GRCh37.p5"), (21, "GRCh37.p9"),
            (22, "GRCh37.p10"), (23, "GRCh37.p11"), (24, "GRCh37.p12"), (25, "GRCh37.p13")
        ]
        manifest.extend([(f"GCF_000001405.{v}", aid) for v, aid in grch37_versions])

        # Human T2T-CHM13 (a separate RefSeq assembly series from GRCh)
        manifest.append(("GCF_009914755.1", "T2T-CHM13v2.0"))

        # Mouse
        manifest.append(("GCF_000001635.27", "GRCm39"))
        manifest.extend([(f"GCF_000001635.{20 + i}", f"GRCm38" if i == 0 else f"GRCm38.p{i}") for i in range(7)])

        # Yeast
        manifest.append(("GCF_000146045.2", "R64"))

        return manifest

    def fetch_assembly_report(self, acc: str, aid: str, force_update: bool = False) -> pd.DataFrame:
        """Retrieves assembly report following priority: Cache -> Bundled -> Live."""
        cache_path = self.cache_dir / f"{acc}.report.txt"

        # 1. Check user cache (~/.genome_cache)
        if cache_path.exists() and not force_update:
            return self._parse_report(cache_path)

        # 2. Check bundled package data (pebbles/data/assemblies/)
        try:
            # Anchor updated to pebbles.data.assemblies
            pkg_data = importlib.resources.files("pebbles.data.assemblies") / f"{acc}.report.txt"
            if pkg_data.exists() and not force_update:
                # We convert to Path for consistent handling in _parse_report
                return self._parse_report(Path(str(pkg_data)))
        except (ImportError, AttributeError, FileNotFoundError):
            # Fallback if package structure is not installed or file is missing
            pass

        # 3. Download from NCBI if not found locally
        return self.refresh_cache(acc, aid)

    def refresh_cache(self, acc: str, aid: str) -> pd.DataFrame:
        """Downloads the report using robust species paths and updates the local cache."""
        # 1. Determine the category and species based on the accession/aid
        # This logic matches the hierarchy used in your build script
        if "000001405" in acc or "009914755" in acc:
            cat, species = "vertebrate_mammalian", "Homo_sapiens"
        elif "000001635" in acc:
            cat, species = "vertebrate_mammalian", "Mus_musculus"
        elif "000146045" in acc:
            cat, species = "fungi", "Saccharomyces_cerevisiae"
        else:
            raise ValueError(f"No NCBI download location is known for {aid} ({acc})")

        base_ftp = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq"
        folder_name = f"{acc}_{aid}"
        url = f"{base_ftp}/{cat}/{species}/all_assembly_versions/{folder_name}/{folder_name}_assembly_report.txt"

        try:
            response = requests.get(url, timeout=30)

            # 2. Handle GRCh37 p0 / base naming inconsistency
            if response.status_code == 404 and "GRCh37" in aid and ".p" not in aid:
                alt_folder = f"{acc}_{aid}.p0"
                url = f"{base_ftp}/{cat}/{species}/all_assembly_versions/{alt_folder}/{alt_folder}_assembly_report.txt"
                response = requests.get(url, timeout=30)

            response.raise_for_status()

            # 3. Save to user cache (~/.genome_cache)
            cache_path = self.cache_dir / f"{acc}.report.txt"
            cache_path.write_text(response.text, encoding='utf-8')
            return self._parse_report(cache_path)

        except requests.exceptions.RequestException as e:
            raise RuntimeError(f"Failed to download assembly report for {aid} ({acc}) from {url}: {e}")

    def _parse_report(self, path: Path) -> pd.DataFrame:
        """Parses NCBI report using fixed column positions to avoid header name shifts."""
        try:
            df = pd.read_csv(
                path,
                sep='\t',
                comment='#',
                header=None,
                low_memory=False
            )

            # NCBI Standard 10-column mapping (stable across species/versions)
            column_map = {
                0: 'name',  # Sequence-Name
                1: 'role',  # Sequence-Role
                2: 'molecule',  # Assigned-Molecule
                3: 'location',  # Location/Type
                4: 'genbank',  # GenBank-Accn
                5: 'relationship',  # Relationship
                6: 'refseq',  # RefSeq-Accn
                7: 'unit',  # Assembly-Unit
                8: 'length',  # Sequence-Length
                9: 'ucsc'  # UCSC-Style-Name
            }

            df = df.rename(columns=column_map)

            # Defensive check: if the file had more or fewer than 10 columns,
            # ensure we didn't crash and at least 'length' is numeric.
            if 'length' in df.columns:
                df['length'] = pd.to_numeric(df['length'], errors='coerce').fillna(0).astype(int)

            return df
        except Exception as e:
            print(f"Error parsing {path}: {e}")
            return pd.DataFrame()

    # Fraction of the header's total bases that must be in sequences whose lengths
    # are in the assembly report for that assembly to be accepted as a match
    MIN_MATCHED_BASE_FRACTION = 0.5

    def identify_assembly(self, bam_path: str) -> Tuple[str, str]:
        """Identify the registered assembly that a SAM/BAM header was aligned to.

        Header sequences are matched to assembly report sequences by length.
        A whole-genome header for an unregistered genome usually contains a few
        scaffolds whose lengths coincide with sequences in a registered assembly,
        so an assembly is only accepted when the matched sequences make up at
        least MIN_MATCHED_BASE_FRACTION of the header's bases and include at
        least one of the assembly's chromosomes (assembled molecules).

        The assembly with the most matched bases (then most matched sequences)
        is returned. Ties within one accession series (eg GRCh38 patch levels)
        go to the newest version; other ties keep the earlier registry entry.

        Raises:
            ValueError: if no registered assembly meets these criteria
        """
        with pysam.AlignmentFile(bam_path, "rb") as sam:
            bam_lengths = [int(length) for length in sam.lengths]
        total_bases = sum(bam_lengths)

        best: Optional[Tuple[Tuple[int, int], str, str]] = None
        for acc, aid in self.registry:
            try:
                df = self.fetch_assembly_report(acc, aid)
            except Exception:
                continue
            if df.empty:
                continue

            report_lengths = set(df['length'].astype(int))
            chromosome_lengths = set(df.loc[df['role'] == 'assembled-molecule', 'length'].astype(int))
            matched = [length for length in bam_lengths if length in report_lengths]
            matched_bases = sum(matched)

            if not any(length in chromosome_lengths for length in matched):
                continue
            if matched_bases < self.MIN_MATCHED_BASE_FRACTION * total_bases:
                continue

            score = (matched_bases, len(matched))
            if best is None or score > best[0] or (score == best[0] and self._is_newer_version(acc, best[1])):
                best = (score, acc, aid)

        if best is None:
            raise ValueError(f"No match for {bam_path}: no registered assembly accounts for "
                             f"at least {self.MIN_MATCHED_BASE_FRACTION:.0%} of the bases in its "
                             f"{len(bam_lengths)} header sequences")

        return best[1], best[2]

    @staticmethod
    def _is_newer_version(acc: str, other: str) -> bool:
        """True if acc is a later version of the same accession as other (eg GCF_000001405.40 vs .26)"""
        base, _, version = acc.partition('.')
        other_base, _, other_version = other.partition('.')
        return base == other_base and int(version) > int(other_version)

    def get_id_translation_table(self, bam_path: str) -> Dict[str, str]:
        """Creates an alias-to-RefSeq translation dictionary.

        Aliases are the sequence name (Ensembl style, eg '1'), GenBank accession
        and UCSC style name (eg 'chr1'). The assigned molecule column is only an
        alias for assembled molecules: scaffolds and patches also record the
        chromosome they belong to there, so using it for them would map '1' to
        a chromosome 1 scaffold. Assembled molecules are added first and are
        never overwritten.
        """
        acc, aid = self.identify_assembly(bam_path)
        df = self.fetch_assembly_report(acc, aid)

        is_assembled = df['role'] == 'assembled-molecule'
        translation_table: Dict[str, str] = {}
        for _, row in pd.concat([df[is_assembled], df[~is_assembled]]).iterrows():
            refseq_id = str(row['refseq'])
            if refseq_id == 'na':
                continue
            aliases = [row['name'], row['genbank'], row['ucsc']]
            if row['role'] == 'assembled-molecule':
                aliases.append(row['molecule'])
            for alias in aliases:
                if pd.notna(alias) and str(alias).lower() != 'na':
                    translation_table.setdefault(str(alias), refseq_id)
        return translation_table

    def create_refseq_header(self, input_bam: str, output_sam: str):
        """Generates a header-only SAM with versioned RefSeq identifiers."""
        acc, aid = self.identify_assembly(input_bam)
        df = self.fetch_assembly_report(acc, aid)

        valid_lengths = {
            str(row['refseq']): int(row['length'])
            for _, row in df.iterrows()
            if str(row['refseq']).lower() != 'na'
        }
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