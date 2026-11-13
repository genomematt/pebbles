import pytest
import pysam
import pathlib
from pebbles.pebbles import *
from pebbles.genome_identifier import GenomeIdentifier

# --- CONFIGURATION & FIXTURES ---

DATA_DIR = pathlib.Path(__file__).parent.joinpath("data")
TEST_SAM = DATA_DIR.joinpath("map.sam")
TEST_BAM = DATA_DIR.joinpath("map.bam")


@pytest.fixture(scope="module")
def gi_tool(tmp_path_factory):
    """Provides a GenomeIdentifier instance with a temporary cache."""
    cache = tmp_path_factory.mktemp("genome_cache")
    return GenomeIdentifier(cache_dir=str(cache))


@pytest.fixture
def mock_bam_factory(tmp_path):
    """Creates BAM files with specific headers for identification tests."""

    def _create_bam(filename: str, sequences: list):
        path = tmp_path / filename
        header = {
            'HD': {'VN': '1.6', 'SO': 'coordinate'},
            'SQ': [{'SN': name, 'LN': length} for name, length in sequences]
        }
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            pass
        return str(path)

    return _create_bam


# --- CORE PEBBLES TESTS ---

@pytest.mark.parametrize("cigar, expected", [
    ('80M5D2M2I10M', 'M' * 80 + 'D' * 5 + 'M' * 2 + 'I' * 2 + 'M' * 10),
    ('2S4M5D2M2I10M', 'SSMMMMDDDDDMMIIMMMMMMMMMM')
])
def test_expand_cigar(cigar, expected):
    assert expand_cigar(cigar) == expected


def test_engap():
    assert engap(seq='MMMMMMIIMMMMMMMMMM', cigar='4M5D2M2I10M') == 'MMMM-----MMIIMMMMMMMMMM'
    assert engap(seq='MMMMDDDDDMMMMMMMMMMMM', cigar='4M5D2M2I10M', is_reference=True) == 'MMMMDDDDDMM--MMMMMMMMMM'


@pytest.mark.parametrize("mdtag, expected", [
    ('7A8', '.......A........'),
    ('2G2A2', '..G..A..'),
    ('G2A', 'G..A'),
    ('7^CAT8', '.......CAT........'),
    ('7', '.' * 7),
    ('7^CAT0G8', '.......CATG........')
])
def test_expand_mdtag(mdtag, expected):
    assert expand_mdtag(mdtag) == expected


def test_call_mutations():
    # Canonical WT check
    assert call_mutations(
        refname='AY286018', pos=0,
        expanded_engapped_md='.' * 110,
        expanded_cigar='M' * 110,
        gapped_read='ATGACACAGGCATGGGACCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGGG'
    ) == []

    # Deletion check
    assert call_mutations(
        refname='AY286018', pos=0,
        expanded_engapped_md='...............GAC' + '.' * 92,
        expanded_cigar='M' * 15 + 'D' * 3 + 'M' * 92,
        gapped_read='ATGACACAGGCATGG---CCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGGG'
    ) == ['AY286018:g.16_18delGAC']


def test_call_mutations_from_pysam():
    expected_mutations = [
        ('WT', []),
        ('16_18delGAC', ['AY286018:g.16_18delGAC']),
        ('18_19insATG', ['AY286018:g.18_19insATG']),
        ('19_20delinsAG', ['AY286018:g.19_20delinsAG']),
        ('19_20delinsAG', ['AY286018:g.19_20delinsAG']),
        ('19_21delinsATG', ['AY286018:g.19_21delinsATG']),
        ('59A>T', ['AY286018:g.59A>T']),
        ('59A>T', ['AY286018:g.59A>T']),
        ('42G>T;59A>T', ['AY286018:g.42G>T', 'AY286018:g.59A>T']),
    ]

    for file_path, mode in [(TEST_SAM, "r"), (TEST_BAM, "rb")]:
        result = list(call_mutations_from_pysam(pysam.AlignmentFile(file_path, mode)))
        assert result == expected_mutations


def test_count_and_dicts():
    expected_dict = {
        'AY286018:g.16_18delGAC': 1,
        'AY286018:g.18_19insATG': 1,
        'AY286018:g.19_20delinsAG': 2,
        'AY286018:g.19_21delinsATG': 1,
        'AY286018:g.59A>T': 2
    }
    assert count_dict(pysam.AlignmentFile(TEST_BAM, "rb")) == expected_dict
    assert count_dict(pysam.AlignmentFile(TEST_BAM, "rb"), row_limit=3) == {
        'AY286018:g.16_18delGAC': 1,
        'AY286018:g.18_19insATG': 1
    }


# --- GENOMEIDENTIFIER TESTS ---

@pytest.mark.parametrize("name, length, expected_acc", [
    ("chr1", 248956422, "GCF_000001405.40"),  # Human p14
    ("1", 249250621, "GCF_000001405.25"),  # Human p13/hg19
    ("chr1", 195154279, "GCF_000001635.27"),  # Mouse mm39
    ("chr1", 248387328, "GCF_009914755.1"),  # Human T2T-CHM13v2.0
])
def test_assembly_identification(gi_tool, mock_bam_factory, name, length, expected_acc):
    """Verifies that the Jaccard-style length matching identifies the correct assembly."""
    bam = mock_bam_factory(f"test_{expected_acc}.bam", [(name, length)])
    acc, aid = gi_tool.identify_assembly(bam)
    assert acc == expected_acc


def test_contaminant_handling(gi_tool, mock_bam_factory, tmp_path):
    """Ensures contaminants (EBV) are preserved while core chromosomes are translated."""
    sequences = [
        ("chr1", 248956422),
        ("EBV", 171823)  # Decoy sequence
    ]
    bam = mock_bam_factory("contam.bam", sequences)
    out_sam = str(tmp_path / "reheaded.sam")

    gi_tool.create_refseq_header(bam, out_sam)

    with pysam.AlignmentFile(out_sam, "r") as sam:
        names = [sq['SN'] for sq in sam.header['SQ']]
        assert "NC_000001.11" in names  # hg38 RefSeq ID
        assert "EBV" in names  # Unchanged

# --- REGRESSION TESTS: ADJACENT AND BOUNDARY EVENTS ---
# Reference used below: ACGTACGTAC (positions 1-10), alignments start at pos 0

def _call(cigar: str, mdtag: str, seq: str, pos: int = 0, refname: str = 'R') -> list:
    return call_mutations(refname, pos,
                          engap(expand_mdtag(mdtag), cigar, is_reference=True),
                          expand_cigar(cigar),
                          engap(seq, cigar))


@pytest.mark.parametrize("cigar, mdtag, seq, expected", [
    # mismatch immediately followed by a deletion
    ('5M2D3M', '4A0^CG3', 'ACGTTTAC', ['R:g.5A>T', 'R:g.6_7delCG']),
    # deletion immediately followed by a mismatch
    ('4M2D4M', '4^AC0G3', 'ACGTCTAC', ['R:g.5_6delAC', 'R:g.7G>C']),
    # mismatch immediately followed by an insertion
    ('5M2I5M', '4A5', 'ACGTTGGCGTAC', ['R:g.5A>T', 'R:g.5_6insGG']),
])
def test_call_mutations_adjacent_events(cigar, mdtag, seq, expected):
    assert _call(cigar, mdtag, seq) == expected


@pytest.mark.parametrize("cigar, mdtag, seq, expected", [
    # mismatch downstream of an insertion
    ('5M2I5M', '8A1', 'ACGTAGGCGTTC', ['R:g.5_6insGG', 'R:g.9A>T']),
    # deletion downstream of an insertion
    ('5M2I1M2D2M', '6^GT2', 'ACGTAGGCAC', ['R:g.5_6insGG', 'R:g.7_8delGT']),
    # insertion immediately followed by a deletion
    ('5M2I2D3M', '5^CG3', 'ACGTAGGTAC', ['R:g.5_6insGG', 'R:g.6_7delCG']),
    # two insertions in one read
    ('2M1I3M1I5M', '10', 'ACTGTATCGTAC', ['R:g.2_3insT', 'R:g.5_6insT']),
])
def test_call_mutations_after_insertion(cigar, mdtag, seq, expected):
    assert _call(cigar, mdtag, seq) == expected


@pytest.mark.parametrize("cigar, mdtag, seq, expected", [
    # insertion as the final alignment operator
    ('10M2I', '10', 'ACGTACGTACGG', ['R:g.10_11insGG']),
    # deletion as the final alignment operator
    ('8M2D', '8^AC', 'ACGTACGT', ['R:g.9_10delAC']),
    # insertion before a trailing soft clip
    ('10M2I2S', '10', 'ACGTACGTACGGTT', ['R:g.10_11insGG']),
])
def test_call_mutations_trailing_events(cigar, mdtag, seq, expected):
    assert _call(cigar, mdtag, seq) == expected


def test_secondary_alignments_ignored(sam_factory):
    sam = sam_factory([
        ['r1', '0', 'R', '1', '60', '10M', '*', '0', '0', 'ACGTTCGTAC', '*', 'MD:Z:4A5'],
        # secondary with a stored sequence must not be counted a second time
        ['r1', '256', 'R', '1', '0', '10M', '*', '0', '0', 'ACGTTCGTAC', '*', 'MD:Z:4A5'],
        # secondary with SEQ of * (the usual case) must not crash
        ['r1', '256', 'R', '1', '0', '10M', '*', '0', '0', '*', '*', 'MD:Z:4A5'],
    ])
    assert list(call_mutations_from_pysam(pysam.AlignmentFile(sam, "r"))) == [('r1', ['R:g.5A>T'])]
    assert count_dict(pysam.AlignmentFile(sam, "r")) == {'R:g.5A>T': 1}


def test_primary_without_sequence_skipped(sam_factory):
    sam = sam_factory([
        ['r1', '0', 'R', '1', '60', '10M', '*', '0', '0', '*', '*', 'MD:Z:4A5'],
        ['r2', '0', 'R', '1', '60', '10M', '*', '0', '0', 'ACGTTCGTAC', '*', 'MD:Z:4A5'],
    ])
    assert list(call_mutations_from_pysam(pysam.AlignmentFile(sam, "r"))) == [('r2', ['R:g.5A>T'])]


@pytest.mark.parametrize("cigar, mdtag, seq, expected", [
    # leading hard clip: clipped bases are absent from SEQ and do not consume reference
    ('2H10M', '4A5', 'ACGTTCGTAC', ['R:g.5A>T']),
    # trailing hard clip
    ('10M2H', '4A5', 'ACGTTCGTAC', ['R:g.5A>T']),
    # hard and soft clips together
    ('3H2S10M', '4A5', 'GGACGTTCGTAC', ['R:g.5A>T']),
    # reference skip: skipped bases are absent from both SEQ and the MD tag
    ('4M3N3M', '5A1', 'ACGTTGC', ['R:g.9A>G']),
])
def test_call_mutations_clip_and_skip_operators(cigar, mdtag, seq, expected):
    assert _call(cigar, mdtag, seq) == expected


def test_engap_non_consuming_operators():
    assert engap(seq='ACGTTGC', cigar='2H4M3N3M') == '--ACGT---TGC'
    assert engap(seq='....T.G', cigar='2H4M3N3M', is_reference=True) == '--....---T.G'


def test_skip_warnings_use_logger(sam_factory):
    class RecordingLogger:
        def __init__(self):
            self.warnings = []

        def warning(self, message, detail=None):
            self.warnings.append(message)

    logger = RecordingLogger()
    sam = sam_factory([['noMD', '0', 'R', '1', '60', '10M', '*', '0', '0', 'ACGTTCGTAC', '*']])
    list(call_mutations_from_pysam(pysam.AlignmentFile(sam, "r"), logger=logger))
    assert logger.warnings == ['skipping noMD as it has no MD tag']


def test_unknown_accession_download_is_a_clear_error(gi_tool):
    with pytest.raises(ValueError, match="No NCBI download location"):
        gi_tool.refresh_cache("GCF_999999999.1", "NotAnAssembly")


def test_every_registry_assembly_is_bundled(gi_tool):
    """Every registry entry ships with the package so identification never needs the network."""
    import importlib.resources
    bundled = importlib.resources.files("pebbles.data.assemblies")
    missing = [acc for acc, aid in gi_tool.registry if not (bundled / f"{acc}.report.txt").is_file()]
    assert missing == []


def test_identification_works_offline(tmp_path, mock_bam_factory, monkeypatch):
    """With an empty cache and no network, bundled reports are enough to identify and reject headers."""
    import requests

    def no_network(*args, **kwargs):
        raise AssertionError(f"unexpected download: {args[0] if args else kwargs}")

    monkeypatch.setattr(requests, "get", no_network)
    gi = GenomeIdentifier(cache_dir=str(tmp_path / "empty_cache"))
    assert gi.identify_assembly(mock_bam_factory("hg38.bam", [("chr1", 248956422)]))[0] == "GCF_000001405.40"
    # an amplicon reference matches nothing, which must not trigger downloads either
    with pytest.raises(ValueError, match="No match"):
        gi.identify_assembly(mock_bam_factory("amplicon.bam", [("AY286018", 1433)]))
    assert list((tmp_path / "empty_cache").iterdir()) == []


@pytest.mark.parametrize("sequences, expected", [
    # GRCh38 Ensembl/NCBI style names, with an unlocalized scaffold by GenBank accession
    ([("1", 248956422), ("X", 156040895), ("MT", 16569), ("KI270706.1", 175055)],
     ["NC_000001.11", "NC_000023.11", "NC_012920.1", "NT_187361.1"]),
    # GRCh37 b37 style names (eg GATK b37 / hs37d5)
    ([("1", 249250621), ("X", 155270560), ("MT", 16569), ("GL000192.1", 547496)],
     ["NC_000001.10", "NC_000023.10", "NC_012920.1", "NT_167207.1"]),
    # GRCm39 Ensembl names, including chromosomes that have unlocalized scaffolds
    ([("1", 195154279), ("4", 156860686), ("X", 169476592)],
     ["NC_000067.7", "NC_000070.7", "NC_000086.8"]),
    # GRCh38 UCSC names, including an unlocalized scaffold
    ([("chr1", 248956422), ("chrX", 156040895), ("chrM", 16569), ("chr1_KI270706v1_random", 175055)],
     ["NC_000001.11", "NC_000023.11", "NC_012920.1", "NT_187361.1"]),
])
def test_refseq_header_translates_name_styles(gi_tool, mock_bam_factory, tmp_path, sequences, expected):
    bam = mock_bam_factory("names.bam", sequences)
    out_sam = str(tmp_path / "reheaded.sam")
    gi_tool.create_refseq_header(bam, out_sam)
    with pysam.AlignmentFile(out_sam, "r") as sam:
        assert [sq['SN'] for sq in sam.header['SQ']] == expected


# Real headers for genomes that are not in the registry. Each includes the
# scaffolds whose lengths coincide with sequences in registered assemblies,
# which previously caused them to be identified as human or mouse.
ZEBRAFISH_GRCZ11 = [
    ('1', 59578282), ('2', 59640629), ('3', 62628489), ('4', 78093715), ('5', 72500376), ('6', 60270059),
    ('7', 74282399), ('8', 54304671), ('9', 56459846), ('10', 45420867), ('11', 45484837), ('12', 49182954),
    ('13', 52186027), ('14', 52660232), ('15', 48040578), ('16', 55266484), ('17', 53461100), ('18', 51023478),
    ('19', 48449771), ('20', 55201332), ('21', 45934066), ('22', 39133080), ('23', 46223584), ('24', 42172926),
    ('25', 37502051), ('MT', 16596),
    ('KN149777.1', 1553), ('KN149846.1', 28664), ('KN149909.1', 40056), ('KN150119.1', 1143),
    ('KZ116038.1', 1428), ('KZ115432.1', 272653),
]
CHICKEN_GRCG7B = [
    ('1', 196449156), ('2', 149539284), ('3', 110642502), ('4', 90861225), ('5', 59506338), ('6', 36220557),
    ('7', 36382834), ('8', 29578256), ('9', 23733309), ('10', 20453248), ('11', 19638187), ('12', 20119077),
    ('13', 17905061), ('14', 15331188), ('15', 12703657), ('16', 2706039), ('17', 11092391), ('18', 11623896),
    ('19', 10455293), ('20', 14265659), ('21', 6970754), ('22', 4686657), ('23', 6253421), ('24', 6478339),
    ('25', 3067737), ('26', 5349051), ('27', 5228753), ('28', 5437364), ('29', 726478), ('30', 755666),
    ('31', 2457334), ('32', 125424), ('33', 3839931), ('34', 3469343), ('35', 554126), ('36', 358375),
    ('37', 157853), ('38', 667312), ('39', 177356), ('W', 9109940), ('Z', 86044486), ('MT', 16784),
    ('JAENSK010000673.1', 2415),
]


@pytest.mark.parametrize("sequences", [ZEBRAFISH_GRCZ11, CHICKEN_GRCG7B], ids=["zebrafish_GRCz11", "chicken_GRCg7b"])
def test_unregistered_genome_with_coincidental_scaffold_lengths_is_not_identified(gi_tool, mock_bam_factory, sequences):
    bam = mock_bam_factory("unregistered.bam", sequences)
    with pytest.raises(ValueError, match="No match"):
        gi_tool.identify_assembly(bam)


def test_human_header_with_many_unlisted_decoys_is_identified(gi_tool, mock_bam_factory):
    """hs38DH-style headers add thousands of decoy and HLA sequences that are not in the assembly report."""
    report = gi_tool.fetch_assembly_report("GCF_000001405.40", "GRCh38.p14")
    chromosomes = [(f"chr{r['name']}", int(r['length'])) for _, r in report.iterrows() if r['role'] == 'assembled-molecule']
    reported_lengths = set(report['length'].astype(int))
    decoy_lengths = [length for length in range(1001, 20000, 7) if length not in reported_lengths][:2900]
    decoys = [(f"chrUn_JTFH0100{i:05d}v1_decoy", length) for i, length in enumerate(decoy_lengths)]
    bam = mock_bam_factory("hs38DH_like.bam", chromosomes + decoys + [("chrEBV", 171823)])
    assert gi_tool.identify_assembly(bam) == ("GCF_000001405.40", "GRCh38.p14")


def test_header_without_any_chromosome_is_not_identified(gi_tool, mock_bam_factory):
    # an unlocalized GRCh38 scaffold on its own is not enough evidence for an assembly
    bam = mock_bam_factory("scaffold_only.bam", [("chr1_KI270706v1_random", 175055)])
    with pytest.raises(ValueError, match="No match"):
        gi_tool.identify_assembly(bam)
