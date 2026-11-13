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
