import unittest
import pathlib
from pebbles import *

TEST_FILE = pathlib.Path(__file__).parent.joinpath("data/map.sam")


class PebblesTestCase(unittest.TestCase):
    def test_expand_cigar(self):
        self.assertEqual(expand_cigar('80M5D2M2I10M'), 'M' * 80 + 'D' * 5 + 'M' * 2 + 'I' * 2 + 'M' * 10)
        self.assertEqual(expand_cigar('2S4M5D2M2I10M'), 'SSMMMMDDDDDMMIIMMMMMMMMMM')
        pass

    def test_engap(self):
        self.assertEqual('MMMM-----MMIIMMMMMMMMMM',
                         engap(seq='MMMMMMIIMMMMMMMMMM', cigar='4M5D2M2I10M'))
        self.assertEqual('MMMMDDDDDMM--MMMMMMMMMM',
                         engap(seq='MMMMDDDDDMMMMMMMMMMMM', cigar='4M5D2M2I10M', is_reference=True))
        pass

    def test_expand_mdtag(self):
        self.assertEqual('.......A........', expand_mdtag('7A8'))
        self.assertEqual('..G..A..', expand_mdtag('2G2A2'))
        self.assertEqual('G..A', expand_mdtag('G2A'))
        self.assertEqual('.......CAT........', expand_mdtag('7^CAT8'))
        self.assertEqual('.' * 7, expand_mdtag('7'))
        self.assertEqual('.......CATG........', expand_mdtag('7^CAT0G8'))
        pass

    def test_call_mutations(self):
        self.assertEqual(call_mutations(
            **{'refname': 'AY286018', 'pos': 0,
               'expanded_engapped_md': '..............................................................................................................',
               'expanded_cigar': "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
               'gapped_read': 'ATGACACAGGCATGGGACCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGGG'
               }), None)
        self.assertEqual(call_mutations(
            **{'refname': 'AY286018', 'pos': 0,
               'expanded_engapped_md': '...............GAC............................................................................................',
               'expanded_cigar': 'MMMMMMMMMMMMMMMDDDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM',
               'gapped_read': 'ATGACACAGGCATGG---CCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGGG'
               }), ['AY286018:g.16_18delGAC'])
        self.assertEqual(call_mutations(
            **{'refname': 'AY286018', 'pos': 0,
               'expanded_engapped_md': '..................---.........................................................................................',
               'expanded_cigar': 'MMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM',
               'gapped_read': 'ATGACACAGGCATGGGACATGCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAA'
               }), ['AY286018:g.18_19insATG'])
        self.assertEqual(call_mutations(
            **{'refname': 'AY286018', 'pos': 0,
               'expanded_engapped_md': '..................CC..........................................................................................',
               'expanded_cigar': 'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM',
               'gapped_read': 'ATGACACAGGCATGGGACAGTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGGG'
               }), ['AY286018:g.19_20delinsAG'])
        self.assertEqual(call_mutations(
            **{'refname': 'AY286018', 'pos': 0,
               'expanded_engapped_md': '-..................CC.........................................................................................',
               'expanded_cigar': 'SMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM',
               'gapped_read': 'GATGACACAGGCATGGGACAGTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGG'
               }), ['AY286018:g.19_20delinsAG'])
        self.assertEqual(call_mutations(
            **{'refname': 'AY286018', 'pos': 1,
               'expanded_engapped_md': '.................CCT..........................................................................................',
               'expanded_cigar': 'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM',
               'gapped_read': 'TGACACAGGCATGGGACATGGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGAGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGGGC'
               }), ['AY286018:g.19_21delinsATG'])
        self.assertEqual(call_mutations(
            **{'refname': 'AY286018', 'pos': 0,
               'expanded_engapped_md': '..........................................................A...................................................',
               'expanded_cigar': 'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM',
               'gapped_read': 'ATGACACAGGCATGGGACCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGTGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGGG'
               }), ['AY286018:g.59A>T'])
        self.assertEqual(call_mutations(
            **{'refname': 'AY286018', 'pos': 5,
               'expanded_engapped_md': '.....................................................A...................................................',
               'expanded_cigar': 'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM',
               'gapped_read': 'ACAGGCATGGGACCCTGCAGGGTTCTTGGCTTGGCGGCGGGACGAGAACGAGGTGACGACTCGGGCAAGCCTTTTTGTTTATACCAACAGCAACAACACAAAGGG'
               }), ['AY286018:g.59A>T'])

    def test_call_mutations_from_pysam(self):
        result = [call for call in call_mutations_from_pysam(pysam.AlignmentFile(TEST_FILE, "r"))]
        self.assertEqual(result,
                         [('WT', None),
                          ('16_18delGAC', ['AY286018:g.16_18delGAC']),
                          ('18_19insATG', ['AY286018:g.18_19insATG']),
                          ('19_20delinsAG', ['AY286018:g.19_20delinsAG']),
                          ('19_20delinsAG', ['AY286018:g.19_20delinsAG']),
                          ('19_21delinsATG', ['AY286018:g.19_21delinsATG']),
                          ('59A>T', ['AY286018:g.59A>T']),
                          ('59A>T', ['AY286018:g.59A>T'])
                          ])


if __name__ == '__main__':
    unittest.main()