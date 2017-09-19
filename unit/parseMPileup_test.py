import unittest

from parseMPileup import extractCigarSeq

class ExtractCigarSeqTestCase(unittest.TestCase):
    """Tests for `parseMPileup.py`."""

    def test_extractCigarSeq_only_ref_allele(self):
        """Test a cigar sequence that has only the reference allele on forward and reverse strand"""
        sequence = ".,,."
        phred = "DF!D"
        mapq = "]]ac"
        info = dict([('chr', "2"), ('pos', "323"), ('ref', "A"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [(35, 60), (37, 60), (0, 64), (35, 66)])
        self.assertTrue(lettersCount['C'] == [])
        self.assertTrue(lettersCount['G'] == [])
        self.assertTrue(lettersCount['T'] == [])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [])
        
if __name__ == '__main__':
    unittest.main()