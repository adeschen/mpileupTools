import unittest
import sys
import StringIO

from parseMPileup import extractCigarSeq
from parseMPileup import extractArguments

class ParseMPileupTestCase(unittest.TestCase):
    """Tests for `parseMPileup.py`."""


    def tearDown(self):
        """ Reset redirected output to default"""
        sys.stdout = sys.__stdout__               
        
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
    
    def test_extractCigarSeq_ref_and_one_alternative_allele(self):
        """Test a cigar sequence that has the reference allele on forward and reverse strand and one alternative allele"""
        sequence = ".,T.,"
        phred = "DF!DG"
        mapq = "]]acb"
        info = dict([('chr', "2"), ('pos', "323"), ('ref', "C"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [])
        self.assertTrue(lettersCount['C'] == [(35, 60), (37, 60), (35, 66), (38, 65)])
        self.assertTrue(lettersCount['G'] == [])
        self.assertTrue(lettersCount['T'] == [(0, 64)])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [])
    
    def test_extractCigarSeq_ref_and_two_alternative_alleles(self):
        """Test a cigar sequence that has the reference allele on forward and reverse strand and two alternative alleles"""
        sequence = ".AG.,"
        phred = "!F!DG"
        mapq = "]hacb"
        info = dict([('chr', "2"), ('pos', "3423"), ('ref', "T"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [(37, 71)])
        self.assertTrue(lettersCount['C'] == [])
        self.assertTrue(lettersCount['G'] == [(0, 64)])
        self.assertTrue(lettersCount['T'] == [(0, 60), (35, 66), (38, 65)])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [])
    
    def test_extractCigarSeq_ref_and_N(self):
        """Test a cigar sequence that has the reference allele on forward and N"""
        sequence = ".N"
        phred = "EI"
        mapq = "[\\"
        info = dict([('chr', "2"), ('pos', "3423"), ('ref', "T"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [])
        self.assertTrue(lettersCount['C'] == [])
        self.assertTrue(lettersCount['G'] == [])
        self.assertTrue(lettersCount['T'] == [(36, 58)])
        self.assertTrue(lettersCount['N'] == [(40, 59)])
        self.assertTrue(lettersCount['O'] == [])
    
    def test_extractCigarSeq_one_insertion(self):
        """Test a cigar sequence that has the reference allele and one insertion"""
        sequence = ".+2AC,"
        phred = "EI"
        mapq = "[\\"
        info = dict([('chr', "2"), ('pos', "3423"), ('ref', "T"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [])
        self.assertTrue(lettersCount['C'] == [])
        self.assertTrue(lettersCount['G'] == [])
        self.assertTrue(lettersCount['T'] == [(36, 58), (40, 59)])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [(-1, -1)])
    
    def test_extractCigarSeq_two_insertions(self):
        """Test a cigar sequence that has the reference allele and two insertions"""
        sequence = ".+2AC+1T,"
        phred = "EI"
        mapq = "[K"
        info = dict([('chr', "2"), ('pos', "3423"), ('ref', "T"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [])
        self.assertTrue(lettersCount['C'] == [])
        self.assertTrue(lettersCount['G'] == [])
        self.assertTrue(lettersCount['T'] == [(36, 58), (40, 42)])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [(-1, -1), (-1, -1)])
    
    def test_extractCigarSeq_one_deletion(self):
        """Test a cigar sequence that has the reference allele and one deletion"""
        sequence = ".-1T,"
        phred = "TI"
        mapq = "[H"
        info = dict([('chr', "2"), ('pos', "3423"), ('ref', "T"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [])
        self.assertTrue(lettersCount['C'] == [])
        self.assertTrue(lettersCount['G'] == [])
        self.assertTrue(lettersCount['T'] == [(51, 58), (40, 39)])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [(-1, -1)])
    
    def test_extractCigarSeq_two_deletions(self):
        """Test a cigar sequence that has the reference allele and two deletions"""
        sequence = ".-1T,-2AT"
        phred = "QI"
        mapq = "!E"
        info = dict([('chr', "2"), ('pos', "3423"), ('ref', "G"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [])
        self.assertTrue(lettersCount['C'] == [])
        self.assertTrue(lettersCount['G'] == [(48, 0), (40, 36)])
        self.assertTrue(lettersCount['T'] == [])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [(-1, -1), (-1, -1)])
    
    def test_extractCigarSeq_end_of_read(self):
        """Test a cigar sequence that has the reference allele and two deletions"""
        sequence = "+2CT-1T$.-2AT$,"
        phred = "UT"
        mapq = "!!"
        info = dict([('chr', "2"), ('pos', "3423"), ('ref', "G"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [])
        self.assertTrue(lettersCount['C'] == [])
        self.assertTrue(lettersCount['G'] == [(52, 0), (51, 0)])
        self.assertTrue(lettersCount['T'] == [])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [(-1, -1), (-1, -1), (-1, -1)])
    
    def test_extractCigarSeq_start_of_read(self):
        """Test a cigar sequence that has reference and alternative alleles, indels, end and start of reads"""
        sequence = "+3CTT-1T$.^D-4ATAT$a"
        phred = "UB"
        mapq = "!["
        info = dict([('chr', "2"), ('pos', "3423"), ('ref', "C"), ('NB', "4")])
        lettersCount = extractCigarSeq(sequence, phred, mapq, info)
        self.assertTrue(len(lettersCount) == 6)
        self.assertTrue(lettersCount['A'] == [(33, 58)])
        self.assertTrue(lettersCount['C'] == [(52, 0)])
        self.assertTrue(lettersCount['G'] == [])
        self.assertTrue(lettersCount['T'] == [])
        self.assertTrue(lettersCount['N'] == [])
        self.assertTrue(lettersCount['O'] == [(-1, -1), (-1, -1), (-1, -1)])
    
    def test_extractArguments_wrong_args(self):
        """Test extraction of arguments when undefined argument passed to function"""
        sys.argv = ["prog", "-d"]
        with self.assertRaises(SystemExit):
                extractArguments()
    
    def test_extractArguments_missing_arg_value(self):
        """Test extraction of arguments when undefined argument passed to function"""
        sys.argv = ["prog", "-i"]
        with self.assertRaises(SystemExit):
                extractArguments()
    
    def test_extractArguments_goods_args(self):
        """Test extraction of arguments when undefined argument passed to function"""
        sys.argv = ["prog", "-i", "toto.txt", "-p", "titi_"]
        (inputA, outputB) = extractArguments()
        self.assertTrue(inputA == "toto.txt")
        self.assertTrue(outputB == "titi_")
        
    def test_extractArguments_help(self):
        """Test extraction of arguments when undefined argument passed to function"""
        capturedOutput = StringIO.StringIO()      # A StringIO object
        sys.stdout = capturedOutput               # Redirect stdout
        sys.argv = ["prog", "-h"]  
        with self.assertRaises(SystemExit):   
            extractArguments()
        self.assertEqual(capturedOutput.getvalue(), "usage: parsePileup.py -i <inputFile> -p <outputPrefix>\n")
        
if __name__ == '__main__':
    unittest.main()