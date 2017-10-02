import unittest
import sys
import StringIO
import os
import re

from parseMPileup import extractArguments
from parseMPileup import extractCigarSeq
from parseMPileup import parsePileup
from tempfile import NamedTemporaryFile

class ParseMPileupTestCase(unittest.TestCase):
    """Tests for `parseMPileup.py`."""

    def createTempPileup01(self):
        tempFile01 = NamedTemporaryFile(delete=False)
        #outfile_path = tempfile.mkstemp()[1]
        #oFile = open(outfile_path, 'w')
        #$oFile.write("chr1\t3153345\tT\t2\t.,+4atta\t^!\t]]\t100,61\n")
        tempFile01.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr1", "3153345", "T", "2", ".,+4atta", "^!", "]]", "100,61"))   
        tempFile01.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr1", "3800923", "A", "3", "Tt^!.", "g!6", "]]!", "52,66,1")) 
        tempFile01.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr2", "4688598", "G", "3", ".Cc", "Cg!", "]]]", "17,10,11"))
        tempFile01.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr4", "5134371", "C", "3", ".-1A,-1a.", "f!C", "]]]", "74,23,7"))
        tempFile01.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr5", "6916649", "G", "2", ".$,$", "K!", "$$", "94,101"))
        tempFile01.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr6", "37108587", "C", "4", ".,.N", "a!D#", "]]]]", "88,89,52,12"))
        tempFile01.close()
        return(tempFile01.name)

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
        (inputA, outputB, separatedFiles) = extractArguments()
        self.assertTrue(inputA == "toto.txt")
        self.assertTrue(outputB == "titi_")
        self.assertFalse(separatedFiles, "test_extractArguments_goods_args: The separatedFiles argument doesn't have the expected value.")
    
    def test_extractArguments_with_s_arg(self):
        """Test extraction of arguments when undefined argument passed to function"""
        sys.argv = ["prog", "-i", "effect.txt", "-p", "oneTest_", "-s"]
        (inputA, outputB, separatedFiles) = extractArguments()
        self.assertTrue(inputA == "effect.txt")
        self.assertTrue(outputB == "oneTest_")
        self.assertTrue(separatedFiles, "test_extractArguments_with_s_arg: The separatedFiles argument doesn't have the expected value.")
    
    def test_extractArguments_with_upper_s_arg(self):
        """Test extraction of arguments when undefined argument passed to function"""
        capturedOutput = StringIO.StringIO()      # A StringIO object
        sys.stdout = capturedOutput               # Redirect stdout
        sys.argv = ["prog", "-p", "test_", "-i", "where.txt", "-S"]
        with self.assertRaises(SystemExit):
            extractArguments()
        self.assertEqual(capturedOutput.getvalue(), "usage: parsePileup.py -i <inputFile> -p <outputPrefix> [-s] [-h]\n")
        
    def test_extractArguments_help(self):
        """Test extraction of arguments when undefined argument passed to function"""
        capturedOutput = StringIO.StringIO()      # A StringIO object
        sys.stdout = capturedOutput               # Redirect stdout
        sys.argv = ["prog", "-h"]  
        with self.assertRaises(SystemExit):   
            extractArguments()
        self.assertEqual(capturedOutput.getvalue(), "usage: parsePileup.py -i <inputFile> -p <outputPrefix> [-s] [-h]\n")
    
    def test_extractArguments_no_arg(self):
        """Test extraction of arguments when no argument passed to function"""
        capturedOutput = StringIO.StringIO()      # A StringIO object
        sys.stdout = capturedOutput               # Redirect stdout
        sys.argv = ["prog"]  
        with self.assertRaises(SystemExit):   
            extractArguments()
        self.assertEqual(capturedOutput.getvalue(), "usage: parsePileup.py -i <inputFile> -p <outputPrefix> [-s] [-h]\n")
    
    def test_parseMPileup_01(self):
        inputfile_path = self.createTempPileup01()
        try:
            parsePileup(inputfile_path, "toto", False)
            contents = open("toto.txt").read()
        finally:
            # NOTE: To retain the tempfile if the test fails, remove
            # the try-finally clauses
            if os.path.exists(inputfile_path):
                os.unlink(inputfile_path)
            files = [ f for f in os.listdir(".") if re.match(r'toto.*\.txt', f)]
            for f in files:
                if os.path.exists(f):
                    os.unlink(f)
            
        expected = "chr1\t3153345\t2\t0\t0\t0\t2\t1\t0\n"
        expected = expected + "chr1\t3800923\t3\t1\t0\t0\t2\t0\t0\n"
        expected = expected + "chr2\t4688598\t3\t0\t2\t1\t0\t0\t0\n"
        expected = expected + "chr4\t5134371\t3\t0\t3\t0\t0\t2\t0\n"
        expected = expected + "chr5\t6916649\t2\t0\t0\t2\t0\t0\t0\n"
        expected = expected + "chr6\t37108587\t4\t0\t3\t0\t0\t0\t1\n"
        self.assertEqual(contents, expected)
        
        
if __name__ == '__main__':
    unittest.main()