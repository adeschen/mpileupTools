#!/usr/bin/python

#############################################################################
## Title: parseMPileup.py
## Description: Parse output of samtools mpileup
##
## Author: Pascal Belleau and Astrid Deschenes
## Creation: 2017-09-14
## License: GPL-3
#############################################################################

##################
## IMPORT
##################

import sys
import getopt
from re import split
from re import match

def extractArguments():
    """
    Extract argument values as input by user.

    Keyword arguments:
    none
    """
    
    inputFile = ''
    outputPrefix = ''
    usage = 'usage: parsePileup.py -i <inputFile> -p <outputPrefix>'
    
    try:
        opts, arg = getopt.getopt(sys.argv[1:], "hi:p:", ["help", "ifile=", "pfile="])
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise

    if len(opts) == 0:
        print usage
        sys.exit(2)
           
    for opt, arg in opts:
        if opt == "-h":
            print usage
            sys.exit(0)
        elif opt in ("-i", "--ifile"):
            inputFile = arg
        elif opt in ("-p", "--pfile"):
            outputPrefix = arg
         
    print 'Input file is "', inputFile, '"'
    print 'Output prefix file is "', outputPrefix, '"'
    return(inputFile, outputPrefix)


def open_read_file(inputfile):
    """
    Open reading file and return file pointer.
    An IOError is raised when the file cannot be opened.
    
    Keyword arguments:
    inputfile -- the name of the input file
    """
    try:    
        input_file = open(inputfile, 'r')
    except IOError:                     
        raise IOError("Cannot open file : %s \n\n" % (inputfile))
    return(input_file)
           
def open_write_file(outputFile):
    """
    Open writing file and return file pointer.
    An IOError is raised when the file cannot be opened.
    
    Keyword arguments:
    outputFile -- the name of the output file
    """
    try:    
        output_file = open(outputFile, 'w')
    except IOError:                     
        raise IOError("Cannot open file : %s \n\n" % (outputFile))
    return(output_file)
      
def extractCigarSeq(sequence, phred, mapq, info):
    """
    Open writing file and return file pointer.
    An IOError is raised when the file cannot be opened.
        
    Keyword arguments:
    sequence -- the read bases aligned at a specific position
    phred    -- the base qualities for the same position
    mapq     -- the map qualities for the same position
    info     -- a dictionary containing the information about the reference base, the position, the chromosome and the number of bases
    """
        
    ref = info['ref']
    
    ## The dictionary that is going to contain all (Phred score, MPAQ) tuple for each base of the sequence
    letters = dict([('A', list()), ('C', list()), ('G', list()), ('T', list()), ('N', list()), ('O', list())])
    lettersKeys = letters.keys()
    
    sequence = sequence.upper()
    
    ## Offset used to obtain the Phred and MPAQ values
    asciiOffset = 33
    
    positionSeq = 0
    positionPhred = 0
    positionMapq = 0
    
    while(positionSeq < len(sequence)):
        currentData = sequence[positionSeq]
        if currentData == "," or currentData == ".":
            phredVal = ord(phred[positionPhred]) - asciiOffset
            mapqVal  = ord(mapq[positionMapq]) - asciiOffset
            letters[ref].append((phredVal, mapqVal))
            positionSeq += 1
            positionPhred += 1
            positionMapq += 1
        elif currentData in lettersKeys:
            phredVal = ord(phred[positionPhred]) - asciiOffset
            mapqVal  = ord(mapq[positionMapq]) - asciiOffset
            letters[currentData].append((phredVal, mapqVal))
            positionSeq += 1
            positionPhred += 1
            positionMapq += 1
        elif currentData == "+" or currentData == "-":
            letters['O'].append((-1, -1))
            res = match(r"[+-](\d+)", sequence[positionSeq:len(sequence)])
            indelLength = res.groups()[0]
            positionSeq = positionSeq + 1 + len(indelLength) + int(indelLength)
            ## Not change in positionPhred because there is not Phred value for an indel
            ## Not change in positionMapq because there is not Mapq value for an indel
        elif currentData == "$":
            positionSeq += 1
            ## Not change in positionPhred because there is not Phred value for the end of a sequence
            ## Not change in positionMapq because there is not Mapq value for the end of a sequence
        elif currentData == "^":
            ## The beginning is composed of 2 characters : "^" followed by the Mapq value
            positionSeq += 2
            ## Not change in positionPhred because there is not Phred value for the beginning of a sequence
            ## Not change in positionMapq because the Mapq value is in the sequence field, not in the Mapq field
        else:
            print("Problem with letter: " + currentData)
            print("in this cigar string: " + sequence)
            sys.exit(2)
    
    return(letters)
    
def parseMPileup(inputFile, outputPrefix):
    """
    Extract argument values as input by user.

    Keyword arguments:
    inputFile -- the name of the input file
    outputPrefix -- the prefix of the output files
    """
    
    iFile = open_read_file(inputFile)
    oFile = open_write_file(outputPrefix + ".txt")
    
    oExtraFile = dict()
    for phred in (15, 20, 30, 25):
        for mapq in (0, 1, 5, 10):
            oExtraFile[str(phred) + str(mapq)] = open_write_file(outputPrefix + "_" + str(phred) + "_" + str(mapq) + ".txt")
    
    for line in iFile:
        data = split("\t", line)
        
        info = dict([('chr', data[0]), ('pos', data[1]), ('ref', data[2]), ('NB', data[3])])

        lettersCount = extractCigarSeq(data[4], data[5], data[6], info);
        
        ## Write information will all aligned data
        oFile.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (info['chr'], info['pos'], info['NB'], len(lettersCount['A']), len(lettersCount['C']), len(lettersCount['G']), len(lettersCount['T']), len(lettersCount['O']), len(lettersCount['N'])))

        for phred in (15, 20, 30, 25):
            for mapq in (0, 1, 5, 10):
                outputF = oExtraFile[str(phred) + str(mapq)] 
                outputF.write("%s\t%s\t%s" % (info['chr'], info['pos'], info['NB']))
                for letterNow in ('A', 'C', 'G', 'T'):
                    lenLetter = len(filter(lambda g: g[0] < phred and g[1] < mapq, lettersCount[letterNow]))
                    outputF.write("\t%d" % (lenLetter))
                outputF.write("\n")
    
    for phred in (15, 20, 30, 25):
        for mapq in (0, 1, 5, 10):
            outputF = oExtraFile[str(phred) + str(mapq)]
            outputF.close()
             
    oFile.close()
    iFile.close()
    
if __name__ == "__main__":
    
    # Extract arguments. Message shown when the number of arguments is not coherent
    (inputFile, outputPrefix) = extractArguments()
    # Parsing pileup file
    parseMPileup(inputFile, outputPrefix)
