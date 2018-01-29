mpileupTools
=====================

[![Build Status](https://travis-ci.org/adeschen/mpileupTools.svg?branch=master)](https://travis-ci.org/adeschen/mpileupTools)
[![codecov](https://codecov.io/gh/adeschen/mpileupTools/branch/master/graph/badge.svg)](https://codecov.io/gh/adeschen/mpileupTools)

This program parses a mpileup file to extract the coverage of each base at each position. 

## Usage 

python parseMPileup.py  -i <inputFile> -p <outputPrefix> [-s] [-h]

* inputFile = Output from samtools pileup
* outputPrefix = Prefix used on output files
* -s = When True, the position file in created in a separated file
* -h = Help

