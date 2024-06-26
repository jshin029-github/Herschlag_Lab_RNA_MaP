#!/usr/bin/env python

'''
Adapted from Anthony Ho, 6/29/2014 by Lauren Chircus
lchircus@stanford.edu
Stanford Univeristy
August 4, 2014

Based on Jason Buenrostro's script pyCompressBBsv2.py

This script takes in a merged fastq file (CPseq formatted or other format with all the data for one cluster per line) that has been sorted by barcode sequence. It merges all the sequences in a barcode block and then has them vote on the consensus sequence.

This script does NOT filter the barcodes at all - filterBarcodesAndAlignToReference does all the filtering in a second step.
'''


''' Import Modules'''
import os, sys
import numpy as np
from collections import Counter
from optparse import OptionParser
import datetime


''' Functions '''

# Writes a single line with a newline to the specified file
# assumes that file is already open
def writeLine(myfile, line):
    myfile.write(line+'\n')

# Prints a line to the console and writes it to the specified file
def printLine(myfile, line):
    print line
    writeLine(myfile, line)


## Phred Q score calculator, default is Phred-33
def qScore(sequence,phred=33):
    if phred == 64:
        return [ord(n)-64 for i,n in enumerate(sequence)]
    else:
        return [ord(n)-33 for i,n in enumerate(sequence)]


## Calculate the average Phred Q score of a given sequence, defualt is Phred-33
def avgQScore(sequence,phred=33):
    if phred == 64:
        l = [ord(n)-64 for i,n in enumerate(sequence)]
    else:
        l = [ord(n)-33 for i,n in enumerate(sequence)]
    return np.mean(l)


## Finding the consensus sequence of a barcode block
## !!! Needs to be refined !!!
def consensusVoting(r1Block,Q1Block,degeneracy):
    ## Find the consensus sequence
    consensus = ""
    charArray = np.array(map(list, r1Block[:]))
    bases = "ACGT"
    for i in range(0,len(charArray[0])):
        ## Base array
        baseArray = charArray[:,i].tolist()
        ## Count bases and vote
        baseCount = (baseArray.count('A'), baseArray.count('C'), baseArray.count('G'), baseArray.count('T'))
        vote = np.argmax(baseCount)
        consensus += bases[vote]
    return consensus


def main():

    ## Get options and arguments from command line
    parser = OptionParser()
    parser.add_option('-i', dest="inFile", help="Merged file to be analyzed")
    parser.add_option('-o', dest="outFolder", help="directory in which to write output")
    parser.add_option("-q", dest="avgQScoreCutoff", type=int, default=28, help="cutoff for average q score of the read")
    parser.add_option("-c", dest='colBarcode', help="column number in which the barcodes appear (1 indexed)", action='store', type='int')
    parser.add_option("-C", dest='colTarget', help="column number in which the target sequences appear (1 indexed)", action='store', type='int')
    options, arguments = parser.parse_args()

    # return usage information if no argvs given
    if len(sys.argv)==1:
        os.system(sys.argv[0]+" --help")
        sys.exit()

    # Making sure all mandatory options appeared.
    mandatories = ['inFile', 'colBarcode', 'colTarget', 'outFolder']
    for m in mandatories:
        if not options.__dict__[m]:
            print "mandatory option is missing\n"
            parser.print_help()
            exit(-1)

    options.colBarcode = options.colBarcode-1
    options.colTarget = options.colTarget-1


    ## Initiate output files
    inFile_name = os.path.splitext(os.path.basename(options.inFile))[0]
    outputFileName =  options.outFolder+inFile_name+"_compressedBarcodes.unq"
    statFileName =  options.outFolder+inFile_name+"_compressedBarcodes.stat"

    logFileName =   options.outFolder+"compressBarcodes_v8.log"
    logFile = open(logFileName,'w')

    ## Write to log
    writeLine(logFile, 'Analysis performed ' + str(datetime.datetime.now()) + ' using compressBarcodes_v8.py')
    writeLine(logFile, 'Barcodes compressed from file: ' + options.inFile)
    writeLine(logFile, 'Average sequence quality score cutoff used: ' + str(options.avgQScoreCutoff))

    writeLine(logFile, '\nCompressed barcodes in output file ' + outputFileName)
    writeLine(logFile, 'Columns:')
    writeLine(logFile, '  (1) consensus sequence')
    writeLine(logFile, '  (2) barcode')
    writeLine(logFile, '  (3) number of sequences with barcode')
    writeLine(logFile, '  (4) voting results (fraction sequences that perfectly match consensus)')
    writeLine(logFile, '  (5) 1 if sequences are all same length, 2 of sequences are not,')
    writeLine(logFile, '  (6) mode target sequence length')
    writeLine(logFile, '  (7) mean barcode quality score')

    writeLine(logFile, 'number of clusters per barcode in ' + statFileName)

    ## Initialization
    numSeqs = 0
    lineCount = 0

    r1Block = []
    Q1Block = []
    aQbBlock = []

    # Stat on reads
    goodSeqs = 0
    crapSeqs = 0

    # Stat of barcodes
    goodBC = 0

    # Stat for merging
    perfectSeqs = 0
    sameLenSeqs = 0
    diffLenSeqs = 0

    ## Initial counting
    with open(options.inFile,"r") as r:
        for line in r:
            numSeqs += 1
            if numSeqs == 1:
                lastBC = line.rstrip().split('\t')[options.colBarcode]

    ## Going through the CPseq file
    with open(options.inFile, "r") as r, open(outputFileName, "w") as wo, open(statFileName, "w") as ws:
        for line in r:

            # Reading line by line
            lineCount += 1
            if lineCount % 10000 == 0:
                print "Processing the "+str(lineCount) +"th sequence"
            seqLine = line.rstrip().split('\t')
            if len(seqLine) > 7:
	        r1 = seqLine[options.colTarget]
                Q1 = seqLine[options.colTarget + 1]
                BC = seqLine[options.colBarcode]
                Qbc = seqLine[options.colBarcode +1]
	    else:
		pass
	    #print(BC)
	    #print(Qbc)
            # At the beginning of each barcode block,
            # and do it for the very last line instead of the very first line
            if BC != lastBC or lineCount == numSeqs:

                # Add in the very last line to the barcode block
                # Append sequence and q-scores to the barcode block
                # Check the sequences before putting into block
                if lineCount == numSeqs:
                    if avgQScore(Q1) <= options.avgQScoreCutoff:
                        crapSeqs += 1
                    else:
                        r1Block.append(r1)
                        Q1Block.append(Q1)
                        aQbBlock.append(avgQScore(Qbc))

                # Analyze barcode block here
                degeneracy = len(r1Block)
		if degeneracy > 0:
                    r1Block = np.array(r1Block)
                    Q1Block = np.array(Q1Block)
                    aQbBlock = np.array(aQbBlock)
                   

		    # Dealing with statistics
		    ws.write(str(degeneracy)+"\n")
                    goodBC += 1
                    goodSeqs += degeneracy

                    # Merging
                    # !!! Needs to be refined !!!
                    # Find perfect matches within a block first:
 		    if r1Block[:].tolist() == [r1Block[0]]*degeneracy:
                        perfectSeqs += 1

                        wo.write(r1Block[0]+"\t"+lastBC+"\t"+str(degeneracy)+"\t1\t1\t"+str(len(r1Block[0]))+'\t'+str(np.mean(aQbBlock))+"\n")

                    # Then deal with sequence blocks of the same length of sequences:
                    elif max(r1Block,key=len) == min(r1Block,key=len):
                        sameLenSeqs += 1
                        consensusSeq = consensusVoting(r1Block,Q1Block,degeneracy)

                        fractionMatchingConsensus = float((len(filter(lambda seq: seq==consensusSeq, r1Block[:].tolist()))))/float(degeneracy)

                        wo.write(consensusSeq+"\t"+lastBC+"\t"+str(degeneracy)+"\t"+str(fractionMatchingConsensus)+"\t1\t"+str(len(r1Block[0]))+'\t'+str(np.mean(aQbBlock))+"\n")
                    # If not, remove the sequences that are not of the same length:
                    else:
                        lenList = [len(x) for x in r1Block]
                        lenCount = Counter(lenList)
                        modeLen = lenCount.most_common(1)[0][0]
                        truncatedR1Block = [x for x in r1Block if len(x) == modeLen]
                        diffLenSeqs += 1
                        consensusSeq = consensusVoting(truncatedR1Block,Q1Block,degeneracy)

                        fractionMatchingConsensus = float((len(filter(lambda seq: seq==consensusSeq, r1Block[:].tolist()))))/float(degeneracy)

                        wo.write(consensusSeq+"\t"+lastBC+"\t"+str(degeneracy)+"\t"+str(fractionMatchingConsensus)+"\t2\t"+str(modeLen)+'\t'+str(np.mean(aQbBlock))+"\n")
                # Initialize for the next barcode block
                r1Block = []
                Q1Block = []
                aQbBlock = []

            # Append sequence and q-scores to the barcode block
            # Check the sequences before putting into block
            if avgQScore(Q1) <= options.avgQScoreCutoff:
                crapSeqs += 1
            else:
                r1Block.append(r1)
                Q1Block.append(Q1)
                aQbBlock.append(avgQScore(Qbc))

            # Make the current barcode the new last barcode
            lastBC = BC

    ## Printing summary
    printLine(logFile, "\nTotal number of sequences analyzed: " + str(numSeqs) + " (100%)")
    printLine(logFile, "Number of sequences passing filter: " + str(goodSeqs) + " ("+str(round(float(goodSeqs)/numSeqs*100,2)) + "%)")
    printLine(logFile, "Number of low quality sequences: " + str(crapSeqs) + " ("+str(round(float(crapSeqs)/numSeqs*100,2)) + "%)")

    printLine(logFile, "\nTotal number of barcodes: " + str(goodBC) + " (100%)")

    ## Testing
    printLine(logFile, "\nNumber of unique sequence with same sequences within a barcode block: " + str(perfectSeqs) + " ("+str(round(float(perfectSeqs)/goodBC*100,2)) + "%)")
    printLine(logFile, "Number of unique sequences with different sequences of all the same length within a barcode block: " + str(sameLenSeqs) + " ("+str(round(float(sameLenSeqs)/goodBC*100,2)) + "%)")
    printLine(logFile, "Number of unique sequences with different sequences of different length: " + str(diffLenSeqs) + " ("+str(round(float(diffLenSeqs)/goodBC*100,2)) + "%)")

    r.close()
    wo.close()
    ws.close()
    logFile.close()

    return 1

if __name__ == "__main__":
    main()

