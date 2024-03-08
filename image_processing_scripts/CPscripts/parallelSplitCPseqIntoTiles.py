#!/usr/bin/env python

#written by Curtis Layton Nov 2013
#splits a CPseq file (or other file where the first token of every line is a cluster id) into an idividual file for every tile
#the cluster id (first field of the tab-delimited file) should be of the format "<machine id>:<run index>:<flowcell id>:<lane #>:<tile id>:<x coord>:<y coord>"
# Updates by John Shin April 2022


import os
import sys
import argparse
import re
from joblib import Parallel, delayed

import CPlibs

def getNumLinesInFile(filename):
    n_lines = 0
    with open(filename) as f:
        for line in f: n_lines += 1
    return n_lines

def checkValidOutputFilename(outputFilename):
    with open(outputFilename,"w") as f: #attempt to open for writing.  will raise error if fails
        pass
    try:
        os.unlink(outputFilename) #remove the file that was temporarily created
    except OSError:
        pass

def getOutputFilename(inFullPath,currTile,currSide,outputPath):
    #split the input filename into parts
    (inPath,inFilename) = os.path.split(inFullPath)
    (inBasename,ext) = os.path.splitext(inFilename)
    
    if currSide == 1:
        currSideString = 'Top'
    elif currSide == 2:
        currSideString = 'Bottom'
    
    outputFilename = os.path.join(outputPath,inBasename + '_tile' + currTile + '_' + currSideString + ext) #append tile designation
    return outputFilename

def helperFcnBoth(currInputLine,inputFilename,outputPath):

        currID = re.split('[\t ]',currInputLine)[0] #get the cluster ID (first token of the line)

        (flowcellSN,currSide,currSwath,currTile) = CPlibs.parseClusterID(currID)

        #NOTE: the swath is ignored in this implementation.

        currTile = '{:03}'.format(int(currTile)) #format currTile in 3-digit string in anticipation of porting the pipeline to sytems with more tiles

        return getOutputFilename(inputFilename,currTile,currSide,outputPath),currInputLine

def helperFcnTop(currInputLine,inputFilename,outputPath):

        currID = re.split('[\t ]',currInputLine)[0] #get the cluster ID (first token of the line)

        (flowcellSN,currSide,currSwath,currTile) = CPlibs.parseClusterID(currID)

        #NOTE: the swath is ignored in this implementation.

        currTile = '{:03}'.format(int(currTile)) #format currTile in 3-digit string in anticipation of porting the pipeline to sytems with more tiles

        if currSide == 1:
            return getOutputFilename(inputFilename,currTile,currSide,outputPath),currInputLine
        else:
            return None,None

def helperFcnBottom(currInputLine,inputFilename,outputPath):

        currID = re.split('[\t ]',currInputLine)[0] #get the cluster ID (first token of the line)

        (flowcellSN,currSide,currSwath,currTile) = CPlibs.parseClusterID(currID)

        #NOTE: the swath is ignored in this implementation.

        currTile = '{:03}'.format(int(currTile)) #format currTile in 3-digit string in anticipation of porting the pipeline to sytems with more tiles

        if currSide == 2:
            return getOutputFilename(inputFilename,currTile,currSide,outputPath),currInputLine
        else:
            return None,None


### MAIN ###
def main(rawArgs, **keywordParams):
    
    #set up command line argument parser
    parser = argparse.ArgumentParser(description="Splits a .CPseq file into individual .CPseq files for each tile. Output filenames are derived automatically from the input filename.")
    
    parser.add_argument('inputFilename', help='input .CPseq file to be split into tile-specific .CPseq files')
    parser.add_argument('-o','--output_dir', help='output directory. Default is the same directory as the input file.')
    parser.add_argument('-s','--side', help='specify if tiles will be output from the top surface (default), bottom surface, or both.', choices=['top','bottom','both'], default='top')
    parser.add_argument('-c','--num_cores', help='number of cores to use for parallelization',type=int,default=-1)
    outputExtension = '.CPseq'
    
    #parse command line arguments
    args = parser.parse_args(rawArgs[1:])

    #get the number of lines in the input file
    numLinesInInputFile = getNumLinesInFile(args.inputFilename)
    
    #get output directory
    if args.output_dir is None:
        (outputPath, inputFilenameStripped) = os.path.split(args.inputFilename) #use the same path as the input file
    else:
        outputPath = args.output_dir
    
    
    #open the input file and begin reading
    print 'Processing ' + str(numLinesInInputFile) + ' entries from input file: ' + args.inputFilename
    with open(args.inputFilename) as inputFile:

        if args.side.lower() == 'both':
            res = (Parallel(n_jobs=args.num_cores,verbose=10)
              (delayed(helperFcnBoth)(currInputLine,args.inputFilename,outputPath) for currInputLine in inputFile))

        elif args.side.lower() == 'top':
            res = (Parallel(n_jobs=args.num_cores,verbose=10)
              (delayed(helperFcnTop)(currInputLine,args.inputFilename,outputPath) for currInputLine in inputFile))

        elif args.side.lower() == 'bottom':
            res = (Parallel(n_jobs=args.num_cores,verbose=10)
              (delayed(helperFcnBottom)(currInputLine,args.inputFilename,outputPath) for currInputLine in inputFile))


    print(res)

    file_dict = {}
    for k,v in res:
        file_dict.setdefault(k, []).append(v)

    # create and fill files
    for filename in file_dict.keys():
        open(filename,'w+').close()

        with open(filename,'a') as f:
            for line in file_dict[file]:
                f.write(line)

    if ('from_command_line' in keywordParams):
        #function was called from the command line
        return 0 #return exit value

if __name__=='__main__':
    sys.exit(main(sys.argv, from_command_line=True))


