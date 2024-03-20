#!/usr/bin/env python

"""
Get the registration offsets for a directory of images relative to a 
directory of CPseqs.

Note: Python 3

Inputs:
   directory of images to register
   directory of CPseq files from which to pick clusters
   

Outputs:
   registration offsets

Ben Ober-Reynolds
"""

import os
import sys
import argparse
import re
import uuid
import subprocess
import time
from collections import OrderedDict
from joblib import Parallel, delayed


##### Gloval vars #####
offset_scale_x = -3.7
offset_scale_y = 3.7


def main():  
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='script for isolating specific \
        clusters from fastq files, based on a set of CPseq files')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-id', '--image_directory', required=True,
        help='directory containing images to register')
    group.add_argument('-sd', '--CPseq_dir', required=True,
        help='directory containing CPseq files')
    group.add_argument('-gv','--global_vars_path', required=True, 
        help='path to the directory in which the "GlobalVars.m" parameter file \
        for the run can be found')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('--nofile', action='store_true',
        help='use flag to prevent output of offsets file. Results printed to STDOUT only.')
    group.add_argument('-f', '--filters', type=str, nargs='+',
        help='Which filter(s) from the CPseq files to register. Default is "FIDUCIAL"')
    group.add_argument('-ds', '--data_scaling', type=str, default='MiSeq_to_TIRFStation1',
        help='Data scaling for registration. Default is "MiSeq_to_TIRFStation1"')
    group.add_argument('-od', '--output_directory',
        help='output directory for registration offsets (default is original \
            image_directory)')
    group.add_argument('-op', '--output_prefix', type=str, default='registration_offsets',
        help='output prefix for registration offsets file (default is "registration_offsets")')
    group.add_argument('-n', '--num_cores', type=int, default=18,
        help='number of cores to use (should be same as number of image \
            files)')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()
    numCores = args.num_cores

    # Pre-defined variables, constants, and settings
    image_extension = 'tif'
    CPseq_extension = 'CPseq'

    # Check directories and output files
    image_dir = args.image_directory
    if not os.path.isdir(image_dir):
        print("Error: invalid image directory selection. Exiting...")
        sys.exit()

    CPseq_dir = args.CPseq_dir
    if not os.path.isdir(CPseq_dir):
        print("Error: invalid CPseq directory selection. Exiting...")
        sys.exit()

    output_dir = args.output_directory
    if not output_dir:
        output_dir = image_dir

    # Check global vars:
    globalVarsFilename = os.path.join(args.global_vars_path, 'GlobalVars.m')
    if not os.path.isfile(globalVarsFilename):
        print('ERROR: GlobalVars.m file not found in directory "' + args.global_vars_path + '".  Aborting')
        sys.exit()

    # Gather image files:
    print("Finding image files in directory {}".format(image_dir))
    image_list = find_files_in_directory(image_dir, 
        extensionList=[image_extension])

    # Gather CPseq files:
    print("Finding CPseq files in directory {}".format(CPseq_dir))
    CPseq_list = find_files_in_directory(CPseq_dir, 
        extensionList=[CPseq_extension])

    # Make tile dict of each tile list
    image_tile_dict = make_tile_dict(image_list)
    CPseq_tile_dict = make_tile_dict(CPseq_list)


    # Pick filters to use
    filter_list = args.filters
    if not filter_list:
        filter_list = ['FIDUCIAL']
    print("Registering images using the following filter(s): {}".format(filter_list))

    # Run registration

    registration_results = []
    if numCores > 1:
        print("Getting registration offsets for {} image files on {} cores...".format(
            len(image_list), numCores))
        registration_results = (Parallel(n_jobs=numCores, verbose=10)\
            (delayed(checkRegistrationOffset)(
                tile, image_tile_dict[tile], CPseq_tile_dict[tile], args.data_scaling, 
                filter_list, args.global_vars_path) for tile in set(image_tile_dict.keys()) & set(CPseq_tile_dict.keys())            )
            )
    else:
        print("Getting registration offsets for {} image files on one core...".format(
            len(image_list)))
        registration_results = [checkRegistrationOffset(
            tile, image_tile_dict[tile], CPseq_tile_dict[tile], args.data_scaling, 
                filter_list, args.global_vars_path) for tile in set(image_tile_dict.keys()) & set(CPseq_tile_dict.keys())]

    # Format and output results:

    # First sort by tile number:
    registration_results.sort(key=lambda x: x[0])

    # Save a file if indicated
    if not args.nofile:
        with open(output_dir + '/' + args.output_prefix + '.txt', 'w') as f:
            f.write('x\ty\n')
            for offset in registration_results:
                f.write("{}\t{}\n".format(round(offset_scale_x*offset[2], 3), round(offset_scale_y*offset[1], 3)))
            # Add zeros at the end of this file, since the imaging stating expects 19 tiles
            f.write("0\t0")
    # print to stdout:
    print("Found offsets:")
    print("\tx\ty")
    for offset in registration_results:
        print("tile {}:\t{}\t{}".format(offset[0], round(offset_scale_x*offset[2], 3), round(offset_scale_y*offset[1], 3)))



def find_files_in_directory(dirPath, extensionList=None, 
                            excludedExtensionList=None):
    """
    Locate files in a given directory path. Optionally, desired files are 
    identified as matching one of the extension types provided in 
    'extensionList'
    Input: 
        dirPath (str) - path to directory
        extensionList (list) - list of acceptable extensions
        excludedExtensionList (list) - list of unacceptable extensions
    Output: 
        fileList (list) - list of found files (with path)
    """
    def extension_match(filename, extensionList=None):
        # from CPlibs
        if extensionList is not None:
            for currExt in extensionList:
                if filename.lower().endswith(currExt.lower()):
                    return True
        return False

    dirList = os.listdir(dirPath)
    fileList = []
    for currFilename in dirList:
        if (extension_match(currFilename, extensionList) 
        and not extension_match(currFilename, excludedExtensionList)): 
            fileList.append(dirPath+currFilename)
    if len(dirList) == 0:
        print('\tNONE FOUND')
    else:
        for filename in fileList:
            print("found:\t\t{}".format(filename))
        return fileList


def get_tile_number_from_filename(inFilename):
    """
    Extract the tile number from a provided filename based on the presence of
    'tile###'
    Input: filename (string)
    Output: three digit tile number (string)
    """
    # from CPlibs
    (path,filename) = os.path.split(inFilename) #split the file into parts
    (root,ext) = os.path.splitext(filename)
    matches = re.findall('tile[0-9]{1,3}',root.lower())
    tileNumber = ''
    if matches != []:
        tileNumber = '{:03}'.format(int(matches[-1][4:]))
    return tileNumber


def make_tile_dict(fileList):
    """
    Make a dictionary of files keyed by tile number.  
    Input: list of files containing tile numbers
    Output: dictionary of file names keyed by tile number
    """
    fileDict = {}
    for f in fileList:
        tile = get_tile_number_from_filename(f)
        if tile == '':
            print("Error: no tile number in file: "+ f)
            sys.exit()
        else:
            if tile in fileDict:
                print("Error: multiple files per tile")
                sys.exit()
            fileDict[tile] = f
    return fileDict


def checkRegistrationOffset(tile, image_file, CPseq_file, data_scaling, filter_list, global_vars_path):
    """
    Run the matlab script 'checkTileRegistration.m'
    Return the raw offsets calculated.
    """
    filter_string = "{{{}}}".format(",".join("'" + x + "'" for x in filter_list))

    matlabFunctionCallString = "checkTileRegistrationV2('{0}','{1}','{2}', {3});".format(
        CPseq_file, image_file, data_scaling, filter_string)

    logstring = spawnMatlabJob(matlabFunctionCallString, global_vars_path)
    # Parse logstring for relevant information:
    center_pos_offsets = ""
    log_lines = logstring.split('\n')
    parens_pat = re.compile(("\(.+?\,.+?\)"))

    for line in log_lines:
        matches = re.findall(parens_pat, line)
        if matches:
            center_pos_offsets = matches[1]

    if center_pos_offsets == "":
        # If registration not found, just return zeros.
        return (int(tile), 0, 0)
    offset_y, offset_x = [float(x) for x in center_pos_offsets[1:-1].split(',')]

    return (int(tile), offset_y, offset_x)



def spawnMatlabJob(matlabFunctionCallString,globalVarsPath):
    """
    Adapted from CPlibs.py 

    """
    try:
        #construct the command-line matlab call 
        functionCallString =                      "try,"
        functionCallString = functionCallString +     "addpath('{0}');".format(globalVarsPath) #placeholder TEMP DEBUG CHANGE
        functionCallString = functionCallString +     matlabFunctionCallString + ';'
        functionCallString = functionCallString + "catch e,"
        functionCallString = functionCallString +     "disp(getReport(e,'extended'));"
        functionCallString = functionCallString + "end,"
        functionCallString = functionCallString + "quit;"
    
        logFilename = 'matlabProcess_' + str(uuid.uuid4()) + str(time.time()) + '.tempLog' #timestamped logfile filename
    
        cmdString ='matlab -nodesktop -nosplash -singleCompThread -r "{0}"'.format(functionCallString)
        cmdString = cmdString + ' 1>> {0}'.format(logFilename)
        cmdString = cmdString + ' 2>> {0}'.format(logFilename)
       
        print('issuing subprocess shell command: ' + cmdString)
       
        returnCode = subprocess.call(cmdString,shell=True) #execute the command in the shell
        returnCode2 = subprocess.call('stty sane',shell=True) #matlab messes up the terminal in a weird way--this fixes it 
    
        #read log file into a string
        try:
            with open(logFilename) as logFilehandle:
                logString = logFilehandle.read()
            # delete logfile
            try:
                os.unlink(logFilename)
            except OSError:
                pass
        except IOError:
            logString = 'Log file not generated for command "' + functionCallString + '".'
    
        # return log
        return logString
    except Exception as e:
        return 'Python exception generated in spawnMatlabJob: ' + e



if __name__ == '__main__':
    main()
