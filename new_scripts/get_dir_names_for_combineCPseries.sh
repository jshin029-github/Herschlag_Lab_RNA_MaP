#!/bin/bash
#
# Bash script that is a helper script for batch_combineCPseries.sh
#
# How to call:
# Usage: sbatch batch_combineCPseries.sh <experiment pattern> <channel color>
#
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments"
    exit 1
fi
###########################


###### USER INPUTS ########
# prefix is the common prefix shared by each set of images (will likely be 'set')

prefix=$1
color=$2
###########################


###### BEGIN SCRIPT #######
dir_list=()
for s in *$prefix*$color*
do
	dir_list+=("$s")
done

echo "${dir_list[@]}"
echo "${#dir_list[@]}"
