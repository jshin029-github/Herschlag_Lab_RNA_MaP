#!/bin/bash
#
# bash script to call gen_dummy_CPfluor.sh for each image condition
#
## Usage: bash gen_dummy_CPfluor.sh <pattern>
## e.g.: bash gen_dummy_CPfluor.sh equilibrium

# raise error for incorrect usage
if [ $# -ne 1 ]
  then
    echo "Incorrect number of arguments"
    exit 1
fi
###########################



###### USER INPUTS ########
# pattern is the common substring shared by each set of images (will likely be 'equilibrium')
pattern=$1
###########################


###### BEGIN SCRIPT #######
for im_dir in *$pattern*/
do
  for file in "$im_dir/"*
    do
    CPfluorFile="CPfluor/"${file%.tif}".CPfluor"
    if [ ! -f $CPfluorFile ]
    then
      touch $CPfluorFile
    fi
  done
done
