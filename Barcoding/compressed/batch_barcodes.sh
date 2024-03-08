#!/bin/bash
#
# bash script to split CPseq file and call barcodesPerLibVariant.sh for each split
#
# Usage: bash barcodesPerLibVariant.sh <CPseq> <libchar> <barcodesPerLibVariant location>
#
#

###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 3 ]
  then
    echo "Incorrect number of arguments: <CPseq> <libchar> <barcodesPerLibVariant location>"
    exit 1
fi

###########################

to_split=$1
lib=$2
path_to_script=$3

###### BEGIN SCRIPT #######

split $to_split CPseq_ -n 10 --additional-suffix _split

for split_file in *split
do
    suffix="${split_file:6:2}"
    if [[ $(head -1 $to_split) != $(head -1 $split_file) ]];
    then
        head -1 $to_split > tmp_file
        cat "$split_file" >> tmp_file
        mv -f tmp_file "$split_file"
    fi
   	sbatch "$path_to_script"/barcodesPerLibVariant.sbatch "$split_file" $lib "barcoding_assignment_"$suffix
done

