#!/bin/bash
#
# bash script to sort CPseq file by barcode
#
# Usage: sbatch sortCPseq.sbatch <input CPseq> <sorted CPseq>
#
#SBATCH --job-name=sortCPseq
#SBATCH --partition biochem,normal
#SBATCH --mem-per-cpu=12G
#SBATCH --time=4:00:00
#SBATCH -e logging/sortCPseq-%j.err
#SBATCH -o logging/sortCPseq-%j.out


###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments: <input CPseq> <sorted CPseq>"
    exit 1
fi

###########################
# Internal variable

column=7 # note this does not work if CPseq has a filepath!

###########################


in=$1
out=$2

###########################

sort -k$column $in > $out
