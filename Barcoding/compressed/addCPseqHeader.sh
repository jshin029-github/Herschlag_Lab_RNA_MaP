#!/bin/bash
#
# bash script to add a header to the .unq file
#
# Usage: bash addCPseqHeader.sh <unq> <CPseq>
#

###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments: <unq> <CPseq>"
    exit 1
fi

###########################

in=$1
out=$2

###########################

( echo -e "sequence\tbarcode\tnumber of sequences with barcode\tvoting results\tsame len flag\tmode target sequence length\tmean barcode quality score"; cat $in ) > $out
