#!/bin/bash
#
# bash script to split CPseq file and call barcodesPerLibVariant.sh for each split
#
## Usage: bash barcodesPerLibVariant.sh <CPseq> <lib> <outfile>



###### BEGIN SCRIPT #######

head -1 barcode_assignment_aa > "all_barcode_assignments.csv"

for split_file in barcode_assignment_*
do
tail -n +2 -q $split_file >> "all_barcode_assignments.csv"
done

