#!/bin/bash
#
# bash script to add barcodes to the input CPseq
#
## Usage: sbatch 2a_genAnyRNA.sbatch <CPseq> <chip name>
#
#
#SBATCH --job-name=genAnyRNA
#SBATCH --partition biochem,normal
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=6
#SBATCH --time=4:00:00
#SBATCH -e logging/genAnyRNA-%j.err
#SBATCH -o logging/genAnyRNA-%j.out

###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments: <CPseq> <chip name>"
    exit 1
fi

CPseq=$1
chip_name=$2

###########################
# Internal variable

outfile=withBCs_"$chip_name".CPseq # note this does not work if CPseq has a filepath!

RNAP="TTATGCTATAATTATT"

###########################

grep $RNAP $CPseq >  anyRNA_"$outfile"

