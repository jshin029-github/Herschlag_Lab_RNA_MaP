#!/bin/bash
#
# bash script to tally the unique barcodes present in the CPseq file
#
# Usage: sbatch compressBarcodes.sbatch <CPseq> <out_dir>
#
#
#SBATCH --job-name=compressBarcodes
#SBATCH --partition biochem,normal
#SBATCH --mem-per-cpu=12G
#SBATCH --time=12:00:00
#SBATCH -e logging/compressBarcodes-%j.err
#SBATCH -o logging/compressBarcodes-%j.out

###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments: <sorted CPseq> <output directory>"
    exit 1
fi

###########################
# Internal variable

barcode_col=7 # 1-indexed
target_col=5 # 1-indexed

###########################


in=$1
out_dir="$2"/

###########################


module load python/2.7.13
source /home/groups/herschla/rna_map/scripts/new_scripts/venv_2_7_13/bin/activate

mkdir -p $out_dir

python /home/groups/herschla/rna_map/scripts/new_scripts/compressBarcodes_v8.py -i $in -o $out_dir -c $barcode_col -C $target_col
