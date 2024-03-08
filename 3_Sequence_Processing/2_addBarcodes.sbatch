#!/bin/bash
#
# bash script to add barcodes to the input CPseq
#
## Usage: sbatch 2_addBarcodes.sbatch <CPseq.gz> <bc>
#
#
#SBATCH --job-name=addBCs
#SBATCH --partition biochem,normal
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=6
#SBATCH --time=4:00:00
#SBATCH -e logging/addBCs-%j.err
#SBATCH -o logging/addBCs-%j.out

###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments: <CPseq.gz> <bc>"
    exit 1
fi

CPseq_gz=$1
barcode_nts=$2

###########################
# Internal variable

CPseq=${CPseq_gz%.gz}
outfile=withBCs_$CPseq # note this does not work if CPseq has a filepath!

###########################
# Call function

module load python/3.6.1
source /home/groups/herschla/rna_map/scripts/new_scripts/venv_3_6_1/bin/activate

gunzip $CPseq_gz

python3 /home/groups/herschla/rna_map/scripts/new_scripts/addBarcodeLines.py -i $CPseq -o $outfile -b $barcode_nts
