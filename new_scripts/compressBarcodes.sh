#!/bin/bash
#################
#set a job name
#SBATCH --job-name=compressBarcodes
#################
#a file for job output, you can check job progress
#SBATCH --output=logging/compressBarcodes.out
#################
# a file for errors from the job
#SBATCH --error=logging/compressBarcodes.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=20:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=biochem,owners,normal
#SBATCH --qos=normal

module load python/2.7.13
source /home/groups/herschla/rna_map/scripts/new_scripts/venv_2_7_13/bin/activate


###### USER INPUTS ########
# CPseq file, should be a sorted, RNA-filtered CPseq such as anyRNA_withBCs_KD8TL.sort.CPseq
CPseq=$1
###########################

barcode_dir=$(pwd)+"Barcode_Output/"

mkdir -p $barcode_dir

python /home/groups/herschla/rna_map/barcode_tools/compressBarcodes_v8.py -i $CPseq -o $barcode_dir -q 28 -c 7 -C 5
