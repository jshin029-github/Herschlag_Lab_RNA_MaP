#!/bin/bash
#
# Quantify array images on sherlock
#
# Ben Ober-Reynolds, boberrey@stanford.edu, 20191105
#
# Usage: bash quantifyImagesFIDOnly.sh <color>
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=quantifyImagesFIDOnly.sh
#################
#a file for job output, you can check job progress
#SBATCH --output=logging/quantifyImagesFIDOnly-%j.out
#SBATCH --error=logging/quantifyImagesFIDOnly-%j.err

#
# raise error for incorrect usage
if [ $# -ne 1 ]
  then
    echo "Incorrect number of arguments"
    exit 1
fi
###########################


######## TO MODIFY ########
# split_tile cpseq files
seq_dir=""
# raise error if left unmodified
if [ -z "$seq_dir" ]
  then
    echo "Please specify a CPseq directory"
    exit 1
fi
###########################


###### OTHER VARIABLES ####
rnamap_scripts="/home/groups/herschla/rna_map/scripts"

roff_dir="./Roff"
fluor_dir="./CPfluor"
gv_path=$IMAGING_DIR
log_dir="./Logfile"
###########################


##### USER INPUTS #########
# color is the common color suffix shared by each set of images ('green' or 'red')
color=$1
###########################


###### BEGIN SCRIPT #######
for s in *$color*/
do
    echo "Current image directory = $s"
    mkdir -p $roff_dir/$s
    mkdir -p $fluor_dir/$s
    mkdir -p $log_dir/$s
    sbatch $IMAGING_DIR/quantify_images_FID.sbatch $s $seq_dir $roff_dir/$s $fluor_dir/$s $gv_path $log_dir/$s
done
