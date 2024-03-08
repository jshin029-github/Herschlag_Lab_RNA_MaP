#!/bin/bash
#
# bash script to call genCPseries.sh for each image condition
#
## Usage: bash batch_genCPseries.sh <pattern>
## e.g.: bash batch_genCPseries.sh equilibrium

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
script_dir="/home/groups/herschla/rna_map/scripts/new_scripts"
###########################


###### USER INPUTS ########
# pattern is the common substring shared by each set of images (will likely be 'equilibrium')
pattern=$1
###########################


###### BEGIN SCRIPT #######
for im_dir in *$pattern*/
do
   	sbatch $script_dir/genCPseries.sh $im_dir $seq_dir $(pwd)
done

