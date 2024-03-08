#!/bin/bash
#
# bash script to call compileImages.py using SLURM scheduler on Sherlock
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=compileImages
#################
#a file for job output, you can check job progress
#SBATCH --output=logging/compileImages.out
#################
# a file for errors from the job
#SBATCH --error=logging/compileImages.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=15:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=biochem,owners,normal
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#tasks to run per node; a "task" is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use "--ntasks-per-node=1 --cpus-per-task=16" instead
# Sherlock nodes have 16 cpus. For some reason, you can request more than that on 'owners' partition, but not on others.
# It can't give you more obviously, but it's interesting that it lets you ask on one partion but not another.
# Note: On Sherlock2, most of the nodes we have access to have 20 cores.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#################

# How to call:
# Usage: sbatch compileImages.sh <pattern> <conditions.txt> <Ntiles> <Image Directory>
#
# raise error for incorrect usage
if [[ $# -ne 2 && $# -ne 4 ]]
  then
    echo "Incorrect number of arguments"
    exit 1
fi
#####################


###### USER INPUTS ########
# pattern is the common pattern shared by each set of images (will likely be 'equilibrium')
pattern=$1
conditions_file=$2
if [[ $# == 2 ]]
  then
    Ntiles=18
	out_dir="Images"
  else
  	Ntiles=$3
  	out_dir=$4
fi

###########################


###### BEGIN SCRIPT #######
# Setup env
module load python/3.6.1
source /home/groups/herschla/rna_map/scripts/new_scripts/venv_3_6_1/bin/activate

# How to call python script:
# python3 compileImages.py <pattern> <conditions.txt> <Ntiles> <Image Directory>

python3 /home/groups/herschla/rna_map/scripts/new_scripts/compileImages.py $pattern $conditions_file $Ntiles $out_dir
