#!/bin/bash
#
# bash script to call plotFluorQCBySublib.py using SLURM scheduler on Sherlock
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=plotFluorQCBySublib
#################
#a file for job output, you can check job progress
#SBATCH --output=plotFluorQCBySublib.out
#################
# a file for errors from the job
#SBATCH --error=plotFluorQCBySublib.err
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
#SBATCH --cpus-per-task=18
#################

# How to call:
# Usage: sbatch plotFluorQCBySublib.sh <experiment pattern>
#
# raise error for incorrect usage
if [ $# -ne 1 ]
  then
    echo "Incorrect number of arguments"
    exit 1
fi
#####################


###### USER INPUTS ########
# pattern is the common pattern shared by each set of images (will likely be 'equilibrium')
pattern=$1
###########################


###### OTHER VARIABLES ####
Libchar="/scratch/groups/herschla/rna_map/scripts/new_scripts/Lib4_Final_Sorted.libchar"
###########################



###### BEGIN SCRIPT #######
# Setup env
module load python/2.7.13
source $py2env/bin/activate

green="$pattern"*green_all.CPseries.gz
red="$pattern"*red_all.CPseries.gz
cond="$pattern"*conditions.txt
outdir="$pattern"_SublibFluorQCPlots

mkdir -p $outdir

# How to call python script:
# Usage: python plotFluorQCBySublib.py [all.annot.green.CPseries] [all.annot.red.CPseries] [conditions.txt] [libchar] [plotdir]

python $rnamap_scripts/new_scripts/plotFluorQCBySublib.py $green $red $cond $Libchar $outdir
