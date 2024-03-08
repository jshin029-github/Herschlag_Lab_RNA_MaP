#!/bin/bash 
#
# bash script to call quantifyTilesDownstream.py using SLURM scheduler on Sherlock
#
# Usage: quantify_images.sbatch image_dir seq_dir roff_dir fluor_dir gv_path
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=makeHeatmap.sh
#SBATCH --mem=32G
#################  
#a file for job output, you can check job progress
#SBATCH --output=makeHeatmap.out
#################
# a file for errors from the job
#SBATCH --error=makeHeatmap.err
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


module load python/2.7.13
source $py2env/bin/activate
####### Change the following before running #######
template="combined_tile%03d.CPseries"           # Template name of CPseries filesc
ntiles=18                                       # No. of subtiles
color="green"                                   # green or red channel images to plot
####################################################

mkdir -p Heatmap_plots
python makeClusterHeatmap.py $template $ntiles $color
