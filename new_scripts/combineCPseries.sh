#!/bin/bash 
#
# bash script to call quantifyTilesDownstream.py using SLURM scheduler on Sherlock
#
# Usage: quantify_images.sbatch image_dir seq_dir roff_dir fluor_dir gv_path
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=combineCPseries
#################  
#a file for job output, you can check job progress
#SBATCH --output=combineCPseries.out
#################
# a file for errors from the job
#SBATCH --error=combineCPseries.err
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
#SBATCH --mem-per-cpu=32G
#################

# Setup env 
module load python/2.7.13
source $py2env/bin/activate

# How to call python script: 
# Usage: python combineCPseries.py [N tiles] [N concs] [all CPseries Dirs] 
#        [imagename template] [AllOutFilename] [anyRNAOutFilename] [CPannot.gz]

# Green Channel:
python $rnamap_scripts/new_scripts/combineCPseriesMissing.py 18 8 13_equilibrium0p91nM_green 14_equilibrium2p7nM_green 15_equilibrium8p2nM_green 16_equilibrium25nM_green 17_equilibrium74nM_green 18_equilibrium222nM_green 19_equilibrium666nM_green 20_equilibrium2000nM_green fid_anyRNA_JL4CY_tile%03d_Bottom_filtered.CPseries green_all.CPseries green_anyRNA.CPseries JL4CY_anyRNA.CPannot.gz

# Red Channel:
python $rnamap_scripts/new_scripts/combineCPseriesMissing.py 18 8 13_equilibrium0p91nM_red 14_equilibrium2p7nM_red 15_equilibrium8p2nM_red 16_equilibrium25nM_red 17_equilibrium74nM_red 18_equilibrium222nM_red 19_equilibrium666nM_red 20_equilibrium2000nM_red fid_anyRNA_JL4CY_tile%03d_Bottom_filtered.CPseries red_all.CPseries red_anyRNA.CPseries JL4CY_anyRNA.CPannot.gz



