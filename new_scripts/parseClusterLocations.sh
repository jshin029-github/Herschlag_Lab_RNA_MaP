#!/bin/bash 
#################
#set a job name  
#SBATCH --job-name=parseClusterLocations
#################  
#a file for job output, you can check job progress
#SBATCH --output=parseClusterLocations.out
#################
# a file for errors from the job
#SBATCH --error=parseClusterLocations.err
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

# Setup env 
module load python/2.7.13
source $py2env/bin/activate

cpseq="JGFNV_anyRNA.CPseq.gz"
outfile="JGFNV_anyRNA.Clus.gz"

python $rnamap_scripts/new_scripts/parseClusterLocations.py $cpseq $outfile

