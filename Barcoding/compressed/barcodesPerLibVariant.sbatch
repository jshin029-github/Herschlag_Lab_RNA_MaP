#!/bin/bash
#################
#set a job name
#SBATCH --job-name=barcodesPerLibVariant
#################
#a file for job output, you can check job progress
#SBATCH --output=barcodesPerLibVariant_%j.out
#################
# a file for errors from the job
#SBATCH --error=barcodesPerLibVariant_%j.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=48:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=biochem,normal
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#################

module load python/2.7.13
source /home/groups/herschla/rna_map/scripts/new_scripts/venv_2_7_13/bin/activate


CPseq=$1
lib=$2
outfile=$3

python /home/groups/herschla/rna_map/scripts/new_scripts/barcodesPerLibVariant.py -i $CPseq -l $lib -o $outfile -c 16
