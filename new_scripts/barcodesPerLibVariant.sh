#!/bin/bash
#################
#set a job name
#SBATCH --job-name=barcodesPerLibVariant
#################
#a file for job output, you can check job progress
#SBATCH --output=barcodesPerLibVariant.out
#################
# a file for errors from the job
#SBATCH --error=barcodesPerLibVariant.err
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

python /home/groups/herschla/rna_map/scripts/new_scripts/barcodesPerLibVariant.py -i "CPseq filename" -l "Libchar filename" -o "output filename"
