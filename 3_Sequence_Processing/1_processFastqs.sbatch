#!/bin/bash
#
# bash script to combine the fastq files and output a CPseq
#
## Usage: sbatch chip_stats.sh <chip name> <R1.fastq.gz> <R2.fastq.gz>
#
#
#SBATCH --job-name=processFastq
#SBATCH --partition biochem,normal
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=6
#SBATCH --time=4:00:00
#SBATCH -e logging/processFastq-%j.err
#SBATCH -o logging/processFastq-%j.out

###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 3 ]
  then
    echo "Incorrect number of arguments"
    exit 1
fi

chip=$1
r1=$2
r2=$3


module load python/2.7.13
source /home/groups/herschla/rna_map/scripts/new_scripts/venv_2_7_13/bin/activate

python /home/groups/herschla/rna_map/scripts/array_fitting_tools/bin/processFastqUpdated.py -r1 $r1 -r2 $r2 -o "$chip".CPseq.gz -s /home/groups/herschla/rna_map/scripts/array_fitting_tools/bin/NUC.4.4
