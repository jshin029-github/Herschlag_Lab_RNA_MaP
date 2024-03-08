#!/bin/bash
#SBATCH --job-nam=genCPannot
#SBATCH -e logging/genCPannot-%j.err
#SBATCH -o logging/genCPannot-%j.out
#SBATCH --partition biochem,owners
#SBATCH --mem-per-cpu=32G
#SBATCH --time=4:00:00

module load python/2.7.13
source /home/groups/herschla/rna_map/scripts/new_scripts/venv_2_7_13/bin/activate

##### To modify ######
cpseq="/scratch/groups/herschla/CVJNL_Experiments/seq/withBCs_CVJNL.CPseq.gz"    # withBCs___CPseq.gz
outfile="/scratch/groups/herschla/CVJNL_Experiments/seq/CVJNL.CPannot.gz"  # CPannot.gz


prebarfile="/scratch/groups/herschla/CVJNL_Experiments/all_barcode_assignments.csv"
######################

## OTHER VARIABLES ###
barcodecol="7"
######################

###### BEGIN SCRIPT #######
python /home/groups/herschla/rna_map/scripts/new_scripts/annotateClusters.py  $cpseq $prebarfile $outfile $barcodecol
