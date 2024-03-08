#!/bin/bash
#SBATCH --job-nam=genCPannot
#SBATCH -o genCPannot.out
#SBATCH --partition biochem,owners
#SBATCH --mem-per-cpu=32G
#SBATCH --time=4:00:00

module load python/2.7.13
source $py2env/bin/activate

##### To modify ######
cpseq=""    # CPseq.gz
outfile=""  # CPannot.gz
######################

prebarfile="/scratch/groups/herschla/rna_map/scripts/new_scripts/barcode_assignments.csv"
barcodecol="7"

python $rnamap_scripts/new_scripts/annotateClusters.py  $cpseq $prebarfile $outfile $barcodecol
