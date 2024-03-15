#!/bin/bash
#SBATCH --job-nam=genCPannot
#SBATCH -e logging/genCPannot-%j.err
#SBATCH -o logging/genCPannot-%j.out
#SBATCH --partition biochem,owners
#SBATCH --mem-per-cpu=32G
#SBATCH --time=4:00:00

module load python/2.7.13
source $py2env/bin/activate

##### To modify ######
cpseq=""    # withBCs___CPseq.gz
outfile=""  # CPannot.gz


prebarfile="" # tsv mapping barcodes to variant number
######################

## OTHER VARIABLES ###
barcodecol="7"
######################

###### BEGIN SCRIPT #######
python $SCRIPT_DIR/annotateClusters.py  $cpseq $prebarfile $outfile $barcodecol
