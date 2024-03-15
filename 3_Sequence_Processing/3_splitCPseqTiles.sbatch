#!/bin/bash
#
# bash script to generate chip statistics
#
# Usage: sbatch 3_splitCPseqTiles.sbatch <withBCs_CPseq> <outdir>
#
#SBATCH --job-name=splitTiles
#SBATCH --partition biochem,normal
#SBATCH --time=12:00:00
#SBATCH -e logging/splitTiles-%j.err
#SBATCH -o logging/splitTiles-%j.out
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#################


# User input variables:
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments: <withBCs_CPseq> <outdir>"
    exit 1
fi

# CPseq=$1
CPseq=$1
outdir=$2

mkdir $outdir

module load python/2.7.13
source $py2env/bin/activate

python $IMAGING_DIR/CPscripts/splitCPseqIntoTiles.py $CPseq -o $outdir -s bottom
