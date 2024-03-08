#!/bin/bash
#
# bash script to generate chip statistics
#
## Usage: bash chip_stats.sh <chip name> <R1.fastq.gz>
#
#
#SBATCH --job-name=chip_stats
#SBATCH --partition biochem,normal
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH -e logging/chip_stats-%j.err
#SBATCH -o logging/chip_stats-%j.out

###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments"
    exit 1
fi

chip_name=$1
fastq_gz=$2


###########################

# Internal variables:

# Output filename
outfile="$chip_name""_stats.txt"

###########################

# Begin Script


echo "Reading R1 data..."

tot=$(zgrep -c @ $fastq_gz)
rnap=$(zgrep -c TTATGCTATAATTATT $fastq_gz)
fid=$(zgrep -c CTTGGGTCCACAGGAC $fastq_gz)

rnap_frac=$(bc <<< "scale=3; $rnap/$tot")
fid_frac=$(bc <<< "scale=3; $fid/$tot")

echo "Making output file..."

echo Read statistics for $chip_name > $outfile
echo Total number of reads: >> $outfile
echo $tot >> $outfile
echo Number of reads containing RNAP initiation sites: >> $outfile
echo $rnap, $rnap_frac >> $outfile
echo Number of reads containing fiduciary marker: >> $outfile
echo $fid, $fid_frac >> $outfile

echo "Done!"
