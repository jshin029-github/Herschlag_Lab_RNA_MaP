#!/bin/bash 
#
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=EFmax
#################  
#a file for job output, you can check job progress
#SBATCH --output=EnforceFmax.out
#################
# a file for errors from the job
#SBATCH --error=EnforceFmax.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=20:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=biochem,owners,normal
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#################

module load python/2.7.13
source $py2env/bin/activate


cpseriesfile="anyRNA_normed_AllRed.CPseries.gz"
cpfittedfile="anyRNA_normed_AllRed.CPfitted.gz"
cpannotfile="/scratch/groups/herschla/roy-test/Exp1_30mM_Mg_Lib4_20210218/seqData/JGFNV_anyRNA_sorted.CPannot.gz"
concentrationsfile="concentrations_corrected.txt"
cpvariantfile="anyRNA_normed_AllRed.CPvariant.gz"
fitfile="anyRNA_normed_AllRed.fmaxdist.p"
finalfile="anyRNA_normed_AllRed_Bootstrapped.CPvariant.gz"


python $rnamap_scripts/array_fitting_tools/bin/enforceFmax.py -b $cpseriesfile -v $cpvariantfile -m $fitfile -a $cpannotfile -x $concentrationsfile -out $finalfile
