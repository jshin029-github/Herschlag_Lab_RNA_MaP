#!/bin/bash
#
# bash script to call combineCPseriesMissing.py using SLURM scheduler on Sherlock
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=combineCPseries
#################
#a file for job output, you can check job progress
#SBATCH --output=logging/combineCPseries-%j.out
#################
# a file for errors from the job
#SBATCH --error=logging/combineCPseries-%j.err
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
#SBATCH --mem-per-cpu=32G
#################

# How to call:
# Usage: sbatch batch_combineCPseries.sh <experiment pattern> <channel color>
#
# raise error for incorrect usage
if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments"
    exit 1
fi
###########################


######## TO MODIFY ########
CPannot=""
# raise error if left unmodified
if [ -z "$CPannot" ]
  then
    echo "Please specify the path to a CPannot.gz file"
    exit 1
fi

imgTemp="tempfilename_tile%03d.CPseries"
if [[ "$imgTemp" == "tempfilename_tile%03d.CPseries" ]]
  then
    echo "Please specify a image name template!"
    exit 1
fi
###########################


###### OTHER VARIABLES ####
N_tiles=18
script_dir=$SCRIPT_DIR
###########################


###### USER INPUTS ########
# pattern is the common pattern shared by each set of images (will likely be 'equilibrium')
# color is the fluorescence channel of the images ('red' or 'green')
pattern=$1
color=$2
###########################


###### BEGIN SCRIPT #######
# Setup env
module load python/2.7.13
source $py2env/bin/activate


dir_list=()
for s in *$pattern**$color*
do
	dir_list+=("$s")
done

N_concs="${#dir_list[@]}"
all_CPseries_Dirs="${dir_list[@]}"


# How to call python script:
# Usage: python combineCPseries.py [N tiles] [N concs] [all CPseries Dirs]
#        [imagename template] [AllOutFilename] [anyRNAOutFilename] [CPannot.gz]

python $script_dir/combineCPseriesMissing.py $N_tiles $N_concs $all_CPseries_Dirs $imgTemp "$pattern"_"$color"_all.CPseries.gz "$pattern"_"$color"_anyRNA.CPseries.gz $CPannot
