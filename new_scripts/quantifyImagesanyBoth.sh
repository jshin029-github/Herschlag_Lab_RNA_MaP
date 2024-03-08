#!/bin/bash
#
# Quantify array images on sherlock
#
# Usage: quantifyImages.sh prefix seq_dir roff_dir fluor_dir gv_path script_dir
#
# Ben Ober-Reynolds, boberrey@stanford.edu, 20191105

######## TO MODIFY ########
seq_dir=""
###########################


# prefix is the common prefix shared by each set of images (will likely be 'set')
prefix=$1
roff_dir="./Roff"
fluor_dir="./CPfluor"
gv_path="$rnamap_scripts/image_processing_scripts"
log_dir="./Logfile"
script_dir="$rnamap_scripts/image_processing_scripts"

for s in *$prefix*/
do
    echo "Current image directory = $s"
    mkdir -p $roff_dir/$s
    mkdir -p $fluor_dir/$s
    mkdir -p $log_dir/$s
    #echo Made dir: $roff_dir/$s
    #echo Made dir: $fluor_dir/$s
    sbatch $rnamap_scripts/image_processing_scripts/quantify_images_BOTH.sbatch $s $seq_dir $roff_dir/$s $fluor_dir/$s $gv_path $log_dir/$s
done
