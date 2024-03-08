#!/bin/bash
#
# Quantify array images on sherlock
#
# Usage: quantifyImages.sh prefix seq_dir roff_dir fluor_dir gv_path script_dir
#
# Ben Ober-Reynolds, boberrey@stanford.edu, 20191105

# prefix is the common prefix shared by each set of images (will likely be 'set')
prefix=$1
seq_dir="/scratch/groups/herschla/roy-test/20210218_30mM_Mg_Lib4_run1_data/seqData/split-tile"
roff_dir="/scratch/groups/herschla/roy-test/20210225_30mMMg+150mMKCl_Lib4/Roff"
fluor_dir="/scratch/groups/herschla/roy-test/20210225_30mMMg+150mMKCl_Lib4/CPfluor"
script_dir="/scratch/groups/herschla/roy-test/image-processing/lena_array_tools_Ben"
gv_path="/scratch/groups/herschla/roy-test/20210225_30mMMg+150mMKCl_Lib4/BGV"
log_dir="/scratch/groups/herschla/roy-test/20210225_30mMMg+150mMKCl_Lib4/Logfile"

for s in *$prefix*/
do
    echo "Current image directory = $s"
    mkdir -p $roff_dir/$s
    mkdir -p $fluor_dir/$s
    mkdir -p $log_dir/$s
    #echo Made dir: $roff_dir/$s
    #echo Made dir: $fluor_dir/$s
    sbatch $script_dir/quantify_images_BOTH.sbatch $s $seq_dir $roff_dir/$s $fluor_dir/$s $gv_path $log_dir/$s
done
