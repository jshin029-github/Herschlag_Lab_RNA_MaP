#!/bin/bash
#
# Quantify array images on sherlock
#
# Usage: quantifyImages.sh prefix seq_dir roff_dir fluor_dir gv_path script_dir
#
# Ben Ober-Reynolds, boberrey@stanford.edu, 20191105

# prefix is the common prefix shared by each set of images (will likely be 'set')
prefix="20_equilibrium2000nM"
seq_dir="CPseq"
roff_dir="roff"
fluor_dir="cpfluor"
gv_path=////GLbalVar.m
script_dir=/scratch/groups/herschla/roy-test/20210218_30mM_Mg_Lib4_run1_data

for s in *$prefix*/
do
    echo "Current image directory = $s"
    #mkdir -p $roff_dir/$s/
    #mkdir -p $fluor_dir/$s/
    echo Made dir: $s/$roff_dir/
    echo Made dir: $s/$fluor_dir/
    #sbatch $script_dir/quantify_images.sbatch $s $seq_dir $roff_dir/$s/ $fluor_dir/$s/ $gv_path
    echo "sbatch $script_dir/quantify_images.sbatch $s $seq_dir $roff_dir/$s/ $fluor_dir/$s/ $gv_path"
done
