#!/bin/bash 
####### Change the following before running #######
template=""           # Template name of CPseries filesc
ntile=""                                      # subtiles no
color=""                                   # green or red channel images to plot
####################################################

mkdir -p Heatmap_plots
for n in{1...$ntile}
do
    sbatch $rnamap_scripts/new_scripts/makeClusterHeatmapSingle.sh $template $n $color
done
