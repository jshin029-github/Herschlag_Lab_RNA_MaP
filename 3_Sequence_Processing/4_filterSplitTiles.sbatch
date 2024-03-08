#!/bin/bash
#
# bash script to generate chip statistics
#
## Usage: sbatch 4_filterSplitTiles.sh <split tile dir>
#
#
#SBATCH --job-name=filterSplitTiles
#SBATCH --partition biochem,normal
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH -e logging/filterSplitTiles-%j.err
#SBATCH -o logging/filterSplitTiles-%j.out

###########################

# User input variables:
# raise error for incorrect usage
if [ $# -ne 1 ]
  then
    echo "Incorrect number of arguments: <split tile dir>"
    exit 1
fi

############################

# Start script, first filter the split tile CPseq

split_tile_dir="$1"

fid_dir="$split_tile_dir/FID/"
anyRNA_dir="$split_tile_dir/anyRNA/"
phiX_dir="$split_tile_dir/phiX/"

echo "Making new directories"

mkdir "$fid_dir"
mkdir "$anyRNA_dir"
mkdir "$phiX_dir"

for full_name in $split_tile_dir/*CPseq
do
  file_name=$(basename "$full_name")

  echo "Filtering $file_name"

  grep CTTGGGTCCACAGGAC $full_name > "$fid_dir/FID_$file_name"
  grep TTATGCTATAATTATT $full_name > "$anyRNA_dir/anyRNA_$file_name"
  grep -E -v 'CTTGGGTCCACAGGAC|TTATGCTATAATTATT' $full_name > "$phiX_dir/phiX_$file_name"
done


# Next, flag the filtered CPseq and merge the split tiles

module load python/3.6.1
source /home/groups/herschla/rna_map/scripts/new_scripts/venv_3_6_1/bin/activate


echo "Flagging and merging FID files"

fid_arr=($(find "$fid_dir" -name "*.CPseq"))

for file_name in "${fid_arr[@]}"
do
  filehandle=$(basename "$file_name")

  python3 /home/groups/herschla/rna_map/scripts/new_scripts/writeSetFlagtoCPseq.py -i "$file_name" -f FID -o "$fid_dir""FID_flagged_$filehandle"

  rm $file_name
done


echo "Flagging and merging anyRNA files"

anyRNA_arr=($(find "$anyRNA_dir" -name "*.CPseq"))

for file_name in "${anyRNA_arr[@]}"
do
  filehandle=$(basename "$file_name")
  python3 /home/groups/herschla/rna_map/scripts/new_scripts/writeSetFlagtoCPseq.py -i "$file_name" -f anyRNA -o "$anyRNA_dir""anyRNA_flagged_$filehandle"
  rm $file_name
done


echo "Flagging and merging phiX files"

phiX_arr=($(find "$phiX_dir" -name "*.CPseq"))

for file_name in "${phiX_arr[@]}"
do
  filehandle=$(basename "$file_name")
  python3 /home/groups/herschla/rna_map/scripts/new_scripts/writeSetFlagtoCPseq.py -i "$file_name" -f phiX -o "$phiX_dir""phiX_flagged_$filehandle"
  rm $file_name
done


echo "Combining FID and anyRNA files"

anyRNA_fid_dir="$split_tile_dir/anyRNA_FID/"
mkdir "$anyRNA_fid_dir"

for i in $(seq -f "%03g" 1 19)
do
  ar=($(find "$anyRNA_dir" "$fid_dir" -name "*$i*.CPseq"))
  file_name=$(basename "${ar[1]}")
  cat "${ar[0]}" "${ar[1]}" > "$anyRNA_fid_dir""anyRNA_""$file_name"
done


echo "Combining all files"

all_dir="$split_tile_dir/all/"
mkdir "$all_dir"

for i in $(seq -f "%03g" 1 19)
do
  ar=($(find "$anyRNA_fid_dir" "$phiX_dir" -name "*$i*.CPseq"))
  file_name=$(basename "${ar[1]}")
  cat "${ar[0]}" "${ar[1]}" > "$all_dir""all_anyRNA_FID_""${ar[1]}"
done


echo "done!"

