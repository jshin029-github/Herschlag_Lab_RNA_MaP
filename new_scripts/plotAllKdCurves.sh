cpvar=$1
vardir=$2
N=10

mkdir -p kd_curves

for((i=1; i<=$N; i++))
do
    sbatch $rnamap_scripts/new_scripts/plotKdCurves.sh $cpvar $vardir/variant_list_${i}.txt
done


