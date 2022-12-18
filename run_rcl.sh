#!/bin/bash


while getopts p:f:o: flag
do
    case "${flag}" in
        p) path=${OPTARG};;
        f) fname=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

#for i in $(seq 1 $nreps)
nreps=0
for i in $fname
do
((nreps++))
while read j; do
        [ "$j" = "X" ] && continue
        [ "$j" = "Y" ] && continue
        cat $path/chr${j}/*${i}.covBga.txt >> $path/rep$nreps.txt
done < $path/chrList.txt
done

echo "Number of reps ${nreps}"
python main.py --datapath $path --n_rep $nreps --modelpath $path/rcl.ckpt

while read j; do
python rcl_score.py --model $path/rcl.ckpt --dpath $path/chr${j}  --preprocess_region $path/bigInputs.txt --id $j --prefix $path
done < $path/chrList.txt

cat $path/rcl_*bed > $path/rcl.bed
rm $path/rcl_*bed 

nreps=0
for i in $fname
do
	((nreps++))
	rm  $path/rep$nreps.txt
done
