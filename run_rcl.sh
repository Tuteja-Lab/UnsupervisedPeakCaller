#!/bin/bash

helpFunction()
{
	   echo ""
	   echo "Usage: $0 -p \"path to preprocessing data\" -f \"indiv1 indiv2\""
	   echo -e "\t-p path to preprocessing data."
	   echo -e "\t-f name of the replicate file (without suffix)"
	   echo -e "\t-e number of epoch (default 25)."
	   echo -e "\t-b batch size (default 256)."
	   exit 1 # Exit script after printing help
}

epoch=25
batch=256

while getopts p:f:e:b flag
do
    case "${flag}" in
        p) path=${OPTARG};;
        f) fname=${OPTARG};;
        e) epoch=${OPTARG};;
        b) batch=${OPTARG};;
	?) helpFunction ;; # Print helpFunction	
    esac
done

echo "Number of epoch ${epoch}"

# Print helpFunction in case parameters are empty
if [ -z "$path" ] || [ -z "$fname" ]  #|| [ -z "$skip" ]
then
	   echo "Some or all required parameters are empty";
	   helpFunction
fi


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

echo "Number of reps ${nreps}, start training"
python main.py --epochs $epoch --batch_size $batch --datapath $path --n_rep $nreps --modelpath $path/rcl.ckpt &> $path/out.err
echo "Finish training, start writing results (if your data is large, please give more memory)"

while read j; do
python rcl_score.py --model $path/rcl.ckpt --dpath $path/chr${j}  --preprocess_region $path/bigInputs.txt --id $j --prefix $path
done < $path/chrList.txt

cat $path/rcl_*bed > $path/rcl.bed
rm $path/rcl_*bed 
rm $path/out.err

nreps=0
for i in $fname
do
	((nreps++))
	rm  $path/rep$nreps.txt
done

echo "Finished!"
