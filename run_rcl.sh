#!/bin/bash
# Run RCL.
# Input: PATH to data files, list of BAM file replicates BFILE1, BFILE2, ...
# Assumes: PATH/chrList.txt, PATH/chr*/BFILE1.covBga.txt, PATH/chr*/BFILE2.covBga.txt, PATH/bigInputs.txt created by preprocessing.bash
# Output:
#	PATH/rcl.ckpt	fitted model
#	PATH/rcl.bed	candidate peaks and scores

path="."
epoch=25
batch=256
ext=.covBga.txt
ref_prefix="chr"
save=0
overwrite=0
debug=""

helpFunction()
{
	echo ""
	echo "Usage: $0 [-d PATH] -b \"BAM_FILE1 BAM_FILE2[ BAM_FILE3...]\""
	echo -e "\t-b STR list of preprocessed replicate BAM files, surrounded by double quotes."
	echo -e "\t       BAM files should have alread been preprocessed by preprocessing.bash."
	echo -e "\t       Example (from tutorial): \"MCF7_chr10_rep1.bam MCF7_chr10_rep2.bam\""
	echo -e "\t-d STR path to Data files (DEFAULT: $path)."
	echo -e "\t-e INT number of Epochs (DEFAULT: $epoch)."
	echo -e "\t-g     turn on debuGging output (DEFAULT: no)."
	echo -e "\t-h INT batcH size (DEFAULT: $batch)."
	echo -e "\t-r STR Reference sequence name prefix (DEFAULT: $ref_prefix)."
	echo -e "\t-s     Save output files (DEFAULT: no)."
	echo -e "\t-w     overWrite existing files (DEFAULT: no)."
	echo -e "\t-x STR eXtension for coverage files (DEFAULT: $ext)."
	echo -e "\t       Chromosome j coverage file for replicate m assumed to be"
	echo -e "\t       PATH/chrj/BASEm_j.covBga.txt, where BASEm is BAM_FILEm (-f) without .bam extension."
	exit 1 # Exit script after printing help
}

while getopts ?gswb:d:e:h:r:x: flag
do
	case "${flag}" in
		b) fname=${OPTARG};;
		d) path=${OPTARG};;
		e) epoch=${OPTARG};;
		g) debug=--debug;;
		r) ref_prefix=${OPTARG};;
		h) batch=${OPTARG};;
		w) overwrite=1;;
		s) save=1;;
		x) ext=${OPTARG};;
		?) helpFunction ;; # Print helpFunction	
	esac
done

path=${path%/}	# remove terminal /

echo "Number of epoch ${epoch}"

# Print helpFunction in case parameters are empty
if [ -z "$path" ] || [ -z "$fname" ]  #|| [ -z "$skip" ]
then
	   echo "Some or all required parameters are empty";
	   helpFunction
fi

rep_names=${fname//.bam/}

#for i in $(seq 1 $nreps)
nreps=0
for rep_name in $rep_names
do
	((nreps++))
	if [ $overwrite -eq 1 -o ! -s "$path"/rep$nreps.txt ]; then
		rm -f "$path"/rep$nreps.txt
		while read chr; do
			[ "$chr" = "X" ] && continue
			[ "$chr" = "Y" ] && continue
			cat "$path"/$ref_prefix"$chr"/"$rep_name""$ext" >> "$path"/rep$nreps.txt
		done < "$path"/chrList.txt
	fi
done

if [ $overwrite -eq 1 -o ! -s "$path"/rcl.ckpt ]; then
	echo "Number of reps ${nreps}, start training"
	python main.py $debug --epochs $epoch --batch_size $batch --datapath "$path" --n_rep $nreps --modelpath "$path"/rcl.ckpt &> "$path"/out.err
	echo "Finish training, start writing results (if your data is large, please give more memory)"
else
	echo "Using existing \"$path/rcl.ckpt\" file"
fi

if [ $overwrite -eq 1 -o ! -s "$path"/rcl.bed ]; then
	while read chr; do
		python rcl_score.py $debug --model "$path"/rcl.ckpt --dpath "$path"/$ref_prefix"$chr" --names $rep_names --preprocess_region "$path"/bigInputs.txt --id $chr --prefix "$path"
	done < "$path"/chrList.txt
	cat "$path"/rcl_*bed > "$path"/rcl.bed
else
	echo "Using existing \"$path/rcl.bed\" file"
fi

if [ $save -eq 0 ]; then
	rm -f "$path"/rcl_*bed 
	rm -f "$path"/out.err

	nreps=0
	for i in $fname
	do
		((nreps++))
		rm -f "$path"/rep$nreps.txt
	done
fi

echo "Finished!"
