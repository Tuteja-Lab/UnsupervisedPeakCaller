#!/bin/bash

module load samtools
module load bedtools2
module load parallel
module load bedops/2.4.35-gl7y6z6
module load gcc/7.3.0-xegsmw4
module load r/4.0.2-py3-icvulwq
module load gsl/2.5-fpqcpxf
module load udunits/2.2.24-yldmp4h
module load gdal/2.4.4-nw2drgf
module load geos/3.8.1-2m7gav4

helpFunction()
{
   echo ""
   echo "Usage: $0 -p \"program directory\" -i \"input directory\" -o \"output directory\" -g hg -c 2 -m \"merged.bam\" -b \"indi1.bam indi2.bam\" -t 12 -n test -L 1000"
   echo -e "\t-p Absolute directory of where the program is installed at."
   echo -e "\t-i Absolute directory of input files."
   echo -e "\t-o Absolute directory of output files."
   echo -e "\t-g Genome that the data is aligned to. Currently support mm10 (Ensembl), hg38 (Ensembl), danRer11 (Ensembl)."
   echo -e "\t-c Cutoff for prefiltering. Either \"median\" or specific number."
   echo -e "\t-m Bam files merged from individual replicates. Only used for preprocessing purpose, not for calling peaks. Must be sorted."
   echo -e "\t-b Individual bam files of every replicates. Must be sorted."
   echo -e "\t-t Number of threads to use."
   echo -e "\t-n File name prefix."
   echo -e "\t-L Length of input segments."

   exit 1 # Exit script after printing help
}


while getopts ":p:i:o:g:c:m:b:t:n:L:" opt
do
    case "$opt" in
        p) pdir="$OPTARG";;
        i) indir="$OPTARG";;
        o) outdir="$OPTARG";;
        g) genome="$OPTARG";;
        c) cutoff="$OPTARG";;
        m) mergedBam="$OPTARG";;
        b) indivBam="$OPTARG";;
        t) thread="$OPTARG";;
        n) prefix="$OPTARG";;
        L) inputLn="$OPTARG";;
        ?) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "$pdir" ] || [ -z "$indir" ] || [ -z "$outdir" ] || [ -z "$genome" ]  || [ -z "$cutoff" ]  || [ -z "$mergedBam" ] || [ -z "$indivBam" ] || [ -z "$thread" ] || [ -z "$prefix" ] || [ -z "$inputLn" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi



#step 1 - get necessary count files
#step 1.0: define chromosome numbers
chrList=`samtools view -H "$indir"/"$mergedBam" | grep "SN:" | cut -f 2 | sed 's/SN://g'`

#step 1.1: get bga genomecov of the data
for i in $chrList
do
	samtools view -hb $indir/$mergedBam $i -@ $thread | bedtools genomecov -ibam - -pc -bga | grep -w "^$i" > $outdir/"$i"_"$prefix"_merged-bga.begGraph
done

for i in $indivBam
do	
	o=`echo $i | sed 's/.bam//g'`
	bedtools genomecov -ibam "$indir"/"$i" -pc -bga > "$outdir"/"$prefix"_"$o"-bga.begGraph
done

#step 1.2: get counts based on the cutoff
if [[ $cutoff == "median" ]]
then
  for i in $indivBam
	do
        	o=`echo $i | sed 's/.bam//g'`
		bash "$pdir"/getMedian.bash "$chrList" "$outdir"/"$prefix"_"$o"-bga.begGraph "$outdir"/"$prefix"_median.txt
	done

  for i in $chrList
	do
		c=`grep -w "^$i" "$outdir"/"$prefix"_median.txt | sort -k2,2n | cut -f 2 | head -n 1`
		mkdir $outdir/chr"$i"
		echo $indivBam | sed 's/ /\n/g' | parallel "bash \"$pdir\"/getAboveThreshold.bash \"$indir\"/{} \"$i\" \"$c\" \"$outdir\"/chr\"$i\"/\"$i\"_{}-temp.txt \"$thread\""
	done
else
  c="$cutoff"
  for i in $chrList
        do
                mkdir $outdir/chr"$i"
                echo $indivBam | sed 's/ /\n/g' | parallel "bash \"$pdir\"/getAboveThreshold.bash \"$indir\"/{} \"$i\" \"$c\" \"$outdir\"/chr\"$i\"/\"$i\"_{}-temp.txt \"$thread\""
        done
fi

#step 2: get positions that pass the threshold in every replicate
for i in $chrList
do
bedops --intersect "$outdir"/chr"$i"/"$i"_* | bedtools merge -d 90 -i - | awk '{if ($3-$2 >= 100) {print}}' | awk '{$(NF+1)="segment"NR}1' | sed 's/ /\t/g' | bedtools intersect -wa -wb -a - -b "$outdir"/"$i"_"$prefix"_merged-bga.begGraph > "$outdir"/chr"$i"/temp"$i".txt
done

#step 3: get $inputLn bp regions
arr=($chrList)
m=${#arr[@]}
for((n=0;n<${#arr[@]};n++)); do
        if (( $(($n % 3 )) == 0 & $(($n < ($m-2) )))); then
                # Run every 3 entries
                echo ${arr[n]} ${arr[n+1]} ${arr[n+2]} | sed 's/ /\n/g' | parallel "Rscript --vanilla \"$pdir\"/getCountFiles.R \"$outdir\"/chr{}/temp{}.txt \"$inputLn\" \"$outdir\"/chr{}/{}regions.txt"
        fi
done

if [[ $(($m % 3 )) == 1 ]]
then
	echo ${arr[$m-1]} | parallel "Rscript --vanilla \"$pdir\"/getCountFiles.R \"$outdir\"/chr{}/temp{}.txt \"$inputLn\" \"$outdir\"/chr{}/{}regions.txt"
elif [[ $(($m % 3 )) == 2 ]]
then
	echo ${arr[$m-2]} ${arr[$m-1]} | sed 's/ /\n/g' | parallel "Rscript --vanilla \"$pdir\"/getCountFiles.R \"$outdir\"/chr{}/temp{}.txt \"$inputLn\" \"$outdir\"/chr{}/{}regions.txt"
fi


#step 3.1: remove blacklist
if [[ $genome == "hg" ]]
then
  echo "Blacklist regions are Ensembl genome. This will not work if the genome used for alignment is from UCSC."
  bl=""$pdir"/hg38-blacklist.v2.ensembl.bed"
elif [[ $genome == "mm" ]]
then
  echo "Blacklist regions are Ensembl genome. This will not work if the genome used for alignment is from UCSC."
  bl=""$pdir"/mm10-blacklist.v2.ensembl.bed"
else
  echo "Only human and mouse genome blacklists are supported for now. Blacklist regions are Ensembl genome. This will not work if the genome used for alignment is from UCSC. If not mouse or human, will skip removing blacklist regions."
fi

if [[ $genome == "hg" ]] | [[ $genome == "mm" ]]
then
	for i in $chrList
	do
		tail -n +2 "$outdir"/chr"$i"/"$i"regions.txt | sed 's/chr//g' | bedtools intersect -v -a - -b "$bl" > "$outdir"/chr"$i"/"$i"regions2.txt
	done
else
        for i in $chrList
        do
                tail -n +2 "$outdir"/chr"$i"/"$i"regions.txt | sed 's/chr//g' > "$outdir"/chr"$i"/"$i"regions2.txt
        done
fi


#step 4: get counts
indivBam2=`echo $indivBam | sed 's/.bam//g'`
for i in $chrList
do
echo $indivBam2 | sed 's/ /\n/g' | parallel "bedtools intersect -wb -a \"$outdir\"/chr\"$i\"/\"$i\"regions2.txt -b \"$outdir\"/\"$prefix\"_{}-bga.begGraph | cut -f1-4,8 > \"$outdir\"/chr\"$i\"/\"$i\".\"$prefix\"_{}.covBga.txt"
done

