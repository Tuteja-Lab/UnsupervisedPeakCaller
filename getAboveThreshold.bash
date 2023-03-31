#!/bin/bash

inputFile=$1
chr=$2
cutoff=$3
out=$4
thread=$5

#module load samtools
#module load bedtools2

samtools view -hb $inputFile $chr -@ $thread | bedtools genomecov -ibam - -pc -bga | grep -w "^$chr" | awk -v c=$cutoff '{if($4 > c){print}}' > $out




