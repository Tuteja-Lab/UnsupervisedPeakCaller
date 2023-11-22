#!/bin/bash

chrList=$1
file=$2
out=$3
threads=$4

for chr in $chrList
do
	# there were bugs here: need -w and quotes
	median=`grep -w "^$chr" "$file" | cut -f 4 | sort -rn -S 50% --parallel=$threads | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'`
	echo "$file" >> "$out"
	echo -e "$chr\t$median" >> "$out"
done


