#!/bin/bash

chrList=$1
file=$2
out=$3

for i in $chrList
do
m=`grep "^$i" $file | cut -f 4 | sort -rn -S 50% --parallel=5 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]'`
echo $file >> $out
echo -e "$i""\t""$m" >> $out
done


