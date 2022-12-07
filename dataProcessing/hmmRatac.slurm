#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=hhvu@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mem=128G # maximum memory per node
#SBATCH --output="./reports/hmmratac-%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="./reports/hmmratac-%j.errors" # job standard error file (%j replaced by job id)


 
module load jdk
module load samtools


cd /work/LAS/geetu-lab-collab/3_GM12878/1_bowtie2/hg19
#sorted-uniqMappedPaired-rmMT-rmDup-sample69.sorted.bam
#samtools merge GM12878-merged6reps-mapq30.bam sorted-uniqMappedPaired-rmMT-rmDup-sample69.sorted.bam sorted-uniqMappedPaired-rmMT-rmDup-sample70.sorted.bam sorted-uniqMappedPaired-rmMT-rmDup-sample71.sorted.bam sorted-uniqMappedPaired-rmMT-rmDup-sample72.sorted.bam sorted-uniqMappedPaired-rmMT-rmDup-sample73.sorted.bam sorted-uniqMappedPaired-rmMT-rmDup-sample74.sorted.bam -@ 16
#samtools index GM12878-merged6reps-mapq30.bam -@ 16
#samtools sort -o sorted.GM12878-merged6reps-mapq30.bam GM12878-merged6reps-mapq30.bam -@ 16
#samtools index sorted.GM12878-merged6reps-mapq30.bam

file="sorted.GM12878-merged6reps-mapq30.bam"
out="GM12878merged6reps-mapq30-hg19"
java -Xmx128G -jar /work/LAS/geetu-lab-collab/HMMRATAC_V1.2.10_exe.jar -b $file -i "$file".bai -g /work/LAS/geetu-lab/hhvu/index_files/homo_sapiens/hg19.chrom.sizes -o /work/LAS/geetu-lab-collab/3_GM12878/5_hmmratac/$out --bedgraph true --bgscore true --modelonly true --window 250000 # -u 40 -l 10



