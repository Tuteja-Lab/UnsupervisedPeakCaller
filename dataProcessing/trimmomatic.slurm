#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=hhvu@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --output="./reports/trim-e9.5placenta-%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="./reports/trim-e9.5placenta-%j.errors" # job standard error file (%j replaced by job id)


cd /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/e9.5placenta/0_rawData/atac-seq

module load trimmomatic
#SRR7912013  SRR7912014  SRR7912015

java -jar /home/hhvu/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 -trimlog SRR7912014.log SRR7912014/SRR7912014_1.fastq.gz SRR7912014/SRR7912014_2.fastq.gz -baseout SRR7912014.trimmed.fastq.gz ILLUMINACLIP:/work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/e9.5placenta/scripts/atacOver.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

java -jar /home/hhvu/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 -trimlog SRR7912015.log SRR7912015/SRR7912015_1.fastq.gz SRR7912015/SRR7912015_2.fastq.gz -baseout SRR7912015.trimmed.fastq.gz ILLUMINACLIP:/work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/e9.5placenta/scripts/atacOver.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
