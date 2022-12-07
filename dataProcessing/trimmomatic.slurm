#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=hhvu@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --output="./reports/trim-A459-%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="./reports/trim-A459-%j.errors" # job standard error file (%j replaced by job id)


cd /work/LAS/geetu-lab-collab/4_A549/0_rawData

module load trimmomatic

#java -jar /home/hhvu/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 -trimlog rep1_1.log ENCFF047YJC.fastq.gz ENCFF048FVL.fastq.gz -baseout rep1_1.trimmed.fastq.gz ILLUMINACLIP:/home/hhvu/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

#java -jar /home/hhvu/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 -trimlog rep1_2.log ENCFF208VTU.fastq.gz ENCFF118CNJ.fastq.gz -baseout rep1_2.trimmed.fastq.gz ILLUMINACLIP:/home/hhvu/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

java -jar /home/hhvu/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 -trimlog rep1_3.log ENCFF123YKW.fastq.gz ENCFF718ANQ.fastq.gz -baseout rep1_3.trimmed.fastq.gz ILLUMINACLIP:/home/hhvu/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

java -jar /home/hhvu/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 -trimlog rep1_4.log ENCFF358RZN.fastq.gz ENCFF736NNJ.fastq.gz -baseout rep1_4.trimmed.fastq.gz ILLUMINACLIP:/home/hhvu/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36




