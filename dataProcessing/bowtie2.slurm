#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mem=128G 
#SBATCH --mail-user=hhvu@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="./reports/bowtie2%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="./reports/bowtie2%j.errors" # job standard error file (%j replaced by job id)




module load bowtie2

cd /work/LAS/geetu-lab-collab/3_GM12878/0_rawData

for i in {68..74}
do
bowtie2 -p 16 --no-discordant -x /work/LAS/geetu-lab/hhvu/index_files/homo_sapiens/bowtie2_index-hg19/hg19 -1 SRR8912"$i"_1.fastq.gz -2 SRR8912"$i"_2.fastq.gz -S /work/LAS/geetu-lab-collab/3_GM12878/1_bowtie2/hg19/sample"$i".bam
done
