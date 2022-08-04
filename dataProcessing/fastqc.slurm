#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mail-user=hhvu@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="reports/fastQC-A549-%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="reports/fastQC-A549-%j.errors" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load fastqc
cd /work/LAS/geetu-lab-collab/4_A549/0_rawData/
for i in *
do
        fastqc *.fastq.gz -o /work/LAS/geetu-lab-collab/4_A549/1_fastQC-beforeTrim/
done

