#!/bin/bash

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mem=128G
#SBATCH --mail-user=hhvu@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --output="./reports/markdup-GM12878hg19-%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="./reports/markdup-GM12878hg19-%j.errors" # job standard error file (%j replaced by job id)



module load samtools
module load picard
cd /work/LAS/geetu-lab-collab/3_GM12878/1_bowtie2/hg19
for i in *bam
	do
		e=`echo $i | sed 's/bam/sorted.bam/g'`
		echo $e
		samtools sort -@ 16 -o $e $i

	done

for i in *sorted.bam
        do 
                e=`echo $i | sed 's/bam/txt/g'`
                echo $e
		picard MarkDuplicates REMOVE_DUPLICATES=true I=$i O=rmDup-$i M=metrics_$e
#picard MarkDuplicates I=/work/LAS/geetu-lab/jain/HIF_ChIP-Project/F11_fastq/F11_fastq.sorted.bam O=/work/LAS/geetu-lab/jain/HIF_ChIP-Project/F11_fastq/F11_fastq.sorted.markDup.bam M=/work/LAS/geetu-lab/jain/HIF_ChIP-Project/F11_fastq/mark_dupl_metrics.txt REMOVE_DUPLICATES=true
	done
