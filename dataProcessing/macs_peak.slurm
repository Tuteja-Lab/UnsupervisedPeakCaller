#!/bin/bash

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 10 processor core(s) per node 
#SBATCH --mail-user=hhvu@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="./reports/macs-GM12878-%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="./reports/macs-GM12878-%j.errors" # job standard error file (%j replaced by job id)


#module load zlib/1.2.11-vhzh5cf
module purge
module load bedtools2
module load py-macs2/2.1.1.20160309-py2-wgwiexf

#cd /work/LAS/geetu-lab-collab/3_GM12878/1_bowtie2/rmDup/
#sorted.GM12878-merged3reps.bam
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.05-nolambda --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.05 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.5-nolambda --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.5 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.1-nolambda --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.1 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.005-nolambda --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.005 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.0005-nolambda --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.0005 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.00005-nolambda --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.00005 --nolambda
#macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample68.bam --name rep1 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.05
#macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample69.bam --name rep2 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.05
#macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample70.bam --name rep3 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.05
#macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample71.bam --name rep4 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.05

#cd /work/LAS/geetu-lab-collab/3_GM12878/1_bowtie2/hg19
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.05-nolambda-hg19 --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.05 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.5-nolambda-hg19 --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.5 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.1-nolambda-hg19 --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.1 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.005-nolambda-hg19 --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.005 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.0005-nolambda-hg19 --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --keep-dup all -q 0.0005 --nolambda
#macs2 callpeak --treatment sorted.GM12878-merged6reps-mapq30.bam --name rep2_7-mapq30.0.00005-nolambda-hg19 --outdir /work/LAS/geetu-lab-collab/3_GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.00005 --nolambda

cd /work/LAS/geetu-lab-collab/3_GM12878/1_bowtie2/uniqPaired/
#GM12878.bam                       sorted.GM12878-merged3reps.bam.bai                    sorted-uniqMappedPaired-rmDup-sample70.bam            sorted-uniqMappedPaired-rmDup-sample73.bam
#GM12878_coorSorted.bam            sorted.GM12878-merged6reps.bam                        sorted-uniqMappedPaired-rmDup-sample70.bam.bai        sorted-uniqMappedPaired-rmDup-sample73.bam.bai
#GM12878_coorSorted.bam.bai        sorted.GM12878-merged6reps.bam.bai                    sorted-uniqMappedPaired-rmDup-sample70_normalized.bw  sorted-uniqMappedPaired-rmDup-sample73_normalized.bw
#GM12878_coorSorted_normalized.bw  sorted-uniqMappedPaired-rmDup-sample68.bam            sorted-uniqMappedPaired-rmDup-sample71.bam            sorted-uniqMappedPaired-rmDup-sample74.bam
#GM12878-merged3reps.bam           sorted-uniqMappedPaired-rmDup-sample68.bam.bai        sorted-uniqMappedPaired-rmDup-sample71.bam.bai        sorted-uniqMappedPaired-rmDup-sample74.bam.bai
#GM12878-merged3reps.bam.bai       sorted-uniqMappedPaired-rmDup-sample68_normalized.bw  sorted-uniqMappedPaired-rmDup-sample71_normalized.bw  sorted-uniqMappedPaired-rmDup-sample74_normalized.bw
#GM12878-merged6reps.bam           sorted-uniqMappedPaired-rmDup-sample69.bam            sorted-uniqMappedPaired-rmDup-sample72.bam
#GM12878-merged6reps.bam.bai       sorted-uniqMappedPaired-rmDup-sample69.bam.bai        sorted-uniqMappedPaired-rmDup-sample72.bam.bai
sorted.GM12878-merged3reps.bam    sorted-uniqMappedPaired-rmDup-sample69_normalized.bw  sorted-uniqMappedPaired-rmDup-sample72_normalized.bw

for i in {68..71}
do
#	macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.05 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.05
#	macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.5 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.5
#	macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.1 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.1
#	macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.005 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.005
#	macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.0005 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.0005
#	macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.00005 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.00005

        macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.01 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.01
#        macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.5 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.5
#        macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.1 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.1
        macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.001 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.001
        macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.0001 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.0001
        macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.00001 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.00001
        macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.00002 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.00002
	macs2 callpeak --treatment sorted-uniqMappedPaired-rmDup-sample"$i".bam --name sample"$i"_q0.000001 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.000001
done


#macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.05 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.05
#macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.5 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.5
#macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.1 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.1
#macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.005 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.005
#macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.0005 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.0005
#macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.00005 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.00005

macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.01 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.01
macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.001 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.001
macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.0001 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.0001
macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.00001 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.00001
macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.00002 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.00002
macs2 callpeak --treatment GM12878_coorSorted.bam --name merged68_71_q0.000001 --outdir /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/GM12878/4_MACS/ -g hs -f BAMPE --bdg --keep-dup all -q 0.000001

