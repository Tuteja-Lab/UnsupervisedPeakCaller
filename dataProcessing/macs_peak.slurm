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


cd /work/LAS/geetu-lab-collab/3_GM12878/1_bowtie2/uniqPaired/
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

