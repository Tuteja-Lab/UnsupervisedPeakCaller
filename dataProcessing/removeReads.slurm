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
#SBATCH --output="./reports/remove-GM12878-%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="./reports/remove-GM12878-%j.errors" # job standard error file (%j replaced by job id)




module load samtools
#cd /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/MCF7/2_bowtie2/uniqPaired/
#for i in {1..22} X Y
#do
#samtools view -hb MCF7_coorSorted.bam $i -@ 16 > "$i".MCF7.bam
#samtools sort -@ 16 "$i".MCF7.bam -o sorted."$i".MCF7.bam
#samtools index sorted."$i".MCF7.bam
#done

cd /work/LAS/geetu-lab-collab/3_GM12878/1_bowtie2/hg19 
for file in rmDup*bam
do
samtools view -h $file | grep -vw 'chrMT' - | samtools view -Sb - -o rmMT-$file
samtools view -c rmMT-$file
samtools view --threads 16 -h -b -q 30 -f 2 rmMT-$file > uniqMappedPaired-rmMT-"$file"
samtools flagstat uniqMappedPaired-rmMT-"$file" --threads 16
samtools sort -@ 16 uniqMappedPaired-rmMT-"$file" -o sorted-uniqMappedPaired-rmMT-"$file"
samtools index sorted-uniqMappedPaired-rmMT-"$file"
done

#cd /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/e9.5placenta/2_bowtie2/H3K27ac/
#for file in rmDup-SRR1796459.sorted.bam rmDup-SRR1796460.sorted.bam rmDup-SRR1796461.sorted.bam rmDup-SRR1796462.sorted.bam
#do
##samtools view -h $file | grep -vw 'MT' - | samtools view -Sb - -o rmMT-$file
#samtools view -c $file -@ 16
#samtools view --threads 16 -h -b -q 20 $file > uniqMappedPaired-"$file"
#samtools flagstat uniqMappedPaired-"$file" --threads 16
#samtools sort -@ 16 uniqMappedPaired-"$file" -o sorted-uniqMappedPaired-"$file"
#samtools index sorted-uniqMappedPaired-"$file"
#done


#cd /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/e9.5placenta/2_bowtie2/H3K4me1
#for file in rmDup*bam
#do
##samtools view -h $file | grep -vw 'MT' - | samtools view -Sb - -o rmMT-$file
#samtools view -c $file -@ 16
#samtools view --threads 16 -h -b -q 20 $file > uniqMappedPaired-"$file"
#samtools flagstat uniqMappedPaired-"$file" --threads 16
#samtools sort -@ 16 uniqMappedPaired-"$file" -o sorted-uniqMappedPaired-"$file"
#samtools index sorted-uniqMappedPaired-"$file"
#done


#cd /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/e9.5placenta/2_bowtie2/H3K4me3
#for file in rmDup*bam
#do
##samtools view -h $file | grep -vw 'MT' - | samtools view -Sb - -o rmMT-$file
#samtools view -c $file -@ 16
#samtools view --threads 16 -h -b -q 20 $file > uniqMappedPaired-"$file"
#samtools flagstat uniqMappedPaired-"$file" --threads 16
#samtools sort -@ 16 uniqMappedPaired-"$file" -o sorted-uniqMappedPaired-"$file"
#samtools index sorted-uniqMappedPaired-"$file"
#done

#cd /work/LAS/geetu-lab/hhvu/project4_ATACpeakCalling/data/e9.5placenta/2_bowtie2/H3K27me3
#for file in rmDup*bam
#do
##samtools view -h $file | grep -vw 'MT' - | samtools view -Sb - -o rmMT-$file
#samtools view -c $file -@ 16
#samtools view --threads 16 -h -b -q 20 $file > uniqMappedPaired-"$file"
#samtools flagstat uniqMappedPaired-"$file" --threads 16
#samtools sort -@ 16 uniqMappedPaired-"$file" -o sorted-uniqMappedPaired-"$file"
#samtools index sorted-uniqMappedPaired-"$file"
#done

