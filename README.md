# Unsupervised Contrastive Peak Caller

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Preprocessing](#preprocessing)
1. [Peak Calling](#peakcalling)
1. [How to Cite](#cite)
1. [Contact](#contact)

# Prerequisites <a name = "prerequisites" />
For input preprocessing steps, the following tools and R libraries are required:
```
samtools (>= 1.10)
bedtools2 (>= 2.27.1)
parallel (>= 20170322)
R (>= 4.0.2)
bedops (>= 2.4.35)

R library dplyr (>= 1.0.7)
R library bedr (>= 1.0.7)
R library doParallel (>= 1.0.16)
```
For the deep learner step, **GPU** is needed. Other packages needed are:
```
Python (>=3.7.10)
PyTorch Lightning (>=1.5.1)
PyTorch (>=1.10.0)
numpy (>=1.21.5)
pandas (>=1.3.5)
argparse (>=1.1)
scikit-learn (>=1.0.1)
```

# Installation
```
git clone https://github.com/Tuteja-Lab/UnsupervisedPeakCaller.git
```

# Preprocessing <a name = "preprocessing" />
```
Usage: preprocessing.bash -p "program directory" -i "input directory" -o "output directory" -g hg -c 2 -m "merged.bam" -b "indi1.bam indi2.bam" -t 12 -n test -L 1000
        -p Absolute directory of where the program is installed at.
        -i Absolute directory of input files.
        -o Absolute directory of output files.
        -g Genome that the data is aligned to. Currently support mm10 (Ensembl) or hg38 (Ensembl).
        -c Cutoff for prefiltering. Either "median" or specific number.
        -m Bam files merged from individual replicates. Only used for preprocessing purpose, not for calling peaks. Must be indexed and sorted.
        -b Individual bam files of every replicate. Must be indexed and sorted.
        -t Number of threads to use.
        -n File name prefix.
        -L Length of input segments.
```
At this step, the script assumes your data has been aligned to mouse or human genome, Ensembl assembly.

## Example
Download example bam files here https://iastate.box.com/s/9uavg2zsy5w0i7v227ei7yaeea7yktr6

```
module load samtools
module load bedtools2
module load parallel
module load bedops/2.4.35-gl7y6z6
module load gcc/7.3.0-xegsmw4
module load r/4.0.2-py3-icvulwq
module load gsl/2.5-fpqcpxf
module load udunits/2.2.24-yldmp4h
module load gdal/2.4.4-nw2drgf
module load geos/3.8.1-2m7gav4

p="/work/LAS/geetu-lab-collab/UnsupervisedPeakCaller"
bash ${p}/preprocessing.bash -p ${p} -i "/work/LAS/geetu-lab-collab/UnsupervisedPeakCaller/example" -o "/work/LAS/geetu-lab-collab/UnsupervisedPeakCaller/example" -g "hg" -c "median" -m "MCF7_chr10_merged.bam" -b "MCF7_chr10_rep1.bam MCF7_chr10_rep2.bam" -t 12 -n "test" -L 1000
```

# Peak Calling <a name = "peakcalling" />

## Example

Train the model and obtain the predictions.

```
bash run_rcl.sh -p example -f "rep1 rep2"
```

## Command-Line Options

```
Input (required):
    --p 
        Path to preprocessing data.
    --f
        Names of the individual BAM files (without suffix). For example, if your BAM files are rep1.bam and rep2.bam, use "rep1 rep2"

Parameters (optional):
    --e  Training epoches.
        default=25
    --b Batch size.
        default=256
```

## Output

The trained model is called `rcl.ckpt` and results are stored in `rcl.bed`. The output will have 

*chromosome name, peak start position, peak end position, peak name, peak score, training region start position, training region end position*, for example
```
10      49829   50258   10segment1      0.18526842      49543   50543
10      73663   74515   10segment2      0.8270205       73589   74589
```

# How to Cite <a name = "cite" />
Vu, H. T., Zhang, Y., Tuteja, G., & Dorman, K. S. (2023). Unsupervised contrastive peak caller for ATAC-seq. Genome Research, gr-277677.   
Link to article: https://genome.cshlp.org/content/33/7/1133.full

# Contact <a name = "contact" />

Ha Vu (hhvu@iastate.edu or vthihong@umich.com), Yudi Zhang (yudiz@iastate.edu)

