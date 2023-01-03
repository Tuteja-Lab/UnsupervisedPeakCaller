# Unsupervised Contrastive PeakCaller

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Preprocessing](#preprocessing)
1. [Peak Calling](#peakcalling)
1. [How to Cite](#cite)
1. [Contact](#contact)

# Prerequisites <a name = "prerequisites" />
For input preprocessing steps, this following tools and R libraries are required:
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
The current pipeline has been tested on a HPC machine with the following modules loaded:
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
```

For the deep learner step, **GPU** is needed. Other packages needed are:
```
Python (>=3.7.10)
PyTorch Lightning (>=1.5.1)
PyTorch (>=1.10.0)
numpy (>=1.21.5)
pandas (>=1.3.5)
argparse (>=1.1)
sklearn (>=1.0.1)
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
        -m Bam files merged from individual replicates. Only used for preprocessing purpose, not for calling peaks. Must be sorted.
        -b Individual bam files of every replicates. Must be sorted.
        -t Number of threads to use.
        -n File name prefix.
        -L Length of input segments.
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
        Name of the replicate file (without suffix), for example rep 1 rep2

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

# Contact <a name = "contact" />

Yudi Zhang (yudiz@iastate.edu), Ha Vu (hhvu@iastate.edu)

