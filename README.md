# Unsupervised Contrastive PeakCaller

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Preprocessing](#preprocessing)
1. [Peak Calling](#peakcalling)
1. [How to Cite](#cite)
1. [Contact](#contact)

# Prerequisites <a name = "prerequisites" />
For input preprocessing steps, here are the modules used:
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
For the deep learner step, you will need Python, PyTorch Lightning, and PyTorch.

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

Train the model.

```
python main.py --datapath rep1.txt rep2.txt --modelpath rcl.ckpt
```

Get the predictions. For each replicate, the predicted scores and labels will be written in a file rcl.txt.

```
python rcl_score.py --model rcl.ckpt --dpath rep1.txt rep2.txt
```

## Command-Line Options

Training with **main.py**: 

```
Input (required):
    --datapath 
        Path to each of the preprocessed replicate files.
    --modelpath
        Path to trained model (default = model.ckpt).

Parameters:
    --epochs  Training epoches.
        default=25
    --lr      Convergence rate.
        default=1e-4
    --batch_size Batch size.
        default=256
    --temperature Temperature parameter.
        default=0.5
```

Obtain predictions with **rcl_score.py**:

```
Input (required):
    --dpath
        Path to each of the preprocessed replicate files.
    --model
        Trained model path.
    --prefix
        Prefix of the output (default = .).
```

# How to Cite <a name = "cite" />

# Contact <a name = "contact" />
