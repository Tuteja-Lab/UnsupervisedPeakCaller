# Unsupervised Contrastive Peak Caller
This page described the Unsupervised contrastive peak caller known as RCL (Replicative Contrastive Learner).
The accompanying publication is available: [10.1101/gr.277677.123](https://doi.org/10.1101/gr.277677.123).

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Installation](#installation)
1. [Input](#input)
	1. [Example](#data_example)
1. [Preprocessing](#preprocessing)
	1. [Example](#preprocessing_example)
	1. [Output](#preprocessing_output)
1. [Peak Calling](#peakcalling)
	1. [Example](#peakcalling_example)
	1. [Command-Line Options](#peakcalling_options)
	1. [Output](#peakcalling_output)
1. [How to Cite](#cite)
1. [Contact](#contact)

# Prerequisites <a name = "prerequisites" />
For input preprocessing steps, the following tools and R libraries are required:
- [bash](https://www.gnu.org/software/bash/) (>= 5.2)
- [coreutils](https://www.gnu.org/software/coreutils/coreutils.html) (>= 9.3)
- [perl](https://www.perl.org/) (>= 5.38)
- [samtools](https://github.com/samtools/samtools) (>= 1.10)
- [bedtools2](https://github.com/arq5x/bedtools2) (>= 2.27.1)
- [parallel](https://www.gnu.org/software/parallel/) (>= 20170322)
- [bedops](https://github.com/bedops/bedops) (>= 2.4.35)
- [R](https://www.r-project.org/) (>= 4.0.2)
- R library [dplyr](https://dplyr.tidyverse.org/) (>= 1.0.7)
- R library [bedr](https://cran.r-project.org/web/packages/bedr/index.html) (>= 1.0.7)
- R library [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) (>= 1.0.16)

For the deep learner step, a **GPU** is needed. Other packages needed are:

- [Python](https://www.python.org/) (>=3.7.10)
- [PyTorch Lightning](https://lightning.ai/docs/pytorch/stable/) (>=1.5.1)
- [PyTorch](https://pytorch.org/) (>=1.10.0)
- [numpy](https://numpy.org/) (>=1.21.5)
- [pandas](https://pandas.pydata.org/) (>=1.3.5)
- [argparse](https://docs.python.org/library/argparse.html) (>=1.1)
- [scikit-learn](https://scikit-learn.org/stable/) (>=1.0.1)


# Installation <a name = "installation" />

After installing the prerequisites, all you have to do is clone RCL:

```
git clone https://github.com/Tuteja-Lab/UnsupervisedPeakCaller.git
cd UnsupervisedPeakCaller
```

# Input <a name = "input" />

The RCL preprocessor requires [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) files for each replicate.
The RCL program requires the [output of the preprocessing step](#preprocessing_output).

## Example <a name = "data_example" />
To demonstrate RCL, we provide the portion of the [MCF-7](https://www.encodeproject.org/search/?type=Experiment&searchTerm=ENCSR422SUG) dataset aligning to human chromosome 10.
The preprocessed data for this example are provided with RCL, so if you want to skip data preprocessing for now, you can go directly to [peak calling](#peakcalling).

Continuing with the preprocessing demonstration, download the necessary [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) files and indices from [https://iastate.box.com/s/9uavg2zsy5w0i7v227ei7yaeea7yktr6](https://iastate.box.com/s/9uavg2zsy5w0i7v227ei7yaeea7yktr6).
If you download the zip file ```RCLexamples.zip``` from the cybox and place it in the root of the RCL git repository, the following commands (executed from the root of the RCL git repository) will place them appropriately:

```
unzip -j RCLexamples.zip -d example
```

# Preprocessing <a name = "preprocessing" />
We have provided a [bash](https://www.gnu.org/software/bash/) preprocessing script to convert input BAM files (see [input](#input)) into the required RCL input.
The script assumes your data have been aligned to the Ensembl assembly of the mouse or human genome.
If not, no blacklist regions will be removed.
For more information about the preprocessing script type ```bash ./preprocessing -?``` from the RCL git repository root.

## Input <a name = "preprocessing_input" />
The preprocessing input is defined in [Input](#input).

## Example <a name = "preprocessing_example" />
After following the example instructions in [input](data_example) to get and place the data, you can run the preprocessing tool on the sample data as demonstrated below.
It is assumed you are at the root of the RCL git repository when you type this command.
Also, you should choose an appropriate number of threads for your system via option -t.

```
bash ./preprocessing.bash -d example -b "MCF7_chr10_rep1.bam MCF7_chr10_rep2.bam" -t 12 -n test
```

## Output of preprocessing <a name = "preprocessing_output" />
The final output of the preprocessing consists of two types of [Bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) files.

	1. [Bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)-formatted (5-column) for each of the replicates and each of the reference sequences containing the coverage for each fixed-length fragment chosen by [preprocessing](#preprocessing) for input to RCL. Specifically, the coverage for replicate REP in reference sequence SEQ is stored in DIR/SEQ/REP.covBga.txt, where REP is the basename (without .bam extension) of either [input](#input) BAM file, SEQ is a reference sequence found in these BAM files, and DIR is the input directory passed ```preprocessing.bash``` via command option ```-d```.
	1. [Bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)-formatted (4-column) containing the candidate peak regions that RCL will score.

# Peak Calling <a name = "peakcalling" />
We provide a [bash](https://www.gnu.org/software/bash/) script ```run_rcl.sh``` that fits RCL and assigns predictive peaks scores to each of the candidate regions prepared by [Preprocessing](#preprocessing).

## Input <a name = "rcl_input" />
The input to peak calling is the [output](#preprocessing_input) of [Preprocessing](#preprocessing).
It is not necessary to use the provided preprocessing script if you want to prepare input in some other way.
Identification of good candidate regions is an important part to the success of RCL.

## Example

Here we demonstrate how to train the model and score the candidate regions for the MCF7, chromosome 10 data.
We assume you run this script from the root of the RCL git repository.

```
bash run_rcl.sh -d example -b "MCF7_chr10_rep1.bam MCF7_chr10_rep2.bam"
```

## Command-Line Options <a name = "peakcalling_options" />

For more information about how to run the script, type ```run_rcl.sh -?```.
Parameters that might be important to adjust for obtaining better fits are the number of epochs (command line option ```-e```) and the batch size (command line option ```-s```).

## Output <a name = "rcl_output" />

The trained model is called `rcl.ckpt` and results are stored in `rcl.bed`, both in the directory passed to the ```-d``` command-line option. The output will have 

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

