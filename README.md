# Unsupervised Contrastive PeakCaller

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Preprocessing](#preprocessing)
1. [Peak Calling](#peakcalling)
1. [How to Cite](#cite)
1. [Contact](#contact)

# Prerequisites <a name = "prerequisites" />

You will need Python, PyTorch Lightning, and PyTorch.

# Preprocessing <a name = "preprocessing" />

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
