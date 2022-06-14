# UnsupervisedPeakCaller

## Example

Train the model 

`python main.py --datapath rep1.txt rep2.txt --modelpath rcl.ckpt`

Get the predictions

`python rcl_score.py --model rcl.ckpt --dpath rep1.txt rep2.txt`

