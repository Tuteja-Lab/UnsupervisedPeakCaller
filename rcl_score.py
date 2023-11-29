import numpy as np
import torch
import torch.nn as nn
import sys
from conclu import *
from main import *
from ios import *
import pickle
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc
import os
import argparse
import pandas as pd
from functools import reduce
import warnings

warnings.filterwarnings("ignore")
class Classify(nn.Module):
    def __init__(self, embeddings_model_path):
        super().__init__()
        
        base_model = ContrastLearn.load_from_checkpoint(embeddings_model_path).model
        
        self.encoder = base_model.encoder
        self.instance_projector = base_model.instance_projector
        self.clustering = base_model.cluster_projector
        
    def forward(self, x, *args):
        for res_block in self.encoder:
            x = res_block(x)
        h = self.instance_projector(x)
        h = h.unsqueeze(1)
        ass = self.clustering(h)    
        return ass

class Embed(nn.Module):
    def __init__(self, embeddings_model_path):
        super().__init__()
        
        base_model = ContrastLearn.load_from_checkpoint(embeddings_model_path).model
        
        self.encoder = base_model.encoder
        self.decoder = base_model.decoder
        self.instance_projector = base_model.instance_projector
        self.unflat = base_model.unflat
        
    def forward(self, x, *args):
        for res_block in self.encoder:
            x = res_block(x)
        h = self.instance_projector(x)
#         x = self.unflat(h)
        decoded = self.decoder(x)
        return h, decoded
    

##
# Oh, please, what is this doing?
# Massage replicated coverage of single segment into the input expected by RCL.
#
# @param pos	rows segment, columns coverage
# @param func	rcl prediction function
def get_conemb(pos, func):
    pos_tensor = []
    for t in pos:
        pos_tensor.append(torch.from_numpy(t).float().unsqueeze(0).unsqueeze(0))
    pos_tensor = torch.cat(pos_tensor)
    
    if len(pos_tensor.size()) == 4:
        pos_tensor = pos_tensor.squeeze(1)
    pos = func(pos_tensor)
    return pos   

##
# Assign peak label to that class most correlated with the highest signal.
#
# @param rep_class    class probabilities predicted by RCL
# @param alllab2      peaks called by sheer nucleotide coverage
def compute_m(rep_class, alllab2, debug=False):
    m = alllab2[0]
    final_lab = []
    for d in rep_class:
        lab1 = torch.argmax(d.squeeze(1), dim=1).tolist()
        f1 = f1_score(m, lab1)

        lab2 = torch.argmin(d.squeeze(1), dim=1).tolist()
        f2 = f1_score(m, lab2)

        if debug:
            print("compute_m:", f1, f2)

        if f2 > f1:
            lab = lab2
        else:
            lab = lab1

#         print("peak number ", lab.count(1))
        final_lab.append(lab)

    return final_lab

##
# Extract chr, start, end, and name of every distinct region in bed file.
def make_region(datapath):
    df = pd.read_csv(datapath, header = None, sep = "\t", names=["chr", "s", "e", "name", "count"])
    a = df.groupby("name", sort=False).s.idxmin()
    b = df.groupby("name", sort=False).e.idxmax()
    d = df.loc[a]                    # entries with minimum start for each segment
    d['e'] =  df.loc[b]['e'].values  # overwrite its end with maximum end for same segment
    d = d[d.columns[1:4]]            # keep only chr, s, e, and name
    d = d.reset_index(drop=True)
    return d


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='metric', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--model", type=str, default="rcl.ckpt", help="fitted RCL model from main.py")
    parser.add_argument("--dpath", type=str, help="directory containing data")
    parser.add_argument("--names", type=str, nargs="+", help="replicate names")
    parser.add_argument("--prefix", type=str, help="directory to place output") # output path
    parser.add_argument("--id", type=str)
    parser.add_argument("--psudo", type=int, default = 1)
    parser.add_argument("--preprocess_region", type=str, default = "None", help="Preprocessing regions in 4-column bed format.")
    parser.add_argument("--threshold", type=int, default = 30000)
    parser.add_argument("--debug", action = "store_true")

    args = parser.parse_args()

    #classification = Classify('chr' + str(args.id) + "/" + str(args.model)) ## this is training seperately
    classification = Classify(str(args.model))    ## train on all (80 or 90 %)
    #datapath = [args.dpath + '/' + f for f in os.listdir(args.dpath) if f.endswith('covBga.txt')]    
    datapath = [args.dpath + '/' + f + ".covBga.txt" for f in args.names]
    if args.debug:
       print("Using input data: " + str(datapath)) 
    dat = []        # tensor: dense coverage per segment per replicate
    dataf = []
    for file in datapath:
        dat.append(read_data_new(file))
        dataf.append(make_region(file))
    if args.debug:
        print("Coverage, data set 1:", dat[0])
        print("Regions:", dataf[0])

    # count nucleotides covering each segment in each replicate
    alllab2 = []
    if args.psudo:
        for d in dat:
            tmp = np.sum(d, axis = 1)
            alllab2.append(tmp > args.threshold) # this is just a rough threshold for determing peak, the label might be flipped
    else:
        alllab2 = pickle.load(open('chr' + str(args.id) + "/chip_nbl.p", "rb" )) # not applicable here

    rep_class = []
    for d in dat:
        rep_class.append(get_conemb(d, classification))

    if args.debug:
        print("rep_class: ", rep_class[0])

    # choose peak as class 0 or 1 based on which one has higher coverage
    #final_lab = compute_m(rep_class, alllab2, args.debug)
    # better yet: primitive peaks (based on nucleotide coverage) from replicate 1 as truth
    y_true = alllab2[0]

    pre_reg = []
    if args.preprocess_region != "None":
        if args.debug:
            print("Using preprocessing regions from " + args.preprocess_region)
        pre_reg = pd.read_csv(args.preprocess_region, sep="\t", names=["chr", "s", "e", "name"])
    
    dicts = []
    i = 1
    for p, o in zip(rep_class, dataf):
        p = p.squeeze(1).detach().cpu().numpy()
        if args.debug:
            print("y_scores:", p[:, 0])
        # decide which RCL label is peak based on correlation with "truth" labels
        y_scores = p[:, 0]
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        auc1 = auc(recall, precision)
        y_scores = p[:, 1]
        precision2, recall2, thresholds = precision_recall_curve(y_true, y_scores)
        auc2 = auc(recall2, precision2)
        if auc2 > auc1:
            if args.debug:
               print("Taking label 2: auc=", auc2)
            y_scores = p[:, 1]
        else:
            if args.debug:
               print("Taking label 1: auc=", auc1)
            y_scores = p[:, 0]

        d = {'chr': str(args.id), 'score': y_scores}
        
        df = pd.DataFrame(data=d)
        df = pd.concat([df, o], axis=1)[["chr", "s", "e", "name", "score"]]
        df = df.rename(columns={'score' : 'score' + str(i)})
        if args.debug:
            print(df)
        dicts.append(df)
        i += 1
        
    df = reduce(lambda df1, df2: pd.merge(df1, df2, on = ["chr", "s", "e", "name"]), dicts)    
    cols = df.columns[~df.columns.isin(["chr", "s", "e", "name"])]
    df["scores"] = df[cols].mean(axis=1)
    df = df[["chr", "s", "e", "name", "scores"]]
    # join with preprocessing regions
    if args.preprocess_region != "None":
        df['name'] = df['name'].str.rsplit(pat='_', n=1).str.get(0) # notice this, make sure the preprocessing use _ to seperate 
        df = df.merge(pre_reg, how ='left', on = "name")
        df.columns = ["chr", "s1", "e1", "name", "score", "chr1", "s", "e"]
        df = df[["chr", "s", "e", "name", "score", "s1", "e1"]]
    df.to_csv(str(args.prefix) + '/rcl_' + str(args.id) + '.bed', index = False, sep = "\t", header = None)



