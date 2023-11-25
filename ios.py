#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset
import math
from random import randrange
device = 'cuda' if torch.cuda.is_available() else 'cpu'



## read coverage in format 
## pos cov region mac_lab
## 3106561	5	region_2709 0
def read_coverage(file, name = 'lab', pos = 3):
  # name can be 'lab', 'cov', 'pos' 'region'
    out = []

    last_name = None
    i = 0
    with open(file, 'r') as f:
        for line in f.readlines():   
            tmp = line.split()
            if name == 'lab':
                v = int(tmp[pos])
            elif name == 'cov':
                v = int(tmp[1])
            elif name == 'pos':
                v = int(tmp[0])
            elif name == 'region':
                v = tmp[2]  
            out.append([v, tmp[2]])
  
    last_name = out[0][1]
    dat = []
    last_i = 0
    i = 0
    for l in out:
        name = l[1]
        # if cnt == 1:
        #     tmp = [j[0] for j in out[0:last_i]]
        #     dat.append(np.array(tmp))
        if name != last_name:
          # print(i, last_i)
            tmp = [j[0] for j in out[last_i:i]]
            dat.append(np.array(tmp))
            last_i = i
            i += 1
            last_name = name
        else:
            i += 1
    tmp = [j[0] for j in out[last_i:i]]
    dat.append(np.array(tmp))
    return np.array(dat)

def read_data_new(file):
    out = []
    last_name = None
    i = 0
    with open(file, 'r') as f:
        for line in f.readlines():   
            tmp = line.split()
            count = int(tmp[2]) - int(tmp[1])
            v = count * [int(tmp[4])]
            out.append([v, tmp[3]])

    last_name = out[0][1]
    dat = []
    last_i = 0
    i = 0
    for l in out:
        name = l[1]
        if name != last_name:
            tmp = [j[0] for j in out[last_i:i]]
            tmp = [item for sublist in tmp for item in sublist]
            dat.append(np.array(tmp))
            last_i = i
            i += 1
            last_name = name
        else:
            i += 1
    tmp = [j[0] for j in out[last_i:i]]
    tmp = [item for sublist in tmp for item in sublist]
    dat.append(np.array(tmp))
    return np.array(dat)

def read_fragment(file, pos = [1]):
  # pos means:
# 2nd col = counts of fragments of length <= 100bp (corresponding to nucleosome free regions), 
# 3rd col = counts of length between 180 and 247 bp (considered to be mononucleosomes), 
# 4th col = frags of length between 315 and 473 bp (considered to be dinucleosomes), 
# 5th col = frags between 558 and 615 bp (considered to be trinucleosomes), 
# 6th col = region name, 
# 7th col = total counts (no matter how long the fragments are)
# can catke multiple fragments
    out = []
    last_name = None
    i = 0
    with open(file, 'r') as f:
        for line in f.readlines():   
            tmp = line.split()
            v = [int(tmp[p]) for p in pos]
            out.append([v, tmp[5]])
    
    last_name = out[0][1]
    dat = []
    last_i = 0
    i = 0
    for l in out:
        name = l[1]
        if name != last_name:
            tmp = [j[0] for j in out[last_i:i]]
            dat.append(np.array(tmp))
            last_i = i
            i += 1
            last_name = name
        else:
            i += 1
    tmp = [j[0] for j in out[last_i:i]]
    dat.append(np.array(tmp))
    return np.array(dat)


# paste reps together and convert to tensor
# def combine_rep_lab(reps, lab, device):
#     out = []
    
#     rep = reps[0]
#     for r in reps[1:]:            
#         rep = np.dstack((rep, r))
#     rep = rep.transpose(0, 2, 1)
        
#     for r1, r2 in zip(rep, lab):
#         r2 = np.array(r2)
#         t1 = torch.from_numpy(r1).float().to(device)
#         t2 = torch.from_numpy(r2).float().to(device)
    
#         pair = (t1, t2.unsqueeze(0))
    
#         out.append(pair)
  
#     return out


def combine_rep_lab(reps, lab, device):
    out = []
    
    rep = reps[0]
    for r in reps[1:]:            
        rep = np.dstack((rep, r))
    rep = rep.transpose(0, 2, 1)

    lab = np.array(lab)
    lab = lab.transpose(1, 0)
    
    for r1, r2 in zip(rep, lab):
        
        t1 = torch.from_numpy(r1).float().to(device)
        t2 = torch.tensor(r2, dtype=torch.long, device=device)
        
        pair = (t1, t2)
    
        out.append(pair)
  
    return out


def combine_reps(reps, device, center = False):
    n_rep = len(reps)
    out = []
    
    if center == True:
        rep2 = reps[0]
        left = list(range(0, n_rep))
    else:
        anchor = randrange(0, n_rep)
        left = list(range(0, n_rep))
        del left[anchor]
        ## to trick the dataloader function, put the first replicate as data,
        ## all the other replicates as labels
        rep2 = reps[left[0]]
        rep1 = reps[anchor]
        
    if n_rep > 2:
        for r in left[1:]:            
            rep2 = np.dstack((rep2, reps[r]))
        rep2 = rep2.transpose(0, 2, 1)
    
    if center == True:
        rep1 = np.mean(rep2, axis = 1)

    if len(rep1.shape) == 3: ## when we add more features other than just coverage
        rep1 = rep1.transpose(0, 2, 1)                    
        
    if n_rep == 2 and len(rep1.shape) == 3:
        rep2 = rep2.transpose(0, 2, 1)
    
    for r1, r2 in zip(rep1, rep2):
        
        if len(rep1.shape) == 3 and n_rep > 2:
           
            dim1 = n_rep - 1
            dim2 = int(r2.shape[0]/dim1)
            r2 = r2.reshape(dim1, dim2, r2.shape[1])
    
        t1 = torch.from_numpy(r1).float().to(device)
        t2 = torch.from_numpy(r2).float().to(device)
        if len(rep1.shape) == 3:
            pair = (t1, t2)
        elif n_rep == 2:
            pair = (t1.unsqueeze(0), t2.unsqueeze(0))
        else:
            pair = (t1.unsqueeze(0), t2)

        out.append(pair)

    return out

