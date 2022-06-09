from torch.utils.data import Dataset, DataLoader, SubsetRandomSampler, SequentialSampler
import numpy as np
import math
import torch
import torch.nn as nn
import random
from torch.multiprocessing import cpu_count
from torch.optim import Adam
import pytorch_lightning as pl
from argparse import Namespace
import argparse
import random
import pickle
# from pytorch_lightning.loggers import WandbLogger

from conclu import *
from ios import read_data_new, read_covearge, read_fragment, combine_reps

device = 'cuda' if torch.cuda.is_available() else 'cpu'

# In[12]:
def align_loss(x, y, alpha=2):
    return (x - y).norm(p=2, dim=1).pow(alpha)

def uniform_loss(x, t=2):
    return torch.pdist(x, p=2).pow(2).mul(-t).exp().mean().log()

class ContrastLearn(pl.LightningModule):
    def __init__(self, hparams):
        hparams = Namespace(**hparams) if isinstance(hparams, dict) else hparams
        super().__init__()
        self.save_hyperparameters(hparams)
#         emb = ResNetEmbedding(self.hparams.input_size, self.hparams.emb_size)
#         emb = RegionEmbedding(self.hparams.hidden_size, self.hparams.emb_size, 
#                               self.hparams.first_kernel_size, self.hparams.dropout_rate,
#                               self.hparams.input_dim)
#         emb = googlenet(self.hparams.hidden_size, self.hparams.emb_size, self.hparams.dropout_rate,
#                               self.hparams.input_dim)
        if self.hparams.modelname == 'ResAE':
            self.model = ResAE(out_channels=self.hparams.hidden_size, embedding_size=self.hparams.emb_size, 
                               cluster_num=self.hparams.class_num, 
                               kernel_size=self.hparams.first_kernel_size, 
                               input_size = self.hparams.input_size)
        if self.hparams.modelname == 'CNNAE':
            emb = DeepAutoencoder(self.hparams.hidden_size, self.hparams.emb_size, 
                                  self.hparams.first_kernel_size, self.hparams.input_size)
            self.model = Network2(emb, self.hparams.fea_dim, self.hparams.class_num)
            
        self.ins_loss = InstanceLoss(self.hparams.ins_temp, self.hparams.n_views)
        self.clu_loss = ClusterLoss(self.hparams.class_num, self.hparams.clu_temp,
                                    self.hparams.n_views)
        self.AEloss = nn.MSELoss()

    def total_steps(self):
        return len(self.train_dataloader()) // self.hparams.epochs
    
    def train_dataloader(self):
        return DataLoader(rep_data,
                          batch_size=self.hparams.batch_size, 
                          sampler=SubsetRandomSampler(list(range(self.hparams.train_size))))
        
    def val_dataloader(self):
        return DataLoader(rep_data,
                      batch_size=self.hparams.batch_size, 
                      shuffle=False,
                      sampler=SequentialSampler(list(range(self.hparams.train_size + 1, 
                              self.hparams.train_size + self.hparams.validation_size))))
        
    def forward(self, X):
        return self.model(X)
    
    def step(self, batch, step_name = "train"):
        all_emb = []
        all_clu = []
        
        X, Y = batch
        loss = 0
        if self.hparams.smooth:
            X = preprocess(X)
        ## split to single replicate
        embX, cluX, decodedX = self.forward(X)
        all_emb.append(embX)
        all_clu.append(cluX)
        
        loss += self.hparams.beta * self.AEloss(X, decodedX)
        #print("shape of input: ", embX.size()) 
        ## compute uniform and align
        uni = 0
        aln = 0
        uni += torch.sum(uniform_loss(normalize(embX.squeeze(1))))

        for i in range(self.hparams.n_rep - 1):
            if self.hparams.n_rep == 2:
                rep = Y
            elif len(Y.size()) == 4: ## with segments (multiple features)
                rep = Y[:, i, :, :].squeeze(1)
            else: ## only coverage
                rep = Y[:, i, :].unsqueeze(1)
                
            if self.hparams.smooth:
                rep = preprocess(rep)
            # flip the tensor to make augmentation
#             X_flipped = torch.flip(X, [1, 2])
#             Y_flipped = torch.flip(rep, [1, 2])
            embY, cluY, decodedY = self.forward(rep)
            loss += self.hparams.beta * self.AEloss(rep, decodedY)
            all_emb.append(embY)
            all_clu.append(cluY)
            uni += torch.sum(uniform_loss(normalize(embY.squeeze(1))))       

        loss /= self.hparams.n_rep
        print('AE loss ', loss)
        
        orders = [(a, b) for idx, a in enumerate(range(self.hparams.n_rep)) for b in range(self.hparams.n_rep)[idx + 1:]]
        
        for ords in orders:
    #             embXf, _ = self.forward(X_flipped)
    #             embYf, _ = self.forward(Y_flipped)
            aln += torch.sum(align_loss(normalize(all_emb[ords[1]].squeeze(1)), normalize(all_emb[ords[0]].squeeze(1))))
            #print(normalize(all_emb[ords[1]].squeeze(1)))
            #print(normalize(all_emb[ords[0]].squeeze(1)))
            loss_instance = self.ins_loss(all_emb[ords[0]], all_emb[ords[1]])
            loss_cluster = self.clu_loss(all_clu[ords[0]], all_clu[ords[1]])
#             loss_aug1 = self.ins_loss(embX, embXf)
#             loss_aug2 = self.ins_loss(embY, embYf)
#             loss += loss_instance + loss_cluster + loss_aug1 + loss_aug2
            loss += loss_instance + loss_cluster
    
        aln /= len(orders)
        loss /= len(orders)    
        loss_key = f"{step_name}_loss"
        tensorboard_logs = {loss_key: loss}
        self.log('align', aln/self.hparams.train_size, on_epoch=True, prog_bar=True)
        self.log('uni', uni/self.hparams.train_size, on_epoch=True, prog_bar=True)

        return { ("loss" if step_name == "train" else loss_key): loss, 'log': tensorboard_logs,
                        "progress_bar": {loss_key: loss}}
    
    def training_step(self, batch, batch_idx):
        return self.step(batch, "train")

        
    def validation_step(self, batch, batch_idx):
        return self.step(batch, "val")
    
    def validation_end(self, outputs):
        if len(outputs) == 0:
            return {"val_loss": torch.tensor(0)}
        else:
            loss = torch.stack([x["val_loss"] for x in outputs]).mean()
            return {"val_loss": loss, "log": {"val_loss": loss}}

    def configure_optimizers(self):
        optimizer = Adam(self.model.parameters(), lr=self.hparams.lr)
        return [optimizer], []
    
class ContrastLearn_lab(pl.LightningModule):
    def __init__(self, hparams):
        hparams = Namespace(**hparams) if isinstance(hparams, dict) else hparams
        super().__init__()
        self.save_hyperparameters(hparams)
        emb = RegionEmbedding(self.hparams.hidden_size, self.hparams.emb_size, 
                              self.hparams.first_kernel_size, self.hparams.dropout_rate,
                              self.hparams.input_dim)
        self.model = Network(emb, self.hparams.fea_dim, self.hparams.class_num)
        # self.model = ImageEmbedding()
        self.ins_loss = InstanceLoss(self.hparams.ins_temp, self.hparams.n_views)
        self.clu_loss = nn.CrossEntropyLoss(weight = self.hparams.class_weights)

    def total_steps(self):
        return len(self.train_dataloader()) // self.hparams.epochs
    
    def train_dataloader(self):
        return DataLoader(rep_data,
                          batch_size=self.hparams.batch_size, 
                          sampler=SubsetRandomSampler(list(range(self.hparams.train_size))))
        
    def val_dataloader(self):
        return DataLoader(rep_data,
                      batch_size=self.hparams.batch_size, 
                      shuffle=False,
                      sampler=SequentialSampler(list(range(self.hparams.train_size + 1, 
                              self.hparams.train_size + self.hparams.validation_size))))
        
    def forward(self, X):
        return self.model(X)
    
    def step(self, batch, step_name = "train"):
        X, Y = batch
        Y.unsqueeze(1);
        loss = 0
        rep = X[:, 0, :].unsqueeze(1)
        embX, cluX = self.forward(rep)
        
        loss_cluster = self.clu_loss(cluX.squeeze(1), Y[:, 0])
        loss = loss_cluster
        ## split to single replicate
        for i in range(self.hparams.n_rep-1):
            rep = X[:, i+1, :].unsqueeze(1)
            embY, cluY = self.forward(rep)
            loss_instance = self.ins_loss(embX, embY)
            loss_cluster += self.clu_loss(cluY.squeeze(1), Y[:, i+1])
            loss += loss_instance + loss_cluster
        
        loss_key = f"{step_name}_loss"
        tensorboard_logs = {loss_key: loss}

        return { ("loss" if step_name == "train" else loss_key): loss, 'log': tensorboard_logs,
                        "progress_bar": {loss_key: loss}}
    
    def training_step(self, batch, batch_idx):
        return self.step(batch, "train")

        
    def validation_step(self, batch, batch_idx):
        return self.step(batch, "val")
    
    def validation_end(self, outputs):
        if len(outputs) == 0:
            return {"val_loss": torch.tensor(0)}
        else:
            loss = torch.stack([x["val_loss"] for x in outputs]).mean()
            return {"val_loss": loss, "log": {"val_loss": loss}}

    def configure_optimizers(self):
        optimizer = Adam(self.model.parameters(), lr=self.hparams.lr)
        return [optimizer], []
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='train', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--datapath', default=[], nargs='*')
    parser.add_argument('--fragpath', default=[], nargs='*')
    parser.add_argument('--modelpath', default='model.ckpt', type=str)
    parser.add_argument('--labpath', default='null', type=str)
    
    parser.add_argument('--model', default='ResAE', type=str) ## model name, can be ResAE, Resnet and CNNAE
    parser.add_argument('--epochs', default=25, type=int)
    parser.add_argument('--lr', default=1e-4, type=float)
    parser.add_argument('--batch_size', default=256, type=int)
    parser.add_argument('--emb_size', default=50, type=int) ## notice the dim after cov is 80
    parser.add_argument('--fea_dim', default=25, type=int)
    parser.add_argument('--first_kernel_size', default=31, type=int)
    parser.add_argument("--dropout_rate", default=0.05, type=float)
    parser.add_argument("--temperature", default=0.5, type=float)
    parser.add_argument("--hidden_size", default=5, type=int) ## num of hidden channels
    parser.add_argument("--gpus", default=1, type=int)
    parser.add_argument("--smooth", default=0, type=int)
    parser.add_argument("--seed", default=0, type=int)
    parser.add_argument("--n_class", default=2, type=int)
    parser.add_argument("--sample", default='null', type=str)

    args = parser.parse_args()
    
    random.seed(args.seed)
    print("Reading coverage files\n")
    n_rep = len(args.datapath)
    d = []
    for i in args.datapath:
        cov = read_data_new(i)
        d.append(cov)
    # test set
    if args.sample is not 'null':
        selected = np.random.choice(d[0].shape[0], int(d[0].shape[0] * 0.85), replace = False)
        pickle.dump(selected, open(str(args.sample) + ".p", "wb"))
        d = np.array(d)
        d = d[:, selected, :]
        d = list(d)

    n_dat = len(d[0])
    n_train = math.ceil(n_dat * 0.8)
    n_val = n_dat - n_train
    input_size = len(d[0][0])
    if args.model == 'ResAE':
        input_size = (1, input_size)
    input_dim = 1

    if len(args.fragpath) > 0:
        input_dim = 2
        print("Reading fragment length files\n")
        for i, f in enumerate(args.fragpath):
            fra = read_fragment(f, [1, 2, 3, 4, 6])
            d[i] = np.dstack((fra, d[i]))
       
    w_pos = 0
    w_neg = 0
    if args.labpath != 'null':
        lab = pickle.load(open(args.labpath, "rb"))
        posn = 0
        negn = 0
        for l in lab:
            negn += l.count(0)
            posn += l.count(1)
        w_pos = (posn + negn) / (2 * posn)
        w_neg = (posn + negn) / (2 * negn)
        rep_data = combine_rep_lab(d, lab, device=device)
    else:
        rep_data = combine_reps(d, device=device)
    class_weights = torch.FloatTensor([w_neg, w_pos]).to(device)
    print("weight ", class_weights)
    print("Finished reading\n")
    hparams = Namespace(lr=args.lr,
                        epochs=args.epochs,
                        batch_size=args.batch_size,
                        train_size=n_train,
                        validation_size=n_val,
                        hidden_size=args.hidden_size,
                        emb_size=args.emb_size,
                        fea_dim=args.fea_dim, ## dim of feature for computing loss
                        input_dim=input_dim,
                        input_size = input_size,
                        class_num=args.n_class,
                        ins_temp=args.temperature,
                        clu_temp=args.temperature,
                        n_views=2,## this is for pairwise comparison
                        n_rep = n_rep,
                        first_kernel_size = args.first_kernel_size, 
                        dropout_rate = args.dropout_rate,
                        device = device,
                        smooth = args.smooth,
                        class_weights = class_weights,
                        modelname = args.model,
                        beta = 1 ## penalty for encoder decoder, lower beta make things worse
                        )

    print("Start training\n")
    if args.labpath != 'null':
        module = ContrastLearn_lab(hparams)
    else:
        module = ContrastLearn(hparams)
    trainer = pl.Trainer(gpus=args.gpus, max_epochs=hparams.epochs)
    trainer.fit(module)
    checkpoint_file = args.modelpath
    trainer.save_checkpoint(checkpoint_file)
    
    
    
    
    
    
    
    
    
    
