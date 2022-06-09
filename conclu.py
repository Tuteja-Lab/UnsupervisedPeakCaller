#!/usr/bin/env python
# coding: utf-8

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset
import math
from random import randrange
import torch.nn.functional as F
from torch.autograd import Variable
from scipy.signal import gaussian
from torch.nn.functional import normalize

device = 'cuda' if torch.cuda.is_available() else 'cpu'


def preprocess(input_dat, gaussian_size=100, smoothing=True, normalize=False, eps=0.00001):
    if smoothing:
        smoothing_filter = gaussian(gaussian_size, 50)/np.sum(gaussian(gaussian_size,50))
        smoothing_filter = torch.FloatTensor(smoothing_filter).cuda()
        smoothing_filter = smoothing_filter.unsqueeze(0).unsqueeze(0)
        cov = nn.Conv1d(1, 1, 
                      kernel_size=gaussian_size, padding='same', bias = False)
        with torch.no_grad():
            cov.weight = torch.nn.Parameter(smoothing_filter)
        out_data = cov(input_dat)
        
    if normalize:
        std, mean = std_mean(input_dat, dim = 2, keepdim = True)
        out_data = (input_dat-mean)/(std+eps)
        
    return out_data

    
class Activation(nn.Module):
    """Configurable activation layer."""

    def __init__(self, afunc='relu'):
        """Initialize layer.
        Args:
            afunc : Type of activation function.
        """
        super(Activation, self).__init__()
        self.act_layer = nn.Identity()
        if afunc == 'relu':
            self.act_layer = nn.ReLU()
        elif afunc == 'prelu':
            self.act_layer = nn.PReLU()
        elif afunc is not None:
            raise NotImplementedError

    def forward(self, x):
        """Execute layer on input.
        Args:
            x : Input data.
        """
        return self.act_layer(x)
    
class L2Norm(nn.Module):
    def forward(self, x):
        return x / x.norm(p=2, dim=2, keepdim=True)
    
class ConvAct1d(nn.Module):
    """1D conv layer with same padding.
    Optional batch normalization and activation layer.
    """

    def __init__(self, in_channels=1, out_channels=5,
                 kernel_size=41, stride=1, dilation=1, bias=False,
                 bn=True, afunc='relu', padding='same', encoder = True):
        super(ConvAct1d, self).__init__()

        if encoder:
            self.conv_layer = nn.Conv1d(
            in_channels, out_channels, kernel_size, stride, padding=padding,
            dilation=dilation, bias=bias)
        else:
            self.conv_layer = nn.ConvTranspose1d(
            in_channels, out_channels, kernel_size, stride, 
            padding=0, dilation=dilation, bias=bias)
        self.bn_layer = nn.BatchNorm1d(out_channels) if bn else None
        self.act_layer = Activation(afunc) if afunc else None

    def forward(self, x):
        """Execute layer on input.
        Args:
            x : Input data.
        """
        
        x = self.conv_layer(x)
        if self.bn_layer:
            x = self.bn_layer(x)
        if self.act_layer:
            x = self.act_layer(x)
        return x
    
class ResBlock(nn.Module):
    """Residual block.
    2 conv/activation layers followed by residual connection
    and third activation.
    """

    def __init__(self, in_channels=1, out_channels=5, kernel_size=41,
                 stride=1, dilation=1, bias=False, bn=True,
                 afunc='relu', conv_input=False, padding='same',
                 downsample = False):
        """Initialize layer.
        Args:
            in_channels : Input channels.
            out_channels : Output channels.
            kernel_size : Filter size.
            stride : Filter stride.
            dilation : Dilation for filter.
            bias : Conv layer bias
            bn : Enable batch norm
            afunc : Activation function
        """
        super(ResBlock, self).__init__()
        
        self.downsample = downsample
        if conv_input:
            self.conv_input = ConvAct1d(in_channels, out_channels, kernel_size=1,
                                        bn=bn, afunc=afunc, padding=padding)
        else:
            self.conv_input = nn.Identity()
        self.conv_act1 = ConvAct1d(
            in_channels, out_channels, kernel_size,
            stride, dilation, bias, bn, afunc, padding=padding)
        self.conv_act2 = ConvAct1d(
            out_channels, out_channels, kernel_size,
            stride, dilation, bias, bn, afunc, padding=padding)
        self.conv_act3 = ConvAct1d(
            out_channels, out_channels, kernel_size,
            stride, dilation, bias, bn, afunc=None, padding=padding)
        self.activation = nn.ReLU()
        self.downs = ConvAct1d(
            out_channels, out_channels, kernel_size,
            stride, dilation, bias, bn, afunc, padding=0)

    def forward(self, input):
        """Execute layer on input.
            input : Input data.
        """
        
        x = self.conv_act1(input)
        x = self.conv_act2(x)
        x = self.conv_act3(x)
        
        resid = self.conv_input(input)
        x = x + resid
        x = self.activation(x)
        
        if self.downsample is True:
            x = self.downs(x)

        return x 

class ResNetEmbedding(nn.Module):
    """Resnet model."""

    def __init__(self, input_size, embedding_size, in_channels=1, out_channels=15,
                 num_blocks=5, kernel_size=31, dilation=8, bn=False, afunc='relu',
                 num_blocks_class=2, kernel_size_class=31, dilation_class=8,
                 out_channels_class=15, padding='same'):
        super(ResNetEmbedding, self).__init__()
        self.rep_dim = embedding_size
        self.res_blocks = nn.ModuleList()
        self.res_blocks_class = nn.ModuleList()

        # Residual blocks for regression
        self.res_blocks.append(
            ResBlock(in_channels, out_channels, kernel_size,
                     dilation=dilation, bn=bn, afunc=afunc,
                     conv_input=True, padding=padding))
        for _ in range(num_blocks - 1):
            self.res_blocks.append(
                ResBlock(out_channels, out_channels,
                         kernel_size,
                         dilation=dilation, bn=bn, afunc=afunc,
                         conv_input=False, padding=padding))
        self.regressor = ConvAct1d(in_channels=out_channels,
                                   out_channels=1, kernel_size=1, dilation=1,
                                   bn=bn, afunc=afunc, padding=padding)
        # Residual blocks for classification
        self.res_blocks_class.append(ResBlock(in_channels=1,
                                              out_channels=out_channels_class,
                                              kernel_size=kernel_size_class,
                                              dilation=dilation_class, bn=bn,
                                              afunc=afunc, conv_input=True,
                                              bias=True, padding=padding))
        for _ in range(num_blocks_class - 1):
            self.res_blocks_class.append(
                ResBlock(out_channels_class, out_channels_class,
                         kernel_size_class, dilation=dilation_class, bn=bn,
                         afunc=afunc, conv_input=False, bias=True, padding=padding))
        self.classifier = ConvAct1d(in_channels=out_channels,
                                    out_channels=1, kernel_size=1, dilation=1,
                                    bn=bn, afunc=None, bias=True, padding=padding)
        self.encoder = nn.Sequential(
                        nn.Flatten(),
                        nn.Linear(input_size, embedding_size))

    def forward(self, x):
        
        for res_block in self.res_blocks:
            x = res_block(x)
        x = self.regressor(x)
        decoded = x
        for res_block in self.res_blocks_class:
            x = res_block(x)
        embedding = self.encoder(self.classifier(x))

        return embedding.unsqueeze(1), decoded
    
class RegionEmbedding(nn.Module):       
        
    def __init__(self, internal_embedding_size, embedding_size, 
                 first_kernel_size, dropout_rate, input_dim):
        super().__init__()
        self.rep_dim = embedding_size
        hidden_kernel_size = 9
        hidden_embedding_size = 10
        pool_size=5
        dilation = 8

        self.embedding = nn.Sequential(
            nn.Conv1d(input_dim, internal_embedding_size, 
                      kernel_size=first_kernel_size, padding='same'),
            nn.BatchNorm1d(internal_embedding_size),
            nn.ReLU(),
            
            
            nn.Conv1d(internal_embedding_size, hidden_embedding_size,
                      kernel_size=hidden_kernel_size, padding='same'),
            nn.BatchNorm1d(hidden_embedding_size),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=pool_size),

            nn.Conv1d(hidden_embedding_size, hidden_embedding_size, 
                      dilation=dilation, kernel_size=hidden_kernel_size, padding='same'),
            nn.BatchNorm1d(hidden_embedding_size),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=pool_size),

            nn.Conv1d(hidden_embedding_size, hidden_embedding_size, 
                      kernel_size=hidden_kernel_size, padding='same'),
            nn.BatchNorm1d(hidden_embedding_size),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=pool_size),

            nn.Flatten(),
            nn.Dropout(dropout_rate),
            nn.LazyLinear(out_features=self.rep_dim),
            nn.ReLU()
          
            # nn.LazyLinear(out_features=embedding_size)
            # L2Norm()
            )

    def forward(self, X):
        
        embedding = self.embedding(X)
        
        return embedding.unsqueeze(1)


# In[10]:
class googlenet(nn.Module):
    def __init__(self, internal_embedding_size, embedding_size, 
                 dropout_rate, input_dim):
        super().__init__()
        self.rep_dim = embedding_size
        pool_size=5
        self.cov1 = nn.Sequential(
            nn.Conv1d(input_dim, internal_embedding_size, 
                      kernel_size=15, padding=1),
            nn.BatchNorm1d(internal_embedding_size),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=pool_size))
        
        self.cov2 = nn.Sequential(
            nn.Conv1d(input_dim, internal_embedding_size, 
                      kernel_size=30, padding=1),
            nn.BatchNorm1d(internal_embedding_size),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=pool_size))
        
        self.cov3 = nn.Sequential(
            nn.Conv1d(input_dim, internal_embedding_size, 
                      kernel_size=50, padding=1),
            nn.BatchNorm1d(internal_embedding_size),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=pool_size))
        
        self.drop = nn.Dropout(dropout_rate)
        self.fc1 = nn.LazyLinear(out_features=2 * self.rep_dim)
        self.fc2 = nn.Linear(2 * self.rep_dim, self.rep_dim)
        
    def forward(self, X):
        
        cov1 = self.cov1(X)
        cov2 = self.cov2(X)
        cov3 = self.cov3(X)
        
        x = torch.cat((cov1, cov2, cov3), 2)
        x = torch.flatten(x, 1)

        x = self.drop(x)
        x = self.fc1(x)
        x = F.relu(x)
        x = self.drop(x)
        embedding = self.fc2(x)

        return embedding.unsqueeze(1)

class PrintLayer(nn.Module):
    def __init__(self):
        super(PrintLayer, self).__init__()

    def forward(self, x):
        for t in x:
            if t.min() < 0:
                print(x.size())
                print('smaller than 0', t.min())

        return x
    
    
class ResDecBlock(nn.Module):

    def __init__(self, in_channels=5, out_channels=5, kernel_size=41,
                 stride=1, dilation=1, bias=False, bn=True,
                 afunc='relu', padding='same'):

        super(ResDecBlock, self).__init__()

        self.deconv = ConvAct1d(
            in_channels, out_channels, kernel_size,
            stride, dilation, bias, bn, afunc, encoder = False)
 
        self.conv_act1 = ConvAct1d(
            in_channels, out_channels, kernel_size,
            stride, dilation, bias, bn, afunc, padding=padding)
        self.conv_act2 = ConvAct1d(
            out_channels, out_channels, kernel_size,
            stride, dilation, bias, bn, afunc, padding=padding)
        self.conv_act3 = ConvAct1d(
            out_channels, out_channels, kernel_size,
            stride, dilation, bias, bn, afunc, padding=padding)
      
    
    def forward(self, input):
        """Execute layer on input.
        Args:
            input : Input data.
        """

        x = self.deconv(input)
        resid = x
        x = self.conv_act1(x)
        x = self.conv_act2(x)
        x = self.conv_act3(x)

        x = x + resid
     
        return x
    
class dblock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size,
                 dilation=8, bias=False, bn=False,
                 afunc='relu', num_blocks=5):
        super().__init__()
        
        # decoder
        self.res_dec = nn.ModuleList()

        for _ in range(num_blocks - 1):
            self.res_dec.append(
                ResDecBlock(in_channels, in_channels,
                         kernel_size,
                         dilation=dilation, bn=bn, afunc=afunc))
        
        self.res_dec.append(ConvAct1d(in_channels,
                            out_channels, kernel_size=kernel_size, dilation=dilation,
                            bn=bn, afunc=afunc))


    def forward(self, x):
        
        for res_block in self.res_dec:
            x = res_block(x)
        return x


class ResAE(nn.Module):
    def __init__(self, out_channels, embedding_size, cluster_num, 
                 in_channels = 1, kernel_size = 31, 
                 dilation=8, bias=False, bn=False, input_size = (1, 1000),
                 afunc='relu', num_blocks=5, bottleneck_kernel = 13, flat = False):
        super().__init__()
        self.rep_dim = embedding_size
        self.cluster_num = cluster_num
        self.input_size = input_size
        self.flat = flat
        # encoder
        self.encoder = nn.ModuleList()
        self.encoder.append(ResBlock(in_channels, out_channels, kernel_size,
                     dilation=dilation, bn=bn, afunc=afunc, conv_input=True))
       
        for _ in range(num_blocks - 1):
            self.encoder.append(
                ResBlock(out_channels, out_channels,
                         kernel_size,
                         dilation=dilation, bn=bn, afunc=afunc,
                         conv_input=False, downsample=True))
        
        # decoder
        self.decoder = dblock(in_channels=out_channels, out_channels=in_channels, 
                              kernel_size=kernel_size, dilation=dilation, bn=bn, afunc=afunc)
        
        self.flat_fts = self.get_flat_fts(self.encoder)
        
        self.unflat = nn.Sequential(
            nn.Linear(self.rep_dim, self.flat_fts),
            nn.Unflatten(1, (out_channels, int(self.flat_fts/out_channels))),
            nn.ReLU())
        
        self.instance_projector = nn.Sequential(
            nn.Flatten(),
            nn.Linear(self.flat_fts, self.rep_dim),
            nn.BatchNorm1d(self.rep_dim),
            nn.ReLU()
        )
        self.cluster_projector = nn.Sequential(
            nn.Linear(self.rep_dim, self.cluster_num),
            nn.Softmax(dim=2)
        )
        
    def get_flat_fts(self, fts):
        tmp = Variable(torch.ones(1, *self.input_size))
        for res_block in self.encoder:
            tmp = res_block(tmp)

        return int(np.prod(tmp.size()[1:]))    
    
    def forward(self, x):
        for res_block in self.encoder:
            x = res_block(x)

        h = self.instance_projector(x)
        if self.flat:
            x = self.unflat(h)
            
        decoded = self.decoder(x)
        # print('encoded', x.size())
        # print('decoded', decoded.size())
        
        h = h.unsqueeze(1)
        # print('h', h.size())
        c = self.cluster_projector(h)
        return h, c, decoded
    
    
class DeepAutoencoder(nn.Module):
    def __init__(self, internal_embedding_size, embedding_size, 
                 first_kernel_size, input_size):
        super().__init__()   
        
        self.rep_dim = embedding_size
        hidden_kernel_size = 15 ## the kernel size - 1 has to be divisible by 2
        hidden_embedding_size = 10
    
        self.encoder = nn.Sequential(
            nn.Conv1d(1, internal_embedding_size, 
                      kernel_size=first_kernel_size, padding='same'),
            nn.BatchNorm1d(internal_embedding_size),
            nn.ReLU(),
            
            nn.Conv1d(internal_embedding_size, hidden_embedding_size,
                      kernel_size=hidden_kernel_size, padding='same'),
            nn.BatchNorm1d(hidden_embedding_size),
            nn.ReLU(),
            nn.Conv1d(hidden_embedding_size, hidden_embedding_size, 
                      kernel_size=hidden_kernel_size, padding='same'),
            nn.BatchNorm1d(hidden_embedding_size),
            nn.ReLU(),
            nn.Conv1d(hidden_embedding_size, hidden_embedding_size, 
                      kernel_size=hidden_kernel_size, padding='same'),
            nn.BatchNorm1d(hidden_embedding_size),
            nn.ReLU(),
            nn.Flatten(),
            nn.Linear(hidden_embedding_size * input_size, self.rep_dim)
            )
        self.decoder = nn.Sequential(
            nn.Linear(self.rep_dim, hidden_embedding_size * input_size),
            nn.Unflatten(1, (hidden_embedding_size, input_size)),
 
            nn.ConvTranspose1d(hidden_embedding_size, hidden_embedding_size, 
                        kernel_size=hidden_kernel_size, padding=(hidden_kernel_size - 1) // 2),
            nn.BatchNorm1d(hidden_embedding_size),
            nn.ReLU(),

            nn.ConvTranspose1d(hidden_embedding_size, hidden_embedding_size, 
                        kernel_size=hidden_kernel_size, padding=(hidden_kernel_size - 1) // 2),
            nn.BatchNorm1d(hidden_embedding_size),
            nn.ReLU(),

            nn.ConvTranspose1d(hidden_embedding_size, internal_embedding_size, 
                        kernel_size=hidden_kernel_size, padding=(hidden_kernel_size - 1) // 2),
            nn.ReLU(),
            nn.BatchNorm1d(internal_embedding_size),

            nn.ConvTranspose1d(internal_embedding_size, 1, 
                        kernel_size=first_kernel_size, padding=(first_kernel_size - 1) // 2),
            nn.ReLU(),

            PrintLayer()
        )
  
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return encoded.unsqueeze(1), decoded

class Network2(nn.Module):
    'deep dembdeeding (encoder and decoder)'
    def __init__(self, embnet, feature_dim, class_num):
        super(Network2, self).__init__()
        self.embnet = embnet
        self.feature_dim = feature_dim
        self.cluster_num = class_num
#         self.instance_projector = nn.Sequential(
#             nn.Linear(self.embnet.rep_dim, self.embnet.rep_dim),
#             nn.ReLU(),
#             nn.Linear(self.embnet.rep_dim, self.feature_dim),
#         )
        self.cluster_projector = nn.Sequential(
            nn.Linear(self.embnet.rep_dim, self.embnet.rep_dim),
            nn.ReLU(),
            nn.Linear(self.embnet.rep_dim, self.cluster_num),
            nn.Softmax(dim=2)
        )

    def forward(self, x):
        h, decoded = self.embnet(x)
#         z = self.instance_projector(h)
        c = self.cluster_projector(h)
        return h, c, decoded


class Network(nn.Module):
    def __init__(self, embnet, feature_dim, class_num):
        super(Network, self).__init__()
        self.embnet = embnet
        self.feature_dim = feature_dim
        self.cluster_num = class_num
        self.instance_projector = nn.Sequential(
            nn.Linear(self.embnet.rep_dim, self.embnet.rep_dim),
            nn.ReLU(),
            nn.Linear(self.embnet.rep_dim, self.feature_dim),
        )
        self.cluster_projector = nn.Sequential(
            nn.Linear(self.embnet.rep_dim, self.embnet.rep_dim),
            nn.ReLU(),
            nn.Linear(self.embnet.rep_dim, self.cluster_num),
            nn.Softmax(dim=2)
        )

    def forward(self, x):
        h = self.embnet(x)
        z = self.instance_projector(h)
        c = self.cluster_projector(h)

        return z, c

#     def forward_cluster(self, x):
#         h = self.embnet(x)
#         c = self.cluster_projector(h)
#         c = torch.argmax(c, dim=2)
#         return c

def poisson_loss(y_true, y_pred):
    y_pred = tf.cast(y_pred, tf.float32)
    y_true = tf.cast(y_true, tf.float32)

    # we can use the Possion PMF from TensorFlow as well
    # dist = tf.contrib.distributions
    # return -tf.reduce_mean(dist.Poisson(y_pred).log_pmf(y_true))

    nelem = _nelem(y_true)
    y_true = _nan2zero(y_true)

    # last term can be avoided since it doesn't depend on y_pred
    # however keeping it gives a nice lower bound to zero
    ret = y_pred - y_true*tf.math.log(y_pred+1e-10) + tf.math.lgamma(y_true+1.0)

    return tf.divide(tf.reduce_sum(ret), nelem)

class InstanceLoss(nn.Module):
    def __init__(self, temperature, n_views):
        super(InstanceLoss, self).__init__()
        self.n_views = n_views
        self.temperature = temperature
            
    def forward(self, emb_i, emb_j):
        """
        emb_i and emb_j are batches of embeddings, where corresponding indices are pairs
        z_i, z_j as per SimCLR paper
        """
        batch_size = emb_i.size()
        batch_size = batch_size[0]
        negatives_mask = (~torch.eye(batch_size * self.n_views, 
                                    batch_size * self.n_views, dtype=bool)).float().to(device)
        emb_i = emb_i.squeeze(1)
        emb_j = emb_j.squeeze(1)
        z_i = F.normalize(emb_i, dim=1)
        z_j = F.normalize(emb_j, dim=1)
        
        representations = torch.cat([z_i, z_j], dim=0)
        similarity_matrix = F.cosine_similarity(representations.unsqueeze(1), representations.unsqueeze(0), dim=2)
        sim_ij = torch.diag(similarity_matrix, batch_size)
        sim_ji = torch.diag(similarity_matrix, -batch_size)
        positives = torch.cat([sim_ij, sim_ji], dim=0)
        
        nominator = torch.exp(positives / self.temperature)
        denominator = negatives_mask * torch.exp(similarity_matrix / self.temperature)
    
        loss_partial = -torch.log(nominator / torch.sum(denominator, dim=1))
        loss = torch.sum(loss_partial) / (self.n_views * batch_size)
        return loss


# In[13]:


class ClusterLoss(nn.Module):
    def __init__(self, class_num, temperature, n_views):
        super(ClusterLoss, self).__init__()
        self.class_num = class_num
        self.n_views = n_views
        self.temperature = temperature
        self.register_buffer("negatives_mask", (~torch.eye(class_num * n_views, 
                                    class_num * n_views, dtype=bool)).float())
          
    def forward(self, emb_i, emb_j):
        
        emb_i = emb_i.squeeze(1)
        emb_j = emb_j.squeeze(1)
        
        # compute entropy of cluster assignment
        p_i = emb_i.sum(0).view(-1)
        p_i /= p_i.sum()
        ne_i = math.log(p_i.size(0)) + (p_i * torch.log(p_i)).sum()
        p_j = emb_j.sum(0).view(-1)
        p_j /= p_j.sum()
        ne_j = math.log(p_j.size(0)) + (p_j * torch.log(p_j)).sum()
        ne_loss = ne_i + ne_j

        # contrsative loss
        emb_i = emb_i.t()
        emb_j = emb_j.t()
        
        representations = torch.cat([emb_i, emb_j], dim=0)
        similarity_matrix = F.cosine_similarity(representations.unsqueeze(1), representations.unsqueeze(0), dim=2)
        sim_ij = torch.diag(similarity_matrix, self.class_num)
        sim_ji = torch.diag(similarity_matrix, -self.class_num)
        positives = torch.cat([sim_ij, sim_ji], dim=0)
        
        nominator = torch.exp(positives / self.temperature)
        denominator = self.negatives_mask * torch.exp(similarity_matrix / self.temperature)
    
        loss_partial = -torch.log(nominator / torch.sum(denominator, dim=1))
        loss = torch.sum(loss_partial) / (self.n_views * self.class_num)
        
        return loss + ne_loss
  



