#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 23:45:23 2018

@author: tomer
"""

#%%
##=========================================##
## Mutational Signatures Deciphering - NMF ##
##=========================================##

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF

#%% Reading the data

mut_cat = pd.read_table('./../../databases/mutation_signatures/Endometrial_mutational_catalog.txt', index_col = 0)
mut_cat_mat = np.asmatrix(mut_cat)


#%% Dimension reduction

n_comp = 4

X = np.copy(mut_cat_mat)
nmf = NMF(n_components=n_comp)

W = nmf.fit_transform(X)
H = nmf.components_


res_labels = W.argmax(axis = 1)
res_clust_mat = X[res_labels.argsort()]
w_clust = W[res_labels.argsort()]
kin_labels = H.argmax(axis = 0)
kin_clust_mat = X[:,kin_labels.argsort()]
h_clust = H[:,kin_labels.argsort()]

kin_labels
