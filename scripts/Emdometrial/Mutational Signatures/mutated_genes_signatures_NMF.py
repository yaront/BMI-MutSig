#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 08:27:36 2018

@author: tomer
"""

#%%
##===========================================##
## Mutated Gene Signatures Deciphering - NMF ##
##===========================================##

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from matplotlib import pyplot as plt
import matplotlib.cm as cm

#%% Reading the data

gene_bmi_mut = pd.read_table('./../../../databases/Endometrial/mutation_gene_bmi/UCEC_bmi_gene_mut.txt', sep = '\t', index_col = 0)
mut_gene_mat = np.asmatrix(gene_bmi_mut.iloc[2:,:-1])


#%% Dimension reduction

#percent_remove = 0
#
#mut_freq = mut_cat.sum(axis=1)/mut_cat.sum().sum()
#mut_freq = mut_freq.sort_values()
#
#low_mut_type = list(mut_freq[np.cumsum(mut_freq) < percent_remove].index)
#mut_cat_high = mut_cat.drop(low_mut_type)
#mut_cat_high_mat = np.asmatrix(mut_cat_high)


#%% Bootstrap



n_comp = 4
row_plt = 4
col_plt = 1

X = np.copy(mut_gene_mat)
#X = mut_cat_high_mat[:,range(50) + range(190,240)]
#X = mut_cat_high_mat.astype(np.float) / mut_cat_high_mat.sum(axis=0)
#X = X[:,range(50) + range(190,240)]
nmf = NMF(n_components=n_comp)

W = nmf.fit_transform(X)
H = nmf.components_

mut_sig = W.argmax(axis = 1)
mut_clust_mat = X[mut_sig.argsort()]
mut_clust = W[mut_sig.argsort()]
pat_sig = H.argmax(axis = 0)
pat_clust_mat = X[:,pat_sig.argsort()]
pat_clust = H[:,pat_sig.argsort()]

plt.figure()
plt.imshow(W,cmap = cm.gray_r)
plt.figure()
plt.imshow(H,cmap = cm.gray_r)
#plt.figure()
#plt.imshow(mut_clust_mat,cmap = cm.gray_r)
#plt.figure()
#plt.imshow(pat_clust_mat,cmap = cm.gray_r)

pat_sig

#%%

plt.figure()
plt.subplot(row_plt,col_plt,1)
plt.bar(range(17740), W[:,0])
for i in range(1,W.shape[1]):
    plt.subplot(row_plt,col_plt,i+1)
    plt.bar(range(17740),W[:,i])
plt.suptitle('Mutational Signatures; n = ' + str(n_comp))


