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
from sklearn.decomposition import PCA,NMF
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

#%% Reading the data

gene_bmi_mut = pd.read_table('./../../../databases/Endometrial/mutation_gene_bmi/UCEC_bmi_gene_mut.txt', sep = '\t', index_col = 0)
mut_gene_mat = np.asmatrix(gene_bmi_mut.iloc[2:,:-1])


#%% Dimension reduction

percent_remove = 0.1

mut_freq = mut_gene_mat.sum(axis=1)/float(mut_gene_mat.shape[1])

gene_bmi_mut_high = gene_bmi_mut.iloc[2:,:-1][mut_freq>percent_remove]
mut_gene_high_mat = np.asmatrix(gene_bmi_mut_high)


#%% Bootstrap



n_comp = 5
row_plt = 3
col_plt = 2

#for n_comp in range(2,6):

#X = np.copy(mut_gene_high_mat)
X = mut_gene_high_mat[:,range(50) + range(188,238)]
nmf = NMF(n_components=n_comp)

W = nmf.fit_transform(X)
H = nmf.components_

mut_sig = W.argmax(axis = 1)
mut_clust_mat = X[mut_sig.argsort()]
mut_clust = W[mut_sig.argsort()]
pat_sig = H.argmax(axis = 0)
pat_clust_mat = X[:,pat_sig.argsort()]
pat_clust = H[:,pat_sig.argsort()]

#    plt.figure()
#    plt.imshow(W,cmap = cm.gray_r)
plt.figure()
plt.imshow(H,cmap = cm.gray_r)
#plt.figure()
#plt.imshow(mut_clust_mat,cmap = cm.gray_r)
#plt.figure()
#plt.imshow(pat_clust_mat,cmap = cm.gray_r)

print gene_bmi_mut_high.index[W[:,0].argsort()[::-1][:5]]
print pat_sig

#%%

plt.figure()
plt.subplot(row_plt,col_plt,1)
plt.bar(range(17740), W[:,0])
for i in range(1,W.shape[1]):
    plt.subplot(row_plt,col_plt,i+1)
    plt.bar(range(17740),W[:,i])
    plt.suptitle('Mutational Signatures; n = ' + str(n_comp))


