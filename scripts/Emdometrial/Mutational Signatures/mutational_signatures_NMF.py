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
from matplotlib import pylab as plt
import matplotlib.cm as cm

#%% Reading the data

mut_cat = pd.read_table('./../../databases/mutation_signatures/Endometrial_mutational_catalog_bmi_sort.txt', index_col = 0)
mut_cat_mat = np.asmatrix(mut_cat)


#%% Dimension reduction

percent_remove = 0.1

mut_freq = mut_cat.sum(axis=1)/mut_cat.sum().sum()
mut_freq = mut_freq.sort_values()

low_mut_type = list(mut_freq[np.cumsum(mut_freq) < percent_remove].index)
mut_cat_high = mut_cat.drop(low_mut_type)
mut_cat_high_mat = np.asmatrix(mut_cat_high)


#%% Bootstrap



n_comp = 10

#X = np.copy(mut_cat_high_mat)
#X = mut_cat_high_mat[:,range(50) + range(190,240)]
X = mut_cat_high_mat.astype(np.float) / mut_cat_high_mat.sum(axis=0)
X = X[:,range(50) + range(190,240)]
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




