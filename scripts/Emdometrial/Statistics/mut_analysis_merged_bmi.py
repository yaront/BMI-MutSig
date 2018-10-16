#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 17:54:44 2018

@author: tomer
"""

#%%
# =================================================
# # Mutation per gene
# =================================================

import numpy as np
import pandas as pd

#%%

mut_data = pd.read_table('./../../../databases/Endometrial/mutations/Endometrial_mutations.txt', sep = '\t')
bmi_data = pd.read_table('./../../../databases/Endometrial/information/TCGA_bmi_data.txt', sep = '\t')
pat_bmi = bmi_data[bmi_data['bmi'] != '--']
pat_bmi = pat_bmi[(18.5 < pd.to_numeric(pat_bmi['bmi'])) & (pd.to_numeric(pat_bmi['bmi']) < 90)]

patients = list(set(np.unique(['-'.join(x.split('-')[0:3]) for x in mut_data['Tumor_Sample_Barcode']])).intersection(list(pat_bmi['submitter_id'].values)))

pat_bmi = pat_bmi[[(x in patients) for x in pat_bmi['submitter_id'].values]].sort_values(by = ['bmi'])
pat_mut = mut_data[[(x in patients) for x in ['-'.join(x.split('-')[0:3]) for x in mut_data['Tumor_Sample_Barcode']]]]
pat_mut = pat_mut[pat_mut['Variant_Classification'].isin(['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Translation_Start_Site'])]

gene_bmi_mut = pd.DataFrame(index = np.unique(pat_bmi['bmi']))

for g in np.unique(pat_mut['Hugo_Symbol']):
    gene_mut = pat_mut[pat_mut['Hugo_Symbol'] == g]
    gene_pat = ['-'.join(x.split('-')[0:3]) for x in gene_mut['Tumor_Sample_Barcode']]

    mut_count = []
    for b in gene_bmi_mut.index:
        mut_count.append(len(list(set(list(pat_bmi[pd.to_numeric(pat_bmi['bmi']) == int(b)]['submitter_id'])).intersection(gene_pat))))
    gene_bmi_mut[g] = mut_count

gene_bmi_mut = gene_bmi_mut.transpose()

#%%

slope = []
for i,j in gene_bmi_mut.iterrows():
    slope.append(np.polyfit(pd.to_numeric(gene_bmi_mut.columns),j.values,1)[0])

gene_bmi_mut['Slope'] = slope
gene_bmi_mut = gene_bmi_mut.sort_values(by = ['Slope'])


#%%

gene_bmi_mut.to_csv('./../output/endometrial_bmi_gene_mut.txt', header = True, index = True, sep = '\t')






