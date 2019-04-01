#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 20:42:43 2018

@author: tomer
"""

#%%
# =================================================
# # Mutation per gene
# =================================================

import numpy as np
import pandas as pd

#%%

#tumor = sys.argv[1]
#tumor = tumor.split('/')[-1].split('.')[0]
#print tumor

tumor = 'UCEC'

#%% Reading data

print("Starting: " + tumor)

mut_data = pd.read_table('./../../../databases/Endometrial/TCGA_MAFs/' + tumor + '.maf', sep = '\t')
bmi_data = pd.read_table('./../../../databases/Endometrial/information/TCGA_bmi_data.txt', sep = '\t')
pat_bmi = bmi_data[bmi_data['bmi'] != '--']
pat_bmi = pat_bmi[(18.5 < pd.to_numeric(pat_bmi['bmi'])) & (pd.to_numeric(pat_bmi['bmi']) < 90)]

patients = list(set(np.unique(['-'.join(x.split('-')[0:3]) for x in mut_data['Tumor_Sample_Barcode']])).intersection(list(pat_bmi['submitter_id'].values)))

pat_bmi = pat_bmi[[(x in patients) for x in pat_bmi['submitter_id'].values]].sort_values(by = ['bmi'])
pat_mut = mut_data[[(x in patients) for x in ['-'.join(x.split('-')[0:3]) for x in mut_data['Tumor_Sample_Barcode']]]]
pat_mut = pat_mut[pat_mut['Variant_Classification'].isin(['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Translation_Start_Site'])]

#%% Creating table of mutations per BMI and mutation burden per patient

gene_bmi_mut = pd.DataFrame(0, columns = ['BMI','Total_Mutations'] + list(np.unique(pat_mut['Hugo_Symbol'])), index = np.sort(pat_bmi[['submitter_id','bmi']])[:,1])
gene_bmi_mut['BMI'] = np.sort(pat_bmi[['submitter_id','bmi']])[:,0]

pat_name_mut = ['-'.join(x.split('-')[0:3]) for x in pat_mut['Tumor_Sample_Barcode']]

for pat in gene_bmi_mut.index:
    gene_bmi_mut.loc[pat,'Total_Mutations'] = pat_name_mut.count(pat)

gene_bmi_mut = gene_bmi_mut[gene_bmi_mut['Total_Mutations'] < 3000]


#%% Assigning mutations per gene per patient

print("Calculating mutations for " + tumor)

for g in np.unique(pat_mut['Hugo_Symbol']):
    gene_mut = pat_mut[pat_mut['Hugo_Symbol'] == g]
    gene_pat = ['-'.join(x.split('-')[0:3]) for x in gene_mut['Tumor_Sample_Barcode']]

    for p in np.unique(gene_pat):
        gene_bmi_mut.loc[p,g] = gene_pat.count(p)

gene_bmi_mut = gene_bmi_mut.transpose()

norm_gene_bmi_mut = []


#%% Finding the slope

print("Calculating slope for " + tumor)

inds = {bmi: ind for ind,bmi in enumerate(set(pd.to_numeric(gene_bmi_mut.loc['BMI',:])))}
bmi_ind = [inds[bmi] for bmi in pd.to_numeric(gene_bmi_mut.loc['BMI',:])]

slope = []
for i,j in gene_bmi_mut.iloc[2:,:].iterrows():
    norm_mut = pd.to_numeric(j) / pd.to_numeric(gene_bmi_mut.loc['Total_Mutations'])
    norm_gene_bmi_mut.append(norm_mut)
    weight_mut = np.bincount(np.array(bmi_ind),weights=list(map(float,norm_mut.values))) / np.bincount(np.array(bmi_ind))
    slope.append(np.polyfit(list(range(len(weight_mut))), weight_mut,1)[0])

norm_gene_bmi_mut = pd.DataFrame(norm_gene_bmi_mut)
norm_gene_bmi_mut = pd.concat([gene_bmi_mut.loc[['BMI','Total_Mutations'],:],norm_gene_bmi_mut])
norm_gene_bmi_mut.index = gene_bmi_mut.index

gene_bmi_mut['Slope'] = [-np.inf,-np.inf] + slope
gene_bmi_mut = gene_bmi_mut.sort_values(by = ['Slope'])
gene_bmi_mut.loc[['BMI','Total_Mutations'],'Slope'] = '-'

norm_gene_bmi_mut['Slope'] = [-np.inf,-np.inf] + slope
norm_gene_bmi_mut = norm_gene_bmi_mut.sort_values(by = ['Slope'])
norm_gene_bmi_mut.loc[['BMI','Total_Mutations'],'Slope'] = '-'


#%% Writing the data

print("Writing " + tumor)

gene_bmi_mut.to_csv('./../output/' + tumor + '_bmi_gene_mut.txt', header = True, index = True, sep = '\t')
norm_gene_bmi_mut.to_csv('./../output/' + tumor + '_bmi_gene_mut_norm.txt', header = True, index = True, sep = '\t')

writer = pd.ExcelWriter('./../output/' + tumor + '_bmi_gene_mut_slope.xlsx', engine='xlsxwriter')
gene_bmi_mut.to_excel(writer, sheet_name = tumor + '_binary')
norm_gene_bmi_mut.to_excel(writer, sheet_name = tumor + '_norm')
writer.save()

print("Done: " + tumor)





