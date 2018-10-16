#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 22:01:00 2018

@author: tomer
"""

#%%
##========================================##
## Creating the Mutational Catalog Matrix ##
##========================================##

import numpy as np
import pandas as pd

#%% Creating the mutational data

print "Reading data..."
mut_data = pd.read_csv('./../../../databases/Endometrial/mutations/Endometrial_SNP_mutations_context.txt', sep = '\t')
bmi_data = pd.read_table('./../../../databases/Endometrial/information/TCGA_bmi_data.txt', sep = '\t')

pat_bmi = bmi_data[bmi_data['bmi'] != '--']
pat_bmi = pat_bmi[(18.5 < pd.to_numeric(pat_bmi['bmi'])) & (pd.to_numeric(pat_bmi['bmi']) < 90)]

patients = list(set(np.unique(['-'.join(x.split('-')[0:3]) for x in mut_data['Tumor_Sample_Barcode']])).intersection(list(pat_bmi['submitter_id'].values)))
patients = np.sort(patients)

pat_bmi = pat_bmi[[(x in patients) for x in pat_bmi['submitter_id'].values]].sort_values(by = ['bmi'])
pat_mut = mut_data[[(x in patients) for x in ['-'.join(x.split('-')[0:3]) for x in mut_data['Tumor_Sample_Barcode']]]]

print "Assigning BMI and context..."
sam_bmi = []
base_sub_type = []
for i,s in pat_mut.iterrows():
    sam_barcode = '-'.join(s['Tumor_Sample_Barcode'].split('-')[:3])
    sam_bmi.append(int(pat_bmi[pat_bmi['submitter_id'] == sam_barcode]['bmi'].values.astype(int)))
    if (((s['Reference_Allele'] == 'C') & (s['Tumor_Seq_Allele2'] == 'A')) | ((s['Reference_Allele'] == 'G') & (s['Tumor_Seq_Allele2'] == 'T'))):
        base_sub_type.append('C:G>A:T')
    elif (((s['Reference_Allele'] == 'C') & (s['Tumor_Seq_Allele2'] == 'G')) | ((s['Reference_Allele'] == 'G') & (s['Tumor_Seq_Allele2'] == 'C'))):
        base_sub_type.append('C:G>G:C')
    elif (((s['Reference_Allele'] == 'C') & (s['Tumor_Seq_Allele2'] == 'T')) | ((s['Reference_Allele'] == 'G') & (s['Tumor_Seq_Allele2'] == 'A'))):
        base_sub_type.append('C:G>T:A')
    elif (((s['Reference_Allele'] == 'T') & (s['Tumor_Seq_Allele2'] == 'A')) | ((s['Reference_Allele'] == 'A') & (s['Tumor_Seq_Allele2'] == 'T'))):
        base_sub_type.append('T:A>A:T')
    elif (((s['Reference_Allele'] == 'T') & (s['Tumor_Seq_Allele2'] == 'C')) | ((s['Reference_Allele'] == 'A') & (s['Tumor_Seq_Allele2'] == 'G'))):
        base_sub_type.append('T:A>C:G')
    elif (((s['Reference_Allele'] == 'T') & (s['Tumor_Seq_Allele2'] == 'G')) | ((s['Reference_Allele'] == 'A') & (s['Tumor_Seq_Allele2'] == 'C'))):
        base_sub_type.append('T:A>G:C')
pat_mut = pat_mut.assign(BMI = sam_bmi)
pat_mut = pat_mut.assign(Mut_Type = base_sub_type)

full_context_type = [''.join(list(x['Context'])[0] + 'p' + x['Mut_Type'] + 'p' + list(x['Context'])[2]) for i,x in pat_mut.iterrows()]
pat_mut = pat_mut.assign(Full_Context_Type = full_context_type)
full_context = [''.join(list(x['Context'])[0] + 'p' + x['Reference_Allele'] + '>' + x['Tumor_Seq_Allele2'] + 'p' + list(x['Context'])[2]) for i,x in pat_mut.iterrows()]
pat_mut = pat_mut.assign(Full_Context = full_context)

mut_type = np.unique(pat_mut['Mut_Type'])
mut_con = np.unique(pat_mut['Full_Context'])
mut_con_type = np.unique(pat_mut['Full_Context_Type'])
pat_mut = pat_mut.sort_values(by=['BMI','Tumor_Sample_Barcode','Reference_Allele','Tumor_Seq_Allele2','Context'])


#%% GEnerating the mutational catalog data

print "Calculating catalog matrix..."
mut_cat = pd.DataFrame(0, index=mut_con_type, columns=pd.unique(pat_mut['Tumor_Sample_Barcode']))
for i,s in pat_mut.iterrows():
    mut_cat.loc[s['Full_Context_Type'],s['Tumor_Sample_Barcode']]+=1

#%% Saving data

print "Writing data..."
mut_cat.to_csv('./../../output/Endometrial_mutational_catalog_bmi_sort.txt', sep = '\t')












