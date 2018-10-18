#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 18:33:50 2018

@author: tomer
"""

#%%
##======================================##
## Merging the RNAseq List and BMI List ##
##======================================##

import numpy as np
import pandas as pd

#%%

clinical_data = pd.read_table('/Users/tomer/Google Drive/PhD/Thesis/Mutation Signature/databases/Endometrial/information/ucec_tcga_clinical_data.txt', sep = '\t', index_col = 0)
tcga_ann =  pd.read_table('/Users/tomer/Google Drive/PhD/Thesis/Mutation Signature/databases/Endometrial/information/TCGA_annotations.txt', sep = '\t')
rna_seq = pd.read_csv('/Users/tomer/Google Drive/PhD/Thesis/Mutation Signature/databases/Endometrial/information/RNASeq-TCGA-UCEC_FPKM_020217.txt.summarise.csv', index_col = 0)

clinical_data = clinical_data[~clinical_data.index.duplicated(keep='first')]


#%%

bmi = []
race = []
age = []
tumor_type = []
stage = []

for r in rna_seq.columns:
    if bool(tcga_ann[tcga_ann['Sample'] == r].any().sum()):
        pat = str(tcga_ann[tcga_ann['Sample'] == r]['sample_name'].values[0])
        if (pat not in clinical_data.index):
            bmi.append('-')
            race.append('-')
            age.append('-')
            tumor_type.append('-')
            stage.append('-')
            continue
        pat_data = clinical_data.loc[pat,:]
        if (pat_data['BMI'] == '#DIV/0!'):
            bmi.append('-')
        else:
            bmi.append(float(pat_data['BMI']))
        race.append(pat_data['Race Category'])
        if (np.isnan(pat_data['Diagnosis Age'])):
            age.append('-')
        else:
            age.append(int(pat_data['Diagnosis Age']))
        tumor_type.append(pat_data['Cancer Type Detailed'])
        stage.append(pat_data['Neoplasm American Joint Committee on Cancer Clinical Group Stage'])
    else:
        bmi.append('-')
        race.append('-')
        age.append('-')
        tumor_type.append('-')
        stage.append('-')


#%%

rna_seq_data = rna_seq.transpose()
rna_seq_data.insert(loc=0, column='STAGE', value=stage)
rna_seq_data.insert(loc=0, column='TUMOR_TYPE', value=tumor_type)
rna_seq_data.insert(loc=0, column='AGE', value=age)
rna_seq_data.insert(loc=0, column='RACE', value=race)
rna_seq_data.insert(loc=0, column='BMI', value=bmi)
rna_seq_data = rna_seq_data.transpose()


#%%

rna_seq_data.to_csv('./../../../output/endometrial_rna_seq_with_clinical_data.txt', sep = '\t')
