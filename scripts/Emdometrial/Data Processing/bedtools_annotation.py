#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 16:13:37 2018

@author: tomer
"""

#%%
# =============================================================================
# # Creating an annotation file for bedtools
# =============================================================================

print("Starting...")

import numpy as np
import pandas as pd

#%%

chr_list = list(range(1,23)) + ['MT','X']
chr_list = [str(x) for x in chr_list]

mut_full = pd.read_table('./../../../databases/Endometrial/mutations/Endometrial_mutations.txt', sep = '\t')
mut_filtered = mut_full[mut_full['Variant_Type'] == 'SNP']
mut_filtered = mut_filtered[mut_filtered['Chromosome'].isin(chr_list)]
mut_filtered['Chromosome'] = mut_filtered['Chromosome'].replace(to_replace='MT', value='M')

bedtools_ann = mut_filtered[['Chromosome','Start_Position','End_Position']]
bedtools_ann['Chromosome'] = 'chr' + bedtools_ann['Chromosome'].astype(str)
bedtools_ann['Start_Position'] = bedtools_ann['Start_Position'] - 2
bedtools_ann['End_Position'] = bedtools_ann['End_Position'] + 1
bedtools_ann['name'] = [(str(x['Tumor_Sample_Barcode']) + '_' + str(x['Hugo_Symbol'])) for i,x in mut_filtered.iterrows()]

bedtools_ann.columns = ['chrom','chromStart','chromEnd','name']

bedtools_ann.to_csv('../../output/endo_mut_bedtools_ann.bed',sep = '\t',index = False, header = False)

print("Done.")


