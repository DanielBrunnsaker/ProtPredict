#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:01:20 2022

@author: danbru
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir('/Volumes/Samsung_T5/ProtPredict/scripts')

#first, read performance
data_AAILP = pd.read_csv('../results/evaluation/AAILP/AAILP_20230311.csv', index_col = 0)
data_AAILP_mean = data_AAILP.groupby(by='Protein').mean() 
data_AAILP_mean = data_AAILP_mean.sort_values(by = 'R-squared', ascending = False)


#read in feature importance information
feature_shap = pd.read_csv('../results/feature_values/AAILP/Proteins_AAILP_SHAP_means_20230311.tsv', sep = '\t', index_col = 0)
feature_shap = feature_shap.iloc[:,np.where(feature_shap.columns.isin(list(data_AAILP_mean[data_AAILP_mean['R-squared'] > 0].index)) == True)[0]] #only pick those with r2 over 0
feature_shap = feature_shap.rename(index={'alanine': 'Ala', 'arginine' : 'Arg', 'asparagine' : 'Asn', 'aspartate' : 'Asp',
                                          'glutamate' : 'Glu', 'glutamine' : 'Gln', 'glycine' : 'Gly', 'histidine' : 'His', 'isoleucine' : 'Ile',
                                          'leucine' : 'Leu', 'lysine' : 'Lys', 'methionine' : 'Met', 'phenylalanine' : 'Phe','proline' : 'Pro', 'serine' : 'Ser',
                                          'threonine' : 'Thr', 'tryptophan' : 'Trp', 'tyrosine' : 'Tyr', 'valine' : 'Val'})

feature_shap = feature_shap.iloc[:19,:].transpose()
feature_columns = list(feature_shap.columns)

#Figure 5, Heatmap
cg = sns.clustermap(data = feature_shap, figsize=(9,10), cmap = 'mako', z_score = 0, dendrogram_ratio = 0.05, row_cluster = True, 
                    col_cluster = True, yticklabels = False, linecolor = 'black')
cg.ax_col_dendrogram.set_visible(True)
cg.ax_row_dendrogram.set_visible(False)
cg.cax.set_visible(False)
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
sns.set_context('paper', font_scale=2)
#cg.savefig('../Figure5.png', dpi=400)
#plt.tight_layout()
plt.show()







