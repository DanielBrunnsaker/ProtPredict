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


#first, read performance?
#data_AAILP = pd.read_csv('intermediateData/scores/GCV2_AAILP_20230119.csv', index_col = 0)
data_AAILP = pd.read_csv('../results/evaluation/AAILP/AAILP_20230311.csv', index_col = 0)
data_AAILP_mean = data_AAILP.groupby(by='Protein').mean() #probably not needed
data_AAILP_mean = data_AAILP_mean.sort_values(by = 'R-squared', ascending = False)
#data_AAILP_mean = data_AAILP_mean[~data_AAILP_mean.index.isin(['P50623', 'P47035'])] #remove for viz purposes, be careful about this tho



feature_shap = pd.read_csv('../results/feature_values/AAILP/Proteins_AAILP_SHAP_means_20230311.tsv', sep = '\t', index_col = 0)
feature_shap = feature_shap.iloc[:,np.where(feature_shap.columns.isin(list(data_AAILP_mean[data_AAILP_mean['R-squared'] > 0].index)) == True)[0]] #only pick those with r2 over 0
feature_shap = feature_shap.rename(index={'alanine': 'Ala', 'arginine' : 'Arg', 'asparagine' : 'Asn', 'aspartate' : 'Asp',
                                          'glutamate' : 'Glu', 'glutamine' : 'Gln', 'glycine' : 'Gly', 'histidine' : 'His', 'isoleucine' : 'Ile',
                                          'leucine' : 'Leu', 'lysine' : 'Lys', 'methionine' : 'Met', 'phenylalanine' : 'Phe','proline' : 'Pro', 'serine' : 'Ser',
                                          'threonine' : 'Thr', 'tryptophan' : 'Trp', 'tyrosine' : 'Tyr', 'valine' : 'Val'})

#Figure 3D, Heatmap
cg = sns.clustermap(data = feature_shap.iloc[:19,:].transpose(), figsize=(12,15), cmap = 'mako', z_score = 0, 
                    dendrogram_ratio = 0.05, row_cluster = True, xticklabels = True, yticklabels = False, 
                    cbar_pos = 'left', linewidths=0.000000001, linecolor = 'black', tree_kws=dict(linewidths=1.5))
cg.ax_col_dendrogram.set_visible(True)
cg.ax_row_dendrogram.set_visible(False)
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
sns.set_context('paper', font_scale=1.5)
#cg.savefig('/Volumes/Samsung_T5/MetPredict/manuscript/prpred/test2png.png', dpi=300)
plt.show()













