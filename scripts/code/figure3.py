#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:01:20 2022

@author: danbru
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

os.chdir('') #change to scripts folder

#load evaluation data
data = pd.read_csv('../results/evaluation/ILP/ILP_20230309.csv', index_col = 0)
data_mean = data.groupby(by='Protein').mean() #average the folds

#Figure 3 - overall feature dependencies for predictions given in 1!
#Only include models with an r2 of above 0!
list_to_include = list(data_mean[data_mean['R-squared'] > 0].index)

#remove features with r2 under 0
feature_shap = pd.read_csv('../results/feature_values/ILP/SHAP/ILP_SHAP_20230310.tsv', sep = '\t', index_col = 0)
feature_shap = feature_shap[feature_shap.columns.intersection(list_to_include)]

#Normalize columns by dividing with total sum and then average
for col in feature_shap.columns:
    feature_shap[col] = feature_shap[col]/feature_shap[col].sum(axis = 0)
fmeans = feature_shap.mean(axis = 1)

#do the same with XGB-values?
features_XGB = pd.DataFrame(index = feature_shap.index)
path = '../results/feature_values/ILP/XGB/'

for files in os.listdir(path):
    
    if files.split('_')[0] in list_to_include:
        
        if '20230310' in files: #just to pick the most recent ones
            xgb_features = pd.DataFrame(pd.read_csv(path+files, sep = '\t', index_col = 0)['gain'])
            xgb_features.columns = [files.split('_')[0]]
            #weigh them?
            features_XGB = pd.merge(features_XGB, xgb_features, left_index = True, right_index = True, how = 'left')
            features_XGB = features_XGB.fillna(0)
            features_XGB[files.split('_')[0]] = features_XGB[files.split('_')[0]]/features_XGB[files.split('_')[0]].sum(axis = 0) 
        
xgbmeans = features_XGB.mean(axis = 1)


#do barplot of features?
#pick top features by mean? 
nr_to_include = 15

fig, axes = plt.subplots(2, 1, figsize=(5, 10))

#shap
sns.barplot(ax = axes[0], data = feature_shap.transpose()[fmeans.nlargest(nr_to_include).index], color = '#3B97B6')#619CFF'
#axes[0].set_title('mean(|SHAP|) (across all proteins)', fontsize = 13)
axes[0].tick_params(labelrotation=45, labelsize = 12)
#axes[0].set_xlabel('Features', fontsize = 13)
axes[0].set_ylabel('Fraction of contribution to model output', fontsize = 13)

#xgb
sns.barplot(ax = axes[1], data = features_XGB.transpose()[xgbmeans.nlargest(nr_to_include).index], color = '#3B97B6')
axes[1].tick_params(labelrotation=45, labelsize = 12)
axes[1].set_ylim([0.00,0.008])
axes[1].set_xlabel('Features', fontsize = 13)
axes[1].set_ylabel('Fraction of contribution to gain', fontsize = 13)
plt.tight_layout()
