#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:01:20 2022

@author: danbru
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#figure 2 should be an overview on predictive performance given
#the relational descriptors. With some global expalantion of what features seem
#to be important. Use this to point towards amino acids as a predictor of 
#proteins.

# Performance metrics for all of the protein models!

#potentially change to this one? e.g with 2048 features and gcv3, plus duplicates removed
data = pd.read_csv('../intermediateData/scores/ILP_20230309.csv', index_col = 0)
data_mean = data.groupby(by='Protein').mean() #probably not needed
#sort by average?
data_mean = data_mean.sort_values(by = 'R-squared', ascending = False)
#data_mean = data_mean[data_mean.index != 'P08536'] #remove for viz purposes, be careful about this tho
#data_mean = data_mean[data_mean.index != 'P38858'] #remove for viz purposes, be careful about this tho
#data_mean = data_mean[data_mean.index != 'Q08438'] #remove for viz purposes, be careful about this tho

fig = plt.figure()
gs = GridSpec(4,4)

ax_joint = fig.add_subplot(gs[0:4,0:3])
ax_marg_y = fig.add_subplot(gs[0:4,3])
ax_joint.scatter(x = data_mean.index,y = data_mean['R-squared'], alpha = 0.5, marker = 's', s = 7, label='ILP', color = '#3B97B6')

ax_marg_y.hist(data_mean['R-squared'],orientation="horizontal", bins = 75, alpha = 0.5, density = True, color = '#3B97B6')

# Turn off tick labels on marginals
plt.setp(ax_marg_y.get_yticklabels(), visible=False)

# Set labels on joint
ax_joint.set_xlabel('Proteins')
ax_joint.set_ylabel('R-squared')

ax_marg_y.spines['top'].set_visible(False)
ax_marg_y.spines['right'].set_visible(False)
ax_marg_y.spines['bottom'].set_visible(False)

ax_joint.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

ax_marg_y.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

ax_joint.legend(['Relational features'])

ax_joint.set_title('Protein model performance')
plt.show()