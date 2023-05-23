#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:01:20 2022

@author: danbru
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os

#figure 2 should be an overview on predictive performance given the relational descriptors. 

os.chdir('') #change to scripts folder

# Performance metrics for all of the protein models!
data = pd.read_csv('../results/evaluation/ILP/ILP_20230309.csv', index_col = 0)
data_mean = data.groupby(by='Protein').mean() #

#calculate STD?
data_sd = pd.DataFrame(data.groupby(by='Protein').std()['R-squared'])
data_sd.columns = ['SD']


#sort by average?
data_mean = data_mean.sort_values(by = 'R-squared', ascending = False)

#merge with sd
data_mean = pd.merge(data_mean, data_sd, left_index = True, right_index = True, how = 'left')


fig = plt.figure()
gs = GridSpec(4,5)

ax_joint = fig.add_subplot(gs[0:4,0:4])
ax_marg_y = fig.add_subplot(gs[0:4,4])


ax_joint.scatter(x = data_mean.index,y = data_mean['R-squared'], alpha = 0.9, marker = 's', s = 7, label='ILP', color = '#3B97B6')

ax_joint.fill_between(data_mean.index, data_mean['R-squared'] - data_mean['SD'], data_mean['R-squared'] + data_mean['SD'], color='#3B97B6', alpha=0.4,
                      label = 'Standard deviation')
ax_marg_y.hist(data_mean['R-squared'],orientation="horizontal", bins = 100, alpha = 0.9, density = True, color = '#3B97B6')
ax_marg_y.set_ylim([-0.6,0.6])
ax_joint.set_ylim([-0.6,0.6])


# Turn off tick labels on marginals
plt.setp(ax_marg_y.get_yticklabels(), visible=False)

# Set labels on joint
ax_joint.set_xlabel('Protein', fontsize = 14)
ax_joint.set_ylabel('$R^2$', fontsize = 14)

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

ax_joint.legend(['Mean $R^2$ (Relational features)', 'SD (Relational features)'], loc = 'lower left', fontsize = 10.5)

ax_joint.set_title('Predictive model performance', fontsize = 14)
plt.show()