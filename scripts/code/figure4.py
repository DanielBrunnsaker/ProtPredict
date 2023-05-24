import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.gridspec import GridSpec

os.chdir('') #change to scripts folder

#Figure 2 should represent the difference in performance between featuresets (AA only vs ILP & AA)

#first, read performance?
data_AA = pd.read_csv('../results/evaluation/AA/AA_20230314.csv', index_col = 0)
data_AA_mean = data_AA.groupby(by='Protein').mean() #probably not needed
data_AA_mean = data_AA_mean.sort_values(by = 'R-squared', ascending = False)
data_AA_SD = pd.DataFrame(data_AA.groupby(by='Protein').std()['R-squared'])
data_AA_SD.columns = ['SD']

#first, read performance?
data_AAILP = pd.read_csv('../results/evaluation/AAILP/AAILP_20230311.csv', index_col = 0)
data_AAILP_mean = data_AAILP.groupby(by='Protein').mean() #probably not needed
data_AAILP_mean = data_AAILP_mean.sort_values(by = 'R-squared', ascending = False)
data_AAILP_SD = pd.DataFrame(data_AAILP.groupby(by='Protein').std()['R-squared'])
data_AAILP_SD.columns = ['SD']

#add in SD
data_AA_mean = pd.merge(data_AA_mean, data_AA_SD, left_index = True, right_index = True, how = 'left')
data_AAILP_mean = pd.merge(data_AAILP_mean, data_AAILP_SD, left_index = True, right_index = True, how = 'left')

#merge the two evaluation-sets for joint plotting
data = pd.merge(data_AA_mean, data_AAILP_mean, left_index = True, right_index = True, how = 'inner')

#Figure plotting
fig = plt.figure()
gs = GridSpec(4,5)

ax_joint = fig.add_subplot(gs[0:4,0:4])
ax_marg_y = fig.add_subplot(gs[0:4,4])
ax_joint.scatter(data.index,data['R-squared_y'], alpha = 0.7, marker = 's', s = 7, label='ILP', color = '#3B97B6')
ax_joint.fill_between(data.index, data['R-squared_y'] - data['SD_y'], data['R-squared_y'] + data['SD_y'], color='#3B97B6', alpha=0.4, label = 'ILP_SD')
ax_joint.scatter(data.index,data['R-squared_x'], alpha = 1, marker = 'o',label='Amino acids', color = '#F8766D', s = 1)

ax_marg_y.hist(data['R-squared_y'],orientation="horizontal", bins = 100, alpha = 0.6, density = True, color = '#3B97B6')
ax_marg_y.hist(data['R-squared_x'],orientation="horizontal", bins = 100, alpha = 0.6, density = True, color = '#F8766D')

ax_marg_y.set_ylim([-0.7,0.7])
ax_joint.set_ylim([-0.7,0.7])

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

ax_joint.legend(['Mean $R^2$ (Amino acids & Relational features)','SD (Amino acids & Relational features)', 'Mean $R^2$ (Amino acids only)'], 
                loc = 'lower left', fontsize = 10.5)

ax_joint.set_title('Predictive model performance', fontsize = 14)
plt.show()






