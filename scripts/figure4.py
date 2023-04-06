import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.gridspec import GridSpec

os.chdir('/Volumes/Samsung_T5/ProtPredict/scripts')

#Figure 2 should represent the difference in performance between featuresets (AA only vs ILP & AA)

#first, read performance?
data_AA = pd.read_csv('../results/evaluation/AA/AA_20230314.csv', index_col = 0)
data_AA_mean = data_AA.groupby(by='Protein').mean() #probably not needed
data_AA_mean = data_AA_mean.sort_values(by = 'R-squared', ascending = False)
#data_AA_mean = data_AA_mean[~data_AA_mean.index.isin(['P50623', 'P47035'])] #remove for viz purposes, be careful about this tho

#first, read performance?
data_AAILP = pd.read_csv('../results/evaluation/AAILP/AAILP_20230311.csv', index_col = 0)
data_AAILP_mean = data_AAILP.groupby(by='Protein').mean() #probably not needed
data_AAILP_mean = data_AAILP_mean.sort_values(by = 'R-squared', ascending = False)
#data_AAILP_mean = data_AAILP_mean[~data_AAILP_mean.index.isin(['P50623', 'P47035'])] #remove for viz purposes, be careful about this tho


#merge the two evaluation-sets for joint plotting
data = pd.merge(data_AA_mean['R-squared'], data_AAILP_mean['R-squared'], left_index = True, right_index = True, how = 'inner')
dfm = data.melt(ignore_index=False).reset_index()
dfm = dfm.replace('R-squared_x','AA')
dfm = dfm.replace('R-squared_y','AAILP')
dfi = dfm[dfm['variable'] == 'AA']
dfi = dfi.reset_index()[['index','Protein']]
dfm_complete = pd.merge(dfm, dfi, left_on = 'Protein', right_on = 'Protein')
dfm_complete['variable'] = dfm_complete['variable'].astype('category')

#Order the sets so that the correct category is plotted underneath
X1 = dfm_complete[dfm_complete.variable == 'AA'].Protein
Y1 = dfm_complete[dfm_complete.variable == 'AA']['value']
X2 = dfm_complete[dfm_complete.variable == 'AAILP'].Protein
Y2 = dfm_complete[dfm_complete.variable == 'AAILP']['value']


#Figure plotting
fig = plt.figure()
gs = GridSpec(4,4)

ax_joint = fig.add_subplot(gs[0:4,0:3])
ax_marg_y = fig.add_subplot(gs[0:4,3])
ax_joint.scatter(X2,Y2, alpha = 0.5, marker = 'o', s = 7, label='Amino acids')
ax_joint.scatter(X1,Y1, alpha = 0.5, marker = 's', s = 7, label='ILP')

ax_marg_y.hist(Y2,orientation="horizontal", bins = 75, alpha = 0.5, density = True)
ax_marg_y.hist(Y1,orientation="horizontal", bins = 75, alpha = 0.5, density = True)

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

ax_joint.legend(['Amino acids & Relational features','Amino acids'])

ax_joint.set_title('Protein model performance')
plt.show()













