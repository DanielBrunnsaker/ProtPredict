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



















dfm = data.melt(ignore_index=False).reset_index()
dfm = dfm.replace('R-squared_x','AA')
dfm = dfm.replace('R-squared_y','AAILP')

dfm = dfm.replace('SD_x','AA_SD')
dfm = dfm.replace('SD_y','AAILP_SD')

dfi = dfm[dfm['variable'] == 'AA']
dfi = dfi.reset_index()[['index','Protein']]
dfm_complete = pd.merge(dfm, dfi, left_on = 'Protein', right_on = 'Protein')
dfm_complete['variable'] = dfm_complete['variable'].astype('category')








#Order the sets so that the correct category is plotted underneath
X1 = dfm_complete[dfm_complete.variable == 'AA'].Protein
Y1 = dfm_complete[dfm_complete.variable == 'AA']['value']
X2 = dfm_complete[dfm_complete.variable == 'AAILP'].Protein
Y2 = dfm_complete[dfm_complete.variable == 'AAILP']['value']

y1SD = dfm_complete[dfm_complete.variable == 'AA_SD']['value']
y2SD = dfm_complete[dfm_complete.variable == 'AAILP_SD']['value']



#Figure plotting
fig = plt.figure()
gs = GridSpec(4,4)

ax_joint = fig.add_subplot(gs[0:4,0:3])
ax_marg_y = fig.add_subplot(gs[0:4,3])
ax_joint.scatter(X2,Y2, alpha = 0.5, marker = 'o', s = 7, label='ILP', color = '#3B97B6')
ax_joint.scatter(X1,Y1, alpha = 0.5, marker = 's', s = 7, label='Amino acids', color = '#F8766D')

ax_joint.fill_between(X2, Y2 - Y2SD, Y2 + Y2SD, color='#3B97B6', alpha=0.4)
ax_joint.fill_between(X1, Y1 - Y1SD, Y1 + Y1SD, color='#3B97B6', alpha=0.4)


ax_marg_y.hist(Y2,orientation="horizontal", bins = 100, alpha = 0.5, density = True)
ax_marg_y.hist(Y1,orientation="horizontal", bins = 100, alpha = 0.5, density = True)

# Turn off tick labels on marginals
plt.setp(ax_marg_y.get_yticklabels(), visible=False)

# Set labels on joint
ax_joint.set_xlabel('Protein', fontsize = 13)
ax_joint.set_ylabel('R-squared', fontsize = 13)

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

ax_joint.legend(['Amino acids & Relational features','Amino acids only'])

ax_joint.set_title('Predictive model performance')
plt.show()




















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
ax_joint.scatter(X2,Y2, alpha = 0.5, marker = 'o', s = 7, label='ILP', color = '#3B97B6')
ax_joint.scatter(X1,Y1, alpha = 0.5, marker = 's', s = 7, label='Amino acids', color = '#F8766D')

ax_joint.fill_between(data_mean.index, data_mean['R-squared'] - data_mean['SD'], data_mean['R-squared'] + data_mean['SD'], color='#3B97B6', alpha=0.4)
ax_joint.fill_between(data_mean.index, data_mean['R-squared'] - data_mean['SD'], data_mean['R-squared'] + data_mean['SD'], color='#3B97B6', alpha=0.4)


ax_marg_y.hist(Y2,orientation="horizontal", bins = 100, alpha = 0.5, density = True)
ax_marg_y.hist(Y1,orientation="horizontal", bins = 100, alpha = 0.5, density = True)

# Turn off tick labels on marginals
plt.setp(ax_marg_y.get_yticklabels(), visible=False)

# Set labels on joint
ax_joint.set_xlabel('Protein', fontsize = 13)
ax_joint.set_ylabel('R-squared', fontsize = 13)

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

ax_joint.legend(['Amino acids & Relational features','Amino acids only'])

ax_joint.set_title('Predictive model performance')
plt.show()













