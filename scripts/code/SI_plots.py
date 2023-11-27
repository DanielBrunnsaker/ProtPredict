# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
os.chdir('')



#### ======================== ILP only methods ======================== ####



# first, read performance for ILP
data_ILP = pd.read_csv('../results/evaluation/ILP/ILP_20230309.csv', index_col = 0)
data_ILP_mean = data_ILP.groupby(by='Protein').mean() #probably not needed
data_ILP_mean = data_ILP_mean.sort_values(by = 'R-squared', ascending = False)

# Dummy regressor
data_dummy = pd.read_csv('../results/evaluation/Comparison/dummymean_results.csv', index_col = 0)
data_dummy_mean = data_dummy.groupby(by='Protein').mean() #probably not needed
data_dummy_mean = data_dummy_mean.reindex(data_ILP_mean.index)

# ElasticNet Regressor
data_ilp_elasticnet = pd.read_csv('/Users/danbru/Desktop/Projects/comparison/elasticnet_ilp_results.csv', index_col = 0)
data_ilp_elasticnet_mean = data_ilp_elasticnet.groupby(by='Protein').mean() #probably not needed
data_ilp_elasticnet_mean = data_ilp_elasticnet_mean.reindex(data_ILP_mean.index)

# SVR Regressor
data_ilp_svr = pd.read_csv('/Users/danbru/Desktop/Projects/comparison/svr_ilp_results.csv', index_col = 0)
data_ilp_svr_mean = data_ilp_svr.groupby(by='Protein').mean() #probably not needed
data_ilp_svr_mean = data_ilp_svr_mean.reindex(data_ILP_mean.index)

# Merge these?
plotdata = pd.merge(data_ILP_mean, data_dummy_mean, left_on = 'Protein', right_on = 'Protein')
plotdata = pd.merge(plotdata, data_ilp_elasticnet_mean, left_on = 'Protein', right_on = 'Protein')
plotdata = plotdata[['R-squared_x','R-squared_y','R-squared']]
plotdata.columns = ['XGBoost','Dummy','ElasticNet']
plotdata = pd.merge(plotdata, data_ilp_svr_mean, left_on = 'Protein', right_on = 'Protein')[['XGBoost','ElasticNet','R-squared','Dummy']]
plotdata.columns = ['XGBoost','ElasticNet','Support vector regression','Naive (mean)']


# Create a figure with subplots for scatter plots and histogram
df = plotdata
# Create a single figure with two subplots
plt.figure(figsize=(10, 4))
plt.title('Predictive model performance comparison (relational features)', fontsize=12)
# Create scatter plots for each column
#scatter_ax.scatter(df.index, df['Naive (mean)'], label='Naive (mean)', alpha = 0.4, s = 4, color = 'tab:red')
plt.scatter(df.index, df['Support vector regression'], label='SVR', alpha = 0.5, s = 8, color = 'tab:green', marker = 'h')
plt.scatter(df.index, df['ElasticNet'], label='ElasticNet', alpha = 0.5, s = 10, color = 'tab:orange', marker = '^')
plt.scatter(df.index, df['XGBoost'], label='XGBoost', alpha = 0.5, s = 2, color = 'tab:blue', marker = 's')
#handles, labels = plt.get_legend_handles_labels()
plt.ylim([-0.6,0.6])
plt.xlabel('Protein index')
plt.ylabel('$R^2$')
plt.legend()

#plt.legend(handles[::-1], labels[::-1])

# Hide x-ticks in the scatterplot
plt.xticks([])

# Adjust layout for better readability
plt.tight_layout()

# Save the plots
plt.savefig("/Users/danbru/Desktop/Projects/SI_figures/benchmark_rel.svg", format="svg", bbox_inches="tight")

# Show the plots
plt.show()

# Save this in vectorized format?





#### ======================== ILP´+AA methods ======================== ####


# first, read performance for ILP
data_AAILP = pd.read_csv('../results/evaluation/AAILP/AAILP_20230311.csv', index_col = 0)
data_AAILP_mean = data_AAILP.groupby(by='Protein').mean() #probably not needed
data_AAILP_mean = data_AAILP_mean.sort_values(by = 'R-squared', ascending = False)

# Dummy regressor
data_dummy = pd.read_csv('../results/evaluation/Comparison/dummymean_results.csv', index_col = 0)
data_dummy_mean = data_dummy.groupby(by='Protein').mean() #probably not needed
data_dummy_mean = data_dummy_mean.reindex(data_AAILP_mean.index)

# ElasticNet Regressor
data_AAILP_elasticnet = pd.read_csv('/Users/danbru/Desktop/Projects/comparison/elasticnet_ilpAA_results.csv', index_col = 0)
data_AAILP_elasticnet_mean = data_AAILP_elasticnet.groupby(by='Protein').mean() #probably not needed
data_AAILP_elasticnet_mean = data_AAILP_elasticnet_mean.reindex(data_AAILP_mean.index)

# SVR Regressor
data_AAILP_svr = pd.read_csv('/Users/danbru/Desktop/Projects/comparison/svr_ilpAA_results.csv', index_col = 0)
data_AAILP_svr_mean = data_AAILP_svr.groupby(by='Protein').mean() #probably not needed
data_AAILP_svr_mean = data_AAILP_svr_mean.reindex(data_AAILP_mean.index)


# Merge these?
plotdata = pd.merge(data_AAILP_mean, data_dummy_mean, left_on = 'Protein', right_on = 'Protein')
plotdata = pd.merge(plotdata, data_AAILP_elasticnet_mean, left_on = 'Protein', right_on = 'Protein')
plotdata = plotdata[['R-squared_x','R-squared_y','R-squared']]
plotdata.columns = ['XGB','Dummy','ElasticNet']
plotdata = pd.merge(plotdata, data_AAILP_svr_mean, left_on = 'Protein', right_on = 'Protein')[['XGB','Dummy','ElasticNet','R-squared']]
plotdata.columns = ['XGBoost','Naive (mean)','ElasticNet','Support vector regression']



# Create a figure with subplots for scatter plots and histogram
df = plotdata
# Create a single figure with two subplots
plt.figure(figsize=(10, 4))
plt.title('Predictive model performance comparison (AA & relational features)', fontsize=12)
# Create scatter plots for each column
    
#scatter_ax.scatter(df.index, df['Naive (mean)'], label='Naive (mean)', alpha = 0.4, s = 4, color = 'tab:red')
plt.scatter(df.index, df['Support vector regression'], label='SVR', alpha = 0.5, s = 8, color = 'tab:green', marker = 'h')
plt.scatter(df.index, df['ElasticNet'], label='ElasticNet', alpha = 0.5, s = 10, color = 'tab:orange', marker = '^')
plt.scatter(df.index, df['XGBoost'], label='XGBoost', alpha = 0.5, s = 2, color = 'tab:blue', marker = 's')
plt.ylim([-0.6,0.6])
plt.xlabel('Protein index')
plt.ylabel('$R^2$')
plt.legend()


# Hide x-ticks in the scatterplot
plt.xticks([])


# Adjust layout for better readability
plt.tight_layout()


# Save figure
plt.savefig("/Users/danbru/Desktop/Projects/SI_figures/benchmark_relAA.svg", format="svg", bbox_inches="tight")

# Show the plots
plt.show()


#### ======================== ILP vs instantiations ======================== ####

# first, read performance for ILP
data_ILP = pd.read_csv('../results/evaluation/ILP/ILP_20230309.csv', index_col = 0)
data_ILP_mean = data_ILP.groupby(by='Protein').mean() #probably not needed
data_ILP_mean = data_ILP_mean.sort_values(by = 'R-squared', ascending = False)

# Dummy regressor
data_instant = pd.read_csv('/Users/danbru/Desktop/Projects/comparison/ElasticNet_instantiation_revised_results.csv', index_col = 0)
data_instant_mean = data_instant.groupby(by='Protein').mean() #probably not needed
data_instant_mean = data_instant_mean.sort_values(by = 'R-squared', ascending = False)


# Dummy regressor
data_elasticnet = pd.read_csv('/Users/danbru/Desktop/Projects/comparison/elasticnet_ilp_results.csv', index_col = 0)
data_elasticnet_mean = data_elasticnet.groupby(by='Protein').mean() #probably not needed
data_elasticnet_mean = data_elasticnet_mean.sort_values(by = 'R-squared', ascending = False)

data_instant_mean = data_instant_mean.reindex(data_ILP_mean.index)
data_elasticnet_mean = data_elasticnet_mean.reindex(data_ILP_mean.index)


plotdata = pd.merge(data_ILP_mean, data_instant_mean, left_on = 'Protein', right_on = 'Protein', how = 'left')[['R-squared_x', 'R-squared_y']]
plotdata = pd.merge(plotdata, data_elasticnet_mean, left_on = 'Protein', right_on = 'Protein', how = 'left')[['R-squared_x', 'R-squared_y','R-squared']]

plotdata.columns = ['Relational features (XGBoost)', 'Instantiated features (ElasticNet)', 'Relational features (ElasticNet)']


# Create a figure with subplots for scatter plots and histogram
df = plotdata
# Create a single figure with two subplots
plt.figure(figsize=(10, 4))
plt.title('Representation comparison', fontsize=12)

plt.scatter(df.index, df['Instantiated features (ElasticNet)'], label='Instantiated features (ElasticNet)', alpha = 0.4, s = 10, color = 'tab:orange')
plt.scatter(df.index, df['Relational features (ElasticNet)'], label='Relational features (ElasticNet)', alpha = 0.4, s = 10, color = 'tab:green', marker = '^')
plt.scatter(df.index, df['Relational features (XGBoost)'], label='Relational features (XGBoost)', alpha = 0.4, s = 4, color = 'tab:blue', marker = 's')



plt.ylim([-0.6,0.6])
plt.xlabel('Protein index')
plt.ylabel('$R^2$')
plt.legend()
# Hide x-ticks in the scatterplot
plt.xticks([])


plt.tight_layout()

# Save figure
#plt.savefig("/Users/danbru/Desktop/Projects/SI_figures/baseline_comp.svg", format="svg", bbox_inches="tight")

# Show the plots
plt.show()





















import shap


shap.summary_plot(his4, max_display = 100)




# P00815
plt.title('His4 (mean)')
raw_shap_values, raw_data_values, df_feature_importance = extract_top_features('/Users/danbru/Desktop/Projects/ProtPredict/results/feature_values/Detailed/P00815_SHAP.pickle', 'asd', 1000)
shap.summary_plot(raw_shap_values, features=raw_data_values, 
                  feature_names = df_feature_importance.feature.values, plot_size = 0.6, color_bar = False, use_log_scale=False)

plt.subplot(2, 2, 4)
plt.title('His4 (max)')
raw_shap_values, raw_data_values, df_feature_importance = extract_top_features('../results/feature_values/Detailed/'+Protein_1+'_SHAP.pickle', 'max', 500)
shap.summary_plot(raw_shap_values, features=raw_data_values, cmap = create_cmap(), 
                  feature_names = df_feature_importance.feature.values, plot_size = 0.6, use_log_scale=False)
plt.tight_layout()



df_feature_importance[df_feature_importance.feature == 'ilp53']
df_feature_importance[df_feature_importance.feature == 'ilp946']
df_feature_importance[df_feature_importance.feature == 'ilp34']

subset_features = df_feature_importance.iloc[np.where(df_feature_importance.feature.isin(['ilp53','ilp946','ilp34']))[0],:]

plt.bar(subset_features.feature, subset_features.importance)
plt.ylabel('mean(|SHAP|) (average impact on model output)')
plt.xlabel('Relational feature')
plt.rcParams.update({'font.size': 11})
plt.tight_layout()







def extract_top_features(shap_path, metric, nr):
    '''
    Extract the top 10 features, according to either mean(abs(shap)) or max model output change (for ilp features)
    
    Inputs:
        shap_path: path to pickle with shap object of interst
        metric: 'max' or 'mean '
    
    Outputs:
        raw_shap_values: shap values taken from the shap-object
        raw_data_values: og. data values taken from the shap-object
        df_feature_importance: Dataframe containing the feature importance values of interest
    
    '''
    import pandas as pd
    import numpy as np
    
    
    shap_values = pd.read_pickle(shap_path)
    df_shap_values = pd.DataFrame(data=shap_values.values,columns=shap_values.feature_names)
    df_feature_importance = pd.DataFrame(columns=['feature','importance'])
    if metric == 'max':
    
        
        for col in df_shap_values.columns:
            importance = df_shap_values[col].abs().max()
            df_feature_importance.loc[len(df_feature_importance)] = [col,importance]
            
        df_feature_importance = df_feature_importance.sort_values('importance',ascending=False)
        df_feature_importance = df_feature_importance[df_feature_importance.feature.str.contains('ilp')][:nr]

    else:
        for col in df_shap_values.columns:
            importance = df_shap_values[col].abs().mean()
            df_feature_importance.loc[len(df_feature_importance)] = [col,importance]
        df_feature_importance = df_feature_importance.sort_values('importance',ascending=False)[:nr]


    ind = df_feature_importance.index.values
    raw_shap_values = np.take(shap_values.values, ind, axis=1)
    raw_data_values = np.take(shap_values.data, ind, axis=1)
    
    return raw_shap_values, raw_data_values, df_feature_importance




