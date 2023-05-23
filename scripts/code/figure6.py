import pandas as pd
import os
from xgboost.sklearn import XGBRegressor
import shap
import pickle
import matplotlib.pyplot as plt

os.chdir('') #change to scripts folder
from other import create_features, return_protein_dataset, extract_top_features, create_cmap

#choose teo proteins to visualize
Protein_1 = 'P00815' #onlyILP
Protein_2 = 'P07342' #AA+ILP

n_parallel = 4

#load proteomics dataset (Messner et al. 2022)
proteomic_dataset = pd.read_csv('../data/proteomics_messner.csv')
proteomic_dataset = proteomic_dataset.groupby('Unnamed: 0').mean() #Average potential duplicate rows

#load the functional metabolomics dataset from Mulleder et al. (2016)
aaDataset = pd.read_excel('../data/functional_metabolomics_mulleder.xls',sheet_name = 'intracellular_concentration_mM')
aaDataset.index = aaDataset['ORF']
aaDataset = aaDataset.iloc[:,2:21] #limit to the measured amino acid concentrations

#Merge the two on the ORF of the deletion
unified_dataset = pd.merge(proteomic_dataset, aaDataset, left_index=True, right_index=True)

#paths to relational features, note that the features derived from Mulleder et al. has been removed prior to generation
path_features = '../intermediateData/generated_features/proteomics_noAA_features.txt' #ordered features
path_labels = '../intermediateData/generated_features/proteomics_noAA.txt' #ordered positive examples

#Reformat and merge the unified dataset with relational features
dataset = create_features(path_features, path_labels, unified_dataset)
X_P1, y_P1, nancount = return_protein_dataset(dataset, Protein_1)

#paths to relational features, note that the features derived from Mulleder et al. has been removed prior to generation
path_features = '../intermediateData/generated_features/proteomics_features.txt' #ordered features
path_labels = '../intermediateData/generated_features/proteomics.txt' #ordered positive examples

#Reformat and merge the unified dataset with relational features
dataset = create_features(path_features, path_labels, proteomic_dataset)
X_P2, y_P2, nancount = return_protein_dataset(dataset, Protein_2)

#these were derived in previous steps (best parameters from the bayesian optimization)
best_params_P1 = pd.read_pickle('../intermediateData/parameters/GCV3_AAILP_20230311.pickle')
best_params_P2 = pd.read_pickle('../intermediateData/parameters/GCV3_protILP_20230309.pickle')


#Protein 1
regressor = XGBRegressor(random_state = 0, n_jobs = n_parallel, **best_params_P1)
regressor.fit(X_P1,y_P1)
explainer = shap.Explainer(regressor, X_P1)
shap_values = explainer(X_P1)
#save as pickle
with open('../results/feature_values/Detailed/'+Protein_1+'_SHAP.pickle', 'wb') as handle:
    pickle.dump(shap_values, handle, protocol=pickle.HIGHEST_PROTOCOL)


#Protein 2
regressor = XGBRegressor(random_state = 0, n_jobs = n_parallel, **best_params_P2)
regressor.fit(X_P2,y_P2)
explainer = shap.Explainer(regressor, X_P2)
shap_values = explainer(X_P2)
#save as pickle
with open('../results/feature_values/Detailed/'+Protein_2+'_SHAP.pickle', 'wb') as handle:
    pickle.dump(shap_values, handle, protocol=pickle.HIGHEST_PROTOCOL)

#Visualize the shap-values
#P07342
fig = plt.figure(figsize = (10,6))
plt.subplot(2, 2, 1)
raw_shap_values, raw_data_values, df_feature_importance = extract_top_features('../results/feature_values/Detailed/'+Protein_2+'_SHAP.pickle', 'asd')
shap.summary_plot(raw_shap_values, features=raw_data_values, cmap = create_cmap(), 
                  feature_names = df_feature_importance.feature.values, plot_size = 0.6, color_bar = False, use_log_scale=False)

plt.subplot(2, 2, 2)
raw_shap_values, raw_data_values, df_feature_importance = extract_top_features('../results/feature_values/Detailed/'+Protein_2+'_SHAP.pickle', 'max')
shap.summary_plot(raw_shap_values, features=raw_data_values, cmap = create_cmap(), 
                  feature_names = df_feature_importance.feature.values, plot_size = 0.6, use_log_scale=False)

#P00815
#fig = plt.figure()
plt.subplot(2, 2, 3)
raw_shap_values, raw_data_values, df_feature_importance = extract_top_features('../results/feature_values/Detailed/'+Protein_1+'_SHAP.pickle', 'asd')
shap.summary_plot(raw_shap_values, features=raw_data_values, cmap = create_cmap(), 
                  feature_names = df_feature_importance.feature.values, plot_size = 0.6, color_bar = False, use_log_scale=False)

plt.subplot(2, 2, 4)
raw_shap_values, raw_data_values, df_feature_importance = extract_top_features('../results/feature_values/Detailed/'+Protein_1+'_SHAP.pickle', 'max')
shap.summary_plot(raw_shap_values, features=raw_data_values, cmap = create_cmap(), 
                  feature_names = df_feature_importance.feature.values, plot_size = 0.6, use_log_scale=False)
plt.tight_layout()
plt.savefig('Figure6.pdf', dpi = 600)