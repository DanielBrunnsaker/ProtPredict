import pandas as pd
import numpy as np
from sklearn.model_selection import cross_validate, KFold
import os
from xgboost.sklearn import XGBRegressor
os.chdir('/Volumes/Samsung_T5/ProtPredict/scripts')
from other import Bayes_CV, extract_means_SHAP, correlation_remover, return_protein_dataset, extract_shap_return_absmeans, extract_xgb_feature_values

#os.chdir('/Users/danbru/Desktop/MetPredict')

n_parallel = 60 #If parallel, define n workers

#define the K in the cross-validation strategies
cvcount_inner = 5 #K-fold for optimization
cvcount_outer = 5 #K-fold for evaluation

#Define the dictionary for the parameter-space used in model optimization
params={'gamma': (0.0,1.0),
        'max_depth': (2, 12),
        'subsample': (0.4, 1.0),
        'learning_rate':(0.005,0.3,'log-uniform'),
        'n_estimators':(10,1000),
        'colsample_bytree':(0.4, 1.0),
        'reg_lambda':(0.0,1.0),
        'reg_alpha':(0.0,0.5)}

#load the proteomics dataset from Messner et al. (2022)
proteomic_dataset = pd.read_csv('../data/proteomics_messner.csv')
proteomic_dataset = proteomic_dataset.groupby('Unnamed: 0').mean() #Average potential duplicate rows

#load the functional metabolomics dataset from Mulleder et al. (2016)
aaDataset = pd.read_excel('../data/functional_metabolomics_mulleder.xls',sheet_name = 'intracellular_concentration_mM')
aaDataset.index = aaDataset['ORF']
aaDataset = aaDataset.iloc[:,2:21] #limit to the measured amino acid concentrations

#Merge the two on the ORF of the deletion
dataset = pd.merge(proteomic_dataset, aaDataset, left_index=True, right_index=True)

#Function to define X and Y
X, y, nancount = return_protein_dataset(dataset, 'A5Z2X5')

#define the bayesian optimization process
best_params = Bayes_CV(params, X, y, cvcount_inner, n_parallel, 42)

#predefine scoring frame
scores = pd.DataFrame(columns = ['Protein','MissingY','R-squared','RMSE'])

for Protein in dataset.columns[:2292]:
    
    X, y, nancount = return_protein_dataset(dataset, Protein)
    
    #define regressor and evaluate model
    regressor = XGBRegressor(random_state = 0, n_jobs = n_parallel, **best_params)
    outer_cv = KFold(n_splits=cvcount_outer, shuffle=True, random_state=0)
    score_all = cross_validate(regressor, X, np.ravel(y), cv=outer_cv, scoring = ('r2', 'neg_root_mean_squared_error'))
    score = score_all['test_r2']
    rmse = score_all['test_neg_root_mean_squared_error']

    #save scores
    temp_df = pd.DataFrame([[Protein, nancount,score[0],rmse[0]], [Protein, nancount,score[1],rmse[1]], [Protein, nancount,score[2],rmse[2]], 
            [Protein, nancount,score[3],rmse[3]], [Protein, nancount,score[4],rmse[4]]], columns=['Protein', 'MissingY', 'R-squared', 'RMSE'])
    scores = scores.append(temp_df)
    scores.to_csv('../intermediateData/protfromAA_results.csv')


shapmeans_proteins = pd.DataFrame(index = dataset.columns[2292:]) #predefine df for feature shap-values
for Protein in dataset.columns[:2292]:
    
    X, y, nancount = return_protein_dataset(dataset, Protein)
    
    #train regressor
    regressor = XGBRegressor(random_state = 0, n_jobs = n_parallel, **best_params)
    regressor.fit(X,y)
    
    #extract shap-values (in the form of abs(average)) for all proteins in the dataset
    shapmeans_proteins = extract_shap_return_absmeans(regressor, X, Protein, shapmeans_proteins)
    shapmeans_proteins.to_csv('../intermediateData/feature_values/AA/Proteins_AA_SHAP_means.tsv', sep = '\t')
    
    #extract all XGB-feature values (gain, total gain, ...) and save protein result each separately
    xgb_feature_values = extract_xgb_feature_values(regressor)
    xgb_feature_values.to_csv('../intermediateData/feature_values/AA/'+Protein+'_AA_XGB.tsv', sep = '\t')






