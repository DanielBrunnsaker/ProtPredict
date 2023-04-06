import pandas as pd
import numpy as np
import os

os.chdir('/Volumes/Samsung_T5/ProtPredict/scripts')

#first, read performance?
data_AAILP = pd.read_csv('../results/evaluation/AAILP/AAILP_20230311.csv', index_col = 0)
data_AAILP_mean = data_AAILP.groupby(by='Protein').mean() #probably not needed

feature_shap = pd.read_csv('../results/feature_values/AAILP/Proteins_AAILP_SHAP_means_20230311.tsv', sep = '\t', index_col = 0)
feature_shap = feature_shap.iloc[:,np.where(feature_shap.columns.isin(list(data_AAILP_mean[data_AAILP_mean['R-squared'] > 0].index)) == True)[0]] #only pick those with r2 over 0

protein_dataset = pd.read_csv('data/proteomics/proteomics.csv', index_col=0)


main_predictors = pd.DataFrame(columns = ['First','Second','Third'], index = protein_dataset.columns[:2292])
for col in feature_shap.columns:
    preds = feature_shap[col].nlargest(3)
    main_predictors.loc[col] = [preds.index[0],preds.index[1],preds.index[2]]    
main_predictors = main_predictors.dropna()


##### FOR THE GO ENRICHMENT
predictors = pd.DataFrame(main_predictors['First']).reset_index()

#convert into wide?
pred_per_aa = predictors.pivot(columns='First', values='index')
pred_per_aa.to_csv('../intermediateData/predictors/AAPredictability.tsv', sep = '\t')
pd.DataFrame(protein_dataset.columns).to_csv('../intermediateData/predictors/ProteinList.tsv', sep = '\t')
pd.DataFrame(feature_shap.columns).to_csv('../intermediateData/predictors/predictableproteins.tsv', sep = '\t')
















