import pandas as pd
import numpy as np
from sklearn.model_selection import cross_validate, KFold
import os
from xgboost.sklearn import XGBRegressor
#os.chdir('../scripts') #change to scripts folder
import pickle
from sklearn.svm import SVR
from sklearn.pipeline import make_pipeline
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer

def Bayes_CV_modified(regressor, params, X, y, cvcount_inner, n_parallel, seed):
    
    from skopt import BayesSearchCV
    
    '''
    Bayesian optimization (via cross-validation)
    
    Inputs:
        params: Dict containing parameter-ranges of choice
        X: training data (independent variables, dataframe)
        y: training data (dependent variables, dataframe)
        cvcount_inner: number of folds used in optimization
        n_parallel: n_workers (int)
        seed: random seed (int)
    
    Outputs:
        best params: Self explanatory (best parameters found during optimization)
    
    '''
    
    inner_cv = KFold(n_splits=cvcount_inner, shuffle=True, random_state=seed) #gets into a local minima at 0
    
    
    regressor_bayes = regressor
    print('Performing bayesian optimization...')
    bayes=BayesSearchCV(regressor_bayes,params, #n_iter = 25,
                        scoring='neg_mean_squared_error',
                        cv=inner_cv,random_state=0, verbose = 1, n_jobs = n_parallel)
    res=bayes.fit(X,y)
    print('Optimal parameters:')
    best_params = res.best_params_
    print(best_params)
    return(best_params)

def correlation_remover(dataset):
    '''
    Remove perfectly correlated features drom a dataframe
    
    Inputs:
        dataset: Dataframe containing your dataset
    
    Outputs:
        dataset: Dataframe with correlated features (except first occurence) removed
    
    '''
    print('Removing perfectly correlated features...')
    dataset = dataset.transpose().drop_duplicates(keep = 'first').transpose()
    return(dataset)


def return_protein_dataset(dataset, protein):
    '''
    Return X and y vectors, given a specific protein
    
    Inputs:
        dataset: Entire dataset, all proteins and all features (dataframe)
        protein: Protein of choice (string)
    
    Outputs:
        X: training data, dataframe
        y: training data, array
        nancount: Number of missing rows for selected protein
    
    '''
    import numpy as np
    import pandas as pd
    
    y = dataset.iloc[:,np.where(dataset.columns == protein)[0]] #define the Y-vector, i.e. the protein of choice
    nancount = pd.isnull(y).sum()[0] #save the amount of missing values for later
    
    X = dataset.iloc[np.where(~np.isnan(y) == True)[0],np.where(dataset.columns != protein)[0]] #remove missing values from y in the corresponding rows in X
    X = X.iloc[:,2291:] #select the metabolites
    y = np.ravel(y.dropna()) #remove nan-valus from the y-vector
    return(X, y, nancount)


#define the K in the cross-validation strategies
cvcount_inner = 5 #K-fold for optimization
cvcount_outer = 5 #K-fold for evaluation

# Instantiations
phenotypedata = pd.read_csv('../data/instantiations_revised.csv',index_col = 0)
phenotypedata = phenotypedata[~phenotypedata.index.duplicated(keep='first')]

# Proteomics data
proteomic_dataset = pd.read_csv('../data/proteomics_messner.csv', index_col = 0)

# Merge
dataset = proteomic_dataset.merge(phenotypedata, left_index = True, right_index = True)
dataset = dataset[~dataset.index.duplicated(keep='first')]

#Function to define X and Y
#X, y, nancount = return_protein_dataset(dataset, 'A5Z2X5')
#X.columns = [f'feature{i}' for i in range(1, (X.shape[1]+1))] # Fixes an error of featurenames in XGB

#regr = ElasticNet(random_state = 0)
print('ElasticNet')
#params={'alpha': (0.01,1.0),
#        'l1_ratio': (0.01, 1.0),
#        'selection': ['random', 'cyclic']}

#best_params = Bayes_CV_modified(regr, params, X, y, cvcount_inner, 5, 0)

scores = pd.DataFrame(columns = ['Protein','MissingY','R-squared','RMSE'])

for Protein in dataset.columns[:2292]:
    
    X, y, nancount = return_protein_dataset(dataset, Protein)

    # define regressor and evaluate model
    #regressor = DummyRegressor(strategy="mean")
    #regressor = Pipeline(steps = [('scaler', StandardScaler()), ('svr', SVR(C = 52.916, epsilon=0.192, gamma = 'auto',kernel = 'rbf', tol = 0.065))])
    #regressor = Pipeline(steps = [('imputer', SimpleImputer(strategy='mean')), ('scaler', StandardScaler()),('elasticnet', ElasticNet(alpha = 1.0, l1_ratio = 0.876, max_iter = 5000, selection = 'cyclic'))]) #Done
    #regressor = Pipeline(steps = [('imputer', SimpleImputer(strategy='mean')), ('xgb', XGBRegressor(random_state = 0, learning_rate = 0.0212 ,colsample_bytree = 0.805, gamma = 0.154, max_depth = 3, n_estimators = 311,reg_alpha = 0.171, reg_lambda = 0.664, subsample = 0.682, n_jobs = 12))])
    regressor = ElasticNet(random_state = 0, alpha = 0.2374, l1_ratio=0.9201, selection='cyclic')
    outer_cv = KFold(n_splits=cvcount_outer, shuffle=True, random_state=0)
    score_all = cross_validate(regressor, X, np.ravel(y), cv=outer_cv, scoring = ('r2', 'neg_root_mean_squared_error'), n_jobs = 10)
    score = score_all['test_r2']
    rmse = score_all['test_neg_root_mean_squared_error']
    print(f"ElasticNet: {Protein} -> {np.mean(score)}")
    # save scores
    temp_df = pd.DataFrame([[Protein, nancount,score[0],rmse[0]], [Protein, nancount,score[1],rmse[1]], [Protein, nancount,score[2],rmse[2]], 
            [Protein, nancount,score[3],rmse[3]], [Protein, nancount,score[4],rmse[4]]], columns=['Protein', 'MissingY', 'R-squared', 'RMSE'])
    # scores = scores.append(temp_df)

    scores = pd.concat([scores, temp_df], ignore_index=True)
    scores.to_csv('../results/evaluation/comparison/ElasticNet_instantiation_revised_results.csv')
