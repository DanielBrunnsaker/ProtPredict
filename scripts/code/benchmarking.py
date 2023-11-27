import pandas as pd
import numpy as np
from sklearn.model_selection import cross_validate, KFold
import os
from xgboost.sklearn import XGBRegressor
#os.chdir('/Users/danbru/Desktop/Projects/ProtPredict/scripts') #change to scripts folder
os.chdir('/Users/danbru/Projects/ProtPredict/scripts') #change to scripts folder
import pickle

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
    X = correlation_remover(X) #remove perfectly correlated features
    y = np.ravel(y.dropna()) #remove nan-valus from the y-vector
    return(X, y, nancount)

def create_features(path_features, path_labels, aaData):
    '''
    Merge the output from aleph with a standard propositional dataset. Merges on the ORF
    
    Inputs:
        path_features: path to the .txt file containing the binary descriptors/features
        path_labels: path to the .txt file containing the row-labels (ORFs)
        aaData = propositional dataset (as a dataframe, index must containg matching values to path_labels)
    
    Outputs:
        dataset: dataframe with our complete dataset (propositional data and data descriptors)
    
    '''
    import pandas as pd
    import numpy as np
    
    ilp_neg_features = pd.read_csv(path_features, sep = ' ', header = None)
    ilp_neg_labels = pd.read_csv(path_labels, header = None)

    start = "networked('"
    end = "')"

    #go through all positive examples and rename
    for i in range(len(ilp_neg_labels)):
        s = ilp_neg_labels.iloc[i][0]
        ilp_neg_labels.iloc[i] = s[s.find(start)+len(start):s.rfind(end)]


    nrfeat = len(ilp_neg_features.columns)
    df = pd.DataFrame({'Feature': ['ilp']*nrfeat, 'Nr': list(map(str, np.arange(1,nrfeat+1)))})
    df['featureName'] = df[['Feature', 'Nr']].apply(lambda x: ''.join(x), axis=1)

    ilp_neg_features.columns = df['featureName']
    ilp_neg_features.index = list(ilp_neg_labels.iloc[:,0])
    dataset = pd.merge(aaData, ilp_neg_features, left_index=True,right_index=True)
    dataset = dataset.iloc[:,:-1]
    return(dataset)

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
    bayes=BayesSearchCV(regressor_bayes,params,#n_iter = 25,
                        scoring='neg_mean_squared_error',
                        cv=inner_cv,random_state=0, verbose = 1, n_jobs = n_parallel)
    res=bayes.fit(X,y)
    print('Optimal parameters:')
    best_params = res.best_params_
    print(best_params)
    return(best_params)

n_parallel = 60 #If parallel, define n workers

#define the K in the cross-validation strategies
cvcount_inner = 5 #K-fold for optimization
cvcount_outer = 5 #K-fold for evaluation

#Define the dictionary for the parameter-space used in model optimization
from sklearn.svm import SVR
from sklearn.pipeline import make_pipeline
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
# Monkey patch for issues with later versons of numpy
import numpy as np
np.int = int

#regr = make_pipeline(steps = [('scaler', StandardScaler()), ('svr', SVR())])
#regr = Pipeline(steps = [('scaler', StandardScaler()), ('svr', SVR())])
#regr = RandomForestRegressor(random_state = 0, n_jobs = 10)

#params={'svr__gamma': ['scale','auto'],
#        'svr__C': (0.1, 100),
#        'svr__epsilon': (0.001, 1.0),
#        'svr__kernel': ['linear','rbf','sigmoid'],
#        'svr__coef0': (-1,1),
#        'svr__degree': (1,7), 
#        'svr__tol': (1e-3,1e-1)} 

#regr = Pipeline(steps = [('scaler', StandardScaler()), ('elasticnet', ElasticNet(random_state = 0, max_iter = 5000))])
#params={'elasticnet__alpha': (0.01,1.0),
#        'elasticnet__l1_ratio': (0.01, 1.0),
#        'elasticnet__selection': ['random', 'cyclic']}

#params={'n_estimators': (50,1000),
#        'max_depth': (5, 20),
#        'min_samples_split':(2,20),
#        'min_samples_leaf':(2,10),
#        'max_features':['auto', 'sqrt','log2']}

#load the proteomics dataset from Messner et al. (2022)
proteomic_dataset = pd.read_csv('../data/proteomics_messner.csv')
proteomic_dataset = proteomic_dataset.groupby('Unnamed: 0').mean() #Average potential duplicate rows

#paths to relational features
path_features = '../intermediateData/generated_features/proteomics_features.txt' #ordered features
path_labels = '../intermediateData/generated_features/proteomics.txt' #ordered positve examples

#reformat the relational features and merge them with the proteomics dataset
dataset = create_features(path_features, path_labels, proteomic_dataset)
scores = pd.DataFrame(columns = ['Protein','MissingY','R-squared','RMSE'])

#define the bayesian optimization process
#Function to define X and Y
#X, y, nancount = return_protein_dataset(dataset, 'A5Z2X5')

# best_params = Bayes_CV(params, X, y, cvcount_inner, n_parallel, 0)
best_params = Bayes_CV_modified(regr, params, X, y, cvcount_inner, 5, 0)

with open('../params/rf_params_ilp.pickle', 'wb') as handle:
    pickle.dump(best_params, handle, protocol=pickle.HIGHEST_PROTOCOL)

'''

for Protein in dataset.columns[:2292]:
    
    X, y, nancount = return_protein_dataset(dataset, Protein)

    # define regressor and evaluate model
    #regressor = DummyRegressor(strategy="mean")
    #regressor = Pipeline(steps = [('scaler', StandardScaler()), ('svr', SVR(C = 36.96, epsilon=1, gamma = 'auto',kernel = 'rbf', tol = 0.1))])
    #regressor = RandomForestRegressor(random_state = 0, max_depth = 12, max_features = 'sqrt', min_samples_leaf = 2, min_samples_split = 17, n_estimators = 50, n_jobs = 5)
    regressor = Pipeline(steps = [('scaler', StandardScaler()), ('elasticnet', ElasticNet(alpha = 1.0, l1_ratio = 0.011, max_iter = 5000, selection = 'cyclic'))]) 
    outer_cv = KFold(n_splits=cvcount_outer, shuffle=True, random_state=0)
    score_all = cross_validate(regressor, X, np.ravel(y), cv=outer_cv, scoring = ('r2', 'neg_root_mean_squared_error'), n_jobs = 5)
    score = score_all['test_r2']
    rmse = score_all['test_neg_root_mean_squared_error']
    print(f"ElasticNet: {Protein} -> {np.mean(score)}")
    # save scores
    temp_df = pd.DataFrame([[Protein, nancount,score[0],rmse[0]], [Protein, nancount,score[1],rmse[1]], [Protein, nancount,score[2],rmse[2]], 
            [Protein, nancount,score[3],rmse[3]], [Protein, nancount,score[4],rmse[4]]], columns=['Protein', 'MissingY', 'R-squared', 'RMSE'])
    # scores = scores.append(temp_df)

    scores = pd.concat([scores, temp_df], ignore_index=True)
    scores.to_csv('../results/evaluation/comparison/elasticnet_ilp_results.csv')



'''


