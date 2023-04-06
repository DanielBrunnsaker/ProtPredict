

#set of assistive functions

def Bayes_CV(params, X, y, cvcount_inner, n_parallel, seed):
    
    from skopt import BayesSearchCV

    #Declare the inner cross-validation strategy
    inner_cv = KFold(n_splits=cvcount_inner, shuffle=True, random_state=seed) #gets into a local minima at 0
    
    
    regressor_bayes = XGBRegressor(random_state = 0, n_jobs = n_parallel)
    print('Performing bayesian optimization...')
    bayes=BayesSearchCV(regressor_bayes,params,#n_iter = 25,
                        scoring='neg_mean_squared_error',
                        cv=inner_cv,random_state=0, verbose = 1)
    res=bayes.fit(X,np.ravel(y))
    print('Optimal parameters:')
    best_params = res.best_params_
    print(best_params)
    return(best_params)


def extract_means_SHAP(shap_values, X, Protein):
    values = abs(shap_values.values)
    avg = pd.DataFrame(values.mean(axis = 0))
    avg.index = X.columns
    avg.columns = [Protein]
    return(avg)


def correlation_remover(dataset):
    print('Removing perfectly correlated features...')
    dataset = dataset.transpose().drop_duplicates(keep = 'first').transpose()
    return(dataset)


def return_protein_dataset(dataset, protein):
    import numpy as np
    import pandas as pd
    
    y = dataset.iloc[:,np.where(dataset.columns == protein)[0]] #define the Y-vector, i.e. the protein of choice
    nancount = pd.isnull(y).sum()[0] #save the amount of missing values for later
    
    X = dataset.iloc[np.where(~np.isnan(y) == True)[0],np.where(dataset.columns != protein)[0]] #remove missing values from y in the corresponding rows in X
    X = X.iloc[:,2291:] #select the metabolites
    X = correlation_remover(X) #remove perfectly correlated features
    y = np.ravel(y.dropna()) #remove nan-valus from the y-vector
    return(X, y, nancount)


def extract_shap_return_absmeans(regressor, X, Protein, predefined_df):
    explainer = shap.Explainer(regressor, X)
    shap_values = explainer(X)
    
    shapmeans = extract_means_SHAP(shap_values, X, Protein)
    predefined_df = pd.merge(predefined_df, shapmeans, left_index = True, right_index = True, how = 'left')
    return(predefined_df)


def extract_xgb_feature_values(regressor):
    xgb_feature_values = pd.DataFrame(index = X.columns)
    types = ['weight', 'gain', 'cover', 'total_gain', 'total_cover']
    for f in types:
        fts = regressor.get_booster().get_score(importance_type= f)
        df = pd.DataFrame(list(fts.items()),columns = ['Feature',f]) 
        xgb_feature_values = pd.merge(xgb_feature_values, df, left_index = True, right_on = 'Feature').set_index('Feature')
    return(xgb_feature_values)

def create_features(path_features, path_labels, aaData):
    
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

def extract_abs_means_SHAP(shap_values, X, name):
    values = abs(shap_values.values)
    #values = shap_values.values
    avg = pd.DataFrame(values.mean(axis = 0))
    avg.index = X.columns
    avg.columns = [name]
    return(avg)

def extract_top_features(shap_path, metric):
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
        df_feature_importance = df_feature_importance[df_feature_importance.feature.str.contains('ilp')][:10]

    else:
        for col in df_shap_values.columns:
            importance = df_shap_values[col].abs().mean()
            df_feature_importance.loc[len(df_feature_importance)] = [col,importance]
        df_feature_importance = df_feature_importance.sort_values('importance',ascending=False)[:10]


    ind = df_feature_importance.index.values
    raw_shap_values = np.take(shap_values.values, ind, axis=1)
    raw_data_values = np.take(shap_values.data, ind, axis=1)
    
    return raw_shap_values, raw_data_values, df_feature_importance




def inter_from_256(x):
    import numpy as np
    return np.interp(x=x,xp=[0,255],fp=[0,1])

def create_cmap():
    
    import numpy as np
    from matplotlib import colors
    from matplotlib import cm

    cdict = {
        'red':((0.0,inter_from_256(64),inter_from_256(64)),
               (1/5*1,inter_from_256(112),inter_from_256(112)),
               (1/5*2,inter_from_256(230),inter_from_256(230)),
               (1/5*3,inter_from_256(253),inter_from_256(253)),
               (1/5*4,inter_from_256(244),inter_from_256(244)),
               (1.0,inter_from_256(169),inter_from_256(169))),
        'green': ((0.0, inter_from_256(57), inter_from_256(57)),
                (1 / 5 * 1, inter_from_256(198), inter_from_256(198)),
                (1 / 5 * 2, inter_from_256(241), inter_from_256(241)),
                (1 / 5 * 3, inter_from_256(219), inter_from_256(219)),
                (1 / 5 * 4, inter_from_256(109), inter_from_256(109)),
                (1.0, inter_from_256(23), inter_from_256(23))),
        'blue': ((0.0, inter_from_256(144), inter_from_256(144)),
                  (1 / 5 * 1, inter_from_256(162), inter_from_256(162)),
                  (1 / 5 * 2, inter_from_256(246), inter_from_256(146)),
                  (1 / 5 * 3, inter_from_256(127), inter_from_256(127)),
                  (1 / 5 * 4, inter_from_256(69), inter_from_256(69)),
                  (1.0, inter_from_256(69), inter_from_256(69))),
    }
        
    new_cmap = colors.LinearSegmentedColormap('new_cmap',segmentdata=cdict)
    
    return(new_cmap)


