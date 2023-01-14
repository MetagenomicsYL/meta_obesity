####################### train the XGBoost model ###########################
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier
from sklearn import model_selection
from sklearn.base import is_classifier, clone
from sklearn.utils import indexable, check_random_state, _safe_indexing
from sklearn.model_selection import check_cv
from sklearn.metrics import check_scoring
from joblib import Parallel, logger
from sklearn.utils.fixes import delayed
from sklearn.utils.metaestimators import _safe_split
from sklearn.utils.validation import _check_fit_params
from utils import train_data, obese
from utils import permutation_test_score_eachfold
from utils import expand_grid
from utils import crossvalidate
from utils import median
from utils import validate
from utils import _shuffle
from utils import cross_validate_with_permutation_test_score
from pykliep import DensityRatioEstimator


#grid search for the best parameter

def xgb_10crossvalidation(all_indiv_df, param_xgb):
    xgb_10cv = pd.DataFrame(columns = ['n_estimators', 'max_depth', 'learning_rate', 'auc'])
    xtrain = all_indiv_df.iloc[:, :len(all_indiv_df.columns)-2]
    ytrain = all_indiv_df["group"]
    for i in range(len(param_xgb)):
        estimators = param_xgb.loc[i, "n_estimators"]
        depth = param_xgb.loc[i, "max_depth"]
        rate = param_xgb.loc[i, "learning_rate"]
        model = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate, n_jobs = -1, random_state = 2021)
        aucscore = crossvalidate(model, xtrain, ytrain,10)
        df = pd.DataFrame({"n_estimators": estimators, "max_depth": depth, "learning_rate": rate, "auc": median(aucscore)}, index = [0])
        xgb_10cv = pd.concat([xgb_10cv, df])
    return xgb_10cv

#use the best paramter for 10-CV

def model_bestpara_10cv(xgb_10cv, all_indiv_df):
    xtrain = all_indiv_df.iloc[:, :len(all_indiv_df.columns)-2]
    ytrain = all_indiv_df["group"]
    best_para = xgb_10cv[xgb_10cv.auc == xgb_10cv.auc.max()]
    best_para = best_para.reset_index()
    estimators = best_para.loc[0, "n_estimators"]
    depth = best_para.loc[0, "max_depth"]
    rate = best_para.loc[0, "learning_rate"]  
    model = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate, n_jobs = -1, random_state = 2021)
    aucscore, Pval = cross_validate_with_permutation_test_score(model, xtrain, ytrain, 10)
    xgb_10cv_best = pd.DataFrame({"Model": "XGBoost", "fold": [i for i in range(1,11,1)], "AUC": aucscore, "Pval": Pval, "Train": "10-CV", "Study": "Combined"})
    return xgb_10cv_best


def gridsearch_xgb(param, X_train, X_test, y_train, y_test):
    xgb_within = pd.DataFrame(columns = ['n_estimators', 'max_depth', 'learning_rate', 'auc'])
    for i in range(len(param)):
        estimators = param.loc[i, "n_estimators"]
        depth = param.loc[i, "max_depth"]
        rate = param.loc[i, "learning_rate"]   
        clf = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate, n_jobs = -1, random_state = 2021)
        auroc1 = validate(clf, np.array(X_train), np.array(X_test), y_train, y_test)
        df = pd.DataFrame({"n_estimators": estimators, "max_depth": depth, "learning_rate": rate, "auc": auroc1}, index = [0])
        xgb_within = pd.concat([xgb_within, df])
    xgb_within = xgb_within.reset_index()
    xgb_within.drop('index', axis = 1, inplace = True)
    return xgb_within

# train the model within each study 

def train_within(all_indiv_df, param_xgb, pvalname):
    alist = all_indiv_df["StudyID"].unique()
    xgb_withsty = pd.DataFrame(columns = ['n_estimators', 'max_depth', 'learning_rate', 'auc', 'StudyID', 'Pval'])
    for s in alist:
        print(s)
        train = all_indiv_df[all_indiv_df["StudyID"] == s]
        train = train.reset_index()
        train.drop("index", axis=1, inplace=True)
        X = train[pvalname]
        Y = train["group"]
        X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Y, test_size = 0.25, random_state = 2021)
        xgb_within = gridsearch_xgb(param_xgb, X_train, X_test, y_train, y_test)
        maxline = xgb_within[xgb_within.auc == xgb_within.auc.max()]
        maxline = maxline.reset_index()
        maxline.drop("index", axis = 1, inplace = True)
        maxline["StudyID"] = s
        maxline = maxline.loc[:0, :]
        estimators = maxline.loc[0, "n_estimators"]
        depth = maxline.loc[0, "max_depth"]
        rate = maxline.loc[0, "learning_rate"]
        trainix = np.array(X_train.index)
        testix = np.array(X_test.index)
        X = np.array(train[pvalname])
        Y = np.array(train["group"])  
        clf = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate, n_jobs = -1, random_state = 2021)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_eachfold(
            clf, X, Y, scoring = "roc_auc", cv = 2, train = trainix, test = testix, n_permutations = 100)
        maxline["Pval"] = pvalue_train
        xgb_withsty = pd.concat([xgb_withsty, maxline])
    return xgb_withsty

# train the model using loov 

def train_loov(all_indiv_df, param_xgb, pvalname):
    alist = all_indiv_df["StudyID"].unique()  
    xgb_loosty = pd.DataFrame(columns = ['n_estimators', 'max_depth', 'learning_rate', 'auc', 'StudyID', 'Pval'])
    for s in alist:
        print(s)
        train = all_indiv_df[all_indiv_df["StudyID"]!= s]
        test = all_indiv_df[all_indiv_df["StudyID"] == s]
        xtrain = train[pvalname]
        ytrain = train["group"]
        xtest = test[pvalname]
        ytest = test["group"]
        xgb_within = gridsearch_xgb(param_xgb, xtrain, xtest, ytrain, ytest)
        maxline = xgb_within[xgb_within.auc == xgb_within.auc.max()]
        maxline = maxline.reset_index()
        maxline.drop("index", axis = 1, inplace = True)
        maxline["StudyID"] = s
        maxline = maxline.loc[:0,:]
        trainix = np.array(train.index)
        testix = np.array(test.index)
        X = np.array(all_indiv_df[pvalname])
        Y = np.array(all_indiv_df["group"])
        estimators = maxline.loc[0,"n_estimators"]
        depth = maxline.loc[0,"max_depth"]
        rate = maxline.loc[0,"learning_rate"]
        clf = XGBClassifier(n_estimators=estimators,max_depth=depth,learning_rate=rate, subsample=0.9, n_jobs = -1, random_state=0)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_eachfold(
            clf, X, Y, scoring="roc_auc", cv = 2, train=trainix, test=testix, n_permutations=100)
        maxline["Pval"] = pvalue_train
        xgb_loosty = pd.concat([xgb_loosty, maxline])
    return xgb_loosty


############## correcting covariate shift for the model training #############

#permutation test for model correcting covariate shift

def _permutation_test_score_covariate(estimator, X, y, groups, cv, train, test, weights, scorer, fit_params):
    fit_params = fit_params if fit_params is not None else {}
    avg_score = []
    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, y_test = _safe_split(estimator, X, y, test, train)
    if len(pd.Series(y_test).unique()) != 2:
        return ("g")
    else:
        fit_params = _check_fit_params(X, fit_params, train)
        estimator.fit(X_train, y_train, sample_weight = weights, **fit_params)
        avg_score.append(scorer(estimator, X_test, y_test))
        return np.mean(avg_score)
 
    
def permutation_test_score_covariate(
    estimator,
    X,
    y,
    *,
    weights,
    groups = None,
    cv = None,
    train, 
    test,
    n_permutations = 100,
    n_jobs = None,
    random_state = 0,
    verbose= 0,
    scoring= None,
    fit_params= None,
):

    X, y, groups = indexable(X, y, groups)

    cv = check_cv(cv, y, classifier=is_classifier(estimator))
    scorer = check_scoring(estimator, scoring=scoring)
    random_state = check_random_state(random_state)
    score = _permutation_test_score_covariate(
        clone(estimator), X, y, groups, cv, train, test, weights, scorer, fit_params=fit_params)
    permutation_scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_permutation_test_score_covariate)(
            clone(estimator),
            X,
            _shuffle(y, groups, random_state),
            groups,
            cv,
            train, 
            test,
            weights,
            scorer,
            fit_params=fit_params,
        )
        for _ in range(n_permutations)
    )
    dupcount = permutation_scores.count("g")
    dup = np.array(permutation_scores)
    ix = np.where(dup == "g")
    for index in reversed(ix[0]):
        permutation_scores.pop(index)
    n_permutations = n_permutations-dupcount
    permutation_scores = np.array(permutation_scores)
    pvalue = (np.sum(permutation_scores >= score) + 1.0) / (n_permutations + 1)
    return score, permutation_scores, pvalue

#grid search for the best parameter for 10-CV with correcting covariate shift

def xgb_10cv_grid_search_covariate_shift(all_indiv_df, param_xgb):
    xgb_10cv = pd.DataFrame(columns = ['para_com', 'n_estimators', 'max_depth', 'learning_rate', 'auc', 'fold'])
    xtrain = all_indiv_df.iloc[:, :len(all_indiv_df.columns)-2]
    ytrain = all_indiv_df["group"]
    cv = 10
    m = 0
    kfold1=StratifiedKFold(n_splits = cv, shuffle = True, random_state = 2021)
    for train, test in kfold1.split(xtrain, ytrain):
        m = m + 1
        X_train_xval = np.array(xtrain)[train, :].astype(float)
        X_test_xval = np.array(xtrain)[test, :].astype(float)
        y_train_xval = np.array(ytrain)[train].astype(float)
        y_test_xval = np.array(ytrain)[test].astype(float)
        kliep = DensityRatioEstimator()
        kliep.fit(X_train_xval, X_test_xval)
        weights = kliep.predict(X_train_xval)
        for i in range(len(param_xgb)):
            estimators = param_xgb.loc[i, "n_estimators"]
            depth = param_xgb.loc[i, "max_depth"]
            rate = param_xgb.loc[i, "learning_rate"]
            model = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate,  n_jobs = -1, random_state = 2021)
            model.fit(X_train_xval, y_train_xval, sample_weight = weights)
            y_pred = model.predict(X_test_xval)
            aucscore = metrics.roc_auc_score(y_test_xval, y_pred)        
            df = pd.DataFrame({"para_com": i, "n_estimators": estimators, "max_depth": depth, "learning_rate": rate, "auc": aucscore, "fold": m},index = [0])
            xgb_10cv = pd.concat([xgb_10cv, df])
            
        xgb_10cv_summ = pd.DataFrame(columns = ['para_com', 'n_estimators', 'max_depth', 'learning_rate', 'auc'])
        for cl in xgb_10cv["para_com"].unique():
            each_para = xgb_10cv[xgb_10cv["para_com"] == cl]
            each_para = each_para.reset_index()
            each_para.drop("index",axis = 1, inplace = True)
            df_summ = pd.DataFrame({"para_com": cl, 
                                    "n_estimators":each_para.loc[0, "n_estimators"], 
                                    "max_depth":each_para.loc[0, "max_depth"], 
                                    "learning_rate":each_para.loc[0, "learning_rate"], 
                                    "auc": median(each_para.loc[:, "auc"].to_list())}, 
                                   index = [0])
            xgb_10cv_summ = pd.concat([xgb_10cv_summ, df_summ])
    return xgb_10cv_summ


# train model using best parameter for 10-CV correcting covariate shift

def CV_covariate_shift(xgb_10cv_summ, traindt, ytarget,cv):   
    maxline = xgb_10cv_summ[xgb_10cv_summ.auc == xgb_10cv_summ.auc.max()]
    maxline = maxline.reset_index()
    maxline.drop("index",axis = 1, inplace = True)
    maxline = maxline.loc[:0, :]
    estimators = maxline.loc[0, "n_estimators"]
    depth = maxline.loc[0, "max_depth"]
    rate = maxline.loc[0, "learning_rate"]
    model = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate,  n_jobs = -1, random_state = 2021)
    aucscore = []
    Pval = []
    kfold1 = StratifiedKFold(n_splits = cv, shuffle = True, random_state = 2021)
    for train,test in kfold1.split(traindt, ytarget):
        X_train_xval = np.array(traindt)[train, :].astype(float)
        X_test_xval = np.array(traindt)[test, :].astype(float)
        y_train_xval = np.array(ytarget)[train].astype(float)
        y_test_xval = np.array(ytarget)[test].astype(float)
        kliep = DensityRatioEstimator()
        kliep.fit(X_train_xval, X_test_xval)
        weights= kliep.predict(X_train_xval)
        X = np.array(traindt)
        Y = np.array(ytarget)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_covariate(
            model, X, Y, scoring = "roc_auc", cv = 2, weights=weights, train = train, test = test, n_permutations = 100)
        model.fit(X_train_xval, y_train_xval, sample_weight = weights)
        y_pred = model.predict(X_test_xval)
        aucscore.append(metrics.roc_auc_score(y_test_xval, y_pred))
        Pval.append(pvalue_train)
    xgb_10cv_covariate_shift = pd.DataFrame({"Model": "XGBoost", "fold": [i for i in range(1,11,1)], "AUC": aucscore, "Pval": Pval, "Train": "10-CV_Covariate_shift", "Study": "Combined"})
    return xgb_10cv_covariate_shift


def gridsearch_xgb_covariate(param, X_train, X_test, y_train, y_test):
    xgb_within = pd.DataFrame(columns = ['n_estimators', 'max_depth', 'learning_rate', 'auc'])
    kliep = DensityRatioEstimator()
    kliep.fit(X_train, X_test)
    weights = kliep.predict(X_train)
    for i in range(len(param)):
        estimators = param.loc[i, "n_estimators"]
        depth = param.loc[i, "max_depth"]
        rate = param.loc[i, "learning_rate"]
        clf = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate, subsample = 0.85, n_jobs = -1, random_state = 2021)
        clf.fit(X_train, y_train, sample_weight = weights)
        y_pred = clf.predict(X_test)
        auroc = metrics.roc_auc_score(y_test, y_pred)
        df = pd.DataFrame({"n_estimators": estimators, "max_depth": depth,"learning_rate": rate, "auc": auroc}, index = [0])
        xgb_within = pd.concat([xgb_within, df])
    xgb_within = xgb_within.reset_index()
    xgb_within.drop('index', axis = 1, inplace = True)
    return xgb_within

# train the model within each study correcting covariate shift

def train_within_covariate_shift(all_indiv_df, param_xgb, pvalname):
    alist = all_indiv_df["StudyID"].unique()
    xgb_withsty = pd.DataFrame(columns = ['n_estimators', 'max_depth', 'learning_rate', 'auc', 'StudyID'])
    for s in alist:
        print(s)
        train = all_indiv_df[all_indiv_df["StudyID"] == s]
        X = np.array(train[pvalname]).astype(float)
        Y = np.array(train["group"]).astype(float)
        X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Y, test_size = 0.25, random_state = 2021)
        xgb_within = gridsearch_xgb_covariate(param_xgb, X_train, X_test, y_train, y_test)
        maxline = xgb_within[xgb_within.auc == xgb_within.auc.max()]
        maxline = maxline.reset_index()
        maxline.drop("index", axis = 1, inplace=True)
        maxline["StudyID"] = s
        maxline = maxline.loc[:0, :]
        xgb_withsty = pd.concat([xgb_withsty, maxline])
    return xgb_withsty


def within_study_covariate_best(all_indiv_df, xgb_withsty, pvalname):
    alist=all_indiv_df["StudyID"].unique()  
    xgb_withsty_cov = pd.DataFrame(columns = ['Pval', 'auc', 'StudyID']) 
    for s in alist:
        train = all_indiv_df[all_indiv_df["StudyID"] == s]
        train = train.reset_index()
        train.drop(columns = ['index'], inplace = True)    
        X = train[pvalname]
        Y = train["group"]
        param = xgb_withsty[xgb_withsty["StudyID"] == s]
        param = param.reset_index()
        estimators = param.loc[0, "n_estimators"]
        depth = param.loc[0, "max_depth"]
        rate = param.loc[0, "learning_rate"]
        X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Y, test_size = 0.25, random_state = 2021)
        trainix = np.array(X_train.index)
        testix = np.array(X_test.index)
        clf = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate, subsample=0.85, n_jobs = -1, random_state = 2021)
        X = np.array(train[pvalname])
        Y = np.array(train["group"])
        X_train = np.array(X_train).astype(float)
        X_test = np.array(X_test).astype(float)
        y_train = np.array(y_train).astype(float)
        y_test = np.array(y_test).astype(float)
        kliep = DensityRatioEstimator()
        kliep.fit(X_train, X_test)
        weights= kliep.predict(X_train)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_covariate(
            clf, X, Y, scoring = "roc_auc", weights = weights, train = trainix, test = testix, n_permutations = 100)
        clf.fit(X_train, y_train, sample_weight = weights)
        y_pred = clf.predict(X_test)
        auroc = metrics.roc_auc_score(y_test, y_pred)
        df = pd.DataFrame({'Pval': pvalue_train, 'auc': auroc, 'StudyID': s}, index = [0])
        xgb_withsty_cov = pd.concat([xgb_withsty_cov, df])
    return xgb_withsty_cov
   
 
# train the model using LOOV correcting covariate shift

def train_loov_covariate_shift(all_indiv_df, param_xgb, pvalname):
    alist = all_indiv_df["StudyID"].unique()
    xgb_loov = pd.DataFrame(columns = ['n_estimators', 'max_depth', 'learning_rate', 'auc', 'StudyID'])
    for s in alist:
        print(s)
        train = all_indiv_df[all_indiv_df["StudyID"] != s]
        train = train.reset_index()
        train.drop("index",axis = 1, inplace = True)
        test = all_indiv_df[all_indiv_df["StudyID"] == s]
        test= test.reset_index()
        test.drop("index", axis = 1, inplace = True)
        xtrain = np.array(train[pvalname]).astype(float)
        ytrain = np.array(train["group"]).astype(float)
        xtest = np.array(test[pvalname]).astype(float)
        ytest = np.array(test["group"]).astype(float)
        xgb_loov_each = gridsearch_xgb_covariate(param_xgb, xtrain, xtest, ytrain, ytest)
        maxline = xgb_loov_each[xgb_loov_each.auc == xgb_loov_each.auc.max()]
        maxline = maxline.reset_index()
        maxline.drop("index", axis = 1, inplace = True)
        maxline["StudyID"] = s
        maxline = maxline.loc[:0, :]
        xgb_loov = pd.concat([xgb_loov, maxline])
    return xgb_loov


def loov_covariate_shift_best(all_indiv_df, bestparameter, pvalname):
    alist = all_indiv_df["StudyID"].unique()
    xgb_loosty = pd.DataFrame(columns = ['Pval', 'auc', 'StudyID'])
    for s in alist:
        print(s)
        param = bestparameter[bestparameter["StudyID"] == s]
        param = param.reset_index()
        estimators = param.loc[0,"n_estimators"]
        depth = param.loc[0,"max_depth"]
        rate = param.loc[0,"learning_rate"]
        train = all_indiv_df[all_indiv_df["StudyID"] != s]
        train = train.reset_index()
        train.drop("index",axis = 1, inplace = True)
        test = all_indiv_df[all_indiv_df["StudyID"] == s]
        test= test.reset_index()
        test.drop("index",axis = 1, inplace = True)
        xtrain = np.array(train[pvalname]).astype(float)
        ytrain = np.array(train["group"]).astype(float)
        xtest = np.array(test[pvalname]).astype(float)
        ytest = np.array(test["group"]).astype(float)
        kliep = DensityRatioEstimator()
        kliep.fit(xtrain, xtest)
        weights= kliep.predict(xtrain)
        trainix = np.array(train.index)
        testix = np.array(test.index)  
        X = np.array(all_indiv_df[pvalname])
        Y = np.array(all_indiv_df["group"])
        clf = XGBClassifier(n_estimators = estimators, max_depth = depth, learning_rate = rate, subsample=0.85, n_jobs = -1, random_state = 2021)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_covariate(
            clf, X, Y, scoring = "roc_auc", weights = weights, train = trainix, test = testix, n_permutations = 100)
        clf.fit(xtrain, ytrain, sample_weight = weights)
        y_pred = clf.predict(xtest)
        auroc = metrics.roc_auc_score(ytest, y_pred)
        df = pd.DataFrame({'Pval': pvalue_train, 'auc': auroc, 'StudyID': s}, index = [0])
        xgb_loosty = pd.concat([xgb_loosty, df])
    return xgb_loosty


def main():
    meta = pd.read_csv("meta_BMI_obesity.csv")
    table = pd.read_csv("genus_combined_rel_rdp_filter.csv")
    pvalname = table["Genus"].str.replace("g__", "")
    
    # train model using all studies
    
    all_indiv_df = train_data(table, meta, pvalname)
    all_indiv_df['group'] = all_indiv_df['group'].apply(lambda x:obese(x))
    params = {'max_depth':range(3, 11), 'n_estimators':range(100, 1100, 100), 'learning_rate':[0.05, 0.1, 0.25, 0.5, 1.0]}
    param_xgb = expand_grid(params)
    all_xgb_10cv = xgb_10crossvalidation(all_indiv_df, param_xgb)
    all_xgb_10cv.to_csv("all_XGB_without_CS_10CV.csv", index=False)
    
    all_xgb_10cv_best = model_bestpara_10cv(all_xgb_10cv, all_indiv_df)
    all_xgb_10cv_best.to_csv("all_XGB_without_CS_best_parameter_10CV.csv", index=False)
    
    all_xgb_withsty = train_within(all_indiv_df, param_xgb, pvalname)
    all_xgb_withsty.to_csv("all_XGB_without_CS_within.csv", index=False)
    
    all_xgb_loosty = train_loov(all_indiv_df, param_xgb, pvalname)
    all_xgb_loosty.to_csv("all_XGB_without_CS_LOOV.csv", index=False)
    
    
    xgb_10cv_summ = xgb_10cv_grid_search_covariate_shift(all_indiv_df, param_xgb)
    traindt = all_indiv_df.copy()
    ytarget = all_indiv_df["group"]
    traindt.drop("group", axis=1, inplace=True)
    traindt.drop("StudyID", axis=1, inplace=True)
    all_10CV_cov_shift = CV_covariate_shift(xgb_10cv_summ, traindt, ytarget, 10)
    all_10CV_cov_shift.to_csv("all_XGB_with_CS_10CV.csv", index=False)
    
    all_grid_cov_shift_within = train_within_covariate_shift(all_indiv_df, param_xgb, pvalname)
    all_within_cov_shift = within_study_covariate_best(all_indiv_df, all_grid_cov_shift_within, pvalname)
    all_within_cov_shift.to_csv("all_XGB_with_CS_within.csv", index=False)
    
    all_grid_cov_shift_loov = train_loov_covariate_shift(all_indiv_df, param_xgb, pvalname)
    all_loov_cov_shift = loov_covariate_shift_best(all_indiv_df, all_grid_cov_shift_loov, pvalname)
    all_loov_cov_shift.to_csv("all_XGB_with_CS_loov.csv", index=False)
    
    
    #train model using obesity association studies
    
    obelist = ["Olsson", "Zupancic", "Ross", "Lippert", "Barengolts", "Chavez", "Gao", "Ahmad"]
    metaobe = meta[meta["StudyID"].isin(obelist)]
    obe_indiv_df = train_data(table, metaobe, pvalname)
    
    obe_indiv_df['group'] = obe_indiv_df['group'].apply(lambda x:obese(x))
    
    obe_xgb_10cv = xgb_10crossvalidation(obe_indiv_df, param_xgb)
    obe_xgb_10cv_best = model_bestpara_10cv(obe_xgb_10cv, obe_indiv_df)
    
    obe_xgb_10cv.to_csv("obe_XGB_without_CS_10CV.csv", index=False)
    obe_xgb_10cv_best.to_csv("obe_XGB_without_CS_best_parameter_10CV.csv", index=False)
    
    obe_xgb_withsty = train_within(obe_indiv_df, param_xgb, pvalname)
    obe_xgb_withsty.to_csv("obe_XGB_without_CS_within.csv", index=False)
    
    obe_xgb_loosty = train_loov(obe_indiv_df, param_xgb, pvalname)
    obe_xgb_loosty.to_csv("obe_XGB_without_CS_LOOV.csv", index=False)
    
    
    obe_10cv_summ = xgb_10cv_grid_search_covariate_shift(obe_indiv_df, param_xgb)
    traindt = obe_indiv_df.copy()
    ytarget = obe_indiv_df["group"]
    traindt.drop("group", axis=1, inplace=True)
    traindt.drop("StudyID", axis=1, inplace=True)
    obe_10CV_cov_shift = CV_covariate_shift(obe_10cv_summ, traindt, ytarget,10)
    obe_10CV_cov_shift.to_csv("obe_XGB_with_CS_10CV.csv", index=False)
    
    obe_grid_cov_shift_within = train_within_covariate_shift(obe_indiv_df, param_xgb, pvalname) 
    obe_within_cov_shift = within_study_covariate_best(obe_indiv_df, obe_grid_cov_shift_within, pvalname)
    obe_within_cov_shift.to_csv("obe_XGB_with_CS_within.csv", index=False)
    
    obe_grid_cov_shift_loov = train_loov_covariate_shift(obe_indiv_df, param_xgb, pvalname)
    obe_loov_cov_shift = loov_covariate_shift_best(obe_indiv_df, obe_grid_cov_shift_loov, pvalname)
    obe_loov_cov_shift.to_csv("obe_XGB_with_CS_loov.csv", index=False)


if __name__ == '__main__':
    main()
    
