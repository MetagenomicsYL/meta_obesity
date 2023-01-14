####################### train the RF model ###################################
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import model_selection
from utils import train_data, obese
from utils import expand_grid
from utils import crossvalidate
from utils import median
from utils import permutation_test_score_eachfold
from utils import validate
from utils import cross_validate_with_permutation_test_score


# grid search for the best combination of parameters

def RF_10crossvalidation(all_indiv_df, param_xgb):
    rf_10cv = pd.DataFrame(columns=['n_estimators', 'max_depth', 'max_features', 'auc'])
    xtrain = all_indiv_df.iloc[:, :len(all_indiv_df.columns) - 2]
    ytrain = all_indiv_df["group"]
    for i in range(len(param_xgb)):
        estimators = param_xgb.loc[i, "n_estimators"]
        depth = param_xgb.loc[i, "max_depth"]
        maxfeature = param_xgb.loc[i, "max_features"]
        model = RandomForestClassifier(n_estimators = estimators, max_depth = depth, 
                                       max_features = maxfeature, n_jobs = -1,
                                       random_state = 2021)
        aucscore = crossvalidate(model, xtrain, ytrain, 10)
        df = pd.DataFrame(
            {"n_estimators": estimators, "max_depth": depth, "max_features": maxfeature, 
             "auc": median(aucscore)}, index=[0])
        rf_10cv = pd.concat([rf_10cv, df])
    return rf_10cv


def model_bestpara_10cv(xgb_10cv, all_indiv_df):
    xtrain = all_indiv_df.iloc[:, :len(all_indiv_df.columns) - 2]
    ytrain = all_indiv_df["group"]
    best_para = xgb_10cv[xgb_10cv.auc == xgb_10cv.auc.max()]
    best_para = best_para.reset_index()
    estimators = best_para.loc[0, "n_estimators"]
    depth = best_para.loc[0, "max_depth"]
    maxfeature = best_para.loc[0, "max_features"]
    model = RandomForestClassifier(n_estimators=estimators, max_depth=depth, 
                                   max_features=maxfeature, n_jobs=-1,
                                   random_state=2021)
    aucscore, Pval = cross_validate_with_permutation_test_score(model, xtrain, 
                                                               ytrain, 10)
    rf_10cv_best = pd.DataFrame(
        {"Model": "RF", "fold": [i for i in range(1, 11, 1)], 
         "AUC": aucscore, "Pval": Pval,  "Train": "10-CV", 
         "Study": "Combined"})
    return rf_10cv_best


def gridsearch_rf(param, X_train, X_test, y_train, y_test):
    rf_within = pd.DataFrame(columns=['n_estimators', 'max_depth', 
                                      'max_features', 'auc'])
    for i in range(len(param)):
        estimators = param.loc[i, "n_estimators"]
        depth = param.loc[i, "max_depth"]
        maxfeature = param.loc[i, "max_features"]
        model2 = RandomForestClassifier(n_estimators = estimators, max_depth = depth,
                                        max_features = maxfeature, n_jobs = -1, 
                                        random_state = 2021)
        auroc = validate(model2, np.array(X_train), np.array(X_test), y_train, y_test)
        df = pd.DataFrame({"n_estimators": estimators, "max_depth": depth, 
                           "max_features": maxfeature, "auc": auroc},
                          index=[0])
        rf_within = pd.concat([rf_within, df])
    rf_within = rf_within.reset_index()
    rf_within.drop('index', axis = 1, inplace = True)
    return rf_within

# train random forest within each study

def train_within(all_indiv_df, param_rf, pvalname):
    alist = all_indiv_df["StudyID"].unique()
    rf_withsty = pd.DataFrame(columns=['n_estimators', 'max_depth', 'max_features', 
                                       'auc', 'StudyID', 'Pval'])
    for s in alist:
        # print(s)
        train = all_indiv_df[all_indiv_df["StudyID"] == s]
        train = train.reset_index()
        train.drop("index", axis=1, inplace=True)
        X = train[pvalname]
        Y = train["group"]
        X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Y, test_size=0.25, random_state=2021)
        rf_within = gridsearch_rf(param_rf, X_train, X_test, y_train, y_test)
        maxline = rf_within[rf_within.auc == rf_within.auc.max()]
        maxline = maxline.reset_index()
        maxline.drop("index", axis=1, inplace=True)
        maxline["StudyID"] = s
        maxline = maxline.loc[:0, :]
        estimators = maxline.loc[0, "n_estimators"]
        depth = maxline.loc[0, "max_depth"]
        maxfeature = maxline.loc[0, "max_features"]
        trainix = np.array(X_train.index)
        testix = np.array(X_test.index)
        X = np.array(train[pvalname])
        Y = np.array(train["group"])
        clf = RandomForestClassifier(n_estimators=estimators, max_depth=depth,
                                     max_features=maxfeature, n_jobs=-1, random_state=2021)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_eachfold(
            clf, X, Y, scoring="roc_auc", cv=2, train=trainix, test=testix, 
            n_permutations=100, n_jobs = -1)
        maxline["Pval"] = pvalue_train
        rf_withsty = pd.concat([rf_withsty, maxline])
    return rf_withsty


# train random forest using LOOV

def train_loov(all_indiv_df, param_rf, pvalname):
    alist = all_indiv_df["StudyID"].unique()
    rf_loosty = pd.DataFrame(columns=['n_estimators', 'max_depth', 'max_features', 'auc', 'StudyID', 'Pval'])
    for s in alist:
        print(s)
        train = all_indiv_df[all_indiv_df["StudyID"] != s]
        test = all_indiv_df[all_indiv_df["StudyID"] == s]
        xtrain = train[pvalname]
        ytrain = train["group"]
        xtest = test[pvalname]
        ytest = test["group"]
        rf_within = gridsearch_rf(param_rf, xtrain, xtest, ytrain, ytest)
        maxline = rf_within[rf_within.auc == rf_within.auc.max()]
        maxline = maxline.reset_index()
        maxline.drop("index", axis=1, inplace=True)
        maxline["StudyID"] = s
        maxline = maxline.loc[:0, :]
        trainix = np.array(train.index)
        testix = np.array(test.index)
        X = np.array(all_indiv_df[pvalname])
        Y = np.array(all_indiv_df["group"])
        estimators = maxline.loc[0, "n_estimators"]
        depth = maxline.loc[0, "max_depth"]
        maxfeature = maxline.loc[0, "max_features"]
        clf = RandomForestClassifier(n_estimators=estimators, max_depth=depth,
                                     max_features=maxfeature, n_jobs=-1, random_state=2021)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_eachfold(
            clf, X, Y, scoring="roc_auc", cv=2, train=trainix, test=testix, 
            n_permutations=100, n_jobs = -1)
        maxline["Pval"] = pvalue_train
        rf_loosty = pd.concat([rf_loosty, maxline])
    return rf_loosty


def main():
    meta = pd.read_csv("meta_BMI_obesity.csv")
    table = pd.read_csv("genus_combined_rel_rdp_filter.csv")
    pvalname = table["Genus"].str.replace("g__", "")
    
    all_indiv_df = train_data(table, meta, pvalname)
    all_indiv_df['group'] = all_indiv_df['group'].apply(lambda x: obese(x))

    params = {'max_depth': range(3, 11, 1), 'n_estimators': range(100, 1100, 100),
              "max_features": {"sqrt", "log2", None}}
    params_rf = expand_grid(params)
    
    all_rf_10cv = RF_10crossvalidation(all_indiv_df, params_rf)
    all_rf_10cv_best = model_bestpara_10cv(all_rf_10cv, all_indiv_df)

    all_rf_10cv.to_csv("all_RF_10CV.csv", index=False)
    all_rf_10cv_best.to_csv("all_RF_10CV_best.csv", index=False)

    rf_withsty = train_within(all_indiv_df, params_rf, pvalname)  
    rf_withsty.to_csv("all_RF_within.csv", index=False)

    rf_loosty = train_loov(all_indiv_df, params_rf, pvalname) 
    rf_loosty.to_csv("all_RF_LOOV.csv", index=False)


if __name__ == '__main__':
    main()
