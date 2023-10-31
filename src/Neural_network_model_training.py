###################### train the NN model ###################################
import pandas as pd
import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn import model_selection
from utils import train_data, obese
from utils import expand_grid
from utils import crossvalidate
from utils import median
from utils import permutation_test_score_eachfold
from utils import validate
from utils import cross_validate_with_permutation_test_score


# grid search for the best combination of parameters

def mlp_10crossvalidation(all_indiv_df, param_mlp):
    mlp_10cv = pd.DataFrame(
        columns=['hidden_layer_size1', "hidden_layer_size2", 'activation', 'solver', 'alpha', 'learning_rate', 'auc'])
    xtrain = all_indiv_df.iloc[:, :len(all_indiv_df.columns)-2]
    ytrain = all_indiv_df["group"]
    for i in range(len(param_mlp)):
        hiddenlayer = param_mlp.loc[i, "hidden_layer_sizes"]
        activate_fun = param_mlp.loc[i, "activation"]
        n_solver = param_mlp.loc[i, "solver"]
        alpha_para = param_mlp.loc[i, "alpha"]
        rate = param_mlp.loc[i, "learning_rate"]
        mlp_gs = MLPClassifier(hidden_layer_sizes = hiddenlayer, activation = activate_fun,
                               learning_rate = rate, alpha = alpha_para, solver = n_solver, 
                               max_iter = 100, random_state = 2021)
        aucscore = crossvalidate(mlp_gs, xtrain, ytrain, 10)
        if len(hiddenlayer) == 2:
            df = pd.DataFrame(
                {"hidden_layer_size1": hiddenlayer[0], "hidden_layer_size2": hiddenlayer[1], "activation": activate_fun,
                 "solver": n_solver, "alpha": alpha_para,
                 "learning_rate": rate, "auc": median(aucscore)}, index=[0])
        else:
            df = pd.DataFrame(
                {"hidden_layer_size1": hiddenlayer[0], "hidden_layer_size2": 0, "activation": activate_fun,
                 "solver": n_solver, "alpha": alpha_para,
                 "learning_rate": rate, "auc": median(aucscore)}, index=[0])
        mlp_10cv = pd.concat([mlp_10cv, df])
    return mlp_10cv


def MLP_bestpara_10cv(param_mlp, all_indiv_df):
    xtrain = all_indiv_df.iloc[:, :len(all_indiv_df.columns)-2]
    ytrain = all_indiv_df["group"]
    best_para = param_mlp[param_mlp.auc == param_mlp.auc.max()]
    best_para = best_para.reset_index()
    hidden_layer_size1 = best_para.loc[0, "hidden_layer_size1"]
    hidden_layer_size2 = best_para.loc[0, "hidden_layer_size2"]
    if hidden_layer_size2 != 0:
        layer_size = (hidden_layer_size1, hidden_layer_size2)
    else:
        layer_size = (hidden_layer_size1,)
    activate_fun = best_para.loc[0, "activation"]
    n_solver = best_para.loc[0, "solver"]
    alpha_para = best_para.loc[0, "alpha"]
    rate = best_para.loc[0, "learning_rate"]
    mlp_gs = MLPClassifier(hidden_layer_sizes = layer_size, activation = activate_fun,
                           learning_rate = rate, alpha = alpha_para, solver = n_solver, 
                           max_iter = 100, random_state = 2021)
    aucscore, Pval = cross_validate_with_permutation_test_score(mlp_gs, xtrain, ytrain, 10)
    mlp_10cv_best = pd.DataFrame(
        {"Model": "MLP", "fold": [i for i in range(1, 11, 1)], "AUC": aucscore, "Pval": Pval, "Train": "10-CV",
         "Study": "Combined"})
    return mlp_10cv_best


def gridsearch_MLP(param, X_train, X_test, y_train, y_test):
    mlp_within = pd.DataFrame(
        columns=['hidden_layer_size1', 'hidden_layer_size2', 'activation', 
                 'solver', 'alpha', 'learning_rate', 'auc'])
    for i in range(len(param)):
        # print(i)
        hiddenlayer = param.loc[i, "hidden_layer_sizes"]
        activate_fun = param.loc[i, "activation"]
        n_solver = param.loc[i, "solver"]
        alpha_para = param.loc[i, "alpha"]
        rate = param.loc[i, "learning_rate"]
        mlp_gs = MLPClassifier(hidden_layer_sizes = hiddenlayer, activation = activate_fun,
                               learning_rate = rate, alpha = alpha_para, solver = n_solver, 
                               max_iter = 100, random_state = 2021)
        auroc1 = validate(mlp_gs, np.array(X_train), np.array(X_test), y_train, y_test)

        if len(hiddenlayer) == 2:
            df = pd.DataFrame(
                {"hidden_layer_size1": hiddenlayer[0], "hidden_layer_size2": hiddenlayer[1], "activation": activate_fun,
                 "solver": n_solver, "alpha": alpha_para,
                 "learning_rate": rate, "auc": auroc1}, index=[0])
        else:
            df = pd.DataFrame(
                {"hidden_layer_size1": hiddenlayer[0], "hidden_layer_size2": 0, "activation": activate_fun,
                 "solver": n_solver, "alpha": alpha_para,
                 "learning_rate": rate, "auc": auroc1}, index=[0])
        mlp_within = pd.concat([mlp_within, df])
    mlp_within = mlp_within.reset_index()
    mlp_within.drop('index', axis=1, inplace=True)
    return mlp_within


def train_within(all_indiv_df, param, pvalname):
    alist = all_indiv_df["StudyID"].unique()
    MLP_withsty = pd.DataFrame(
        columns=['hidden_layer_size1', 'hidden_layer_size2', 'activation', 
                 'solver', 'alpha', 'learning_rate', 'auc', 'StudyID', 'Pval'])
    for s in alist:
        print(s)
        train = all_indiv_df[all_indiv_df["StudyID"] == s]
        train = train.reset_index()
        train.drop("index", axis=1, inplace=True)
        X = train[pvalname]
        Y = train["group"]
        X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Y, test_size=0.25, random_state=2021)
        MLP_within = gridsearch_MLP(param, X_train, X_test, y_train, y_test)
        maxline = MLP_within[MLP_within.auc == MLP_within.auc.max()]
        maxline = maxline.reset_index()
        maxline.drop("index", axis=1, inplace=True)
        maxline["StudyID"] = s
        maxline = maxline.loc[:0, :]
        hidden_layer_size1 = maxline.loc[0, "hidden_layer_size1"]
        hidden_layer_size2 = maxline.loc[0, "hidden_layer_size2"]
        if hidden_layer_size2 != 0:
            layer_size = (hidden_layer_size1, hidden_layer_size2)
        else:
            layer_size = (hidden_layer_size1,)
        activate_fun = maxline.loc[0, "activation"]
        n_solver = maxline.loc[0, "solver"]
        alpha_para = maxline.loc[0, "alpha"]
        rate = maxline.loc[0, "learning_rate"]
        trainix = np.array(X_train.index)
        testix = np.array(X_test.index)
        X = np.array(train[pvalname])
        Y = np.array(train["group"])
        clf = MLPClassifier(hidden_layer_sizes=layer_size, activation=activate_fun,
                            learning_rate=rate, alpha=alpha_para, solver=n_solver, 
                            max_iter=100, random_state=2021)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_eachfold(
            clf, X, Y, scoring="roc_auc", cv=2, train=trainix, test=testix, 
            n_permutations=100, n_jobs=-1)
        maxline["Pval"] = pvalue_train
        MLP_withsty = pd.concat([MLP_withsty, maxline])
    return MLP_withsty


def train_loov(all_indiv_df, param, pvalname):
    alist = all_indiv_df["StudyID"].unique()
    MLP_loosty = pd.DataFrame(
        columns=['hidden_layer_size1', 'hidden_layer_size2','activation', 
                 'solver', 'alpha', 'learning_rate', 'auc', 'StudyID', 'Pval'])
    for s in alist:
        print(s)
        train = all_indiv_df[all_indiv_df["StudyID"] != s]
        test = all_indiv_df[all_indiv_df["StudyID"] == s]
        xtrain = train[pvalname]
        ytrain = train["group"]
        xtest = test[pvalname]
        ytest = test["group"]
        MLP_within = gridsearch_MLP(param, xtrain, xtest, ytrain, ytest)
        maxline = MLP_within[MLP_within.auc == MLP_within.auc.max()]
        maxline = maxline.reset_index()
        maxline.drop("index", axis=1, inplace=True)
        maxline["StudyID"] = s
        maxline = maxline.loc[:0, :]
        trainix = np.array(train.index)
        testix = np.array(test.index)
        X = np.array(all_indiv_df[pvalname])
        Y = np.array(all_indiv_df["group"])
        hidden_layer_size1 = maxline.loc[0, "hidden_layer_size1"]
        hidden_layer_size2 = maxline.loc[0, "hidden_layer_size2"]
        if hidden_layer_size2 != 0:
            layer_size = (hidden_layer_size1, hidden_layer_size2)
        else:
            layer_size = (hidden_layer_size1,)
        activate_fun = maxline.loc[0, "activation"]
        n_solver = maxline.loc[0, "solver"]
        alpha_para = maxline.loc[0, "alpha"]
        rate = maxline.loc[0, "learning_rate"]
        clf = MLPClassifier(hidden_layer_sizes = layer_size, activation = activate_fun,
                            learning_rate = rate, alpha = alpha_para, solver = n_solver, 
                            max_iter = 100, random_state = 2021)
        score_train, perm_scores_train, pvalue_train = permutation_test_score_eachfold(
            clf, X, Y, scoring = "roc_auc", cv = 2, train = trainix, test = testix, 
            n_permutations = 100, n_jobs = -1)
        maxline["Pval"] = pvalue_train
        MLP_loosty = pd.concat([MLP_loosty, maxline])
    return MLP_loosty


def main():
    meta = pd.read_csv("meta_BMI_obesity.csv")
    table = pd.read_csv("genus_combined_rel_rdp_filter.csv")
    pvalname = table["Genus"].str.replace("g__", "")

    all_indiv_df = train_data(table, meta, pvalname)
    all_indiv_df['group'] = all_indiv_df['group'].apply(lambda x: obese(x))

    parameter_space = {
        'hidden_layer_sizes': [(258, 86), (126,), (94,), (18,)],
        'activation': ['tanh', 'relu'],
        'solver': ['sgd', 'adam'],
        'alpha': [0.0001, 0.001, 0.01],
        'learning_rate': ['constant', 'adaptive'],
    }

    param_mlp = expand_grid(parameter_space)

    mlp_10cv = mlp_10crossvalidation(all_indiv_df, param_mlp)
    mlp_10cv.to_csv("all_NN_10cv.csv", index=False)
    
    MLP_10cv_best = MLP_bestpara_10cv(mlp_10cv, all_indiv_df)
    MLP_10cv_best.to_csv("all_NN_10CV_best.csv", index=False)

    #train MLP within each study

    MLP_withsty = train_within(all_indiv_df, param_mlp, pvalname)
    MLP_withsty.to_csv("all_NN_within.csv", index=False)

    #train MLP using LOOV

    MLP_loosty = train_loov(all_indiv_df, param_mlp, pvalname)
    MLP_loosty.to_csv("all_NN_LOOV.csv", index=False)


if __name__ == '__main__':
    main()
