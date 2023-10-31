"""

# some codes referred permutation_test_score function from the sklearn package

"""


from sklearn.base import is_classifier, clone
from sklearn.utils import indexable, check_random_state, _safe_indexing
from sklearn.model_selection import check_cv
from sklearn.metrics import check_scoring
from joblib import Parallel, logger
from sklearn.utils.fixes import delayed
from sklearn.utils.metaestimators import _safe_split
from sklearn.utils.validation import _check_fit_params
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn import metrics
from itertools import product


def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns = dictionary.keys())


def train_data(table, meta, pvalname):
    table["Genus"] = table["Genus"].str.replace("g__", "")
    genus = table.T
    genus = genus.reset_index()
    genus.columns = genus.iloc[0]
    genus = genus.drop(genus.index[0])
    genus.rename(columns={"Genus": "Run"}, inplace=True)
    outcome = meta[["Run", "group", "StudyID"]]
    obegenus = pd.merge(left=genus, right=outcome, how="inner", on="Run")
    colname = pvalname.to_list()
    colname = colname + ["StudyID", "group"]
    nostyid = obegenus[colname]
    return nostyid


def obese(x):
    if x.startswith("Obese"):
        return 1
    else:
        return 0


def median(data):
    data.sort()
    half = len(data) // 2
    return (data[half] + data[~half]) / 2


def _shuffle(y, groups, random_state):
    """Return a shuffled copy of y eventually shuffle among same groups."""
    if groups is None:
        indices = random_state.permutation(len(y))
    else:
        indices = np.arange(len(groups))
        for group in np.unique(groups):
            this_mask = groups == group
            indices[this_mask] = random_state.permutation(indices[this_mask])
    return _safe_indexing(y, indices)


def permutation_test_score_eachfold(
        estimator,
        X,
        y,
        *,
        groups = None,
        cv = None,
        train,
        test,
        n_permutations = 100,
        n_jobs = None,
        random_state = 0,
        verbose = 0,
        scoring = None,
        fit_params = None,
):
    X, y, groups = indexable(X, y, groups)

    cv = check_cv(cv, y, classifier=is_classifier(estimator))
    scorer = check_scoring(estimator, scoring=scoring)
    random_state = check_random_state(random_state)

    score = _permutation_test_score_eachfold(
        clone(estimator), X, y, groups, cv, train, test, scorer, fit_params=fit_params)
    permutation_scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_permutation_test_score_eachfold)(
            clone(estimator),
            X,
            _shuffle(y, groups, random_state),
            groups,
            cv,
            train,
            test,
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
    n_permutations = n_permutations - dupcount
    permutation_scores = np.array(permutation_scores)
    pvalue = (np.sum(permutation_scores >= score) + 1.0) / (n_permutations + 1)
    return score, permutation_scores, pvalue


def _permutation_test_score_eachfold(estimator, X, y, groups, cv, train, test, scorer, fit_params):
    fit_params = fit_params if fit_params is not None else {}
    avg_score = []
    # The train and test label comes from the data that have been splitted
    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, y_test = _safe_split(estimator, X, y, test, train)
    X_train = np.array(X_train)
    X_test = np.array(X_test)
    y_train = np.array(y_train)
    y_test = np.array(y_test)
    if len(pd.Series(y_test).unique()) != 2:
        return ("g")
    else:
        fit_params = _check_fit_params(X, fit_params, train)
        estimator.fit(X_train, y_train, **fit_params)
        avg_score.append(scorer(estimator, X_test, y_test))
        return np.mean(avg_score)


def crossvalidate(model, X_resampled, y_resampled, cv):
    auccore = []
    kfold1 = StratifiedKFold(n_splits=cv, shuffle=True, random_state=2021)
    for train, test in kfold1.split(X_resampled, y_resampled):
        X_train_xval = np.array(X_resampled)[train, :]
        X_test_xval = np.array(X_resampled)[test, :]
        y_train_xval = np.array(y_resampled)[train]
        y_test_xval = np.array(y_resampled)[test]
        model.fit(X_train_xval, y_train_xval)
        y_pred = model.predict(X_test_xval)
        auccore.append(metrics.roc_auc_score(y_test_xval, y_pred))
    return auccore


def validate(model, X_train, X_test, y_train, y_test):
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    auroc = metrics.roc_auc_score(y_test, y_pred)
    return auroc

def cross_validate_with_permutation_test_score(model, X_resampled, y_resampled, cv):
    aucscore = []
    kfold1 = StratifiedKFold(n_splits = cv, shuffle = True, random_state = 2021)
    Pval = []
    X = np.array(X_resampled)
    Y = np.array(y_resampled)
    for train, test in kfold1.split(X_resampled, y_resampled):
        X_train_xval = np.array(X_resampled)[train, :]
        X_test_xval = np.array(X_resampled)[test, :]
        y_train_xval = np.array(y_resampled)[train]
        y_test_xval = np.array(y_resampled)[test]
        score_train, perm_scores_train, pvalue_train = permutation_test_score_eachfold(
            model, X, Y, scoring="roc_auc", cv = 2, train = train, 
            test = test, n_permutations = 100, n_jobs = -1)
        Pval.append(pvalue_train)
        model.fit(X_train_xval, y_train_xval)
        y_pred = model.predict(X_test_xval)
        aucscore.append(metrics.roc_auc_score(y_test_xval, y_pred))
    return aucscore, Pval

