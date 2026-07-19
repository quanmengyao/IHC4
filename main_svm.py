import numpy as np
import pandas as pd
import os
import joblib

from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVC
from sklearn.metrics import (
    roc_curve,
    auc,
    accuracy_score,
    confusion_matrix
)

from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline


RANDOM_STATE = 42


def Index_AUC_ACC_SEN_SPE_CM(
    y_true,
    y_score,
    cutpoint=0.5,
    pos_label=1,
    sample_weight=None,
    drop_intermediate=False
):
    fpr, tpr, thresholds = roc_curve(
        y_true,
        y_score,
        pos_label=pos_label,
        sample_weight=sample_weight,
        drop_intermediate=drop_intermediate
    )

    Auc = auc(fpr, tpr)

    y_class = np.where(y_score >= cutpoint, 1, 0)

    Acc = accuracy_score(y_true, y_class)

    tn, fp, fn, tp = confusion_matrix(
        y_true,
        y_class,
        labels=[0, 1]
    ).ravel()

    Sen = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    Spe = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    YI = Sen + Spe - 1

    return (
        Auc,
        Acc,
        Sen,
        Spe,
        YI,
        (tn, fp, fn, tp),
        y_class
    )


def get_youden_cutpoint(y_true, y_score):
    fpr, tpr, thresholds = roc_curve(
        y_true,
        y_score,
        pos_label=1,
        drop_intermediate=False
    )

    youden_index = tpr - fpr
    idx_max = np.argmax(youden_index)

    return thresholds[idx_max]


feature0 = "feature_0.xlsx"
feature1 = "feature_1.xlsx"

col_list = np.arange(8)

data0 = pd.read_excel(
    feature0,
    usecols=col_list
)

data1 = pd.read_excel(
    feature1,
    usecols=col_list
)

X0 = data0.to_numpy(dtype=float)
X1 = data1.to_numpy(dtype=float)

label0 = np.zeros(
    len(data0),
    dtype=int
)

label1 = np.ones(
    len(data1),
    dtype=int
)

X_total = np.concatenate(
    (X0, X1),
    axis=0
)

Y_total = np.concatenate(
    (label0, label1),
    axis=0
)


def SVM(X_train, y_train, X_test):

    pipeline = Pipeline(
        steps=[
            (
                "scaler",
                MinMaxScaler()
            ),
            (
                "smote",
                SMOTE(
                    random_state=RANDOM_STATE
                )
            ),
            (
                "svc",
                SVC(
                    probability=True,
                    random_state=RANDOM_STATE
                )
            )
        ]
    )

    parameters = [
        {
            "svc__kernel": ["rbf"],
            "svc__C": [
                2 ** x for x in range(-10, 11)
            ],
            "svc__gamma": [
                2 ** x for x in range(-10, 11)
            ]
        },
        {
            "svc__kernel": ["linear"],
            "svc__C": [
                2 ** x for x in range(-10, 11)
            ]
        }
    ]

    inner_cv = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=RANDOM_STATE
    )

    clf = GridSearchCV(
        estimator=pipeline,
        param_grid=parameters,
        scoring="roc_auc",
        cv=inner_cv,
        n_jobs=-1,
        refit=True,
        return_train_score=False
    )

    clf.fit(
        X_train,
        y_train
    )

    best_model = clf.best_estimator_

    score_test = best_model.predict_proba(
        X_test
    )[:, 1]

    score_train = best_model.predict_proba(
        X_train
    )[:, 1]

    return (
        best_model,
        score_test,
        score_train,
        clf.best_params_,
        clf.best_score_
    )


skf = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=RANDOM_STATE
)

i = 0

for train_index, test_index in skf.split(
    X_total,
    Y_total
):
    i += 1

    Xtrain = X_total[train_index]
    Ytrain = Y_total[train_index]

    Xtest = X_total[test_index]
    Ytest = Y_total[test_index]

    (
        model,
        y_score,
        y_score_train,
        best_params,
        best_cv_auc
    ) = SVM(
        Xtrain,
        Ytrain,
        Xtest
    )

    train_cutpoint = get_youden_cutpoint(
        Ytrain,
        y_score_train
    )

    res = Index_AUC_ACC_SEN_SPE_CM(
        Ytest,
        y_score,
        cutpoint=train_cutpoint
    )

    print(
        "----------- Fold",
        i,
        "results ------------"
    )

    print("Best parameters :", best_params)
    print("Best inner CV AUC:", best_cv_auc)
    print("Training cutpoint:", train_cutpoint)

    print("AUC : ", res[0])
    print("ACC : ", res[1])
    print("SEN : ", res[2])
    print("SPC : ", res[3])
    print("YI  : ", res[4])
    print("CM  : ", res[5])

    os.makedirs(
        "saved_models",
        exist_ok=True
    )

    joblib.dump(
        model,
        f"saved_models/svm_fold_{i}.pkl"
    )
