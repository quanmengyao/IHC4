import numpy as np
import pandas as pd
import os
import joblib

from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC

from imblearn.over_sampling import SMOTE
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix

def Index_AUC_ACC_SEN_SPE_CM(y_true, y_score, pos_label=1, sample_weight=None, drop_intermediate=False):
    fpr,tpr,thresholds = roc_curve(y_true, y_score,
                                    pos_label=pos_label,
                                    sample_weight=sample_weight,
                                    drop_intermediate=drop_intermediate)
    Auc = auc(fpr,tpr)
    YI = tpr - fpr
    idxMax = np.argmax(YI)
    cutpoint = thresholds[idxMax]
    Sen = tpr[idxMax]
    Spe = 1-fpr[idxMax]
    y_class = np.where(y_score >= cutpoint, 1, 0)
    Acc = accuracy_score(y_true,y_class)
    tn, fp, fn, tp = confusion_matrix(y_true,y_class).ravel()
    return Auc, Acc, Sen, Spe, np.max(YI), (tn, fp, fn, tp), y_class

RANDOM_STATE = 42

# 特征路径
feature0 = 'feature_0.xlsx'
feature1 = 'feature_1.xlsx'

# 读取特征
col_list = np.arange(8)
data0 = pd.read_excel(feature0,usecols=col_list)
data1 = pd.read_excel(feature1,usecols=col_list)

# 设置标签
label0 = np.zeros(len(data0),dtype=int)
label1 = np.ones(len(data1),dtype=int)

# 合并数据
X_total = np.concatenate((data0, data1), axis=0)
Y_total = np.concatenate((label0, label1), axis=0)


# 特征归一化
scaler = MinMaxScaler()
X_total = scaler.fit_transform(X_total)


# SVM分类
def SVM(X_train, y_train, X_test):
    parameters = [{'C': [2 ** x for x in range(-10, 11)],
                   'gamma': [2 ** x for x in range(-10, 11)],
                   'kernel': ['rbf', 'linear']}]
    clf = GridSearchCV(estimator=SVC(probability=True), param_grid=parameters, scoring = 'roc_auc', cv=5,
                       n_jobs=-1)
    clf.fit(X_train, y_train)
    best_model = clf.best_estimator_
    score_test = best_model.predict_proba(X_test)[:, 1]
    score_train = best_model.predict_proba(X_train)[:, 1]
    return best_model, score_test, score_train


# 交叉验证
skf = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=RANDOM_STATE
)


# SMOTE处理
smote = SMOTE(random_state=RANDOM_STATE)


# 训练循环
i = 0
for train_index, test_index in skf.split(X_total, Y_total):
    i += 1

    Xtrain = X_total[train_index]
    Ytrain = Y_total[train_index]

    Xtest = X_total[test_index]
    Ytest = Y_total[test_index]


    # SMOTE过采样
    Xtrain, Ytrain = smote.fit_resample(Xtrain, Ytrain)

    # 模型训练
    model, y_score, y_score_train = SVM(Xtrain, Ytrain, Xtest)

    # 计算指标
    res = Index_AUC_ACC_SEN_SPE_CM(Ytest, y_score)
    print('----------- Fold', i, 'results ------------')
    print('AUC : ', res[0])
    print('ACC : ', res[1])
    print('SEN : ', res[2])
    print('SPC : ', res[3])
    print('YI  : ', res[4])

    # 保存模型
    os.makedirs("saved_models", exist_ok=True)

    joblib.dump(
        model,
        f"saved_models/svm_fold_{i}.pkl"
    )