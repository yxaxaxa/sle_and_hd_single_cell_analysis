#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
import numpy as np
import glob,os
from glob import iglob
#import scanpy as sc
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import RocCurveDisplay
from sklearn.datasets import load_wine
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
import joblib
import time
import random
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42


# # RA PBMC data for machine learning

# In[6]:


### training data import
ra=pd.read_csv('../RNA_seq_for_autoimmune_disease/RA_bulk/GSE90081/GSE90081_ra_part.csv',index_col=0)
hd=pd.read_csv('../RNA_seq_for_autoimmune_disease/RA_bulk/GSE90081/GSE90081_hd_part.csv',index_col=0)
hd1=pd.read_csv('../RNA_seq_for_autoimmune_disease/health_bulk/GSE183204_HC_fpkm.csv',sep=',',index_col=0)


# In[7]:


### feature import
features=pd.read_csv('../script4paper2/combined_gene_for_machine_learning.csv',index_col=1).index.values
features=np.append(features,'patient')

features=[i for i in features if i in ra.index.values]
features=[i for i in features if i in hd1.index.values ]


# # remove unwanted gene

# In[8]:


### remove unwanted gene from validation data

hd1=hd1.loc[features,:].T
ra_part=ra.loc[features,:].T
hd_part=hd.loc[features,:].T


# # label data

# In[9]:


### label training data
ra_part['patient']=1
hd_part['patient']=0
hd1['patient']=0


# # machine learning data training

# In[39]:


### merge training data
df=pd.concat([ra_part,hd_part,hd1],axis=0)

### get data labels
label=df.patient.values

### split data with ratio 30% for test and 70% for training
Xtrain, Xtest, Ytrain, Ytest = train_test_split(df.drop(columns=['patient']),label,test_size=0.3)

### rf model initialization
rfc = RandomForestClassifier(random_state=43,class_weight='balanced',oob_score=True)
rfc = rfc.fit(Xtrain,Ytrain)

### document model score
score_r = rfc.score(Xtest,Ytest)

### save feature importance
ra_pbmc=pd.DataFrame(rfc.feature_importances_)
ra_pbmc['feature_importance']=features
ra_pbmc.to_csv('./model/ra_pbmc_feature_importance_bulk.csv')

### print F score and Out of bag score
print("Random Forest:{}".format(score_r))
print("OOB score:",rfc.oob_score_)





# # Figure 7A

# In[40]:


### Generating ROC curve
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
rfc_disp = RocCurveDisplay.from_estimator(rfc, Xtest, Ytest, ax=ax, alpha=0.8)
plt.legend(loc=4,prop={'size': 10})
plt.xlabel('False Positive Rate', fontsize=18)
plt.ylabel('True Positive Rate', fontsize=16)
ax.plot([0, 1], [0, 1], ls="--", c=".3")
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
plt.savefig('./figure6_and_7/7a_ra_pbmc_bulk_auc.pdf',width=4,height=5)


# # save/load best performance model

# In[24]:


### save the best performance model
#joblib.dump(rfc, './model/ra_synovial_bulk_best.model')

### load model
#rfc=joblib.load('./model/sle_best.model')


# In[19]:


### 10-fold cross validation
print(cross_val_score(rfc,df.drop(columns=['patient']),label,cv=10).mean())
print(cross_val_score(rfc,df.drop(columns=['patient']),label,cv=10).var())


# # Figure 7D

# In[42]:


ra_feature=pd.read_csv('./model/ra_pbmc_feature_importance_bulk.csv')
fig, ax = plt.subplots(figsize=(15, 5))
ax.bar(x=ra_feature['feature_importance'], height=ra_feature['0'])
ax.set_title("Feature importance for RA bulk RNA PBMC model", fontsize=15)
plt.xticks(rotation = 90)
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
plt.savefig('./figure6_and_7/7d_ra_pbmc_bulk.pdf',width=15,height=5)


# # Hyper-parameter adjust

# In[795]:


data=df.drop(columns=['patient'])
label=df.patient.values
start=time.time()
scorel = []
for i in range(0,200,10): # loop for 0-200 decision trees
    rfc = RandomForestClassifier(n_estimators=i+1,n_jobs=-1,random_state=0)
    score = cross_val_score(rfc,data,label,cv=10).mean()
    scorel.append(score)
print(max(scorel),(scorel.index(max(scorel))*10)+1)
end=time.time()
print('Running time: %s Seconds'%(end-start))
plt.figure(figsize=[20,5])
plt.plot(range(1,201,10),scorel)
plt.show()


# In[801]:


scorel = []
for i in range(185,205):
  rfc = RandomForestClassifier(n_estimators=i+1,n_jobs=-1,random_state=0)
  score = cross_val_score(rfc,data,label,cv=10).mean()
  scorel.append(score)
print(max(scorel),([*range(185,205)][scorel.index(max(scorel))]))
plt.figure(figsize=[20,5])
plt.plot(range(185,205),scorel)
plt.show()


# In[802]:


start=time.time()
param_grid = {'max_depth':np.arange(1, 90,2)} 
alg = RandomForestClassifier(n_estimators=190,random_state=0) 
GS = GridSearchCV(alg,param_grid,cv=10) 
GS.fit(data,label)
print(GS.best_params_)
print(GS.best_score_)
end=time.time()
print('Running time: %s Seconds'%(end-start))


# In[803]:


start=time.time()
param_grid = {'max_features':np.arange(5,80,1)}
rfc = RandomForestClassifier(n_estimators=190,random_state=0)
GS = GridSearchCV(rfc,param_grid,cv=10)
GS.fit(data,label)
print(GS.best_params_)
print(GS.best_score_)
end=time.time()
print('Running time: %s Seconds'%(end-start))


# # 100 loop of 10-fold cross validation

# In[35]:


df_n=df.drop(columns=['patient'])
rfc_l = []
fpr_l=[]
tpr_l=[]
acc_l=[]
skf =StratifiedKFold(n_splits=10)
for i in range(100):
    for train_index, test_index in skf.split(df_n,label):
        rfc = RandomForestClassifier(random_state=0,class_weight="balanced",oob_score=True)
        rfc = rfc.fit(df_n.iloc[train_index],label[train_index])
        rfc_l.append(roc_auc_score(label[test_index], rfc.predict_proba(df_n.iloc[test_index])[:, 1]))
        acc_l.append(accuracy_score(label[test_index], rfc.predict(df_n.iloc[test_index])))
    
    


# In[36]:


### average AUC and its standard deviation error
print(np.mean(rfc_l))
print(np.std(rfc_l))

