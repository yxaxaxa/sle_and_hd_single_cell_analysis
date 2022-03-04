#!/usr/bin/env python
# coding: utf-8

# In[11]:


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
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
import os
import matplotlib as mpl
#os.environ["KMP_DUPLICATE_LIB_OK"]  =  "FALSE"
import time
import random
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42


# # single cell sle part

# In[2]:


features=pd.read_csv('./combined_gene_for_machine_learning.csv',index_col=1)
features=np.append(features.index.values,'patient')
features=np.delete(features,[3,7,16,17,18,76,78,79])


# In[3]:


path = '../GSE135779_SLE/test/aHD/'
file = glob.glob(os.path.join(path, "*.csv.gz"))
hd = []
for f in file:
    hd.append(pd.read_csv(f,index_col=0).loc[features[1:80],:].T)
for i in range(len(hd)):
    hd[i]['patient']=0


# In[4]:


path = '../GSE135779_SLE/test/cHD/'
file = glob.glob(os.path.join(path, "*.csv.gz"))
chd = []
for f in file:
    chd.append(pd.read_csv(f,index_col=0).loc[features[1:80],:].T)
for i in range(len(chd)):
    chd[i]['patient']=0


# In[5]:



hd_m=hd[0]
for i in range(1,len(hd)):
    hd_m=pd.concat([hd_m,hd[i]],axis=0)


# In[6]:



chd_m=chd[0]
for i in range(1,len(chd)):
    chd_m=pd.concat([chd_m,chd[i]],axis=0)


# In[7]:


path = '../GSE135779_SLE/test/aSLE/'
file = glob.glob(os.path.join(path, "*.csv.gz"))
asle = []
for f in file:
    asle.append(pd.read_csv(f,index_col=0).loc[features[1:80],:].T)
for i in range(len(asle)):
    asle[i]['patient']=1
asle_m=asle[0]
for i in range(1,len(asle)):
    asle_m=pd.concat([asle_m,asle[i]],axis=0)


# In[8]:


path = '../GSE135779_SLE/test/cSLE/'
file = glob.glob(os.path.join(path, "*.csv.gz"))
csle = []
for f in file:
    csle.append(pd.read_csv(f,index_col=0).loc[features[1:80],:].T)
for i in range(len(csle)):
    csle[i]['patient']=1
csle_m=csle[0]
for i in range(1,len(csle)):
    csle_m=pd.concat([csle_m,csle[i]],axis=0)


# In[9]:


df=pd.concat([hd_m,chd_m,asle_m,csle_m],axis=0)


# In[20]:


scorel = []
for i in range(185,205):
  rfc = RandomForestClassifier(n_estimators=i+1,n_jobs=-1,random_state=0)
  score = cross_val_score(rfc,data,label,cv=10).mean()
  scorel.append(score)
print(max(scorel),([*range(185,200)][scorel.index(max(scorel))]))
plt.figure(figsize=[20,5])
plt.plot(range(185,205),scorel)
plt.show()


# In[19]:


data=df.drop(columns=['patient'])
start=time.time()
scorel = []
for i in range(0,200,10): # 迭代建立包含0-200棵决策树的RF模型进行对比
  rfc = RandomForestClassifier(n_estimators=i+1,n_jobs=-1,random_state=0)
  score = cross_val_score(rfc,data,label,cv=10).mean()
  scorel.append(score)
print(max(scorel),(scorel.index(max(scorel))*10)+1)
end=time.time()
print('Running time: %s Seconds'%(end-start))
plt.figure(figsize=[20,5])
plt.plot(range(1,201,10),scorel)
plt.show()


# In[14]:


label=df.patient.values
Xtrain, Xtest, Ytrain, Ytest = train_test_split(df.drop(columns=['patient']),label,test_size=0.3)
rfc = RandomForestClassifier(random_state=0,class_weight='balanced',n_estimators=199)
rfc = rfc.fit(Xtrain,Ytrain)
score_r = rfc.score(Xtest,Ytest)
c=pd.DataFrame(rfc.feature_importances_)
#a.index=df.columns.values
print("Random Forest:{}".format(score_r))
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
rfc_disp = RocCurveDisplay.from_estimator(rfc, Xtest, Ytest, ax=ax, alpha=0.8)
plt.legend(loc=4,prop={'size': 10})
plt.xlabel('False Positive Rate', fontsize=18)
plt.ylabel('True Positive Rate', fontsize=16)
ax.plot([0, 1], [0, 1], ls="--", c=".3")
#plt.savefig('./figure6_and_7/auc_result_of_sc_sle_model.pdf',width=4,height=5)


# # reproduction of Figure6A

# In[16]:


pre_auc=pd.DataFrame(skf.split(df_n,label))
test_index=pre_auc.iloc[1,1]
train_index=pre_auc.iloc[1,0]
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
rfc = rfc.fit(df_n.iloc[train_index],label[train_index])
rfc_disp = RocCurveDisplay.from_estimator(rfc,df_n.iloc[test_index] ,label[test_index],
 ax=ax, alpha=0.8)
plt.legend(loc=4,prop={'size': 10})
plt.xlabel('False Positive Rate', fontsize=18)
plt.ylabel('True Positive Rate', fontsize=16)
ax.plot([0, 1], [0, 1], ls="--", c=".3")
#plt.savefig('./figure6_and_7/auc_result_of_sc_sle_model.pdf',width=4,height=5)


# In[ ]:


df_n=df.drop(columns=['patient'])
rfc_l = []
fpr_l=[]
tpr_l=[]
acc_l=[]
skf =StratifiedKFold(n_splits=10)
#for i in range(1000):
for train_index, test_index in skf.split(df_n,label):
    rfc = RandomForestClassifier(random_state=0,class_weight="balanced",oob_score=True)
    rfc = rfc.fit(df_n.iloc[train_index],label[train_index])
    rfc_l.append(roc_auc_score(label[test_index], rfc.predict_proba(df_n.iloc[test_index])[:, 1]))
        #fpr_l.append(metrics.roc_curve(label[test_index],rfc.predict(df_n.iloc[test_index]), pos_label=1))
    acc_l.append(accuracy_score(label[test_index], rfc.predict(df_n.iloc[test_index])))
    
print(np.mean(acc_l))
print(np.std(acc_l))    


# In[ ]:


print(cross_val_score(rfc,df.drop(columns=['patient']),label,cv=10).mean())
print(cross_val_score(rfc,df.drop(columns=['patient']),label,cv=10).var())
print(np.mean(rfc_l))
print(np.std(rfc_l))


# In[13]:


print(np.mean(rfc_l))
print(np.std(rfc_l))  


# In[43]:


c['feature_importance']=features[1:80]
fig, ax = plt.subplots(figsize=(15, 5))
ax.bar(x=c['feature_importance'], height=c[0])
ax.set_title("Feature importance for SLE single cell model", fontsize=15)
plt.xticks(rotation = 90)
plt.savefig('./figure6_and_7/sc_model_sle.pdf',width=15,height=5)

