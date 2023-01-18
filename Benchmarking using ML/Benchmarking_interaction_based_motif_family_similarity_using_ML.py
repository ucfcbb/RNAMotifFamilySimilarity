#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import seaborn as sns
import statistics
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
from scipy.stats import norm, skew

from sklearn.metrics import classification_report, make_scorer, accuracy_score, precision_score, recall_score, f1_score
from sklearn.model_selection import cross_val_score, KFold, RepeatedKFold, train_test_split, GridSearchCV, cross_validate, RepeatedStratifiedKFold, cross_val_predict
from sklearn.preprocessing import PowerTransformer
from sklearn.naive_bayes import GaussianNB

pd.options.mode.chained_assignment = None


TEST_PERCEN = 0.2
RANDOM_STATE = 0
REPEATS = 5


# In[2]:


get_ipython().run_cell_magic('javascript', '', 'IPython.OutputArea.prototype._should_scroll = function(lines) {\n    return false;\n}')


# In[3]:


### Reading CSV File
def load_data():
    Data = pd.read_csv('Feature_sets/Features_ScanX.csv')
    pd.set_option('display.max_columns', None)
    return Data


# In[4]:


def labeling(Data):
    Data.columns = Data.columns.str.replace(' ', '')
    Data['motif_family'] = Data['motif_family'].str.replace(' ', '')   
    family = {'Tandem-shear': 1, 'E-loop': 2, 'Sarcin-ricin': 3, 'Kink-turn': 4, 'Hook-turn': 5, 'C-loop': 6, 'Rope-sling': 7, 'Tetraloop-receptor': 8, 'reverse-Kink-turn': 9, 'L1-complex': 10, 'T-loop': 11}    
    Data.motif_family = [family[item] for item in Data.motif_family]
    Data = Data.drop(['motif_str'], axis = 1)
    return Data


# In[5]:


def check_class_bias(Data):
    print('------------- Distribution of Class in the Label Column -------------')
    
    mapp = {1: 'TS', 2: 'EL', 3: 'SR', 4: 'KT', 5: 'HT', 6: 'CL', 7: 'RS', 8: 'TR', 9: 'rKT', 10: 'L1C', 11: 'TL'}
    Data['Family'] = Data["motif_family"].copy()
    Data.Family = [mapp[item] for item in Data.Family]
    item_counts = Data["Family"].value_counts()
    print(item_counts)
    print ('----------------------------')
    sns.set(font_scale=1.5)
    f, ax = plt.subplots(figsize=(12, 9))
    ax = sns.countplot(x="Family", data=Data, label="Label Count")
    sns.despine(bottom=True)
    Data = Data.drop('Family', axis = 1)
    return Data


# In[6]:


Data = load_data()
Data = labeling(Data)
# Separating Features and Label Column
Data = check_class_bias(Data)


# ## Selected motif families based on interaction-based similarity
# <p> Families with high similarity: Tandem-shear (TS), E-loop (EL), Kink-turn (KT) <br>
# Families with low similarity: Sarcin-ricin(SR) </p>

# In[7]:


def Hyper_Tuning(X, Y, x_train, x_test, y_train, y_test, SPLIT):
    ####################################################################
    cv_method = RepeatedStratifiedKFold(n_splits=SPLIT, 
                                        n_repeats=REPEATS, 
                                        random_state=RANDOM_STATE)

    params_NB = {'var_smoothing': np.logspace(0,-9, num=100)}

    gs_NB = GridSearchCV(estimator=GaussianNB(), 
                         param_grid=params_NB, 
                         cv=cv_method,
                         verbose=1, 
                         scoring='accuracy')
        
    #########################################
    Data_transformed_X = PowerTransformer().fit_transform(X)

    ######################## Score #################
    scorer = {
    'accuracy': make_scorer(accuracy_score),
    'sensitivity': make_scorer(recall_score),
    'specificity': make_scorer(recall_score,pos_label=0)
    }
    scores = cross_validate(gs_NB, X, Y, scoring=scorer, cv=cv_method, n_jobs=-1)
    print('Accuracy: %.3f, Sensitivity: %.3f, Specificity: %.3f' % (np.mean(scores['test_accuracy']), np.mean(scores['test_sensitivity']), np.mean(scores['test_specificity'])))  


# ## GNB Binary Classifier for EL
# <p> Class 1: EL <br>
# Class 0: TS, KT, SR </p>

# In[8]:


### Drop class 0
def Data_Preprocess_EL(Data):
    Data = Data[Data["motif_family"] != 11]
    Data = Data[Data["motif_family"] != 10]
    Data = Data[Data["motif_family"] != 9]
    Data = Data[Data["motif_family"] != 8]
    Data = Data[Data["motif_family"] != 7]
    Data = Data[Data["motif_family"] != 6]
    Data_Filtered = Data[Data["motif_family"] != 5]
   
    ##### E-loop #####
    EL_map = {2:1, 1:0, 3:0, 4:0}
    Data_Filtered.motif_family = [EL_map[item] for item in Data_Filtered.motif_family]
    
    Balance_Data = Data_Filtered.groupby('motif_family')    
    Data_Filtered = pd.DataFrame(Balance_Data.apply(lambda x: x.sample(Balance_Data.size().min()).reset_index(drop=True)))
    
    X = Data_Filtered.drop('motif_family', axis = 1)
    Y = Data_Filtered['motif_family']
        
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size = TEST_PERCEN, random_state = 0)
    return X, Y, x_train, x_test, y_train, y_test


# In[29]:


Data = load_data()
Data = labeling(Data)
X, Y, x_train, x_test, y_train, y_test = Data_Preprocess_EL(Data)
print("EL performance evaluation:")
Hyper_Tuning(X, Y, x_train, x_test, y_train, y_test, 3)


# ## GNB Binary Classifier for TS
# <p> Class 1: TS <br>
# Class 0: SR, KT, EL <p>

# In[10]:


### Drop class 0
def Data_Preprocess_TS(Data):
    Data = Data[Data["motif_family"] != 11]
    Data = Data[Data["motif_family"] != 10]
    Data = Data[Data["motif_family"] != 9]
    Data = Data[Data["motif_family"] != 8]
    Data = Data[Data["motif_family"] != 7]
    Data = Data[Data["motif_family"] != 6]
    Data_Filtered = Data[Data["motif_family"] != 5]
   
    ##### Tandem-shear #####
    TS_map = {1:1, 2:0, 3:0, 4:0}
    Data_Filtered.motif_family = [TS_map[item] for item in Data_Filtered.motif_family]
    
    Balance_Data = Data_Filtered.groupby('motif_family')    
    Data_Filtered = pd.DataFrame(Balance_Data.apply(lambda x: x.sample(Balance_Data.size().min()).reset_index(drop=True)))
      
    X = Data_Filtered.drop('motif_family', axis = 1)
    Y = Data_Filtered['motif_family'] 
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size = TEST_PERCEN, random_state = 0)
    return X, Y, x_train, x_test, y_train, y_test


# In[61]:


Data = load_data()
Data = labeling(Data)
X, Y, x_train, x_test, y_train, y_test = Data_Preprocess_TS(Data)
print("TS performance evaluation:")
Hyper_Tuning(X, Y, x_train, x_test, y_train, y_test, 3) 


# ## GNB Binary Classifier for KT
# <p> Class 1: KT <br>
# Class 0: SR, TS, EL </p>

# In[12]:


### Drop class 0
def Data_Preprocess_KT(Data):
    Data = Data[Data["motif_family"] != 11]
    Data = Data[Data["motif_family"] != 10]
    Data = Data[Data["motif_family"] != 9]
    Data = Data[Data["motif_family"] != 8]
    Data = Data[Data["motif_family"] != 7]
    Data = Data[Data["motif_family"] != 6]
    Data_Filtered = Data[Data["motif_family"] != 5]
  
    ##### Tandem-shear #####
    KT_map = {4:1, 1:0, 2:0, 3:0}
    Data_Filtered.motif_family = [KT_map[item] for item in Data_Filtered.motif_family]
    
    Balance_Data = Data_Filtered.groupby('motif_family')    
    Data_Filtered = pd.DataFrame(Balance_Data.apply(lambda x: x.sample(Balance_Data.size().min()).reset_index(drop=True)))
      
    X = Data_Filtered.drop('motif_family', axis = 1)
    Y = Data_Filtered['motif_family']
        
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size = TEST_PERCEN, random_state = 0)
    return X, Y, x_train, x_test, y_train, y_test


# In[43]:


Data = load_data()
Data = labeling(Data)
X, Y, x_train, x_test, y_train, y_test = Data_Preprocess_KT(Data)
print("KT performance evaluation:")
Hyper_Tuning(X, Y, x_train, x_test, y_train, y_test, 3)


# ## GNB Binary Classifier for SR
# <p> Class 1: SR <br>
# Class 0: KT, TS, EL </p>

# In[14]:


### Drop class 0
def Data_Preprocess_SR(Data):
    Data = Data[Data["motif_family"] != 11]
    Data = Data[Data["motif_family"] != 10]
    Data = Data[Data["motif_family"] != 9]
    Data = Data[Data["motif_family"] != 8]
    Data = Data[Data["motif_family"] != 7]
    Data = Data[Data["motif_family"] != 6]
    Data_Filtered = Data[Data["motif_family"] != 5]
    
    ##### Tandem-shear #####
    SR_map = {3:1, 1:0, 2:0, 4:0}
    Data_Filtered.motif_family = [SR_map[item] for item in Data_Filtered.motif_family]
    
    Balance_Data = Data_Filtered.groupby('motif_family')    
    Data_Filtered = pd.DataFrame(Balance_Data.apply(lambda x: x.sample(Balance_Data.size().min()).reset_index(drop=True)))
      
    X = Data_Filtered.drop('motif_family', axis = 1)
    Y = Data_Filtered['motif_family']
        
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size = TEST_PERCEN, random_state = 0)
    return X, Y, x_train, x_test, y_train, y_test


# In[66]:


Data = load_data()
Data = labeling(Data)
X, Y, x_train, x_test, y_train, y_test = Data_Preprocess_SR(Data)
print("SR performance evaluation:")
Hyper_Tuning(X, Y, x_train, x_test, y_train, y_test, 3)


# ## GNB Multi-class classifier

# In[16]:


### Drop class 0
def Data_Preprocess_4(Data):
    Data = Data[Data["motif_family"] != 11]
    Data = Data[Data["motif_family"] != 10]
    Data = Data[Data["motif_family"] != 9]
    Data = Data[Data["motif_family"] != 8]
    Data = Data[Data["motif_family"] != 7]
    Data = Data[Data["motif_family"] != 6]
    Data_Filtered = Data[Data["motif_family"] != 5]

    Balance_Data = Data_Filtered.groupby('motif_family')
    Data_Filtered = pd.DataFrame(Balance_Data.apply(lambda x: x.sample(Balance_Data.size().min()).reset_index(drop=True)))

    X = Data_Filtered.drop('motif_family', axis = 1)
    Y = Data_Filtered['motif_family']
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size = TEST_PERCEN, random_state = 0)
    return X, Y, x_train, x_test, y_train, y_test


# In[17]:


def Hyper_Tuning_Multiclass(X, Y, x_train, x_test, y_train, y_test):
    ####################################################################
    from sklearn.model_selection import RepeatedStratifiedKFold
    from sklearn.preprocessing import PowerTransformer
    from sklearn.metrics import precision_recall_fscore_support as score
    
    cv_method = RepeatedStratifiedKFold(n_splits=3, 
                                        n_repeats=5, 
                                        random_state=991)

    
    params_NB = {'var_smoothing': np.logspace(0,-9, num=100)}

    gs_NB = GridSearchCV(estimator=GaussianNB(), 
                         param_grid=params_NB, 
                         cv=cv_method,
                         verbose=1, 
                         scoring='accuracy')
    
    #########################################
    Data_transformed_X = PowerTransformer().fit_transform(X)

    ######################## Score #################
    scorer = {
    'accuracy': make_scorer(accuracy_score),
    'sensitivity': make_scorer(recall_score, average='weighted'),
    'specificity': make_scorer(recall_score, pos_label=0, average='weighted')
    }

    scores = cross_validate(gs_NB, X, Y, scoring=scorer, cv=cv_method, n_jobs=-1)
    print('Accuracy: %.3f' % (np.mean(scores['test_accuracy'])))  
    
#     return predict_prob


# In[18]:


X, Y, x_train, x_test, y_train, y_test = Data_Preprocess_4(Data)
Hyper_Tuning_Multiclass(X, Y, x_train, x_test, y_train, y_test)

