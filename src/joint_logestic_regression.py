#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 08:26:06 2021

This code is modelling the input matrix after the annotated pipeline

@author: Tianyun Zhang
"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score,roc_auc_score,precision_score,roc_curve
from imblearn.over_sampling import RandomOverSampler


class scRVARefine(object):
    def __init__(self, label, seq_cols, qual_cols, test_size):
        self._label_y = label
        self._seq_cols = seq_cols
        self._qual_cols = qual_cols
        self._test_size_ = test_size
        self._x_seq = []
        self._qual_dataframe_columns = []
        self._true_y = []
        self._pred_y = []
        self._pred_prob_y = []
        
    def OverSampler(self,dataset):        
          dataset_t = dataset[self._label_y]
          dataset_v = dataset.drop([self._label_y], axis = 1)
          ros = RandomOverSampler(random_state = 10)
          values, target = ros.fit_resample(dataset_v, dataset_t)
          dataset = pd.concat([values, target],axis=1)
          return dataset
      
    def dataset_split(self,dataset):
        train_dataset,test_dataset = train_test_split(dataset,
                                             test_size = self._test_size_,
                                             random_state = 0)
        return train_dataset,test_dataset
    
    def train_seq(self,train_dataset):
        X = pd.get_dummies(train_dataset[self._seq_cols])
        X = X.sort_index(axis = 1)
        self._x_seq = X.columns
        y = train_dataset[self._label_y]
        train_X,test_X, train_y, test_y = train_test_split(X,
                                                   y,
                                                   test_size = self._test_size_,
                                                   random_state=0)
        clf = LogisticRegression(C = 0.5, penalty = 'l2', tol = 0.01,
                               solver = 'liblinear',class_weight='balanced')
        clf = clf.fit(train_X,train_y)
        return clf,train_X
    
    def train_qual(self,train_dataset):
        list_cols = [self._label_y]
        for cols in self._qual_cols:
            list_cols.append(cols)
        self._qual_dataframe_columns = list_cols
        qual_data = train_dataset[self._qual_dataframe_columns]
        qual_data = qual_data.dropna()
        M_data = qual_data[self._qual_cols]
        n = qual_data[self._label_y]
        self._ss_ = StandardScaler().fit(M_data)
        M = self._ss_ .transform(M_data)
        train_M,test_M, train_n, test_n = train_test_split(M,
                                                   n,
                                                   test_size = self._test_size_,
                                                   random_state=0)
        clf_qual = LogisticRegression(C = 1, penalty = 'l1',solver = 'liblinear',
                                    tol = 0.01,class_weight='balanced')
        clf_qual = clf_qual.fit(train_M,train_n)
        return clf_qual
    
    def train_test(self,test_dataset,model1,model2):
        seq_dataset = pd.get_dummies(test_dataset[self._seq_cols])
        diff_cols = list(set(self._x_seq).difference(set(seq_dataset)))
        for i in diff_cols:
            seq_dataset.insert(seq_dataset.shape[1],i,0) 
        seq_dataset = seq_dataset[self._x_seq]
        seq_pre = model1.predict_proba(seq_dataset)
        seq_pre = pd.DataFrame(seq_pre,columns =['seq_neg','seq_pos'])
        qual_data = test_dataset[self._qual_cols]
        qual_data = self._ss_ .transform(qual_data)
        qual_pre = model2.predict_proba(qual_data)
        qual_pre = pd.DataFrame(qual_pre,columns = ['qual_neg','qual_pos'])
        input_pre = pd.concat([seq_pre,qual_pre],axis=1)
        prob_1_list = []
        pred_y=[]
        for n in range(input_pre.shape[0]):
            if input_pre['seq_pos'][n] >= 0.5:
                a = input_pre['seq_pos'][n] * 1
            else:
                a = input_pre['seq_pos'][n] * 0
            if input_pre['qual_pos'][n] >= 0.5:
                b = input_pre['qual_pos'][n] * 1
            else:
                b = input_pre['qual_pos'][n] * 0
            prob = (a + b)/2
            prob_1_list.append(prob)
        for m in range(len(prob_1_list)):
            if prob_1_list[m] >= 0.5:
                pred_y.append(1)
            else:
                pred_y.append(0)
        true_y = test_data[self._label_y].tolist()
        self._true_y = true_y
        self._pred_y = pred_y
        self._pred_prob_y = prob_1_list
            
    def test_accuracy(self):
        accuracy = accuracy_score(self._true_y, self._pred_y)
        return accuracy
    
    def test_auc_score(self):
        auc_score = roc_auc_score(self._true_y, self._pred_prob_y)
        return auc_score
    
    def test_precision(self):
        precision = precision_score(self._true_y, self._pred_y)
        return precision
    
    def test_roc_curve(self):
        fpr,tpr,threshold = roc_curve(self._true_y, self._pred_prob_y)
        return fpr,tpr,threshold
    
    def predict(self,unsure_dataset,model1,model2):
        unsure_seq = pd.get_dummies(unsure_dataset[self._seq_cols])
        uns_cols = unsure_seq.columns
        unsure_seq = unsure_seq.dropna()
        diff_cols = list(set(self._x_seq).difference(set(uns_cols)))
        for i in diff_cols:
            unsure_seq.insert(unsure_seq.shape[1],i,0) 
        unsure_seq = unsure_seq[self._x_seq]
        predict_porb = pd.DataFrame(model1.predict_proba(unsure_seq),
                                   columns = ['seq_neg','seq_pos'])
        unsure_out = pd.concat([unsure_dataset[['index','CHROM','POS']],predict_porb],
                               axis = 1)
        unsure_qual = unsure_dataset[self._qual_cols].dropna()
        unsure_qual = self._ss_.transform(unsure_qual)
        predict_porb_qual = pd.DataFrame(model2.predict_proba(unsure_qual),
                                         columns = ['qual_neg','qual_pos'])
        unsure_out = pd.concat([unsure_out,predict_porb_qual],axis=1)
        list_chrom = []
        list_pos = []
        prob_list=[]  
        list_index=[]
        for i in range(unsure_out.shape[0]): 
            if unsure_out['seq_pos'][i] >= 0.5:
                a = unsure_out['seq_pos'][i] * 1
            else:
                a = unsure_out['seq_pos'][i] * 0
            if unsure_out['qual_pos'][i] >= 0.5:
                b = unsure_out['qual_pos'][i] * 1
            else:
                b = unsure_out['qual_pos'][i] * 0
            prob = (a + b)/2  
            prob_list.append(prob)
        for j in range(len(prob_list)):
            if prob_list[j] >= 0.5:
                list_chrom.append(unsure_out['CHROM'][j])
                list_pos.append(unsure_out['POS'][j])
                list_index.append(unsure_out['index'][j])
        dict_out = {'index':list_index,
                    'CHROM':list_chrom,
                    'POS':list_pos}
        result_out = pd.DataFrame(dict_out)
        return result_out 

    def mut_spec(self,unsure_dataset,model1,model2):
        unsure_seq = pd.get_dummies(unsure_dataset[self._seq_cols])
        uns_cols = unsure_seq.columns
        unsure_seq = unsure_seq.dropna()
        diff_cols = list(set(self._x_seq).difference(set(uns_cols)))
        for i in diff_cols:
            unsure_seq.insert(unsure_seq.shape[1],i,0) 
        unsure_seq = unsure_seq[self._x_seq]
        predict_porb = pd.DataFrame(model1.predict_proba(unsure_seq),
                                   columns = ['seq_neg','seq_pos'])
        unsure_out = pd.concat([unsure_dataset[['CHROM','POS','Mut_spec','CELL']],predict_porb],
                               axis = 1)
        unsure_qual = unsure_dataset[self._qual_cols].dropna()
        unsure_qual = self._ss_.transform(unsure_qual)
        predict_porb_qual = pd.DataFrame(model2.predict_proba(unsure_qual),
                                         columns = ['qual_neg','qual_pos'])
        unsure_out = pd.concat([unsure_out,predict_porb_qual],axis=1)
        list_chrom = []
        list_pos = []
        list_cell=[]
        mut_spec=[]
        prob_list=[]
        for i in range(unsure_out.shape[0]): 
            if unsure_out['seq_pos'][i] >= 0.5:
                a = unsure_out['seq_pos'][i] * 1
            else:
                a = unsure_out['seq_pos'][i] * 0
            if unsure_out['qual_pos'][i] >= 0.5:
                b = unsure_out['qual_pos'][i] * 1
            else:
                b = unsure_out['qual_pos'][i] * 0
            prob = (a + b)/2  
            prob_list.append(prob)
        for j in range(len(prob_list)):
            if prob_list[j] >= 0.5:
                list_chrom.append(unsure_out['CHROM'][j])
                list_pos.append(unsure_out['POS'][j])
                list_cell.append(unsure_out['CELL'][j])
                mut_spec.append(unsure_out['Mut_spec'][j])
        dict_out = {'CHROM':list_chrom,
                    'POS':list_pos,
                    'CELL':list_cell,
                    'Mut_spec':mut_spec}                
        result_out = pd.DataFrame(dict_out)
        return result_out
 
