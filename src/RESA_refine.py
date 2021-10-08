#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 03:30:42 2021

@author: Tianyun Zhang
"""


import sys
import numpy as np
import pandas as pd
import os
import argparse
import rpy2.robjects as robjects
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score,roc_auc_score,precision_score,roc_curve
from imblearn.over_sampling import RandomOverSampler


dir_sh = os.path.dirname(sys.argv[0])
parser = argparse.ArgumentParser()
parser.add_argument('-P', type=str,help='Positive file')
parser.add_argument('-N', type=str,help='Negative file')
parser.add_argument('-U', type=str, help='Undefined file')
parser.add_argument('-O', type=str, help='Output dir')
parser.add_argument('-S', type=str, help = 'Statistical summary',
                    default = 'True', 
                    choices=['True','False'])
args = parser.parse_args()

def fileConvert_PreR(input_file):
    with open (input_file,'r') as f:
        file=pd.read_csv(f,sep='\t')
        if len(file.columns) != 16:
            print(input_file)
            sys.exit("Error: The matrix does not match the correct input file format.")
        else:
            with open (input_file,'r') as f:
                file=pd.read_csv(f,names=['CELL','CHROM','POS','POS1','REF','ALT',
                                  'CHROM2','POS2','STRAND','REF1','ALT1',
                                  'QUAL','ID1','INFO','FORMAT',
                                  'GT:AD:DP:GQ:PL'],sep='\t')
                file['GT'] = file['GT:AD:DP:GQ:PL'].map(lambda x:x.split(':')[0])
                file['AD'] = file['GT:AD:DP:GQ:PL'].map(lambda x:x.split(':')[1])
                file['DP1'] = file['GT:AD:DP:GQ:PL'].map(lambda x:x.split(':')[2])
                file['GQ'] = file['GT:AD:DP:GQ:PL'].map(lambda x:x.split(':')[3])
                file['PL'] = file['GT:AD:DP:GQ:PL'].map(lambda x:x.split(':')[4])
                file=file.drop(['GT:AD:DP:GQ:PL'],axis=1)
            result = file.drop(['INFO','FORMAT','POS2',
                      'REF1','ALT1','ID1','CHROM2'],axis=1)
        data_df = file['INFO'].str.split(';|=', expand=False)
        store_list = []
        for n in range(0,len(data_df)):
            x = data_df[n]
            store_list.append({x[i]:x[i+1]for i in range(0,len(x)-1,2)})
        info_data = pd.DataFrame.from_dict(store_list).applymap(str).replace(r'\s+',
                                      np.nan,regex=True).replace('',np.nan)
        result = pd.concat([result, info_data], axis=1)
        annotation = result['ANN'].str.split('|', expand=True,n=14)
        annotation = annotation.iloc[:,1:14]
        cols = [2,3]
        annotation = annotation.drop(annotation.columns[cols], 
                                 axis=1)
        annotation.columns = ['Annotation','Putative_impact','Type',
                        'NM','Biotype','Rank','HGVS.c','HGVS.p','A','B','D']
        result = pd.concat([result,annotation],axis=1)
        #result = result.drop(['ANN'],axis=1)
        result['SNP'] = result['REF'] + result['ALT']
        result = result[['CHROM','POS','STRAND','REF','ALT','QUAL']]
        result.columns = ['#CHROM','POS','STRAND','REF','ALT','QUAL']
        return result


def fileConvert(input_file,mutsig_file,label_number):
    with open(mutsig_file,'r') as m:
        mut=pd.read_table(m,)
    with open (input_file,'r') as f:
        file=pd.read_csv(f,names=['CELL','CHROM','POS','POS1','REF','ALT','CHROM2',
                                  'POS2','ID','REF1','ALT1','QUAL','ID1','INFO',
                                  'FORMAT','GT:AD:DP:GQ:PL'],sep='\t')
        file.insert(file.shape[1], 'LABEL',label_number) 
        for i,label in enumerate(['GT','AD','DP','GQ','PL']):
            file[label]=file['GT:AD:DP:GQ:PL'].map(lambda x:x.split(':')[i])
        for j,label_PL in enumerate(['PL1','PL2','PL3']):
            file[label_PL]=file['PL'].map(lambda x:x.split(',')[j])
        for k,label_AD in enumerate(['AD1','AD2']):
            file[label_AD]=file['AD'].map(lambda x:x.split(',')[k])
        file=file.drop(['GT:AD:DP:GQ:PL'],axis=1)
    result=pd.concat([mut,file], axis=1)
    result=result.drop(['INFO','FORMAT','POS1','POS2','ID','REF1','ALT1','ID1','CHROM2'],axis=1)
    data_df = file['INFO'].str.split(';|=', expand=False)
    store_list=[]
    for n in range(0,len(data_df)):
        x=data_df[n]
        store_list.append({x[i]:x[i+1]for i in range(0,len(x)-1,2)})
    info_data=pd.DataFrame.from_dict(store_list).applymap(str).replace(r'\s+',np.nan,regex=True).replace('',np.nan)
    info_data=info_data.drop(['DP'],axis=1)
    result=pd.concat([result, info_data], axis=1)
    return result

def fileConvert_Und(input_file,mutsig_file):
    with open(mutsig_file,'r') as m:
        mut=pd.read_table(m)
    with open (input_file,'r') as f:
        file=pd.read_csv(f,names=['CELL','CHROM','POS','POS1','REF','ALT',
                                  'CHROM2','POS2','ID','REF1','ALT1','QUAL',
                                  'ID1','INFO','FORMAT','GT:AD:DP:GQ:PL'],
                            sep='\t')
        for i,label in enumerate(['GT','AD','DP','GQ','PL']):
            file[label]=file['GT:AD:DP:GQ:PL'].map(lambda x:x.split(':')[i])
        for j,label_PL in enumerate(['PL1','PL2','PL3']):
            file[label_PL]=file['PL'].map(lambda x:x.split(',')[j])
        for k,label_AD in enumerate(['AD1','AD2']):
            file[label_AD]=file['AD'].map(lambda x:x.split(',')[k])
        file=file.drop(['GT:AD:DP:GQ:PL'],axis=1)
    result=pd.concat([mut,file], axis=1)
    result=result.drop(['INFO','FORMAT','POS1','POS2','ID','REF1','ALT1','ID1','CHROM2'],axis=1)
    data_df = file['INFO'].str.split(';|=', expand=False)
    store_list=[]
    for n in range(0,len(data_df)):
        x=data_df[n]
        store_list.append({x[i]:x[i+1]for i in range(0,len(x)-1,2)})
    info_data=pd.DataFrame.from_dict(store_list).applymap(str).replace(r'\s+',np.nan,regex=True).replace('',np.nan)
    info_data=info_data.drop(['DP'],axis=1)
    result=pd.concat([result, info_data], axis=1)
    return result

class joint_Model(object):
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
        true_y = test_dataset[self._label_y].tolist()
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

    def Predict(self,unsure_dataset,model1,model2):
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


N = fileConvert_PreR(args.N)
pd.DataFrame(['##fileformat=VCFv4.2']).to_csv(args.O+'negative.vcf',
                  sep='\t',mode='a',index=False,header=False)
N.to_csv(args.O+'negative.vcf',sep='\t',mode='a',index=False)
P = fileConvert_PreR(args.P)
pd.DataFrame(['##fileformat=VCFv4.2']).to_csv(args.O+'postive.vcf',
                  sep='\t',mode='a',index=False,header=False)
P.to_csv(args.O+'postive.vcf',sep='\t',mode='a',index=False)
U = fileConvert_PreR(args.U)
pd.DataFrame(['##fileformat=VCFv4.2']).to_csv(args.O+'undefined.vcf',
                  sep='\t',mode='a',index=False,header=False)
U.to_csv(args.O+'undefined.vcf',sep='\t',mode='a',index=False)

#Rscript
robjects.r.source(dir_sh+'/get_mutation_type_seq_context.R')
robjects.r.anno(args.O+'/negative.vcf','negative',args.O)
robjects.r.anno(args.O+'/postive.vcf','postive',args.O)
robjects.r.anno(args.O+'/undefined.vcf','undefined',args.O)

N_file=fileConvert(args.N,args.O+'negative',0)
P_file=fileConvert(args.P,args.O+'postive',1)
input_matrix=pd.concat([P_file,N_file],axis=0)
input_matrix.to_csv(args.O+'input_train.csv')

Und_file=fileConvert_Und(args.U,args.O+'/undefined')
Und_file.to_csv(args.O+'undefined.csv')

data = pd.read_csv(args.O+'input_train.csv')

label = 'LABEL'
seq_cols = ['Mut_type','Seq_context','Mut_spec']
qual_cols = ['QUAL','DP','VAF','PL1','PL2','PL3','AD1','AD2']

data = data.drop(data[data['VAF'].isna()==True].index)

if args.S == 'True':
    sc = joint_Model(label,seq_cols,qual_cols,1/3)
    train_data,test_data = sc.dataset_split(data)
    train_data = sc.OverSampler(train_data)
    m1,train_x = sc.train_seq(train_data)
    m2 = sc.train_qual(train_data)
    sc.train_test(test_data,m1,m2)
    inf = {'Accuracy': [sc.test_accuracy()],
            'Auc score': [sc.test_auc_score()],
            'Precision':[sc.test_precision()]}
    summary = pd.DataFrame(inf)
    summary.to_csv(args.O+'summary_info.txt',sep='\t',index=False)
    und_test=pd.read_csv(args.O+'undefined.csv',index_col=0)
    und_test = und_test.reset_index()
    und_B = und_test[['CHROM','POS','CELL','Mut_spec']]
    sc.Predict(und_test,m1,m2).to_csv(args.O+'prdicted_value.csv')
else: 
    sc = joint_Model(label,seq_cols,qual_cols,1/3)
    train_data,test_data = sc.dataset_split(data)
    train_data = sc.OverSampler(train_data)
    m1,train_x = sc.train_seq(train_data) 
    m2 = sc.train_qual(train_data)
    sc.train_test(test_data,m1,m2)
    und_test=pd.read_csv(args.O+'undefined.csv',index_col=0)
    und_test = und_test.reset_index()
    und_B = und_test[['CHROM','POS','CELL','Mut_spec']]
    sc.Predict(und_test,m1,m2).to_csv(args.O+'prdicted_value.csv')




