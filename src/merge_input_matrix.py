#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 08:34:57 2021

This code is used to produce the input matrix of the joint 
logistic regression model.

@author: Tianyun Zhang
"""



import numpy as np
import pandas as pd

def fileConvert(input_file,mutsig_file,label_number):
    with open(mutsig_file,'r') as m:
        mut=pd.read_table(m,)
        #mut_oh=pd.get_dummies(mut[['Mut_type','Seq_context','Mut_spec']])
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
    annotation=result['ANN'].str.split('|', expand=True,n=14)
    annotation=annotation.iloc[:,1:14]
    cols = [2,3]
    annotation=annotation.drop(annotation.columns[cols], axis=1)
    annotation.columns=['Annotation','Putative_impact','Type','NM',
                        'Biotype','Rank','HGVS.c','HGVS.p','cDNA','CDS','AA']
    result=pd.concat([result,annotation],axis=1)
    result=result.drop(['ANN'],axis=1)
    return result

