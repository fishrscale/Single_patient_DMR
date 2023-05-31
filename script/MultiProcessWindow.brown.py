#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 09:12:55 2020

@author: robingrolaux
"""
import sys
import pandas as pd
import numpy as np
import logging
from multiprocessing import Pool, cpu_count

from scipy.stats import chi2, combine_pvalues

from EmpiricalBrownsMethod import *

def QueueGenChunks(window, case_pval, pop_Beta, annot_file):
    
    annot = annot_file
    
    case = case_pval.reindex(annot.index.tolist())
    pop_Beta = pop_Beta.reindex(annot.index.tolist())
    
    for i in range(1,23):
        
        sub = annot.loc[annot['chr'] == i]
        
        yield(window.loc[window['chr'] == i], case.reindex(sub.index), pop_Beta.reindex(sub.index))
        
    
def ProcessingDF(IndependentData):
    
    window, case_pval, pop_Beta = IndependentData
    
    for i in window.index:

        print(window.loc[i].name)
        print("################")
        cgs = window.loc[i]['cgs']#['IlmnID']

        bvals = pop_Beta.reindex(cgs)#[j]
        bvals = np.array(bvals)[~np.isnan(np.array(bvals)).any(axis=1)]
        covar_matrix = CalculateCovariances(bvals)

        for j in case_fdr.columns.values:
            
                    pvals = case_fdr.reindex(cgs)[j]#case_fdr.loc[cgs][j]
                    pvals = np.array(pvals)[~np.isnan(np.array(pvals))]
                    if not len(pvals):
                        pvals = float("nan")
                    else:
                        if len(pvals) > 1:
                            pvals = CombinePValues(covar_matrix, pvals, extra_info = False)

                    window.at[i,j] = pvals
        
    return(window)
    
    
    
def MultiProcessing(window, case_pval, annot_file, cores_number, pop_Beta):
    
    df = DataFrame()
    
    if cores_number:
        pool = Pool(cores_number)
        
    else:
        pool = Pool()
    
    for result in pool.imap_unordered(ProcessingDF, QueueGenChunks(window, case, caserrBeta, annot_file)):
        
        df = pd.concat([df, result], axis = 0)

    pool.close()
    pool.join()
    
    return(df)
# Path to pd.Dataframe() containing adjusted pvalues  
case_path=sys.argv[1]
# Path to pd.Dataframe() containing the genomic location of regions to combine
annot_window_path=sys.argv[2]
output_path = sys.argv[3]
cores_number = int(sys.argv[4])
# Annotation file containing the CpGs position
annot_path=sys.argv[5]

# Path to pd.Dataframe() containing the Beta values from the control population
pop_pathBeta=sys.argv[6]


case_pval = pd.read_csv(case_path,
                index_col = 0)

window = pd.read_csv(annot_window_path,
                  index_col = 0,
                  #As read_csv reads list as str, transformation of the first columns into a list
                  converters={4: lambda x: x[1:-1].replace('\'','').split(', ')})


case_pval = caserr.dropna(axis=0, how='any')


pop_Beta = pd.read_csv(pop_pathBeta,
                index_col = 0)



annot_file = pd.read_csv(annot_path, index_col = 0)
    
output = MultiProcessing(window, case_pval, annot_file, cores_number, pop_Beta)


output.to_csv(output_path+"/pvalPerWindowBrown.tsv", sep = "\t")
