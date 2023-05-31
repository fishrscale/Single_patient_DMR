#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 17:47:42 2021

@author: robingrolaux
"""
from difmet_functions import *
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import fdrcorrection


def DifMeth(patient, control, score_method='zscore', df=None, corr_method='bonferroni'):
    """
    Compute pvalue of differential methylation analysis
    
    Parameters
    ----------
    patient : pd.DataFrame
        Dataframe containing the patient m-values to evaluate for differential methylation       
    
    control : pd.DataFrame
        Dataframe containing the population statistics with column 'mean' and 'std'
        See z_score()
        
    score_method : str
        method to compute pvalue of differential methylation analysis
        value: ['zscore', 'CH']
        default: 'zscore'
        
    df : int
        Degrees of Freedom. Only necessary when using  Crawford Howell scoring_method.
        value: len(size control population) - 1
        default: None
        
    corr_method : str
        method to use for multiple testing correction
        value: see statsmodels.stats.multitest.multipletests()
        default: 'bonferroni' / Bonferroni multiple testing correction method
        
    Returns
    -------

    pval_corrected: pd.Dataframe
        pvalues of differential methylation analysis corrected for multiple testing
    
    """
    #Ensure matching index between patient and control population
    if not patient.index.equals(control.index):
        
        patient = patient.reindex(control.index)
        patient = patient.dropna()
        control = control.reindex(patient.index)
    
    #Scoring Method
    ##Crawford Howell
    if score_method == 'CH':
        
        score = CH_df(patient, control, df)
        
        pvals = TstatToPval(score, df)
     
    ##Z_score
    else:
        
        score = z_score_df(patient, control)
        
        pvals = ZscoreToPval(score)
     
    #Multiple testing correction method
    if corr_method:
        
        pval_corrected = fdr_correction(pvals, corr_method)
        
    return(pval_corrected)

# Path to pd.Dataframe() containing samples methylation values (M-values recommended)
case_path=sys.argv[1]
# Path to pd.Dataframe() containing the summary statistics (see sum_stat file) of the control population (M-values recommended)
sum_stat_path=sys.argv[2]
output_path = sys.argv[3]
score_method= sys.argv[4]
df= sys.argv[5]
corr_method= sys.argv[6]
cores_number = sys.argv[7]

case = pd.read_csv(case_path, index_col=0)
sum_stat = pd.read_csv(sum_stat_path, index_col=0)

output = DifMeth(case, sum_stat, score_method, df, corr_method)

output.to_csv(output_path+'difmeth_adjusted_pvals.tsv', sep='\t')
    
    
