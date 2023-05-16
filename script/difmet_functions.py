#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:20:30 2020

@author: robingrolaux
"""
from statsmodels.stats.multitest import multipletests
import pandas as pd
import scipy.stats as st
import numpy as np

def z_score(mean, sd, mu):
    """
    Compute z-score from sample and population parameters
    
    Parameters
    ----------
    mean : int
        mean of the population
        
    sd : int
        standard deviation of the population
        
    mu : int
        sample value
        

    Returns
    -------

    z-score: int
    
    """
    return(abs((mu - mean)/sd))



def z_score_df(sample, stat_df):
    """
    Apply z_score() to DataFrame
    
    Parameters
    ----------
    sample : pd.Series
        sample to calculate z-score from
        
    stat_df : pd.DataFrame
        Dataframe containing the population statistics with parameters 'mean' and 'std' of the population

        
    Returns
    -------
    z-score : pd.Series
        z-score of all the sample variables
    """
        
    return(sample.apply(lambda x: z_score(stat_df['mean'], stat_df['std'], x)))



def ZscoreToPval(z_score_matrix):
    """
    Transform matrix of absolute zscores into corresponding pvalue
    
    Parameters
    ----------
    z_score_matrix : pd.DataFrame
        DataFrame containing z-score
        

    Returns
    -------
    
    P-value DF : pd.DataFrame 
        Dataframe containing pvalues derived from z-score
    """
    
    return(pd.DataFrame(st.norm.sf(z_score_matrix)*2, index=z_score_matrix.index,
                     columns = z_score_matrix.columns.values))


def CrawfordHowell(sample, pop_mean, pop_std, df):
    """
    Compute Crawford Howell T statistic from sample and population parameters
    
    Parameters
    ----------
    sample: pd.DataFrame
        Dataframe containing the patient values
    
    pop_mean : int
        mean of the population
        
    pop_sd : int
        standard deviation of the population
        
    df : int
        degrees of freedom = len(population)-1
        

    Returns
    -------

    tval: int
        T statistic
    
    """     
    tval = (sample-pop_mean)/(pop_std*(np.sqrt((df+1)/df)))
    
    #pval = scipy.stats.t.sf(abs(tval), df=df)*2
    
    return(tval)


def CH_df(sample, stat_df, df):
    """
    Apply CrawfordHowell() to DataFrame
    
    Parameters
    ----------
    sample : pd.Series
        sample to calculate T statistics from
        
    stat_df : pd.DataFrame
        Dataframe containing the population statistics with column 'mean' and 'std'
        See CrawfordHowell()
        
    df : int
        Degrees of Freedom = (number of sample) - 1
        
    Returns
    -------
    T statistics DF : pd.DataFrame
        T statistics of all the sample variables
    """  
        
    return(sample.apply(lambda x: CrawfordHowell(x, stat_df['mean'], stat_df['std'], df)))


def TstatToPval(T_stat_matrix, df):
    """
    Transform matrix of T statistics into corresponding pvalue
    
    Parameters
    ----------
    T_score_matrix : pd.DataFrame
        DataFrame containing T statistics
        

    Returns
    -------
    P-value DF : pd.DataFrame 
        Dataframe containing pvalues derived from T statistic
    """  
    return(pd.DataFrame(st.t.sf(abs(T_stat_matrix), df=df)*2, index = T_stat_matrix.index,
                     columns = T_stat_matrix.columns.values))

def fdr_correction(Samples, method):
    """
    Apply multiple testing correction on a matrix of pvalues
    
    Parameters
    ----------
    Samples : pd.DataFrame
        Dataframe containing pvalues
        
    method : str
        method to use for multiple testing correction
        see statsmodels.stats.multitest.multipletests()
    
    Returns
    -------
    P-value DF : pd.DataFrame 
        Dataframe containing the corrected pvalues
    """
    
    df = pd.DataFrame(index=Samples.index)
    
    for i in Samples.columns.values:
        i_wo_na = Samples[i].dropna()
        corrected = multipletests(i_wo_na, method=method)
        df = pd.concat([df, pd.DataFrame(corrected[1], columns=[i], index=i_wo_na.index).reindex(Samples.index)],
                       axis=1,
                       ignore_index=False)
    return(df)

def to_mval(beta):
    """
    Transform Beta values into M values

    Parameters
    ----------
    beta : pd.DataFrame
        Dataframe containing Beta values

    Returns
    -------
    mval : pd.DataFrame
        Dataframe containing the corresponding M values
    """
    beta = beta.replace(1, 0.999)
    beta = beta.replace(0, 0.001)

    mval = np.log2(beta / (1 - beta))

    return(mval)