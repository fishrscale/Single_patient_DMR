#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 13:01:47 2020

@author: robingrolaux
"""

"""
Module 2

Creation of population from different GSE files
Formatting this population to Infinium 450k or 800k
Computing stats from the population (mean, std, percentiles)

"""

import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import entropy

#Combine GSE like files into a single DataFrame reindexed based on *index*
def combine_data(index, files):
    
    
    file_list = []
    sep = ','
    
    for file in files:
        

        
        #Read .txt 
        if '.txt' in file:
            sep = '\t'
            continue
    
        file_list.append(pd.read_csv(file, index_col = 0, sep = sep).reindex(index))
    
        sep = ','
    
    return(pd.concat(file_list, axis = 1))

#Generate row-wise summary statistics about a population df
def make_stats_df(data, percentile_list = [1, 5, 25, 75, 95, 99]):

    infos = {}
    infos['mean'] = data.apply(np.nanmean, axis = 1)
    infos['std'] = data.apply(np.nanstd, axis = 1)
    infos['max'] = data.apply(np.nanmax, axis = 1)
    infos['min'] = data.apply(np.nanmin, axis = 1)
    #infos['entropy'] = data.apply(lambda x: entropy(np.histogram(x, range=(0,1))[0])/np.log(10),
    #                            axis = 1)
    #infos['var'] = data.apply(np.nanvar, axis = 1)


    for percentile in percentile_list:

        infos[str(percentile)] = data.apply(lambda x: np.nanpercentile(x, percentile), axis = 1)

    return(pd.DataFrame(infos))




args = sys.argv

i = 1
while i <= len(args)-1:
    
    #1st argument: File to reindex the data with either Infinium 450k or 800k
    #First column should contain CG index
    if args[i] == '-i':
        
        index = pd.read_csv(sys.argv[i+1], index_col = 0).index
        i += 2
    
    #2nd argument: out_folder
    elif args[i] == '-o':
        
        out_folder = args[i+1]
        i += 2
    
    
    #3rd argument: -all to combine all files in the current directory
    #None if files path given directly
    elif sys.argv[i] == '-all':
        
        files = os.listdir(args[i+1])
    
    else:
        files = args[i:]


def to_mval(data):
    data = data.replace(1, 0.999)
    data = data.replace(0, 0.001)

    return(np.log2(data/(1-data)))


"""
Ctrl_pop = combine_data(index, files)

ctrls = pd.read_csv('//Volumes/Data/ctrl_epic/ctrl_519_epic.csv.gz',
                         index_col = 0)

Sum_stat =  make_stats_df(Ctrl_pop)

Sum_stat.to_csv(out_folder+'/sum_stat.csv')

Sum_stat_m =  make_stats_df(to_mval(Ctrl_pop))

Sum_stat_m.to_csv(out_folder+'/sum_stat_m.csv')
"""
