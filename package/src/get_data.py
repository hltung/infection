# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 00:02:20 2021

@author: janes
"""
import pickle

filename = 'road_outputs'

with open(filename, 'rb') as f:
    test_type,outlist,succ,times = pickle.load(f)

n_inf = test_type['n_inf']

for i in range(len(outlist)):
    saved_dict = outlist[i]
    true_order = saved_dict['true_order']
    freq = saved_dict['freq']
    cred_sizes = saved_dict['cred_sizes']


print(succ)
print(times)