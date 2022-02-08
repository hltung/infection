#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues Feb  7 09:42:49 2022

@author: hltung
"""


import igraph
from infect_tools_noisy_conv import *
import scipy.stats
import numpy as np
from random import *

foo = Graph.Lattice(dim=[100, 100], circular=False)

n = len(foo.vs)
m = len(foo.es)


n_inf = 100
q = 0.75
start = 5050

# epsilon levels for the credible sets we are creating
eps_ls = [0.3, 0.2, 0.1, 0.05, 0.01] # must be decreasing

mcmc_params = {"M_burn" : 200,
               "k_root" : 15,
               "k" : 50,
               "M_pass" : 1,
               "step_ratio" : 0.4,
               "M_rootsamp" : 10,
               "acc_block" : 30,
               "acc_cut" : 0.075,
               "k_decr" : 5}


n_trials = 10
in_set = 0

out_list = []

# simulate an infection process on the graph foo
true_order = simulateInfection(foo, start, n_inf, q)    

# calculate the priors for each node given the graph and detected infected nodes
freq = inferInfection(foo, q, min_iters=6000, max_iters=40000, conv_thr=0.1, **mcmc_params)

ordered_freq = [0] * n_inf

print('results:')

# print a list of node indices and calculated prior for the node in the order that each node was infected  
for ii in range(n_inf):
    print((true_order[ii], freq[true_order[ii]]))
    ordered_freq[ii] = freq[true_order[ii]]
ordered_freq = np.array(ordered_freq)
    
tot = 0
sorted_inds = np.flip(np.argsort(freq))

# print a list of node indices and calculated prior for the node from largest to smallest prior
print('sorted priors:')
cred_sizes = []

# create credible sets for the patient zero based on the calculated priors 
cur_eps_ix = 0
for i in range(n_inf):
    k = sorted_inds[i]
    tot = tot + freq[k]
    print((k, freq[k]))
    if tot > 1 - eps_ls[cur_eps_ix]:
        cred_sizes.append(i)
        cur_eps_ix = cur_eps_ix + 1
        
    if (cur_eps_ix >= len(eps_ls)):
        break

print("credible set sizes:")
print(cred_sizes)

root_post_ix = np.where(sorted_inds == true_order[0])[0][0]


