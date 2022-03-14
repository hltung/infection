#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 15:01:55 2022

@author: minx
"""



import igraph
from infect_tools_noisy_conv import *
import scipy.stats
import numpy as np
from random import *

import time

import pickle

fname = "lattice_vary_q.pkl"
init = False

n_inf = 250

q_ls = [1, 0.7, 0.4]

eps_ls = [0.5, 0.3, 0.2, 0.1, 0.05, 0.01] # must be decreasing

foo = Graph.Lattice(dim=[100, 100], circular=False)

n = len(foo.vs)
m = len(foo.es)

ntrial = 200


mcmc_params = {"M_burn" : 200,
               "k_root" : 10,
               "k" : 80,
               "M_pass" : 1,
               "step_ratio" : 0.4,
               "M_rootsamp" : 15,
               "acc_block" : 30,
               "acc_cut" : 0.1,
               "k_decr" : 5}

calctimes = np.zeros([ntrial, len(q_ls)])
all_covered = np.zeros([ntrial, len(q_ls), len(eps_ls)])
all_cred_sizes = np.zeros([ntrial, len(q_ls), len(eps_ls)])

for i in range(len(q_ls)):

    for j in range(ntrial):
    
        print((i, j))    
    
        if (init and i <= iprev and j <= jprev):
            continue
    
        q = q_ls[i]
        
        start = randrange(n)
        true_order = simulateInfection(foo, start, n_inf, q)    

        stime = time.time()
        freq = inferInfection(foo, q, min_iters=5000, max_iters=60000, conv_thr=0.1, **mcmc_params)
        etime = time.time()
        
        elapsed = etime - stime
        
        ordered_freq = [0] * n_inf
        
        # create credible sets for the patient zero based on the calculated priors 
        sorted_inds = np.flip(np.argsort(freq))

        covered = np.array([0 for x in eps_ls])
        cred_sizes = []
        tot = 0
        cur_eps_ix = 0
        
        for ii in range(n_inf):
            k = sorted_inds[ii]
            tot = tot + freq[k]

            if (k == start):
                covered[cur_eps_ix:] = 1

            if tot > 1 - eps_ls[cur_eps_ix]:
                cred_sizes.append(ii)
                cur_eps_ix = cur_eps_ix + 1
        
            if (cur_eps_ix >= len(eps_ls)):
                break

        calctimes[j, i] = elapsed
        all_covered[j, i, :] = covered
        all_cred_sizes[j, i, :] = cred_sizes
        
        with open(fname, "wb") as f:
            pickle.dump([calctimes, all_covered, all_cred_sizes, i, j], f)

        #with open(fname, "rb") as f:
        #    calctimes, all_covered, all_cred_sizes, iprev, jprev = pickle.load(f)
