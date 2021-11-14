#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 13:13:37 2020

@author: minx
"""

from tree_tools import *
from random import *
import numpy as np
import math
import time
import copy
import matplotlib.pyplot as plt

"""

INPUT:  graf, igraph object
        start, integer index of initial node
        n_inf, num of infected nodes
        q, probability of correctly detecting an infected node


OUTPUT: n_inf -- dim vector of true infection ordering
    
EFFECT: adds "infected" binary node attribute
        adds "inf_tree" binary edge attribute
        adds "detected" binary edge attribute
"""

def simulateInfection(graf, start, n_inf, q):
    
    ## iterate and sample
    ## mark infected nodes
    graf.vs["detected"] = False
    graf.vs["infected"] = False
    graf.es["inf_tree"] = False
    
    det_list = []
    
    edge_ixs = graf.incident(start)
    
    graf.vs[start]["infected"] = True
    
    ## Needed?
    if (random() < q):
        graf.vs[start]["detected"] = True
    
    ii = 0
    
    true_order = [start]
    
    while (ii < n_inf-1):
        
        rand_ix = choices(edge_ixs)[0]
        edge_ixs.remove(rand_ix)
        
        cur_e = graf.es[rand_ix]
        
        if (graf.vs[cur_e.source]["infected"]):
            u_next = cur_e.target
        else:
            u_next = cur_e.source
            
        if (not graf.vs[u_next]["infected"]):
            ii = ii + 1
            cur_e["inf_tree"] = True
            graf.vs[u_next]["infected"] = True
            edge_ixs = edge_ixs + graf.incident(u_next)
            
            if (random() < q):
                graf.vs[u_next]["detected"] = True
                det_list.append(u_next)

            true_order.append(u_next)
    #print(det_list)
            
    return(true_order)
    

""" 
INPUT:  graf, igraph object
        M_outer, number of outer Gibbs iterations
        M_pass, number of passes per outer iteration

OUTPUT:     dict mapping node_index -> posterior root probability

REQUIRE: some nodes of graf has "infected" (binary) attribute

"""

def inferInfection(graf, q, **mcmc_params):
    
    
    ## generates an initial tree and initial sequence from the tree
    graf1 = graf.copy()
    graf2 = graf.copy()
    
    guess_inf1 = generateInfectionTree(graf1)
    perm1 = generateSeqFromTree(graf1, guess_inf1)
    outward1 = computeOutDegreeFromSeq(graf1, perm1)
    
    
    guess_inf2 = generateInfectionTree(graf2)
    perm2 = generateSeqFromTree(graf2, guess_inf2)
    outward2 = computeOutDegreeFromSeq(graf2, perm2)
    
    adjustSubtreeSizes(graf1, perm1, perm1[0])
    adjustSubtreeSizes(graf2, perm2, perm2[0])
    
    n_inf1 = len(perm1)
    n_inf2 = len(perm2)
    #print(n_inf)
    
    
    n = len(graf1.vs)
    freq1 = np.zeros(n)
    freq2 = np.zeros(n)
    
    ## Outer LOOP
    ii = 0
    done = False
    dist = 0
    prop_accs = []
    
    conv_thr = mcmc_params["conv_thr"]
    min_iters = mcmc_params["min_iters"]
    max_iters = mcmc_params["max_iters"]
    
    M_burn = mcmc_params["M_burn"]
    k = mcmc_params["k"]
    step_ratio = mcmc_params["step_ratio"]
    
    acc_block = mcmc_params["acc_block"]
    acc_cut = mcmc_params["acc_cut"]
    k_decr = mcmc_params["k_decr"]
    
    while not done and ii < max_iters:
        if ii%1000 == 0:
            print('loop:', ii)
        
        burn_in = ii < M_burn
        mcmc_params["burn_in"] = burn_in
        
        step = int(np.ceil(k * step_ratio))
        if ii % acc_block == acc_block-1:
            cur_acc = prop_accs[(ii-acc_block-1):(ii+1)]
            
            if (np.mean(cur_acc) < acc_cut):
                
                k = max(k-k_decr, 3)
                step = int(np.ceil(k * step_ratio))
                print((ii, k, step, np.mean(cur_acc)))
        
        mcmc_params["step"] = step
        mcmc_params["k"] = k
        
        perm1, n_inf1, freq1, outward1, prop_acc1 = updatePerm(graf1, perm1, q, n_inf1, freq1, outward1, **mcmc_params)
        perm2, n_inf2, freq2, outward2, prop_acc2 = updatePerm(graf2, perm2, q, n_inf2, freq2, outward2, **mcmc_params)
        
        freq_sum1 = np.sum(freq1)
        freq_sum2 = np.sum(freq2)
        
        if ii > M_burn and freq_sum1 > 0 and freq_sum2 > 0:
            norm_freq1 = freq1 / freq_sum1
            norm_freq2 = freq2 / freq_sum2
            dist = np.sum(np.abs(norm_freq1 - norm_freq2)) / 2
            # dist = np.linalg.norm(norm_freq1 - norm_freq2)
            if ii > min_iters + M_burn:
                done = (dist < conv_thr)
                if done:
                    print('total loops:', ii)
                    break
                # elif ii > 2000 + M_burn and tv > 0.1:
                #     graf1 = graf.copy()
                #     graf2 = graf.copy()
                    
                #     guess_inf1 = generateInfectionTree(graf1)
                #     perm1 = generateSeqFromTree(graf1, guess_inf1)
                #     outward1 = computeOutDegreeFromSeq(graf1, perm1)
                    
                #     guess_inf2 = generateInfectionTree(graf2)
                #     perm2 = generateSeqFromTree(graf2, guess_inf2)
                #     outward2 = computeOutDegreeFromSeq(graf2, perm2)
                    
                #     adjustSubtreeSizes(graf1, perm1, perm1[0])
                #     adjustSubtreeSizes(graf2, perm2, perm2[0])
                    
                #     n_inf1 = len(perm1)
                #     n_inf2 = len(perm2)
                #     #print(n_inf)
                    
                #     freq1 = np.zeros(n)
                #     freq2 = np.zeros(n)
                    
                #     ii = 0
                
            if ii % 1000 == 0:
                print(dist)
                print(np.mean(prop_accs))
        # if prop_acc1 < 0.85 and prop_acc1 > 0 and k_mid1 > 5:
        #     k_mid1 = k_mid1 - 1
        #     print('loop:', ii)
        #     print(dist)
        #     print('k ', k_mid1)
        
        # if prop_acc2 < 0.85 and prop_acc2 > 0 and k_mid2 > 5:
        #     k_mid2 = k_mid2 - 1
            
        prop_accs.append(prop_acc1)
        #print(prop_acc1)
        ii = ii + 1
    
    plt.scatter(list(range(len(prop_accs))), prop_accs)
    plt.show()

    #print(prop_accs)

    print("done:", done)
    distr1 = freq1 / np.sum(freq1)
    distr2 = freq2 / np.sum(freq2)
    
    return((distr1 + distr2) / 2)
    

def updatePerm(graf, perm, q, n_inf, freq, outward, **mcmc_params):
    tot_acc = 0
    
    burn_in = mcmc_params["burn_in"]
    M_pass = mcmc_params["M_pass"]
    k = mcmc_params["k"]
    step = mcmc_params["step"]
    
    if random() < 0.5:
        ## Inner transposition loop, swapping        
        h_weight = countAllHist(graf, perm[0], False)[0]
        
        for jj in range(M_pass):
            perm, outward, w, acc = nodesSwap(graf, n_inf, perm, outward, h_weight, **mcmc_params)
            tot_acc = acc + tot_acc
            if not burn_in:
                freq = w + freq
                
        ## re-orient edges, pick tree in a way that sequence from perm preserved
        graf.es["tree"] = False
        #tree_start = time.time()
        
        for kk in range(1,n_inf):
            cur_vix = perm[kk]
            
            cur_edges = graf.incident(cur_vix)
            
            valid_edges = [eix for eix in cur_edges if 
                           otherNode(graf.es[eix], cur_vix) in perm[0:kk]]
            assert(len(valid_edges) > 0)
            my_edge = choices(valid_edges)[0]
            graf.es[my_edge]["tree"] = True
            
            graf.vs[cur_vix]["pa"] = otherNode(graf.es[my_edge], cur_vix)
            
        countSubtreeSizes(graf, perm[0])
        #tree_end = time.time()
        #print('remake tree:', tree_end - tree_start)
        tot_acc = tot_acc / (M_pass * (1 + np.ceil( (n_inf - k)/step )) )
    else:
        # change_len_start = time.time()
        perm, outward, acc = changeLength(graf, n_inf, perm, outward, q)
        
        n_inf = len(perm)
        # change_len_end = time.time()            
        # print('change_len:', change_len_end - change_len_start)
    return perm, n_inf, freq, outward, tot_acc

"""
Propose lengthening or shortening ordering 

EFFECT:     creates "tree" binary edge attribute

"""   
def changeLength(graf, n_inf, perm, outward, q):
    n = len(graf.vs)
    acc = 0
    #assert len(perm) == n_inf
    if random() < 0.5:
        # propose increasing ordering
        if random() < 1-q and n_inf < n:
            acc = 1
            e_list = []
            
            for i in range(n_inf):
                pot_edges = graf.incident(perm[i])
                valid_edges = [eix for eix in pot_edges if 
                               otherNode(graf.es[eix], perm[i]) not in perm]
                e_list = valid_edges + e_list 
            
            my_edge = choices(e_list)[0]
            graf.es[my_edge]["tree"] = True
            head = graf.es[my_edge].source
            tail = graf.es[my_edge].target
            leaf = tail
            if tail in perm:
                leaf = head
                graf.vs[leaf]["pa"] = tail
            else:
                graf.vs[leaf]["pa"] = head
            perm.append(leaf)
            ii_nbs = graf.neighbors(leaf)
            graf.vs[leaf]["subtree_size"] = 0
            ancs = getAncestors(graf, leaf)
            for v in ancs:
                graf.vs[v]["subtree_size"] = graf.vs[v]["subtree_size"] + 1
            
            prev = outward[-1]
            num_backward = len([vix for vix in ii_nbs if vix in perm])
            outward.append(prev - num_backward + (len(ii_nbs) - num_backward))
            n_inf = n_inf + 1
            
    else:
        # propose decreasing ordering
        if not graf.vs[perm[-1]]["detected"]:
            acc = 1
            ancs = getAncestors(graf, perm[-1])
            for v in ancs:
                graf.vs[v]["subtree_size"] = graf.vs[v]["subtree_size"] - 1
            e = graf.get_eid(perm[-1], graf.vs[perm[-1]]["pa"])
            graf.es[e]["tree"] = False
            graf.vs[perm[-1]]["pa"] = -1
            n_inf = n_inf - 1
            perm = perm[:-1]
            outward = outward[:-1]
    
    # orig_out = computeOutDegreeFromSeq(graf, perm)
    # if outward != orig_out:
    #     print(perm)
    #     print(orig_out)
    #     print(outward)
    #     assert False
    return perm, outward, acc


"""
Potentially swap nodes in ordering

EFFECT:     creates "tree" binary edge attribute

"""    
def nodesSwap(graf, n_inf, perm, outward, all_weight, **mcmc_params):
    
    M_rootsamp = mcmc_params["M_rootsamp"]
    step = mcmc_params["step"]
    k = mcmc_params["k"]
    k_root = mcmc_params["k_root"]
    
    acc = 0

    w = np.zeros(len(graf.vs))
    
    starts = []
    for i in range(0, n_inf):
        cur_start = step * i 
        if (cur_start >= n_inf - k):
            starts.append(n_inf - k)
            break
        else:
            starts.append(cur_start)
            
    #starts = []
    #for i in range(0, n_inf-k_mid+1):
    #    starts.append(n_inf-k_mid - (i*step % (n_inf-k_mid+1)))
    #starts.append(0)
    
    for i in starts:
        cur_pos = i
        if (cur_pos == 0):
            #start_block_start = time.time()
            #print('switch block 0 to k')
            ## deal with root separately
            
            h_weight = [0] * k
            for i in range(k):
                h_weight[i] = all_weight[perm[i]]
                
            
            new_perm, root_dict = switchStart(graf, perm, k, h_weight, {})
            
            
            h_weight = [0] * k_root
            for i in range(k_root):
                h_weight[i] = all_weight[perm[i]]
            
            
            for i in range(M_rootsamp):
                p, root_dict = switchStart(graf, perm, k_root, h_weight, root_dict)
                cur_out = computeOutDegreeFromSeq(graf, p)
                
                denom1 = np.sum(np.log(outward[1:k]))
                denom2 = np.sum(np.log(cur_out[1:k]))
                thr = denom1 - denom2
            
                if random() < np.exp(min(0, thr)):
                    w[p[0]] = w[p[0]] + 1
                
            
            #w = w / np.sum(w)
            
            out_new = computeOutDegreeFromSeq(graf, new_perm)
            
            #thr = np.prod(np.divide(outward[1:k], out_new[1:k]))
            denom1 = np.sum(np.log(outward[1:k]))
            denom2 = np.sum(np.log(out_new[1:k]))
            thr = denom1 - denom2
            
            acc = acc + np.exp(min(0, thr))
            
            if random() < np.exp(min(0, thr)):
                perm[0:k] = new_perm
                outward[0:k] = out_new
            adjustSubtreeSizes(graf, perm[0:k], perm[0])
            
            #start_block_end = time.time()
            #print('start blck:', start_block_end - start_block_start)
            
        else:
            #print('switch block start to start + k')
    
            ## regular swapping
            #mid_switch_start = time.time()
            
            pot_perm = switchMiddle(graf, perm, cur_pos, k)
            new_out_subseq = computeOutDegreeSubseq(graf, pot_perm, outward[cur_pos - 1], cur_pos, k)
            
            #thr = np.prod(np.divide(outward[cur_pos:cur_pos + k_mid], new_out_subseq))
            denom1 = np.sum(np.log(outward[cur_pos:cur_pos + k - 1]))
            denom2 = np.sum(np.log(new_out_subseq[:-1]))
            thr = denom1 - denom2
            
            
            
            # if thr < 0:
            #     print(outward)
            #     print(outward[cur_pos: cur_pos + k_mid])
            #     print(perm)
            #     print(new_out_subseq)
            #     print(pot_perm)
            #     assert False
            # if (random() < 0.1):
            #     print(thr)
            
            acc = acc + np.exp(min(0, thr))
            
            if random() < np.exp(min(0, thr)):
                perm = pot_perm
                outward[cur_pos: cur_pos + k] = new_out_subseq
            #mid_switch_end = time.time()
            #print('mid block', mid_switch_end - mid_switch_start)
    # print(perm)
    # # print('\n')
    # # print(n_inf)
    # for ind in range(n_inf):
    #     print(ind)
    #     v = perm[ind]
    #     g = getAncestors(graf, v)
    #     if g[-1] != perm[0]:
    #         print(v)
    #         print(getAncestors(graf, v))
    #         assert False
    #     for anc in g:
    #         if anc not in perm[0:ind+1]:
    #             print(v)
    #             print(getAncestors(graf, v))
    #             assert False
    return perm, outward, w, acc



"""
Used for initialization

EFFECT:     creates "tree" binary edge attribute

"""    
def generateInfectionTree(graf):
    n = len(graf.vs)
    
    ## generate an initial random tree
    graf.es["tree"] = False
    det_vtxs = [ix for ix in range(n) if graf.vs[ix]["detected"]]
    n_det = len(det_vtxs)
    tree_vtxs = det_vtxs.copy()
    b = graf.subgraph([graf.vs[ix] for ix in det_vtxs])
    
    paths = b.get_shortest_paths(0, to=list(range(1, n_det)))
    
            
    for i in range(1, n_det):
        if len(paths[i-1]) < 1:
            p = graf.get_shortest_paths(det_vtxs[i], to=det_vtxs[0])[0]
            for v in p:
                if v not in tree_vtxs:
                    tree_vtxs.append(v)
    tree_vtxs.sort()
    
    bar = graf.subgraph([graf.vs[ix] for ix in tree_vtxs]) 
    
    wilsonTree(bar, 0)
    tree_eixs = [ix for ix in range(len(bar.es)) if bar.es[ix]["tree"]]
    
    for eix in tree_eixs:        
        head = bar.es[eix].source
        tail = bar.es[eix].target
        
        old_eix = graf.get_eid(tree_vtxs[head], tree_vtxs[tail])
        graf.es[old_eix]["tree"] = True    
    
    return tree_vtxs
    
"""
Modify input seq so that it outputs an ordering
that is valid topological sorting of the tree.


REQUIRE: nodes have "parent" attribute

"""
def straightenSeq(graf, seq):
    
    perm = seq.copy()
    
    n_inf = len(perm)
    
    perm_pos = {}
    for ii in range(n_inf):
        perm_pos[perm[ii]] = ii
                
        
    inf_mark = {}
    inf_mark[perm[0]] = True
    
    for ii in range(1, n_inf):
        cur_ix = perm[ii]
        
        pa_ix = graf.vs[cur_ix]["parent"] ## parent of current node
        
        if (pa_ix in inf_mark):
            inf_mark[cur_ix] = True
            continue
                
        while (True):
            pa_ix2 = graf.vs[pa_ix]["parent"]
            
            if (pa_ix2 in inf_mark):
                break
            else:
                pa_ix = pa_ix2
        
        
        ## pa_ix is most ancestral node that violates
        ## constraint
        
        ## swap the position of cur_ix and pa_ix
        cur_pos = ii
        pa_pos = perm_pos[pa_ix]
        
        perm[cur_pos] = pa_ix
        perm[pa_pos] = cur_ix
        
        inf_mark[pa_ix] = True
        
        perm_pos[pa_ix] = cur_pos
        perm_pos[cur_ix] = pa_pos
        
    
    return(perm)    
    
    
    
"""
REQUIRE: edges have "tree" attribute

EFFECT: creates "parent" node attribute

"""
def generateSeqFromTree(graf, guess_inf):
    normalized_h = countAllHist(graf, guess_inf[0])[0]
    n_inf = len(guess_inf)
    h_weight = [0] * n_inf
    for i in range(n_inf):
        h_weight[i] = normalized_h[guess_inf[i]]
    
    perm = [0] * n_inf
    perm[0] = choices(guess_inf, h_weight)[0]
    
    remain_nodes = [i for i in guess_inf if i != perm[0]]
    perm[1:] = np.random.permutation(remain_nodes)
    
    root = perm[0]
    
    graf.es["marked"] = False
    graf.vs[root]["parent"] = -1
    
    cur_set = [root]
    while (len(cur_set) > 0):
        vix = cur_set.pop(0)
        
        for eix in graf.incident(vix):
            e = graf.es[eix]
            if (not e["tree"] or e["marked"]):
                continue
            
            new_vix = otherNode(e, vix)
            cur_set.append(new_vix)
            
            graf.vs[new_vix]["parent"] = vix
            e["marked"] = True
    
            
    perm = straightenSeq(graf, perm)
    
    return(perm)

"""
REQUIRE: edges have "tree" attribute

EFFECT: creates "parent" node attribute

"""
def generateAltSeq(graf, guess_inf, other_perm):
    perm = []
    for elt in reversed(other_perm):
        if elt in guess_inf:
            perm.append(elt)
    
    root = perm[0]
    
    graf.es["marked"] = False
    graf.vs[root]["parent"] = -1
    
    cur_set = [root]
    while (len(cur_set) > 0):
        vix = cur_set.pop(0)
        
        for eix in graf.incident(vix):
            e = graf.es[eix]
            if (not e["tree"] or e["marked"]):
                continue
            
            new_vix = otherNode(e, vix)
            cur_set.append(new_vix)
            
            graf.vs[new_vix]["parent"] = vix
            e["marked"] = True
    perm = straightenSeq(graf, perm)
    
    return(perm)
    

def computeOutDegreeSubseq(graf, perm, old_out, start, k):
    ## pre-compute the number of outward edges
    outward = [0] * k
    for ii in range(0, k):
        if (ii == 0):
            prev = old_out

            ii_nbs = graf.neighbors(perm[ii + start])
            num_backward = len([vix for vix in ii_nbs if vix in perm[0:start+ii]])
        
            outward[ii] = prev - num_backward + (len(ii_nbs) - num_backward)
        else:
            
            prev = outward[ii - 1]
        
            ii_nbs = graf.neighbors(perm[ii + start])
            num_backward = len([vix for vix in ii_nbs if vix in perm[0:start+ii]])
        
            outward[ii] = prev - num_backward + (len(ii_nbs) - num_backward)
        
    return(outward)        
    
def computeOutDegreeFromSeq(graf, perm, end=-1):
    """
    

    Parameters
    ----------
    graf : igraph object
        Input graph.
    perm : array
        A permutation of a subset of nodes, 
        represented as integers.
    end : TYPE, optional
        DESCRIPTION. The default is -1.

    Returns
    -------
    outward : an array of integers
        outward[i] is the number of edges from 
        nodes perm[1:i] to other nodes in the graph
    

    """
    if (end == -1):
        end = len(perm)
    
    n = len(graf.vs)
    #n_inf = len(inf_nodes)
    
    ## pre-compute the number of outward edges
    outward = [0] * end
    for ii in range(end):
        if (ii == 0):
            outward[ii] = len(graf.incident(perm[0]))
        else:
            
            prev = outward[ii - 1]
        
            ii_nbs = graf.neighbors(perm[ii])
            num_backward = len([vix for vix in ii_nbs if vix in perm[0:ii]])
        
            outward[ii] = prev - num_backward + (len(ii_nbs) - num_backward)
        
    return(outward)


def getConfidenceSet(freq, eps):
    tot = 0
    sorted_inds = np.flip(np.argsort(freq))
    cred_set = []
    for k in sorted_inds:
        cred_set.append(k)
        tot = tot + freq[k]
        if tot > 1 - eps:
            return cred_set
    


def infnb(graf, vix):    
    return([v for v in graf.neighbors(vix) if graf.vs[v]["infected"]])


