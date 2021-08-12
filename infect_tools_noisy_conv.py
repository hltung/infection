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
        M_trans, number of transpositions per outer iteration

OUTPUT:     dict mapping node_index -> posterior root probability

REQUIRE: some nodes of graf has "infected" (binary) attribute

"""

def inferInfection(graf, q, min_iters=500, max_iters=10000, M_trans=1000, M_burn=50, k=4, conv_thr=0.05):
    
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
    
    while not done and ii < max_iters:
        if ii%200 == 0:
            print('loop:', ii)
        
        burn_in = ii < M_burn
        perm1, n_inf1, freq1, outward1 = updatePerm(graf1, perm1, q, n_inf1, freq1, outward1, burn_in, k, M_trans)
        perm2, n_inf2, freq2, outward2 = updatePerm(graf2, perm2, q, n_inf2, freq2, outward2, burn_in, k, M_trans)
        
        freq_sum1 = np.sum(freq1)
        freq_sum2 = np.sum(freq2)
        if ii > M_burn and freq_sum1 > 0 and freq_sum2 > 0:
            norm_freq1 = freq1 / freq_sum1
            norm_freq2 = freq2 / freq_sum2
            tv = np.sum(np.abs(norm_freq1 - norm_freq2)) / 2
            if ii > min_iters:
                done = (tv < conv_thr)
            if ii % 200 == 0:
                print(tv)
        ii = ii + 1
    
    print("done:", done)
    return(freq1 / np.sum(freq1))
    

def updatePerm(graf, perm, q, n_inf, freq, outward, burn_in, k, M_trans):
    if random() < 0.5:
        ## Inner transposition loop, swapping
        
        h_weight = countAllHist(graf, perm[0])[0]
        for jj in range(M_trans):
            perm, outward, w = nodesSwap(graf, n_inf, perm, outward, h_weight, k)
            
            if not burn_in:
                for idx in range(k):
                    freq = w + freq
        ## re-orient edges, pick tree in a way that sequence from perm preserved
        graf.es["tree"] = False
        
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
    else:
        perm, outward = changeLength(graf, n_inf, perm, outward, q)
        n_inf = len(perm)
            
        
    return perm, n_inf, freq, outward 

"""
Propose lengthening or shortening ordering 

EFFECT:     creates "tree" binary edge attribute

"""   
def changeLength(graf, n_inf, perm, outward, q):
    n = len(graf.vs)
    if random() < 0.5:
        # propose increasing ordering
        if random() < 1-q and n_inf < n:
            n_inf = n_inf + 1
            e_list = []
            
            for i in range(n_inf-1):
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
            ii_nbs = graf.neighbors(tail)
            graf.vs[leaf]["subtree_size"] = 0
            ancs = getAncestors(graf, leaf)
            for v in ancs:
                graf.vs[v]["subtree_size"] = graf.vs[v]["subtree_size"] + 1
            
            prev = outward[n_inf - 2]
            num_backward = len([vix for vix in ii_nbs if vix in perm])
            outward.append(prev - num_backward + (len(ii_nbs) - num_backward))
    else:
        # propose decreasing ordering
        if not graf.vs[perm[-1]]["detected"]:   
            ancs = getAncestors(graf, perm[-1])
            for v in ancs:
                graf.vs[v]["subtree_size"] = graf.vs[v]["subtree_size"] - 1
            e = graf.get_eid(perm[-1], graf.vs[perm[-1]]["pa"])
            graf.es[e]["tree"] = False
            graf.vs[perm[-1]]["pa"] = -1
            n_inf = n_inf - 1
            perm = perm[:-1]
            outward = outward[:-1]
            
    return perm, outward


"""
Potentially swap nodes in ordering

EFFECT:     creates "tree" binary edge attribute

"""    
def nodesSwap(graf, n_inf, perm, outward, all_weight, k):
    cur_pos = choices(list(range(n_inf-k+1)))[0]
    M_0 = math.factorial(k) 
    w = np.zeros(len(graf.vs))
   
    if (cur_pos == 0):
        # print('switch block 0 to k')
        ## deal with root separately
        h_weight = [0] * k
        for i in range(k):
            h_weight[i] = all_weight[perm[i]]
        new_perm, root_dict = switchStart(graf, perm, k, h_weight, {})
        for i in range(M_0):
            p, root_dict = switchStart(graf, perm, k, h_weight, root_dict)
            out = computeOutDegreeFromSeq(graf, p)
            w[p[0]] = w[p[0]] + 1/np.prod(out[1:])
        
        w = w / np.sum(w)
        
        out_new = computeOutDegreeFromSeq(graf, new_perm)
        denom1 = np.prod(outward[1:k])
        denom2 = np.prod(out_new[1:])
        
        if random() < min(1, denom1/denom2):
            perm[0:k] = new_perm
            outward[0:k] = out_new
        adjustSubtreeSizes(graf, perm[0:k], perm[0])
    
            
    else:
        # print('switch block start to start + k')

        ## regular swapping
        ## cur_vix NOT rot
        pot_perm = switchMiddle(graf, perm, cur_pos, k)
        new_out_subseq = computeOutDegreeSubseq(graf, pot_perm, outward[cur_pos - 1], cur_pos, k)
        denom1 = np.prod(outward[cur_pos:cur_pos + k])
        denom2 = np.prod(new_out_subseq)
        
        if (random() < min(1, denom1/denom2)):
            perm = pot_perm
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
    return perm, outward, w



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
            num_backward = len([vix for vix in ii_nbs if vix in perm[start:start+ii]])
        
            outward[ii] = prev - num_backward + (len(ii_nbs) - num_backward)
        else:
            
            prev = outward[ii - 1]
        
            ii_nbs = graf.neighbors(perm[ii + start])
            num_backward = len([vix for vix in ii_nbs if vix in perm[start:start+ii]])
        
            outward[ii] = prev - num_backward + (len(ii_nbs) - num_backward)
        
    return(outward)        
    
def computeOutDegreeFromSeq(graf, perm, end=-1):
    
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



def infnb(graf, vix):    
    return([v for v in graf.neighbors(vix) if graf.vs[v]["infected"]])


