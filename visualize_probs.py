#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 09:42:49 2020

@author: minx
"""


import igraph
from tree_tools import *
from infect_tools_noisy2 import *
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from random import *

foo = Graph.Lattice(dim=[10, 10], circular=False)
#foo = Graph.Erdos_Renyi(n=50, m=125)

# foo = igraph.read("data/global-net.dat")
# foo.delete_vertices(0)
# foo.to_undirected()
# foo.simplify()
# #foo = foo.clusters().giant()


# degs = np.array(foo.degree())
# plt.hist(degs, bins='auto')
# #plt.xscale('log')
# plt.yscale('log')
# plt.show()


n = len(foo.vs)
m = len(foo.es)

n_inf = 30
q = 0.9
eps = 0.1


## wuhan 3426
## new york 2229
## saskatoon 3769


#generateInfectionTree(foo)
#tmp = generateSeqFromTree(foo)
#outward = computeOutDegreeFromSeq(foo, tmp)

first = choices(list(range(n)), foo.degree())[0]
true_order = simulateInfection(foo, first, n_inf, q)
#freq = inferInfection(foo, M_outer=70000, M_trans=2000)
freq = inferInfection(foo, q, M_outer=600, M_trans=200)

ordered_freq = [0] * n_inf

print('results:')

for ii in range(n_inf):
    if (true_order[ii] not in freq):
        freq[true_order[ii]] = 0
    print((true_order[ii], freq[true_order[ii]]))
    ordered_freq[ii] = freq[true_order[ii]]
ordered_freq = np.array(ordered_freq)
    
inf_deg = [0] * n
freqls = [0] * n
for ix in range(n):
    inf_deg[ix] = len(infnb(foo, ix))
    if (ix in freq):
        freqls[ix] = freq[ix]

print('epsilon:', eps)
print('credible set')
tot = 0
d = {k: v for k, v in sorted(freq.items(), key=lambda item: item[1], reverse = True)} 
cred_set = []
for k in d:
    cred_set.append(k)
    tot = tot + d[k]
    print((k, d[k]))
    if tot > 1 - eps:
        break
print(cred_set)


foo.es["color"] = "gray"
foo.vs["label_dist"] = 2
for i in range(n):
    foo.vs[i]["size"] = freqls[i] * 50 + 4

for mye in foo.es:
    if mye["inf_tree"]:
        mye["color"] = "black"

for j in range(n_inf):
    foo.vs[true_order[j]]["label"] = j + 1

# for idx in cred_set:
#     foo.vs[idx]["color"] = "blue"

det_set = [v for v in foo.vs if v["detected"]]

for v in det_set:
    v["color"] = "blue"
# foo.vs[first]["shape"] = "rectangle"

out_fig = "tmp.eps"
graftree = foo.copy()
graftree.delete_edges([mye for mye in graftree.es if not mye["inf_tree"]])

mylayout = foo.layout_fruchterman_reingold(niter=500)
plot(foo, out_fig, layout=mylayout)


# inf_ixs = [vix for vix in range(n) if foo.vs[vix]["infected"]]
# bar = foo.subgraph([foo.vs[i] for i in inf_ixs])
# bar.vs["inf_label"] = inf_ixs


# label_to_subix = {}
# for i in range(n_inf):
#     label_to_subix[bar.vs[i]["inf_label"]] = i


# # out_fig = "airport.eps"
# bar.vs["size"] = 5
# #bar.vs["label"] = inf_ixs
# bar.vs["label"] = ""
# bar.vs["color"] = "darkgrey"

# if 1175 in label_to_subix:
#     bar.vs[label_to_subix[1175]]["label"] = "Hong Kong"
#     bar.vs[label_to_subix[1175]]["color"] = "orange"

# if 325 in label_to_subix:
#     bar.vs[label_to_subix[325]]["label"] = "Beijing"
#     bar.vs[label_to_subix[325]]["color"] = "orange"

# if 487 in label_to_subix:
#     bar.vs[label_to_subix[487]]["label"] = "Guangzhou"
#     bar.vs[label_to_subix[487]]["color"] = "orange"

# if 3426 in label_to_subix:
#     bar.vs[label_to_subix[3426]]["label"] = "Wuhan"
#     bar.vs[label_to_subix[3426]]["color"] = "orange"

# if 1946 in label_to_subix:
#     bar.vs[label_to_subix[1946]]["label"] = "Milan"
#     bar.vs[label_to_subix[1946]]["color"] = "orange"

# if 2773 in label_to_subix:
#     bar.vs[label_to_subix[2773]]["label"] = "San Francisco"
#     bar.vs[label_to_subix[2773]]["color"] = "orange"

# if 811 in label_to_subix:
#     bar.vs[label_to_subix[811]]["label"] = "Dubai"
#     bar.vs[label_to_subix[811]]["color"] = "orange"

# if 330 in label_to_subix:
#     bar.vs[label_to_subix[330]]["label"] = "Bangkok"
#     bar.vs[label_to_subix[330]]["color"] = "orange"

# bar.vs["label_size"] = 14
# bar.vs["label_dist"] = 3
# bar.es["color"] = "gray"
# my_layout = bar.layout_fruchterman_reingold(niter=1000)
# plot(bar, out_fig, layout = my_layout)


## add names to the plot





# eign = bar.eigenvector_centrality()
# eign = np.array(eign)
# #lbls = bar.vs["inf_label"]
# #np.c_[eign, lbls]


# label_to_order = {}
# for i in range(n_inf):
#     label_to_order[true_order[i]] = i

# ## reorder eign so that it matches with true_order

# ordered_eign = [0] * n_inf
# for ix in range(n_inf):
#     ordered_eign[label_to_order[bar.vs[ix]["inf_label"]]] = eign[ix]
# ordered_eign = np.array(ordered_eign)


scipy.stats.kendalltau(-ordered_freq, range(n_inf))
# scipy.stats.kendalltau(-ordered_eign, range(n_inf))
## output the label of the top 3 







