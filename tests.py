# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 14:37:29 2021

@author: janes
"""

import igraph
from tree_tools import *
from infect_tools_noisy_conv import *
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from random import *


def outSubseqTest():
    foo = Graph.Lattice(dim=[10, 10], circular=False)    
    perm = [65, 64, 74, 84, 63, 53, 52, 51, 62, 55, 83, 72, 42, 61, 75, 93, 82, 54, 60]
    out = computeOutDegreeFromSeq(foo, perm)
    print(out)
    # perm = [0, 1, 5, 2, 4]
    # out = computeOutDegreeSubseq(foo, perm, out[1], 2, 3)
    
def smallLatticeTest():
    foo = Graph.Lattice(dim=[4, 4], circular=False)    
    n_inf = 10
    q = 1    
    n = len(foo.vs)
    
    first = choices(list(range(n)), foo.degree())[0]
    true_order = simulateInfection(foo, first, n_inf, q)
    
    freq = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.abs(np.sum(freq) - 1) < 1e-3
    freq2 = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.sum(np.abs(freq - freq2))/2 < 0.1
    
def smallNoisyLatticeTest():
    foo = Graph.Lattice(dim=[4, 4], circular=False)    
    n_inf = 10
    q = 0.3    
    n = len(foo.vs)
    
    first = choices(list(range(n)), foo.degree())[0]
    true_order = simulateInfection(foo, first, n_inf, q)
    
    freq = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.abs(np.sum(freq) - 1) < 1e-3
    freq2 = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.sum(np.abs(freq - freq2))/2 < 0.2
    
def largeLatticeTest():
    foo = Graph.Lattice(dim=[4, 4], circular=False)    
    n_inf = 10
    q = 0.95   
    n = len(foo.vs)
    
    first = choices(list(range(n)), foo.degree())[0]
    true_order = simulateInfection(foo, first, n_inf, q)
    
    freq = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.abs(np.sum(freq) - 1) < 1e-3
    freq2 = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.sum(np.abs(freq - freq2))/2 < 0.1
    
def largeNoisyLatticeTest():
    foo = Graph.Lattice(dim=[40, 40], circular=False)    
    n_inf = 100
    q = 0.7
    n = len(foo.vs)
    
    first = choices(list(range(n)), foo.degree())[0]
    true_order = simulateInfection(foo, first, n_inf, q)
    
    freq = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.abs(np.sum(freq) - 1) < 1e-3
    freq2 = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.sum(np.abs(freq - freq2))/2 < 0.1

def erdosRenyiTest():
    foo = Graph.Erdos_Renyi(n=50, m=250)
    n_inf = 20
    q = 0.95   
    n = len(foo.vs)
    
    first = choices(list(range(n)), foo.degree())[0]
    true_order = simulateInfection(foo, first, n_inf, q)
    
    freq = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.abs(np.sum(freq) - 1) < 1e-3
    freq2 = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.sum(np.abs(freq - freq2))/2 < 0.15
    
def airportTest():
    foo = igraph.read("data/global-net.dat")
    foo.delete_vertices(0)
    foo.to_undirected()
    foo.simplify()
    
    n_inf = 50
    q = 0.9   
    n = len(foo.vs)
    
    first = choices(list(range(n)), foo.degree())[0]
    true_order = simulateInfection(foo, first, n_inf, q)
    
    freq = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.abs(np.sum(freq) - 1) < 1e-3
    freq2 = inferInfection(foo, q, max_iters=1000000, M_trans=200, M_burn=100, k=4)
    assert np.sum(np.abs(freq - freq2))/2 < 0.2

outSubseqTest()