# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 15:23:17 2021

@author: gualandi
"""

import networkx as nx
from random import seed

# Versione della libreria
__version__ = '0.1.0'

# Lista delle funzioni che voglio esportare
__all__ = [
    'SafeFloor', 'Cubic', 'BuildRandomGraph', 'PlotGraph', 'PlotMatching', 
    'MaxMatching'
]


from math import floor
def SafeFloor(x, eps=0.01):
    """ Return a safe floor of 'x' with precision 'eps' """
    if abs(int(round(x))-x) < eps:
        return int(round(x))
    else:
        return int(floor(x))


def Cubic(n, _seed=17):
    """ Return a cubic graph on 'n' nodes """
    seed(_seed)
    return nx.generators.random_graphs.random_regular_graph(3, n)


def BuildRandomGraph(n, d=0.5, s=13):
    """ Return a random graph with 'n' nodes and expected density 'd' """
    return nx.erdos_renyi_graph(n, d, seed=s)


def PlotGraph(G):
    """ Plot graph 'G' with a circular layout """
    nx.draw(G, pos=nx.circular_layout(G))


def PlotMatching(G, Vs, Es):
    """ Plot in red the total matching given by 'Vs' and 'Es' on graph 'G' """
    H = nx.Graph()

    for u,v in G.edges():
        if (u,v) in Es or (v,u) in Es:
            H.add_edge(u, v, color='r', weight=4)
        else:
            H.add_edge(u, v, color='b', weight=1)

    vcols = []

    for u in H.nodes():
        if u in Vs:
            vcols.append('r')
        else:
            vcols.append('b')

    colors = [H[u][v]['color'] for u,v in H.edges()]
    weights = [H[u][v]['weight'] for u,v in H.edges()]

    nx.draw(H, pos=nx.circular_layout(G), node_color=vcols,
            edge_color=colors, width=weights, with_labels=True)


def ParseGraph(filename, path=''):
    """ Parse a graph from the House of Graph 
        (https://hog.grinvin.org/ViewGraphInfo.action?id=6728) 
    """
    fh = open(path+filename, 'r', encoding="utf-8")
    
    G = nx.Graph()
    
    for row in fh:
        if len(row) > 1:
            line = list(map(int, row.replace('\n','').replace(':','').split(' ')))        
            i = line[0]
            for j in line[1:]:
                G.add_edge(i-1, j-1)        
        
    return G, filename.split('.')[0]


def MaxMatching(G):
    """ Return a maximum matching in 'G' """
    return nx.algorithms.matching.max_weight_matching(G, maxcardinality=True) 
