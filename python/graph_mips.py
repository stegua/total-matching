# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 11:32:40 2021

@author: gualandi
"""

import networkx as nx
from random import seed

# Versione della libreria
__version__ = '0.1.0'

# Lista delle funzioni che voglio esportare
__all__ = [
    'MIP_StableSet', 'MIP_TotalMatching'
]

from gurobipy import Model, quicksum, GRB

INT_TOL = 1e-06
    

def MIP_StableSet(G):
    """ Solve Max Stable Set on grapb 'G' with basic formulation """
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)

    # Create variables
    x = {}
    for i in G.nodes:
        x[i] = mod.addVar(obj=1, vtype=GRB.BINARY)

    mod.update()

    # Add constraints
    for i,j in G.edges():
        mod.addConstr(x[i] + x[j] <= 1)
        
    # Solve problem
    mod.optimize()
    
    # Check the status
    if mod.Status != GRB.OPTIMAL:
        return 'NaN'
    
    obj = mod.getAttr(GRB.Attr.ObjVal)
    xbar = [v for v in x if x[v].X > 0.5]
    
    return obj, xbar
    

def MIP_TotalMatching(G):
    """ Solve Total Matching IP basic formulation """
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)

    # Create variables
    x = {}
    for i in G.nodes:
        x[i] = mod.addVar(obj=1, vtype=GRB.BINARY)

    y = {}
    for e in G.edges:            
        i,j = e
        if j < i:
            e = j,i
        y[e] = mod.addVar(obj=1, vtype=GRB.BINARY)
        
    mod.update()

    # Add constraints
    for v in G.nodes:
        A = []
        for e in G.edges(v):
            i,j = e
            if j < i:
                A.append((j,i))
            else:
                A.append((i,j))
        mod.addConstr(quicksum(y[e] for e in A) + x[v] <= 1)
        
    for e in y:
        i, j = e
        mod.addConstr(x[i] + x[j] + y[e] <= 1)
        
    
    # Solve problem
    mod.optimize()
    
    # Check the status
    if mod.Status != GRB.OPTIMAL:
        return 'NaN'
    
    obj = mod.getAttr(GRB.Attr.ObjVal)
    xbar = [v for v in x if x[v].X > 0.5]
    ybar = [e for e in y if y[e].X > 0.5]
    
    return obj, xbar, ybar
