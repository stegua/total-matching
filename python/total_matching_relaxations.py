# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 18:30:01 2021

@author: gualandi
"""

from gurobipy import Model, quicksum, GRB
from graph_tools import SafeFloor, BuildRandomGraph, PlotGraph, Cubic, MaxMatching
from graph_separators import SepClique, OddClique, Sep2k3Cycle, ConflictGraph
from graph_mips import MIP_StableSet, MIP_TotalMatching

import logging
import networkx as nx

logging.basicConfig(filename='match.log', 
                    level=logging.DEBUG,
                    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

INT_TOL = 1e-06


def TotalMatchingRel(G):
    """ Solve Total Matching IP basic formulation """
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)

    # Create variables
    x = {}
    for i in G.nodes:
        x[i] = mod.addVar(obj=1, vtype=GRB.CONTINUOUS)

    y = {}
    for e in G.edges:
        i, j = e
        if j < i:
            e = j, i
        y[e] = mod.addVar(obj=1, vtype=GRB.CONTINUOUS)

    mod.update()

    # Add constraints
    for v in G.nodes:
        A = []
        for e in G.edges(v):
            i, j = e
            if j < i:
                A.append((j, i))
            else:
                A.append((i, j))
        mod.addConstr(quicksum(y[e] for e in A) + x[v] <= 1)

    for e in y:
        i, j = e
        mod.addConstr(x[i] + x[j] + y[e] <= 1)

    SeparatorClique = SepClique(G)
    
    SeparatorOddClique = OddClique(G)
        
    Separator2k3Cycle = Sep2k3Cycle(G)


    it = 0
    UB = len(x) + len(y)
    while it <= 1000:
        it += 1
        mod.optimize()
        print(it, "LB: ", mod.getAttr(GRB.Attr.ObjVal))

        if mod.Status != GRB.OPTIMAL:
            return 'NaN'

        obj = mod.getAttr(GRB.Attr.ObjVal)
        xbar = dict((v, x[v].X) for v in x)
        ybar = dict((e, y[e].X) for e in y)

        UBk = int(SafeFloor(sum(x[v].X for v in x) + sum(y[e].X for e in y)))
        if False and UBk < UB:
            UB = UBk
            mod.addConstr(
                quicksum(x[v] for v in x) + quicksum(y[e] for e in y) <= UB)
        else:
            sep, clique = SeparatorClique.solve(xbar)
            if sep > 1.01:
                # print('clique', sep)
                mod.addConstr(quicksum(x[v] for v in clique) <= 1)

                E = []
                for i in clique:
                    for j in clique:
                        if i < j:
                            E.append((i, j))
                if len(clique) % 2 == 1:
                    mod.addConstr(
                        quicksum(y[e] for e in E) <= len(clique) // 2)
                else:
                    mod.addConstr(
                        quicksum(x[v] for v in clique) +
                        quicksum(y[e] for e in E) <= len(clique) // 2)
            else:
                sep, xb, yb, zb = SeparatorOddClique.solve(xbar, ybar)
                if sep > 0.01:
                    # print('oddclique', sep)
                    mod.addConstr(
                        quicksum(x[v]
                                 for v in xb) + quicksum(y[e]
                                                         for e in yb) <= zb)
                else:
                    ObjC, xc, yc = Separator2k3Cycle.solve(xbar, ybar)
                    Cx = [i for i in xc if xc[i] > 0.5]
                    Cy = [j for j in yc if yc[j] > 0.5]
                    k = len(Cx)
                    if ObjC > 0.01 and k > 0:
                        # print('2k3-cycle')
                        mod.addConstr(
                            quicksum(x[c] for c in Cx) +
                            quicksum(y[c]
                                     for c in Cy) <= SafeFloor((2 * k) / 3))
                    else:
                        break

    return obj, it


def TotalMatchingOneAtTime(G, method):
    """ Solve Total Matching IP basic formulation """
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)

    # Create variables
    x = {}
    for i in G.nodes:
        x[i] = mod.addVar(obj=1, vtype=GRB.CONTINUOUS)

    y = {}
    for e in G.edges:
        i, j = e
        if j < i:
            e = j, i
        y[e] = mod.addVar(obj=1, vtype=GRB.CONTINUOUS)

    mod.update()

    # Add constraints
    for v in G.nodes:
        A = []
        for e in G.edges(v):
            i, j = e
            if j < i:
                A.append((j, i))
            else:
                A.append((i, j))
        mod.addConstr(quicksum(y[e] for e in A) + x[v] <= 1)

    for e in y:
        i, j = e
        mod.addConstr(x[i] + x[j] + y[e] <= 1)

    # mod.addConstr(quicksum(x[v] for v in x) <= 7)

    if method == 'clique':
        Separator = SepClique(G)
    
    if method == 'conflict':
        Separator = ConflictGraph(G)
    
    if method == 'odd-clique':
        Separator = OddClique(G)
        
    if method == '2k3-cycle':
        Separator = Sep2k3Cycle(G)
        
    it = 0
    UB = len(x) + len(y)
    while it <= 1000:
        it += 1
        mod.optimize()
        print(it, "LB: ", mod.getAttr(GRB.Attr.ObjVal))

        if mod.Status != GRB.OPTIMAL:
            return 'NaN'

        obj = mod.getAttr(GRB.Attr.ObjVal)
        xbar = dict((v, x[v].X) for v in x)
        ybar = dict((e, y[e].X) for e in y)

        UBk = int(SafeFloor(sum(x[v].X for v in x) + sum(y[e].X for e in y)))
        if False and UBk < UB:
            UB = UBk
            mod.addConstr(
                quicksum(x[v] for v in x) + quicksum(y[e] for e in y) <= UB)
        else:

            if method == 'clique':
                sep, clique = Separator.solve(xbar)
                if sep > 1.01:
                    # print('clique', sep)
                    mod.addConstr(quicksum(x[v] for v in clique) <= 1)
                else:
                    break
                
            if method == 'odd-clique':
                sep, xb, yb, zb = Separator.solve(xbar, ybar)
                if sep > 0.01:
                    # print('oddclique', sep)
                    mod.addConstr(
                        quicksum(x[v]
                                 for v in xb) + quicksum(y[e]
                                                         for e in yb) <= zb)
                else:
                    break
                
            if method == '2k3-cycle':
                ObjC, xc, yc = Separator.solve(xbar, ybar)
                Cx = [i for i in xc if xc[i] > 0.5]
                Cy = [j for j in yc if yc[j] > 0.5]
                k = len(Cx)
                if ObjC > 0.01 and k > 0:
                    # print('2k3-cycle')
                    mod.addConstr(
                        quicksum(x[c] for c in Cx) +
                        quicksum(y[c] for c in Cy) <= SafeFloor((2 * k) / 3))
                else:
                    break

            if method == 'conflict':
                sep, clique = Separator.solve(xbar, ybar)
                if sep > 1.01:
                    # print('clique', sep)
                    I = [v for v in clique if type(v) == int]
                    J = [v for v in clique if type(v) != int]
                    mod.addConstr(quicksum(x[v] for v in I) + quicksum(y[e] for e in J) <= 1)
                else:
                    break

    return obj, it


            
#-----------------------------------------------
# MAIN function
#-----------------------------------------------
if __name__ == "__main__":

    if False:
        # Test su grafi cubici
        for n in [80, 100]:
            for s in [67, 71, 73, 79, 83, 101, 103, 107, 109, 113]:
                G = Cubic(n, s)
                mu = MaxMatching(G)
                al = MIP_StableSet(G)
                mt = MIP_TotalMatching(G)
                nu1, it1 = TotalMatchingOneAtTime(G, 'clique')
                nu2, it2 = TotalMatchingOneAtTime(G, '2k3-cycle')
                nu3, it3 = TotalMatchingOneAtTime(G, 'odd-clique')
                nu4, it4 = TotalMatchingRel(G)
                logging.info(" cubic s {} n {} m {} v(G) = {}, alpha(G) = {}, vt(G) = {}, UB(G) = {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {}".format(
                    s, len(G.nodes()), len(G.edges()), len(mu), al[0], mt[0], nu1, nu2, nu3, nu4, it1, it2, it3, it4))
                
        # Test su grafi random
        for n in [80, 100]:
            for d in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
                for s in [67, 71, 73, 79, 83, 101, 103, 107, 109, 113]:
                    G = BuildRandomGraph(n, d, s)
                    mu = MaxMatching(G)
                    al = MIP_StableSet(G)
                    mt = MIP_TotalMatching(G)
                    nu1, it1 = TotalMatchingOneAtTime(G, 'clique')
                    nu2, it2 = TotalMatchingOneAtTime(G, '2k3-cycle')
                    nu3, it3 = TotalMatchingOneAtTime(G, 'odd-clique')
                    nu4, it4 = TotalMatchingRel(G)
                    logging.info(" random s {} n {} m {} v(G) = {}, alpha(G) = {}, vt(G) = {}, UB(G) = {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {}".format(
                        s, len(G.nodes()), len(G.edges()), len(mu), al[0], mt[0], nu1, nu2, nu3, nu4, it1, it2, it3, it4))
         
    if True:
        G = Cubic(50, 17)
        # G = BuildRandomGraph(35, 0.5)
        # G = nx.Graph()
        # G.add_edge(0,1)
        # G.add_edge(1,2)
        # G.add_edge(0,2)
        
        # PlotGraph(G)
    
        mu = MaxMatching(G)
        al = MIP_StableSet(G)
        mt = MIP_TotalMatching(G)
        nu1, it1 = TotalMatchingOneAtTime(G, 'clique')
        nu2, it2 = TotalMatchingOneAtTime(G, '2k3-cycle')
        nu3, it3 = TotalMatchingOneAtTime(G, 'odd-clique')
        nu4, it4 = TotalMatchingRel(G)
        logging.info(" n {} m {} v(G) = {}, alpha(G) = {}, vt(G) = {}, UB(G) = {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {}".format(
            len(G.nodes()), len(G.edges()), len(mu), al[0], mt[0], nu1, nu2, nu3, nu4, it1, it2, it3, it4))

    logging.shutdown()
