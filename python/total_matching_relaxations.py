# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 18:30:01 2021

@author: gualandi
"""

from gurobipy import Model, quicksum, GRB
from graph_tools import SafeFloor, BuildRandomGraph, Cubic, MaxMatching, ParseGraph
from graph_separators import SepClique, OddClique, Sep2k3Cycle, ConflictGraph
from graph_mips import MIP_StableSet, MIP_TotalMatching

from time import perf_counter

import logging

INT_TOL = 1e-06


def TotalMatchingRel(G, timelimit=3600):
    """ Solve Total Matching IP basic formulation """
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)
    mod.setParam(GRB.Param.TimeLimit, timelimit)

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
    t0 = perf_counter()
    t1 = timelimit
    while t1 >= 0:
        mod.optimize()
        print(it, "LB: ", mod.getAttr(GRB.Attr.ObjVal))

        if mod.Status == GRB.TIME_LIMIT:
            break
        
        if mod.Status != GRB.OPTIMAL:
            return 'NaN'

        obj = mod.getAttr(GRB.Attr.ObjVal)
        xbar = dict((v, x[v].X) for v in x)
        ybar = dict((e, y[e].X) for e in y)

        t1 = 3600 - (perf_counter() - t0)
        if t1 < 0:
            break
        
        UBk = int(SafeFloor(sum(x[v].X for v in x) + sum(y[e].X for e in y)))
        if False and UBk < UB:
            UB = UBk
            mod.addConstr(
                quicksum(x[v] for v in x) + quicksum(y[e] for e in y) <= UB)
        else:
            sep, clique = SeparatorClique.solve(xbar, t1)
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
                sep, xb, yb, zb = SeparatorOddClique.solve(xbar, ybar, t1)
                if sep > 0.01:
                    # print('oddclique', sep)
                    mod.addConstr(
                        quicksum(x[v]
                                 for v in xb) + quicksum(y[e]
                                                         for e in yb) <= zb)
                else:
                    ObjC, xc, yc = Separator2k3Cycle.solve(xbar, ybar, t1)
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
            
        it += 1


    return obj, it, perf_counter() - t0


def TotalMatchingOneAtTime(G, method, timelimit=3600):
    """ Solve Total Matching IP basic formulation """
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)
    mod.setParam(GRB.Param.TimeLimit, timelimit)

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
    t0 = perf_counter()
    t1 = timelimit
    while t1 >= 0:
        mod.optimize()
        print(it, "LB: ", mod.getAttr(GRB.Attr.ObjVal))

        if mod.Status == GRB.TIME_LIMIT:
            break

        if mod.Status != GRB.OPTIMAL:
            return 'NaN'

        obj = mod.getAttr(GRB.Attr.ObjVal)
        xbar = dict((v, x[v].X) for v in x)
        ybar = dict((e, y[e].X) for e in y)

        t1 = 3600 - (perf_counter() - t0)
        if t1 < 0:
            break
        
        UBk = int(SafeFloor(sum(x[v].X for v in x) + sum(y[e].X for e in y)))
        if False and UBk < UB:
            UB = UBk
            mod.addConstr(
                quicksum(x[v] for v in x) + quicksum(y[e] for e in y) <= UB)
        else:

            if method == 'clique':
                sep, clique = Separator.solve(xbar, t1)
                if sep > 1.01:
                    # print('clique', sep)
                    mod.addConstr(quicksum(x[v] for v in clique) <= 1)
                else:
                    break

            if method == 'odd-clique':
                sep, xb, yb, zb = Separator.solve(xbar, ybar, t1)
                if sep > 0.01:
                    # print('oddclique', sep)
                    mod.addConstr(
                        quicksum(x[v]
                                 for v in xb) + quicksum(y[e]
                                                         for e in yb) <= zb)
                else:
                    break

            if method == '2k3-cycle':
                ObjC, xc, yc = Separator.solve(xbar, ybar, t1)
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
                sep, clique = Separator.solve(xbar, ybar, t1)
                if sep > 1.01:
                    # print('clique', sep)
                    I = [v for v in clique if type(v) == int]
                    J = [v for v in clique if type(v) != int]
                    mod.addConstr(
                        quicksum(x[v] for v in I) + quicksum(y[e]
                                                             for e in J) <= 1)
                else:
                    break
        it += 1

    return obj, it, perf_counter() - t0


def TotalMatchingLB(G, timelimit=3600):
    """ Solve Total Matching IP basic formulation """
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)
    mod.setParam(GRB.Param.TimeLimit, timelimit)

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
        
    mod.optimize()

    obj = mod.getAttr(GRB.Attr.ObjVal)

    return obj


#-----------------------------------------------
# MAIN function
#-----------------------------------------------
if __name__ == "__main__":

    if False:
        logging.basicConfig(
            filename='cubicUB.log',
            level=logging.DEBUG,
            format=
            '%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        # Test su grafi cubici
        for n in [50, 60, 70, 80, 90, 100]:
            for s in [67, 71, 73, 79, 83, 101, 103, 107, 109, 113]:
                G = Cubic(n, s)
                mu = MaxMatching(G)
                al = MIP_StableSet(G)
                mt = MIP_TotalMatching(G)
                lb0 = TotalMatchingLB(G)
                nu1, it1, t1 = TotalMatchingOneAtTime(G, 'clique')
                nu2, it2, t2 = TotalMatchingOneAtTime(G, '2k3-cycle')
                nu3, it3, t3 = TotalMatchingOneAtTime(G, 'odd-clique')
                nu4, it4, t4 = TotalMatchingRel(G)
                nu5, it5, t5 = TotalMatchingOneAtTime(G, 'conflict')
                logging.info(
                    " cubic s {} n {} m {} v(G) = {}, alpha(G) = {}, vt(G) = {}, UB(G) = {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {} {} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}"
                    .format(s, len(G.nodes()), len(G.edges()), len(mu), al[0],
                            mt[0], lb0, nu1, nu2, nu3, nu4, nu5, it1, it2, it3, it4, it5, t1, t2, t3, t4, t5))

    if False:
        logging.basicConfig(
            filename='randomUB.log',
            level=logging.DEBUG,
            format=
            '%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        # Test su grafi random
        for n in [80]:
            for d in [0.05, 0.1, 0.15, 0.2, 0.25]:
                for s in [67, 71, 73, 79, 83, 101, 103, 107, 109, 113]:
                    G = BuildRandomGraph(n, d, s)
                    mu = MaxMatching(G)
                    al = MIP_StableSet(G)
                    mt = MIP_TotalMatching(G)
                    lb0 = TotalMatchingLB(G)
                    nu1, it1, t1 = TotalMatchingOneAtTime(G, 'clique')
                    nu2, it2, t2 = TotalMatchingOneAtTime(G, '2k3-cycle')
                    nu3, it3, t3 = TotalMatchingOneAtTime(G, 'odd-clique')
                    nu4, it4, t4 = TotalMatchingRel(G)
                    nu5, it5, t5 = TotalMatchingOneAtTime(G, 'conflict')
                    logging.info(
                        " random s {} n {} m {} v(G) = {}, alpha(G) = {}, vt(G) = {}, UB(G) = {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {} {} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}"
                        .format(s, len(G.nodes()), len(G.edges()), len(mu),
                                al[0], mt[0], lb0, nu1, nu2, nu3, nu4, nu5, it1, it2, it3, it4, it5, t1, t2, t3, t4, t5))

    if True:
        logging.basicConfig(
            filename='housegraphs.log',
            level=logging.DEBUG,
            format=
            '%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')

        # House of Graphs
        HoG = ['graph_6630.lst', 'graph_6710.lst', 'graph_6714.lst', 
                'graph_6720.lst', 'graph_6724.lst', 'graph_6728.lst', 
                'graph_6708.lst', 'graph_6712.lst', 'graph_6718.lst', 
                'graph_6722.lst', 'graph_6726.lst']
    
        Gs = list(map(lambda x: ParseGraph(x, "..\\data\\"), HoG))
    
        # Total matching
        for G, name in Gs:
            mu = MaxMatching(G)
            al = MIP_StableSet(G)
            mt = MIP_TotalMatching(G)
            lb0 = TotalMatchingLB(G)
            nu1, it1, t1 = TotalMatchingOneAtTime(G, 'clique')
            nu2, it2, t2 = TotalMatchingOneAtTime(G, '2k3-cycle')
            nu3, it3, t3 = TotalMatchingOneAtTime(G, 'odd-clique')
            nu4, it4, t4 = TotalMatchingRel(G)
            nu5, it5, t5 = TotalMatchingOneAtTime(G, 'conflict')
            logging.info(
                " {} n {} m {} v(G) {} alpha(G) {} vt(G) {} UB(G) {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {} {} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}"
                .format(name, len(G.nodes()), len(G.edges()), len(mu),
                        al[0], mt[0], lb0, nu1, nu2, nu3, nu4, nu5, it1, it2, it3, it4, it5, t1, t2, t3, t4, t5))
            
            
            
    if False:
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
        logging.info(
            " n {} m {} v(G) = {}, alpha(G) = {}, vt(G) = {}, UB(G) = {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {}"
            .format(len(G.nodes()), len(G.edges()), len(mu), al[0], mt[0], nu1,
                    nu2, nu3, nu4, it1, it2, it3, it4))

    logging.shutdown()
