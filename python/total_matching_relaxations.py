# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 18:30:01 2021

@author: gualandi
"""

from gurobipy import Model, quicksum, GRB
from graph_tools import SafeFloor, BuildRandomGraph, PlotGraph, Cubic
from graph_separators import SepClique, OddClique, Sep2k3Cycle, ConflictGraph

import logging
import networkx as nx

logging.basicConfig(filename='match.log', level=logging.DEBUG)

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

    # Max stable set constraint
    # mod.addConstr(quicksum(x[v] for v in x) <= 7)

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
            sep, clique = SepClique(G, xbar)
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
                sep, xb, yb, zb = OddClique(G, xbar, ybar)
                if sep > 0.01:
                    # print('oddclique', sep)
                    mod.addConstr(
                        quicksum(x[v]
                                 for v in xb) + quicksum(y[e]
                                                         for e in yb) <= zb)
                else:
                    ObjC, xc, yc = Sep2k3Cycle(G, xbar, ybar)
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

    print(sum(x[v].X for v in x), sum(y[e].X for e in y))
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

    print(sum(x[v].X for v in x), sum(y[e].X for e in y))
    return obj, it


    def __init__(self, G):
        # Build conflict graph
        H = nx.Graph()
    
        idx = 0
        V = {}
        W = {}
        for v in G.nodes():
            H.add_node(idx)
            V[v] = idx
            W[idx] = v
            idx += 1
            
        for e in G.edges():
            H.add_node(idx)
            V[e] = idx
            V[e[1], e[0]] = idx
            W[idx] = e
            idx += 1
            
        for e in G.edges():
            f = e
            i, j = e
            if j < i:
                f = j, i
            H.add_edge(V[i], V[j])
            H.add_edge(V[f], V[j])
            H.add_edge(V[f], V[i])
    
        for v in G.nodes():
            for e in G.edges(v):
                for f in G.edges(v):
                    if e < f:
                        H.add_edge(V[e], V[f])
    
        # Maximal cliques on the conflict graph
        mod = Model()
        mod.setParam(GRB.Param.OutputFlag, 0)
    
        mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
        mod.setParam(GRB.Param.Method, 1)
        mod.setParam(GRB.Param.IntFeasTol, INT_TOL)
    
        # Create variables
        x = {}
        for i in H.nodes():
            # x[i] = mod.addVar(obj=H.nodes[i]['weight'] + 0.0000001, vtype=GRB.BINARY)
            x[i] = mod.addVar(obj=0.0, vtype=GRB.BINARY)
    
        mod.update()
    
        # Add constraints
        for i in H.nodes():
            for j in H.nodes():
                if i < j and (i, j) not in H.edges():
                    mod.addConstr(x[i] + x[j] <= 1)

        self.mod = mod
        self.x = x
        self.W = W
        self.V = V
        
        
    def solve(self, xbar, ybar):
        """ Solve single separation problem """
        
        for v in xbar:
            self.x[self.V[v]].setAttr(GRB.Attr.Obj, max(0.0, xbar[v]))
            
        for e in ybar:
            self.x[self.V[e]].setAttr(GRB.Attr.Obj, max(0.0, ybar[e]))
            
        self.mod.optimize()

        if self.mod.Status != GRB.OPTIMAL:
            return 0, []
   
        xbar = [self.W[v] for v in self.x if self.x[v].X > 0.5]
        
        return self.mod.getAttr(GRB.Attr.ObjVal), xbar
            
#-----------------------------------------------
# MAIN function
#-----------------------------------------------
if __name__ == "__main__":

    # G = Cubic(50)
    G = BuildRandomGraph(35, 0.5)
    # G = nx.Graph()
    # G.add_edge(0,1)
    # G.add_edge(1,2)
    # G.add_edge(0,2)
    
    # PlotGraph(G)
    # nu5, it5 = TotalMatchingOneAtTime(G, 'conflict')

    # nu1, it1 = TotalMatchingOneAtTime(G, 'clique')
    nu2, it2 = TotalMatchingOneAtTime(G, '2k3-cycle')
    # nu3, it3 = TotalMatchingOneAtTime(G, 'odd-clique')
    # nu4, it4 = TotalMatchingRel(G)
    # logging.info(" v(G) = {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {}".format(
    #     nu1, nu2, nu3, nu4, it1, it2, it3, it4))
