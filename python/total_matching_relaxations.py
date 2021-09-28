# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 18:30:01 2021

@author: gualandi
"""

from gurobipy import Model, quicksum, GRB
from graph_tools import SafeFloor, BuildRandomGraph, PlotGraph, Cubic

import logging

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
    return obj, xbar


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
                sep, clique = SepClique(G, xbar)
                if sep > 1.01:
                    # print('clique', sep)
                    mod.addConstr(quicksum(x[v] for v in clique) <= 1)
                else:
                    break
            if method == 'odd-clique':
                sep, xb, yb, zb = OddClique(G, xbar, ybar)
                if sep > 0.01:
                    # print('oddclique', sep)
                    mod.addConstr(
                        quicksum(x[v]
                                 for v in xb) + quicksum(y[e]
                                                         for e in yb) <= zb)
                else:
                    break
            if method == '2k3-cycle':
                ObjC, xc, yc = Sep2k3Cycle(G, xbar, ybar)
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

    print(sum(x[v].X for v in x), sum(y[e].X for e in y))
    return obj, xbar


#-----------------------------------------------------------------------------
def SepClique(G, zbar):
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)

    # Create variables
    x = {}
    for i in zbar:
        x[i] = mod.addVar(obj=zbar[i] + 0.0000001, vtype=GRB.BINARY)

    mod.update()

    # Add constraints
    for i in G.nodes():
        for j in G.nodes():
            if i < j and (i, j) not in G.edges():
                mod.addConstr(x[i] + x[j] <= 1)

    mod.optimize()

    if mod.Status != GRB.OPTIMAL:
        return 0, []

    obj = mod.getAttr(GRB.Attr.ObjVal)
    xbar = [v for v in x if x[v].X > 0.5]

    return obj, xbar


#-----------------------------------------------------------------------------
def OddClique(G, xbar, ybar):
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)

    # Create variables
    x = {}
    for i in xbar:
        x[i] = mod.addVar(obj=xbar[i] + 0.000001, vtype=GRB.BINARY)

    y = {}
    for e in ybar:
        y[e] = mod.addVar(obj=ybar[e] + 0.000001, vtype=GRB.BINARY)

    z = mod.addVar(obj=-1, lb=0, ub=len(xbar) // 2, vtype=GRB.INTEGER)

    mod.update()

    # Add parity constraints
    mod.addConstr(quicksum(x[v] for v in xbar) == 2 * z)

    # Add clique constraints
    for i in G.nodes():
        for j in G.nodes():
            if i < j and (i, j) not in G.edges():
                mod.addConstr(x[i] + x[j] <= 1)

    for e in ybar:
        i, j = e
        mod.addConstr(x[i] >= y[e])
        mod.addConstr(x[j] >= y[e])
        mod.addConstr(x[i] + x[j] <= 1 + y[e])

    mod.optimize()

    if mod.Status != GRB.OPTIMAL:
        return 0, []

    obj = mod.getAttr(GRB.Attr.ObjVal)
    xbar = [v for v in x if x[v].X > 0.5]
    ybar = [e for e in y if y[e].X > 0.5]

    return obj, xbar, ybar, z.X


#-----------------------------------------------------------------------------
def Sep2k3Cycle(G, xbar, ybar):
    mod = Model()
    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.OutputFlag, 0)
    n = len(G.nodes)
    m = len(G.edges)

    x = {}
    for v in G.nodes:
        x[v] = mod.addVar(obj=xbar[v] + 0.000001, vtype=GRB.BINARY)

    y = {}
    for e in G.edges:
        i, j = e
        if j < i:
            e = j, i
        y[e] = mod.addVar(obj=ybar[e] + 0.000001, vtype=GRB.BINARY)

    # Cardinality constraints
    w = mod.addVar(ub=2, lb=1, vtype=GRB.INTEGER)
    z = mod.addVar(ub=n + m, lb=1, vtype=GRB.INTEGER)
    k = mod.addVar(obj=-1, ub=n + m, lb=1, vtype=GRB.INTEGER)

    # Flow variables
    f = {}
    for e in G.edges:
        i, j = e
        f[(i, j)] = mod.addVar(ub=n, lb=-n, vtype=GRB.CONTINUOUS)
        f[(j, i)] = mod.addVar(ub=n, lb=-n, vtype=GRB.CONTINUOUS)

    # Vertex-origin-amount flow variables
    s = {}
    u = {}
    for v in G.nodes:
        s[v] = mod.addVar(ub=1, lb=0, vtype=GRB.BINARY)
        u[v] = mod.addVar(ub=n, lb=0, vtype=GRB.INTEGER)

    mod.update()

    # Length of the cycle
    mod.addConstr(
        quicksum(x[v] for v in x) + quicksum(y[e] for e in y) == 6 * z + 2 * w)

    # Floor modeled in the objective function
    mod.addConstr(k <= 2 * z + w * (2 / 3))
    mod.addConstr(k + 1 >= 2 * z + w * (2 / 3) + INT_TOL)

    # Choice of the origin for the external flow
    mod.addConstr(quicksum(s[i] for i in G.nodes) == 1)
    for i in G.nodes:
        mod.addConstr(u[i] <= n * s[i])

    # Capacity constraints
    for e in y:
        i, j = e
        mod.addConstr(f[(i, j)] <= n * y[e])
        mod.addConstr(f[(j, i)] <= n * y[e])
        mod.addConstr(f[(i, j)] >= -n * y[e])
        mod.addConstr(f[(j, i)] >= -n * y[e])

    # Flow preservation principle \sum(out)-\sum(in)+x_v = 0
    for v in G.nodes:
        A_in = []
        A_out = []
        for e in G.edges(v):
            i, j = e
            A_in.append((j, i))
            A_out.append((i, j))
        mod.addConstr(u[v] + quicksum(f[e] for e in A_in) -
                      quicksum(f[e] for e in A_out) - x[v] == 0)

    # Degree constraints
    for v in G.nodes:
        A = []
        for e in G.edges(v):
            i, j = e
            if j < i:
                A.append((j, i))
            else:
                A.append((i, j))
        mod.addConstr(quicksum(y[e] for e in A) == 2 * x[v])

    # Solve
    mod.optimize()

    xbar = mod.getAttr('x', x)
    ybar = mod.getAttr('x', y)

    return mod.getAttr(GRB.Attr.ObjVal), xbar, ybar


#-----------------------------------------------------------------------------
def ConflictGraph(G, xbar, ybar):
    # Build conflict graph
    H = nx.Graph()

    for v in xbar:
        H.add_node(v, weight=xbar[v])

    for e in ybar:
        H.add_node(e, weight=ybar[e])

    for e in G.edges():
        i, j = e
        if j < i:
            i, j = j, i
        H.add_edge(i, j)
        H.add_edge(e, j)
        H.add_edge(e, i)

    for v in G.nodes():
        for e in G.edges(v):
            for f in G.edges(v):
                if e < f:
                    H.add_edge(e, f)

    # Maximal cliques on the conflict graph
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag, 0)

    mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
    mod.setParam(GRB.Param.Method, 1)
    mod.setParam(GRB.Param.IntFeasTol, INT_TOL)

    # Create variables
    x = {}
    for i in H.nodes():
        x[i] = mod.addVar(obj=H[i]['weight'] + 0.0000001, vtype=GRB.BINARY)

    mod.update()

    # Add constraints
    for i in H.nodes():
        for j in H.nodes():
            if i < j and (i, j) not in H.edges():
                mod.addConstr(x[i] + x[j] <= 1)

    mod.optimize()

    if mod.Status != GRB.OPTIMAL:
        return 0, []

    obj = mod.getAttr(GRB.Attr.ObjVal)
    xbar = [v for v in x if x[v].X > 0.5]
    print(xbar)

    return obj, xbar


#-----------------------------------------------
# MAIN function
#-----------------------------------------------
if __name__ == "__main__":

    G = Cubic(76)
    # G = BuildRandomGraph(50, 0.1)

    # PlotGraph(G)

    nu1, it1 = TotalMatchingOneAtTime(G, 'clique')
    nu2, it2 = TotalMatchingOneAtTime(G, '2k3-cycle')
    nu3, it3 = TotalMatchingOneAtTime(G, 'odd-clique')
    nu4, it4 = TotalMatchingRel(G)
    logging.info(" v(G) = {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} {}".format(
        nu1, nu2, nu3, nu4, it1, it2, it3, it4))
