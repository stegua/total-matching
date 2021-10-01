# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 18:30:01 2021

@author: gualandi
"""

# Versione della libreria
__version__ = '0.1.0'

# Lista delle funzioni che voglio esportare
__all__ = [
    'SepClique', 'OddClique', 'Sep2k3Cycle', 'ConflictGraph'
]


from gurobipy import Model, quicksum, GRB

INT_TOL = 1e-06

#-----------------------------------------------------------------------------
class SepClique(object):
    def __init__(self, G):
        """ C'tor for maximal clique separator """
        mod = Model()
        mod.setParam(GRB.Param.OutputFlag, 0)
    
        mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
        mod.setParam(GRB.Param.Method, 1)
        mod.setParam(GRB.Param.IntFeasTol, INT_TOL)
    
        # Create variables
        x = {}
        for i in G.nodes():
            x[i] = mod.addVar(obj=0.0, vtype=GRB.BINARY)
    
        mod.update()
    
        # Add constraints
        for i in G.nodes():
            for j in G.nodes():
                if i < j and (i, j) not in G.edges():
                    mod.addConstr(x[i] + x[j] <= 1)

        self.mod = mod
        self.x = x
        
    def solve(self, xbar):
        """ Solve single separation problem """
        for v in xbar:
            self.x[v].setAttr(GRB.Attr.Obj, max(0.0, xbar[v]))
                        
        self.mod.optimize()

        if self.mod.Status != GRB.OPTIMAL:
            return 0, []
   
        xbar = [v for v in self.x if self.x[v].X > 0.5]
        
        return self.mod.getAttr(GRB.Attr.ObjVal), xbar


#-----------------------------------------------------------------------------
class OddClique(object):
    def __init__(self, G):
        """ C'tor for odd clique separator """
        mod = Model()
        mod.setParam(GRB.Param.OutputFlag, 0)
    
        mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
        mod.setParam(GRB.Param.Method, 1)
        mod.setParam(GRB.Param.IntFeasTol, INT_TOL)
    
        # Create variables
        x = {}
        for i in G.nodes():
            x[i] = mod.addVar(obj=0.0, vtype=GRB.BINARY)
    
        y = {}
        for e in G.edges():
            y[e] = mod.addVar(obj=0.0, vtype=GRB.BINARY)
    
        z = mod.addVar(obj=-1, lb=0, ub=len(G.nodes()) // 2, vtype=GRB.INTEGER)
    
        mod.update()
    
        # Add parity constraints
        mod.addConstr(quicksum(x[v] for v in x) == 2 * z)
    
        # Add clique constraints
        for i in G.nodes():
            for j in G.nodes():
                if i < j and (i, j) not in G.edges():
                    mod.addConstr(x[i] + x[j] <= 1)
    
        for e in y:
            i, j = e
            mod.addConstr(x[i] >= y[e])
            mod.addConstr(x[j] >= y[e])
            mod.addConstr(x[i] + x[j] <= 1 + y[e])

        self.mod = mod
        self.x = x
        self.y = y
        self.z = z
        
        
    def solve(self, xbar, ybar):
        """ Solve single separation problem """
        for v in xbar:
            self.x[v].setAttr(GRB.Attr.Obj, max(0.0, xbar[v]))
        
        for e in ybar:
            self.y[e].setAttr(GRB.Attr.Obj, max(0.0, ybar[e]))
                        
        self.mod.optimize()

        if self.mod.Status != GRB.OPTIMAL:
            return 0, []
   
        xbar = [v for v in self.x if self.x[v].X > 0.5]
        ybar = [e for e in self.y if self.y[e].X > 0.5]
        
        return self.mod.getAttr(GRB.Attr.ObjVal), xbar, ybar, self.z.X


#-----------------------------------------------------------------------------
class Sep2k3Cycle(object):
    def __init__(self, G):
        """ C'tor for odd clique separator """
        mod = Model()
        mod.setAttr(GRB.Attr.ModelSense, GRB.MAXIMIZE)
        mod.setParam(GRB.Param.OutputFlag, 0)
        n = len(G.nodes)
        m = len(G.edges)
    
        x = {}
        for v in G.nodes():
            x[v] = mod.addVar(obj=0.0, vtype=GRB.BINARY)
    
        y = {}
        for e in G.edges():
            i, j = e
            if j < i:
                e = j, i
            y[e] = mod.addVar(obj=0.0, vtype=GRB.BINARY)
    
        # Cardinality constraints
        w = mod.addVar(ub=2, lb=1, vtype=GRB.INTEGER)
        z = mod.addVar(ub=n + m, lb=1, vtype=GRB.INTEGER)
        k = mod.addVar(obj=-1, ub=n + m, lb=1, vtype=GRB.INTEGER)
    
        # Flow variables
        f = {}
        for e in G.edges():
            i, j = e
            f[(i, j)] = mod.addVar(ub=n, lb=-n, vtype=GRB.CONTINUOUS)
            f[(j, i)] = mod.addVar(ub=n, lb=-n, vtype=GRB.CONTINUOUS)
    
        # Vertex-origin-amount flow variables
        s = {}
        u = {}
        for v in G.nodes():
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
        for i in G.nodes():
            mod.addConstr(u[i] <= n * s[i])
    
        # Capacity constraints
        for e in y:
            i, j = e
            mod.addConstr(f[(i, j)] <= n * y[e])
            mod.addConstr(f[(j, i)] <= n * y[e])
            mod.addConstr(f[(i, j)] >= -n * y[e])
            mod.addConstr(f[(j, i)] >= -n * y[e])
    
        # Flow preservation principle \sum(out)-\sum(in)+x_v = 0
        for v in G.nodes():
            A_in = []
            A_out = []
            for e in G.edges(v):
                i, j = e
                A_in.append((j, i))
                A_out.append((i, j))
            mod.addConstr(u[v] + quicksum(f[e] for e in A_in) -
                          quicksum(f[e] for e in A_out) - x[v] == 0)
    
        # Degree constraints
        for v in G.nodes():
            A = []
            for e in G.edges(v):
                i, j = e
                if j < i:
                    A.append((j, i))
                else:
                    A.append((i, j))
            mod.addConstr(quicksum(y[e] for e in A) == 2 * x[v])

        self.mod = mod
        self.x = x
        self.y = y
        
        
    def solve(self, xbar, ybar):
        """ Solve single separation problem """
        for v in xbar:
            self.x[v].setAttr(GRB.Attr.Obj, max(0.0, xbar[v]))
        
        for e in ybar:
            self.y[e].setAttr(GRB.Attr.Obj, max(0.0, ybar[e]))
                        
        self.mod.optimize()

        if self.mod.Status != GRB.OPTIMAL:
            return 0, []
   
        xbar = self.mod.getAttr('x', self.x)
        ybar = self.mod.getAttr('x', self.y)
        
        return self.mod.getAttr(GRB.Attr.ObjVal), xbar, ybar


#-----------------------------------------------------------------------------
class ConflictGraph(object):
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
