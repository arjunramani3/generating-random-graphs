#This file contains methods to count various motifs of graphs from the networkx package. It also contains methods to plot properties of 
#graphs like the degree distribution.


#!/usr/bin/python
import sys
import itertools
import numpy
import pylab as plt
from collections import Counter
from math import factorial
from functools import reduce
import numpy
import math
import timeit
import networkx as nx
from scipy.optimize import curve_fit
from collections import defaultdict

#all graphs are from the networkx package

#input: a graph
#output: number of hairpins
def hairpin(graph):
    hairpins = 0;    
    for e in graph.degree().values():
        hairpins += int(e*(e-1)/2);
    return hairpins;

#input: a graph
#output: number of triangles
def triangle(graph):
    triangles = 0;    
    for e in nx.triangles(graph).values():
        triangles += e;
    return triangles/3;    

#graph.number_of_edges()
#returns number of edges

#input: a graph
#output: number of tripins
def tripin(graph):
    tripins = 0;    
    for e in graph.degree().values():
        tripins += int(e*(e-1)*(e-2)/6);
    return tripins;

#input: graph G, and nodes u and v
#output: boolean true or false based on whether u and v are neighbors
def is_connected(G, u, v):
    return u in G.neighbors(v)

#input: a graph
#output: the number of simple four-cycles (no diagonals)
def four_cycle(graph):
    cycles= 0;
    bob = defaultdict(set);
    for i in graph.nodes(): #iterate through all the nodes
        for j in graph.neighbors(i): #iterate through all the hairpins
            for k in graph.neighbors(i):
                if j<k and not is_connected(graph,j,k):
                    bob[(j,k)].add(i);
    for set_of_nodes in bob.values():
        for a in set_of_nodes:
            for b in set_of_nodes:
                if a<b and not is_connected(graph,a,b):
                    cycles +=1;       
    return cycles/2;

#Find max eccentricity
def diameter(G):
    Gc  = max(nx.connected_component_subgraphs(G),key=len);
    max_ecc = max(nx.eccentricity(Gc).values());
    return max_ecc

#generate graph
array=[]
g=nx.Graph();
for i in range(k):
    array.append(0)
while(array[0]!=-1):
    n=count_permutations(array)
    #Use Geometric random variable on (n,kron(array))
    for thing in get_geo_edges(n,array):
        #print(unrank_array(thing))
        if unrank_x(thing)>unrank_y(thing):
            g.add_edge(unrank_x(thing),unrank_y(thing));
    array = update(array);


stop = timeit.default_timer();
print('time='+ str(stop-start))




#plot degree distribution
plt.rcParams.update({'font.size': 18})
degs = {}
for n in g.nodes():
    deg = g.degree(n)
    if deg not in degs:
        degs[deg]=0;
    degs[deg]+=1;
items = sorted(degs.items())
plt.loglog([k for (k,v) in items], [v for (k,v) in items],marker='o', linestyle='-')
plt.ylabel('Number of Nodes');
plt.xlabel('Degree');
plt.title('Degree Distribution')
plt.xlim(0,)


#power law model
def func(x, a, b):
    return a * pow(x,-b);

xdata = numpy.array([k for (k,v) in items])
ydata = numpy.array([v for (k,v) in items])

popt, pcov = curve_fit(func, xdata, ydata)
print(popt[0],popt[1])

#r_squared calculation

y_avg = sum(ydata)/len(ydata);
ss_res = sum(pow(ydata - func(xdata,popt[0],popt[1]),2));
ss_tot = sum(pow((ydata - y_avg),2));
r_squared = 1-ss_res/ss_tot;

#plt.plot(xdata,func(xdata,popt[0],popt[1]), 'b-', label = 'fit',color='green')

print('r_square = ' + str(r_squared))
plt.show();
#nx.draw(g);
#plt.show();
