# 8/10/16
# edited 12/21/17

import sys
import pylab as plt
import networkx as nx
import random
from scipy.special import comb

# unrank takes as input
#   I: the rank of the graph
#   max_val: the number of potential edges
#   length: the number of edges
# and returns the Ith non-decreasing sequence of edges
# of size 'length' and maximum edge of value (max_val-1)
def unrank(I, max_val, length):
  array = [];
  # finds the next edge value
  for i in range(max_val):
    if length == 0:
      break;
    if length == 1:
      array.append(int(i + I));
      length=length-1;
    # checks if next edge value is i
    elif(comb(max_val-i-1, length-1) > I):
      array.append(int(i));
      length = length - 1; 
    else:
      I = I - comb(max_val-i-1, length-1);
  return array;

# fixed_er takes as input:
#   n: the number of nodes
#   m: the number of edges
# And returns a list of edges for a graph with n nodes and m edges

def fixed_er(n,m):
  graph = [];
  #computes the number of possible edges
  possible_edges = int(n*(n-1)/2);
  #computes the number of possible graphs
  possible_graphs = int(comb(possible_edges,m));

  #randomly generates a "rank" or number of the graph
  I = int(random.random()*possible_graphs)

  #finds the m edges that make up the graph
  edges = unrank(I, possible_edges, m);

  #finds the nodes of edges: a set of edges.
  for edge in edges:
    nodes = unrank(edge, n, 2);
    graph.append([nodes[0], nodes[1]]);
  return graph;
      
    
g = nx.Graph();
edges = erdos_renyi(8,11);

for e in edges:
  g.add_edge(e[0],e[1]);
nx.draw(g);
plt.show();
    
