# !/usr/bin/python
# 1/23/18
import numpy as np
import math
from collections import defaultdict


# ball_drop_kronecker takes as input
#     P : an n-by-n initiator matrix stored as a 2-D numpy array
#     k : the Kronecker power     
# and returns a set of edges that correspond with a realization
# of a Kronecker graph by using a ball-dropping algorithm after generating
# a fixed number of edges
def ball_drop_kronecker(P, k, num_edges = -1):
    if num_edges < 0:
        num_edges = expected_edges(P,k);
    edges = set();
    while (len(edges) < num_edges):
        loc = ball_drop_edge(norm(P), k);
        edges.add(loc);
    return edges;

# ball_drop_edge takes as input
#   P: an normalized initiator matrix (probability values sum to 1) stores as a
#       2-D numpy array
#   k: the kronecker power
# and returns a single edge generated using the ball-dropping method
def ball_drop_edge(P, k):
    n = len(P);
    n2 = n*n;
    P = P.flatten(); #row major order
    x=0; y=0;
    offset = 1;
    for i in range(k):
        ind = np.random.choice(n2, 1, p=P)[0];
        # update row and col value based on row major order
        row = ind // n; col = ind % n;
        x += row*offset
        y += col*offset;
        offset *= n;
    return (int(x),int(y));
               
# expected_edges calculates the expected number of edges for a kronecker graph
def expected_edges(P, k):   
    return int(math.pow(np.sum(P), k));

# norm(P) normalizes a 2-D numpy array  
def norm(P):
    total = np.sum(P)
    for i in range(P.shape[0]):
        for j in range (P.shape[1]):
            P[i,j] = (P[i,j])/total;
    return P;


N = 10000
A = np.zeros((8,8))
K = np.array([[0.99,0.5],[0.5,0.2]])
for t in range(N):
    for e in ball_drop_kronecker(K, 3):
        A[e[0],e[1]] += 1
#print(A/N - np.kron(np.kron(np.array(K),np.array(K)),np.array(K)))
B = np.kron(np.kron(K,K),K);
for i in range(8):
    A[i,i] = 0.0
    B[i,i] = 0.0
print(np.max(abs(A/N - B)))
    
