# !/usr/bin/python
# tested and works as of july
import sys
import itertools
import numpy
#import pylab as plt
from collections import Counter
from math import factorial
from functools import reduce
import numpy
import math
import timeit
import networkx as nx
from collections import defaultdict

# sample_kronecker_fast takes as input
#     K : an n-by-n initiator matrix stored as a list-of-lists
#            P[i][j] is the i, jth element. 
#     k : the Kronecker power     
# and returns a set of edges that correspond with a realization
# of a Kronecker graph
def sample_kronecker_fast(K,k):
    array=[]
    # creates the first Erdos-Renyi subregion,
    # an array with k 0's
    n=len(K);
    v = [K[i][j] for j in range(n) for i in range(n)] # vectorize by cols
    edges_mult = []
    for r in regions(v,k):            # for each region of the mult. table
        edges_mult.extend(grass_hop_region(r, v)) # get edges in mult. table
    edges_kron = []
    for e in edges_mult:              # map edges from mult. table to kron
        edges_kron.append(map_mult_to_kron(e, n))
    return edges_kron;

# update takes as input:
#    n: the dimension of the n x n initiator matrix
#    array: the sequence to be updated
# and returns the next non-decreasing sequence
def update(n, array):
    max = math.pow(n,2)-1;
    if(len(array)==1):
        if(array[0]==max):
            array[0]=-1
            return array
        array[0]=array[0]+1
        return array
    elif(array[len(array)-1]==max):
            place=array[0:len(array)-1]
            place=update(n, place)
            last=place[len(place)-1]
            place.append(last)
            return place
    else:
        array[len(array)-1]=array[len(array)-1]+1
        return array;

# regions takes as input
#     v : a vectorized initiator matrix
#     k : the Kronecker power     
# and returns a list of lists representing all ER subregions
def regions(v,k):
    subregions = [];                     
    array=[]
    for i in range(k):
        array.append(0)
    while(array[0]!=-1):
        subregions.append(list(array));
        array = update(int(math.sqrt(len(v))),array)
    return subregions;
                      

# grass_hop_region takes as input
#    r: the subregion to be sampled represented by an array
#        in which the numbers correspond to letters of the subregion                      
#    v: an n^2 x 1 inititator matrix of probability values represented as a
#        column vector
# and returns a list of edges represented by multisets of numbers
# corresponding to letters.
def grass_hop_region(r,v):
    n=count_permutations(r);
    collection = []
    # p is the common probability value of the subregion
    p=kron(r,v)
    i=-1
    # a is the geometric random variable or length of the hop
    a=numpy.random.geometric(p)
    while(i<=n-a-1):
        i=i+a
        # finds ith multiset permutation
        thearray=unrank(r,i)
        collection.append(thearray);
        # gets the next geometric random variable or hop length
        a=numpy.random.geometric(p);
    return collection;


# kron takes as input:
#    v: an n^2 x 1 iniator matrix of probability values stored as a single
#       column-wise vector
#    r: an array with k elements specifying a cell in the k-space
# and returns the value at the specified location
def kron(r,v):
    n = int(math.sqrt(len(v)));
    final = 1;
    for val in r:
        final = final * v[val]
    return final;
    

# count_permutations takes as input:
#    counter1: any iterable list or multiset
# and returns the number of permutations of the multiset
def count_permutations(counter1):
    counter1=Counter(counter1)
    values = counter1.values()
    return (
        factorial(sum(values))/reduce(lambda a, v: a * factorial(v), values,1)
    )

# unrank takes as input:
#    C: a multiset represented by a list
#    n: the lexicographic rank or index to be found
# and returns the nth permutation of C in lexicographic order.

# Precondition: C must be in non-decreasing order
# This function is implemented recursively, so it is not
# appropriate for extremely long sequences

# Examples:
# unrank([0,1,1,3], 0) returns [0,1,1,3]
# unrank([0,1,1,3], 1) returns [0,1,3,1]
# unrank([0,1,1,3], 2) returns [0,3,1,1]

def unrank(C, n):
    counter = 0;
    # base case of recursion
    if n == 0:
        return numpy.array(C);
    
    # find the element at the start of the unranked permutation
    for i,s in enumerate(C):
        # we want to test if element s starts the
        # new permutation

        # checks for repeated elements to prevent overcounting
        if(i!=0 and s==C[i-1]):
            continue;
        # creates placeholder list and removes s at position i
        place = numpy.delete(C,i)
        # checks if s starts the new permutation. Then start
        # the recursive step to unrank the rest of the list
        if(count_permutations(place) > n-counter):
            return numpy.append(s, unrank(place, n-counter));
        # updates the counter of the total number of permutations
        counter += count_permutations(place);

# map_mult_to_kron takes as input
#    e: an array (probability sequence) consisting of a permutation of letters
#        represented as numbers 0,1,2...(n^2 - 1)
#    n: the dimension of the n x n initiator matrix
# and returns the corresponding row and column index of the probability sequence
# in the kronecker product graph matrix

def map_mult_to_kron(e,n):
    k=len(e);
    I = multindex_to_linear(e,pow(n,2));
    return morton_decode(I,n,k);  

# morton_decode takes as input:
#     I: the linear index of a multiset permutation
#     n: the dimension of the n x n initiator matrix
#     k: the Kronecker power
# and returns the row and column indices of the multiset permutation in the
# Kronecker matrix as an array with 2 elements
def morton_decode(I, n, k):
    #convert I into a base n number
    num = change_base(I,n);
    row=[]
    col=[]
    for i in range(len(num)):
        if i % 2 ==0:
            row.insert(0,num[len(num)-i-1]);
        else:
            col.insert(0,num[len(num)-i-1]);
    r = multindex_to_linear(row, n);
    c = multindex_to_linear(col, n);
    return ([r,c]);

# multindex_to_linear takes as input
#     multind: an array representing a multiset permutation in the
#              multiplication table
#     N: the number of values, n^2, in the n x n initiator matrix
# and returns the linear index of the multiset permutation
# note: this algorithm takes a base N string and converts it to base 10
def multindex_to_linear(multind, N):
    rank = 0;
    for i in range(len(multind)):
        #finds the lexicographic rank of the multiset by converting it from
        #base N to base 10
        rank += multind[i] * pow(N,len(multind)-i-1);
    return rank;

# change_base takes as input:
#     I: the linear index of a multiset permutation
#     n: the dimension of the n x n initiator matrix
# and returns an array with I expressed in base n
def change_base(I,n):
    num=[]
    while I>= n:
        num.insert(0,I%n);
        I=int(I/n);
    num.insert(0,I);
    return num;


#monte_carlo simulation
# test grass_hop_kron
# After 1000000 samples (takse a few mins)
#we'd expect about 0.001 difference based on CLT.
def monte_carlo(K,k,N):
    A = numpy.zeros((8,8))
    for t in range(N):
        for e in sample_kronecker_fast(K,k):
            A[e[0],e[1]] += 1
#print(A/N - np.kron(np.kron(np.array(K),np.array(K)),np.array(K)))
    return numpy.max(abs(A/N - numpy.kron(numpy.kron(numpy.array(K),numpy.array(K)),numpy.array(K))))

print(monte_carlo([[.99,.5],[.5,.2]],3,10000))

#print(sample_kronecker_fast([.4,.7,.2,.6],2))
#nx.draw(g);
#plt.show();


    

