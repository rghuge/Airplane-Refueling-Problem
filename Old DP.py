#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 20:55:59 2020

@author: binyao
"""

from math import *

from gurobipy import *
from itertools import *
import numpy as np
import time

soln = []
tm = []
vec = []
num = 1 # number of experiments

for king in range(num-1,num):
    ######################### setup instance ###################
    n = 20
    np.random.seed(king) 
    V = np.random.uniform(1000,10000,n) 
    C = np.random.uniform(1,100,n) 
    start_time = time.time()
    eps = 2**(1/2)-1
    r = ceil(log(max(C),1+eps)) # number of classes
    R = [(1+eps)**i for i in range(r+1)]
    Class = [] #records the type of each plane
    N = [] # number of planes for each class
    NI = [[] for i in range(r)] # for each class provide the index of planes in decreasing volume
    for j in range(n):
        Class.append(floor(log(C[j],1+eps))) #plane -> class
        NI[floor(log(C[j],1+eps))].append(j) #class ->plane
    for i in range(r):
        N.append(Class.count(i))
    for i in range(r):
        lt = NI[i]
        l = len(lt)
        for k in range(l-1,0,-1): # Bubble Sort
            for j in range(0,k):
                if V[lt[j]]<V[lt[j+1]]:
                    lt[j],lt[j+1] = lt[j+1],lt[j]
        NI[i] = lt
    
    ######################### Preparation ######################
    M = list(range(N[0]+1)) #M stores all the vectors less than or equals to N 
    for i in range(1,r):
        l = list(range(N[i]+1))
        M = list(product(M,l))
        for j in range(len(M)):
            if i == 1:
                tp1 = (M[j][0],)
            else:
                tp1 = M[j][0] 
            tp2 = (M[j][1],)
            M[j] = tp1+tp2   
            
    ######################### Phase 1: Original DP ######################      
    MV = [0 for i in range(len(M))] #MV stores the DP value for each vector
    MPI = [0 for i in range(len(M))] #MPI stores the index of previous vector at DP
    U = [0 for i in range(len(M))] #U:upperbound
    for i in range(len(M)):
        for j in range(r):
            if M[i][j] > 0:
                MM = [x for x in M[i]]
                MM[j] -= 1
                P = tuple(MM)
                pi = M.index(P)
                if V[NI[j][P[j]]]/(U[pi]+(1+eps)**(j+1))+MV[pi] > MV[i]:
                    MV[i] = V[NI[j][P[j]]]/(U[pi]+(1+eps)**(j+1))+MV[pi]
                    MPI[i] = pi
                    U[i] = U[pi]+(1+eps)**(j+1)
    p = -1 
    MPV = [] #MPV traces all the previous vector
    while p != 0:
        MPV.append(M[p])
        #print(M[p],MV[p])
        p = MPI[p]
    MPV.append(tuple([0 for i in range(r)]))
    MPV.reverse()
    
    apx_perm = [] #permutation given by the DP
    for i in range(1,len(MPV)):
        idx = [list(MPV[i])[j]-list(MPV[i-1])[j] for j in range(r)].index(1)
        apx_perm.append(NI[idx][MPV[i-1][idx]])
    end_time = time.time() 
    
    tm.append(end_time-start_time)
    soln.append(MV[-1])  
    vec.append(len(M))
       
print("avg_time: ",sum(tm)/num,'\n')
print('avg_num_vec:',sum(vec)/num,'\n')
print('soln:', soln,'\n')

'''
######################### find exact solution ##########
start_time = time.time()
l = list(permutations(range(n))) 

opt_dist = 0
opt_perm = []
for perm in l:
    #dist = sum(V[perm[j]]/sum((1+eps)**(1+Class[perm[k]]) for k in range(j+1)) for j in range(len(perm)))
    dist = sum(V[perm[j]]/sum(C[perm[k]] for k in range(j+1)) for j in range(len(perm)))
    if dist > opt_dist:
        opt_dist = dist
        opt_perm = perm
print("exact opt permutation:", opt_perm,"\n")
print("exact opt distance:", opt_dist,"\n")
end_time = time.time() 
print("enumeration time taken = "+str(end_time-start_time),"\n") 
#print("ratio:", obj/opt_dist) # empirically at least 95% / theoretically at least 1/(1+eps)^2=82.6%
'''
