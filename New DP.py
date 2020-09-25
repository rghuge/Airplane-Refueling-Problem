 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 20:57:02 2020

@author: binyao
"""

from math import *

from gurobipy import *
from itertools import *
import numpy as np
import time

soln2 = []
tm2 = []
vec2 = []
num = 3 # number of experiments

for king in range(num-1,num):
    ######################### setup instance ###################
    n = 100
    np.random.seed(king) 
    V = np.random.uniform(1000,10000,n) 
    C = np.random.uniform(1,10,n) 
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
    #print(N)
    print(eps*sum(C)/r,(1+eps)**r) # large n, small upperbound for c
    '''
    ############################# approximation ############################
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
    apx_M = [tuple([0 for i in range(r)])]
    for i in range(1,len(M)):
        tcr = sum((1+eps)**j*M[i][j] for j in range(r))
        C_star = 2**(floor(log(tcr,2)))
        T = eps*C_star/r
        apx_m = [min(floor(ceil((1+eps)**(j+1)*M[i][j]/T)*T/(1+eps)**(j+1)),N[j]) for j in range(r)]
        if apx_m not in apx_M:
            apx_M.append(tuple(apx_m))
    print(len(M),len(apx_M))      
    '''
    '''
    k = ceil(log(sum(C),1+eps))
    B = [r/eps,(2+eps)*r/eps]
    apx_M = [] #apx_M stores all the approximate vectors
    for stg in range(floor(k/2)+1):
        C_star = 2**stg
        alpha = eps*C_star/r
        cddt = [] #cddt is the candidate of all k
        for i in range(r):
            cdt = []
            for j in range(N[i]+1):
                if ceil(j*(1+eps)**(i+1)/alpha)<=B[1]:
                    cdt.append(ceil(j*(1+eps)**(i+1)/alpha))
            if len(cdt)>cdt[-1]+1:
                cdt = list(range(cdt[-1]+1))   # e.g. cdt = [0,1,1,2,2,3] -> cdt = [0,1,2,3]
            cddt.append(cdt)
        K = cddt[0]  
        for i in range(1,r):
            K = list(product(K,cddt[i]))
            count = 0
            while count < len(K):
                if i == 1:
                    tp1 = (K[count][0],)
                else:
                    tp1 = K[count][0]   
                tp2 = (K[count][1],)
                K[count] = tp1+tp2 
                count += 1
            for v in K:
                if i == r-1:
                    if sum(v)<B[0] or sum(v)>B[1]:
                        K.remove(v)
                else:
                    if sum(v)>B[1]:
                        K.remove(v)  
        if stg == floor(k/2):
            print(cddt,len(K))
        for v in K:
            apx_m = [min(floor(v[i]*alpha/(1+eps)**(i+1)),N[i]) for i in range(r)]
            if apx_m not in apx_M:
                apx_M.append(apx_m)
    print(len(apx_M))
    '''
    '''        
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