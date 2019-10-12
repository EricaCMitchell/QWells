#!/usr/bin/

#Packages
import sys
import zeros
import math
import numpy as np
from scipy import special, optimize
import matplotlib.pyplot as plt
from time import process_time

def main():
    lmax = 10
    nmax = 10
    Bmax = bMax(lmax,nmax) 
    Bzeros = bZeros(Bmax,lmax,nmax)
    print(Bzeros)
    sumIndex = orderPair(Bzeros)
    print(sumIndex)

def bMax(lmax,nmax):
    sph_jn_max = []
    sph_jn_max.append(zeros.main(lmax,nmax))
    sph_jn_max = np.array(sph_jn_max)
    return sph_jn_max[0][-1]

def bZeros(Bmax,lmax,nmax):
    sph_jn = []
    l = lmax
    n = nmax
    while l > -1 :
        temp = []
        for i in range(len(zeros.main(l,n))):
            if zeros.main(l,n)[i] <= Bmax:
                temp.append(zeros.main(l,n)[i])
        temp = np.array(temp)
        sph_jn.append(temp)
        l = l - 1
        n = n + 2
    sph_jn = np.array(sph_jn)
    return sph_jn

def orderPair(Bzeros):
    ln = []
    for i in range(len(Bzeros)-1,0,-1):
        for j in range(0,len(Bzeros[i])):
            ln.append([abs(i-11),j+1])
    ln = np.array(ln)
    return ln

main()
