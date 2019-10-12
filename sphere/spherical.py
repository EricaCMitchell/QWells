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
    Bmax = BMax(lmax,nmax) 

def BMax(lmax,nmax):
    sph_jn_max = []
    sph_jn_max.append(zeros.main(lmax,nmax))
    sph_jn_max = np.array(sph_jn_max)
    return sph_jn_max[-1][-1]

def BZeros(Bmax,lmax):
    sph_jn = []

main()
