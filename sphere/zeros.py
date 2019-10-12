#! /usr/bin/env python
# Code from https://scipy-cookbook.readthedocs.io/items/SphericalBesselZeros.html

### recursive method: computes zeros ranges of Jn(r,n) from zeros of Jn(r,n-1)
### (also for zeros of (rJn(r,n))')
### pros : you are certain to find the right zeros values;
### cons : all zeros of the n-1 previous Jn have to be computed;
### note : Jn(r,0) = sin(r)/r

from scipy import arange, pi, sqrt, zeros
from scipy.special import jv, jvp
from scipy.optimize import brentq
from sys import argv
from pylab import *

def Jn(r,n):
  return (sqrt(pi/(2*r))*jv(n+0.5,r))
def Jn_zeros(n,nt):
  zerosj = zeros((n+1, nt))
  zerosj[0] = arange(1,nt+1)*pi
  points = arange(1,nt+n+1)*pi
  racines = zeros(nt+n)
  for i in range(1,n+1):
    for j in range(nt+n-i):
      foo = brentq(Jn, points[j], points[j+1], (i,))
      racines[j] = foo
    points = racines
    zerosj[i][:nt] = racines[:nt]
  return (zerosj)

def main(lmax,nmax):
    n = int(lmax)  # n'th spherical bessel function
    nt = int(nmax) # number of zeros to be computed

    jnz = Jn_zeros(n,nt)[n]
    return jnz

