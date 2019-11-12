#!/usr/bin/

#Packages
import sys
import math
import numpy as np
from scipy import special, optimize
import matplotlib.pyplot as plt
from time import process_time

def main():
    mmax = int(input("Enter highest order Bessel function to calculate: "))     # Maximum Bessel function order
    nmax = int(input("Enter highest order Bessel zero to calculate: "))     # Maximum Bessel zero taken at mmax
    Bzeros = bZeros(mmax,nmax) 
    list1 = orderPair(Bzeros)
    print(list1)

# Creates the original array of Bessel zeros 
def bZeros(mmax,nmax):
    temp = []
    temp.append(special.jn_zeros(mmax,nmax))
    Bmax = temp[0][-1]
    sph_jn = []
    m = mmax
    n = nmax + 1
    j = 0
    while m > -1 :   # Loop calculating which elements from the Bessel zeros array to include
        temp = []
        k = special.jn_zeros(m,n)
        for i in range(len(k)):
            if k[i] <= Bmax:     # Condition to assure only values less than or equal to Bmax are included
                temp2 = [m,k[i]]
                temp.append(np.array(temp2))
        temp = np.array(temp)
        sph_jn.append(temp)
        m = m - 1
        j = j + 1
        if j==2:
            n = n + 1   # Tested nmax increase to ensure no state is left out
            j = 0
    m = mmax + 1
    j = 0
    while n > 0:   # Loop calculating which elements from the Bessel zeros array to include
        temp = []
        l = special.jn_zeros(m,n)
        for i in range(len(l)):
            if l[i] <= Bmax :     # Condition to assure only values less than or equal to Bmax are included
                temp2 =[m,l[i]]
                temp.append(np.array(temp2))
        temp = np.array(temp)
        if len(temp) > 0 :
            sph_jn.insert(0,temp)
        m = m + 1
        j = j + 1
        if j==2:
            n = n - 1
            j = 0
    sph_jn = np.array(sph_jn)
    return sph_jn   # 1D Array of Bessel function zeros from largest function to smallest function

# Creates an array with the Bessel zero, function number, and zero number
def orderPair(Bzeros):
    Bmn = []
    for i in range(len(Bzeros)-1,-1,-1):   # Decreasing loop going through the 1D array of bessel functions
        for j in range(0,len(Bzeros[i])):
            Bmn.append(np.array([Bzeros[i][j,-1],abs(i-(len(Bzeros)-1)),j+1]))
    Bmn = np.array(Bmn)
    Bmn = Bmn[Bmn[:,0].argsort()]   # Sorted list based on Bessel zero values from smallest to largest
    FS = filledStates(Bmn)   # Call to the filledStates function
    FE = fermiEnergy(Bmn)   # Call to the fermiEnergy function
    return [FS]   # 2D Array of the filled state with its associated Fermi energy

# Calculates the number of filled states based on 2(Sum_lmax(Sum_nmax(2l+1)))
def filledStates(mn):
    fs = []
    for i in range(len(mn)):   # Loop over the m values
        if len(fs) == False:
            if mn[i][1]==0:
                fs.append(2)
            else:
                fs.append(4)
        else:
            if mn[i][1]==0:
                fs.append(2 + fs[-1])   # Adds all states below it
            else:
                fs.append(4 + fs[-1])
    fs = np.array(fs)
    return fs   # 1D array of total number of filled states

# Creates array of the Fermi energy
def fermiEnergy(Z):
    fermi = []
    for i in range(len(Z)):   # Loop over the number of Bessel zeros
        for j in range(nz):
            energy = float( (((np.pi ** 2)*(a ** 2) / (3 * (L ** 2))) ** (2/3)) * (Z[i][0] ** 2 + ((np.pi * a * j) / L) ** 2))   # Calculated dimensionless Fermi energy
        fermi.append(energy)  
    fermi = np.array(fermi)
    return fermi   # 1D Array of the Fermi energy

main()
