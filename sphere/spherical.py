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
#    t0=process_time()
    lmax = int(input("Enter highest order Bessel function to calculate: "))     # Maximum Bessel function order
    nmax = int(input("Enter highest order Bessel zero to calculate: "))     # Maximum Bessel zero taken at lmax
    Bmax = bMax(lmax,nmax)      # Bessel zero at maximum lmax and nmax
    print("Calculating and creating arrays...")
    Bzeros = bZeros(Bmax,lmax,nmax)     # Array of Bessel function zeros
    list1 = orderPair(Bzeros)
    list2 = arrayAnlyt(int(list1[0][-1]))
    scatterLog(list1[0],list1[1],list2[0],list2[1])  # Creates logarthmic plots of the numerical and analytical solutions
    plt.show()
#    plt.show(block=False) 
#    t1=process_time()
#    print("Time elapsed:", t1-t0, "seconds")
#    print("Job has completed successfully... Fare thee well")

# Creates initial array of Bessel zeros of the lmax function up to the nmax zero and returns the largest value, Bmax
def bMax(lmax,nmax):
    sph_jn_max = []     
    sph_jn_max.append(zeros.main(lmax,nmax))
    sph_jn_max = np.array(sph_jn_max)
    return sph_jn_max[0][-1]   # Largest Bessel function zero

# Creates the original array of Bessel zeros 
def bZeros(Bmax,lmax,nmax):
    sph_jn = []
    l = lmax 
    n = nmax
    while l > -1 :   # Loop calculating which elements from the Bessel zeros array to include
        temp = []
        for i in range(len(zeros.main(l,n))):
            if zeros.main(l,n)[i] <= Bmax:     # Condition to assure only values less than or equal to Bmax are included
                temp.append(zeros.main(l,n)[i])
        temp = np.array(temp)
        sph_jn.append(temp)
        l = l - 1
        n = n + 5   # Tested nmax increase to ensure no state is left out
    sph_jn = np.array(sph_jn)
    return sph_jn   # 1D Array of Bessel function zeros from largest function to smallest function

# Creates an array with the Bessel zero, function number, and zero number
def orderPair(Bzeros):
    Bln = []
    for i in range(len(Bzeros)-1,-1,-1):   # Decreasing loop going through the 1D array of bessel functions
        for j in range(0,len(Bzeros[i])):
            Bln.append(np.array([Bzeros[i][j],abs(i-(len(Bzeros)-1)),j+1]))
    Bln = np.array(Bln)
    Bln = Bln[Bln[:,0].argsort()]   # Sorted list based on Bessel zero values from smallest to largest
    FS = filledStates(Bln)   # Call to the filledStates function
    FE = fermiEnergy(Bln)   # Call to the fermiEnergy function
    return [FS,FE]   # 2D Array of the filled state with its associated Fermi energy

# Calculates the number of filled states based on 2(Sum_lmax(Sum_nmax(2l+1)))
def filledStates(ln):
    fs = []
    for i in range(len(ln)):   # Loop over the l and n values
        l = (2*ln[i][1])+1   # Application of the degeneracies
        n = l * ln[i][2]   # Accounting for spin-1/2 particles
        if len(fs) == False:
            fs.append(int(2*n))
        else:
            fs.append((2*n) + fs[-1])   # Adds all states below it
    fs = np.array(fs)
    return fs   # 1D array of total number of filled states

# Creates array of the Fermi energy with its number of filled states
def fermiEnergy(Z):
    fermi = []
    for i in range(len(Z)):   # Loop over the number of Bessel zeros
        energy = float( ((4/(9*np.pi))**(2/3)) * (Z[i][0]**2))   # Calculated dimensionless Fermi energy
        fermi.append(energy)  
    fermi = np.array(fermi)
    return fermi   # 1D Array of the Fermi energy

# Analytical analysis 
def arrayAnlyt(lim):
    analytic = []
    xvals = []
    for i in range(1,lim):   # Loop to create the analytical results of N^(2/3)
        analytic.append(i**(2/3))
        xvals.append(i)
    return [xvals,analytic]  

# Logarithmic plots
def scatterLog(x1,y1,x2,y2):
    plot2 = plt.figure(1,figsize=(6,6))
    plt.xlabel("N")
    plt.ylabel("Dimensionless Fermi Energy")
    plt.scatter(x1,y1,s=0.5)
    plt.plot(x2,y2,'--')
    plt.xscale('log')
    plt.yscale('log')
    return plot2   # Returns the info to make a logarithmic plot of numerical results

main()
