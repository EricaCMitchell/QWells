#!/usr/bin/

#Packages
import sys
import math
import numpy as np
from scipy import special, optimize
import matplotlib.pyplot as plt
from time import process_time
np.set_printoptions(threshold=sys.maxsize)

def main():
#    t0=process_time()
    mmax = int(input("Enter highest order Bessel function to calculate: "))     # Maximum Bessel function order
    nmax = int(input("Enter highest order Bessel zero to calculate: "))     # Maximum Bessel zero taken at mmax
    nzmax = int(input("Enter highest value of n to calculate: "))
    a = int(input("Enter radius of the cylinder: "))     # Radius of the cylinder
    L = int(input("Enter length of the cylinder: "))     # Length of the cylinder
    Bzeros = bZeros(mmax,nmax,nzmax,a,L) 
    list1 = orderPair(Bzeros,a,L,mmax,nmax,nzmax)
    list2 = arrayAnlyt(list1[0][-1])
    scatterLog(list1[0],list1[1],list2[0],list2[1])
    plt.show()
#    plt.show(block=False) 
#    t1=process_time()
#    print("Time elapsed:", t1-t0, "seconds")
#    print("Job has completed successfully... Fare thee well")
    
# Creates the original array of Bessel zeros 
def bZeros(mmax,nmax,nzmax,a,L):
    temp = []
    temp.append(special.jn_zeros(mmax,nmax))
    Bmax = temp[0][-1]
    Fmax = float( ((Bmax ** 2)/((np.pi * (a/L)) ** 2)) + (nzmax ** 2) ) 
    sph_jn = []
    m = mmax
    n = nmax + 2
    while m > -1 :   # Loop calculating which elements from the Bessel zeros array to include
        temp = []
        k = special.jn_zeros(m,n) 
        for i in range(len(k)):
            if ((k[i] ** 2) / ((np.pi * (a/L)) ** 2) + 1) <= Fmax:     # Condition to assure only values less than or equal to Bmax are included
                temp2 = [m,k[i]]     # Appends both the m value of the Bessel function and the value of the Bessel zero
                temp.append(np.array(temp2))
        temp = np.array(temp)
        if m > 0:
            m = m - 1
            sph_jn.append(temp)
        else:
            if ((k[-1] ** 2) / ((np.pi * (a/L)) ** 2) + 1)  >= Fmax:
                np.delete(temp,-1,0)
                sph_jn.append(temp)
                break
        n = n + 1 # Tested nmax increase to ensure no state is left out
    m = mmax + 1
    j = 0
    while n > 0:   # Loop calculating which elements from the Bessel zeros array to include
        temp = []
        l = special.jn_zeros(m,n)
        for i in range(len(l)):
            if ((l[i] ** 2) / ((np.pi * (a/L)) ** 2) + 1) <= Fmax :     # Condition to assure only values less than or equal to Bmax are included
                temp2 =[m,l[i]]      # Appends both the m value of the Bessel function and the value of the Bessel zero
                temp.append(np.array(temp2))
        temp = np.array(temp)
        if len(temp) > 0 :
            sph_jn.insert(0,temp)
        m = m + 1
        j = j + 1
        if j==3:
            n = n - 1
            j = 0
    sph_jn = np.array(sph_jn)
    return sph_jn   # 1D Array of Bessel function zeros from largest function to smallest function

# Creates an array with the Bessel zero, function number, and zero number
def orderPair(Bzeros,a,L,mmax,nmax,nzmax):
    Bmn = []
    for i in range(len(Bzeros)-1,-1,-1):   # Decreasing loop going through the 1D array of bessel functions
        for j in range(0,len(Bzeros[i])):
            Bmn.append(np.array([Bzeros[i][j,-1],abs(i-(len(Bzeros)-1)),j+1]))
    Bmn = np.array(Bmn)
    Bmn = Bmn[Bmn[:,0].argsort()]   # Sorted list based on Bessel zero values from smallest to largest
    FE = fermiEnergy(Bmn,a,L,mmax,nmax,nzmax)   # Call to the fermiEnergy function
    FS = filledStates(FE)   # Call to the filledStates function
    Fermi = []
    for i in range(len(FE)):
        Fermi.append(FE[i][0])
    Fermi = np.array(Fermi)
    return [FS,Fermi]   # 2D Array of the filled state with its associated Fermi energy

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
def fermiEnergy(Z,a,L,mmax,nmax,nzmax):
    temp = []
    temp.append(special.jn_zeros(mmax,nmax))
    Bmax = temp[0][-1]
    fermi = []
    lim = float( ((Bmax ** 2)/((np.pi * (a/L)) ** 2)) + (nzmax ** 2) ) 
    for i in range(len(Z)):   # Loop over the number of Bessel zeros
        nz = 1
        while float( ((Z[i][0] ** 2)/((np.pi * (a/L)) ** 2)) + (nz ** 2) ) <=  lim: 
            energy = float( ((((np.pi ** 2) * (a **2)) / (3 * (L ** 2))) ** (2/3)) * (((Z[i][0] ** 2)/((np.pi * (a/L)) ** 2)) + (nz ** 2)) )     # Calculated dimensionless Fermi energy
            fermi.append(np.array([energy,Z[i][1],Z[i][2],nz]))  
            nz = nz + 1
    fermi = np.array(fermi)
    fermi = fermi[fermi[:,0].argsort()]
    return fermi   # 1D Array of the Fermi energy

# Analytical analysis 
def arrayAnlyt(lim):
    analytic = []
    xvals = []
    for i in range(1,lim):   # Loop to create the analytical results of N^(2/3)
        analytic.append(i**(2/3))
        xvals.append(i)
    analytic = np.array(analytic)
    xvals = np.array(xvals)
    return [xvals,analytic]  

# Logarithmic plots
def scatterLog(x1,y1,x2,y2):
    plot1 = plt.figure(1)
    plt.xlabel("N")
    plt.ylabel("Dimensionless Fermi Energy")
    plt.scatter(x1,y1,s=0.5)
    plt.plot(x2,y2,'--')
    plt.xscale('log')
    plt.yscale('log')
    return plot1   # Returns the info to make a logarithmic plot of numerical results

main()
