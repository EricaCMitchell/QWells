#!/usr/bin/
# Python code to run calculations of the Fermi energy for an infinite cubic well

# Packages
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from time import process_time

# Main block to change parameters
def main():
    #t0=process_time()
    lim=60000     # Maximum number of filled Fermi states included
    list1 = arrayNumer(lim)
    print(list1)
    list2 = arrayAnlyt(lim)
    scatterReg(list1[0],list1[1],list2[0],list2[1])   # Creates Cartesian plots of the numerical and analytical solutions
    scatterLog(list1[0],list1[1],list2[0],list2[1])   # Creates logaritmic plots of the numerical and analytical solutions
    plt.show()
    #plt.show(block=False)
    #t1=process_time()
    #print("Time elapsed:", t1-t0, "seconds")
    #print("Job has completed successfully... Fare thee well")

# Numerical analysis
def arrayNumer(lim):
    stopper = float( ((3/np.pi) * lim) ** (2/3) )

    allval = []
    for nx in range(1, int(math.sqrt(stopper)+1), 1):   # Loop over nx, ny, and nz to create a list of Fermi energies
        for ny in range(1, int(math.sqrt(stopper - (nx ** 2))+1) , 1):
            for nz in range(1, int(math.sqrt(stopper - (nx ** 2) - (ny ** 2))+1) , 1):
                allval.append(((np.pi/3)**(2/3))*((nx ** 2) + (ny ** 2) + (nz ** 2)))
    allval = np.asarray(allval)   # List to array
    allval, val = np.unique(allval,return_counts=True) # Unique Fermi energies and number of unique states

    newval = []
    for i in range(len(val)):   # Loop over the number of unique Fermi energies to get the number of filled state under each Fermi energy
        if i > 0:
            newval.append(2*(val[i])+newval[i-1])
        else:
            newval.append(2*(val[i]))
    newval = np.asarray(newval)
    return [newval,allval]   # Returns x-values as number of filled states and y-values as the unique Fermi energies

# Analytical analysis
def arrayAnlyt(lim):
    analytic = []
    xvals = []
    for i in range(1,lim):   # Loop to create the analytical results of N^(2/3)
        analytic.append(i**(2/3))
        xvals.append(i)
    return [xvals,analytic]

# Cartesian plots
def scatterReg(x1,y1,x2,y2):
    plot1 = plt.figure(1,figsize=(12,4))
    plt.subplot(1,2,1)
    plt.xlabel("N")
    plt.ylabel("Dimensionless Fermi Energy")
    plt.scatter(x1,y1,s=0.5)
    plt.plot(x2,y2,'--')
    return plot1    # Returns the info to make Cartesian plots of numerical results

# Logarithmic plots
def scatterLog(x1,y1,x2,y2):
    plot2 = plt.figure(1,figsize=(12,4))
    plt.subplot(1,2,2)
    plt.xlabel("N")
    plt.ylabel("Dimensionless Fermi Energy")
    plt.scatter(x1,y1,s=0.5)
    plt.plot(x2,y2,'--')
    plt.xscale('log')
    plt.yscale('log')
    return plot2     # Returns the info to make a logarithmic plot of numerical results

main()
