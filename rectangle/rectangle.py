#!/user/bin
# Python code to run calculations of the Fermi energy for an infinite rectangular well

# Packages
import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from time import process_time

# Main block to change parameters
def main():
    #t0=process_time()
    lim= 60000     # Maximum number of filled Fermi states included
    a= int(input("Enter length of side Lx: "))       # Changes side Lx 
    b= int(input("Enter length of side Ly: "))       # Changes side Ly
    c= int(input("Enter length of side Lz: "))     # Changes side Lz
    list1 = arrayNumer(lim,a,b,c)     
    list2 = arrayAnlyt(lim)   
    scatterLog(list1[0],list1[1],list2[0],list2[1])  # Creates logarthmic plots of the numerical and analytical solutions
    plt.show()
    #plt.show(block=False) 
    #t1=process_time()
    #print("Time elapsed:", t1-t0, "seconds")
    #print("Job has completed successfully... Fare thee well")

# Numerical analysis
def arrayNumer(lim,a,b,c):
    stopper = float( ( ( 3 / ( np.pi * a * b * c ) ) * lim ) ** (2/3) )

    allval = []
    for nx in range(1, int( ( math.sqrt( stopper ) * a ) + 1 ), 1):   # Loop over nx, ny, and nz to create a list of Fermi energies
        for ny in range(1, int( ( math.sqrt( stopper - ( ( nx / a ) ** 2 ) ) * b ) + 1 ) , 1):
            for nz in range(1, int( ( math.sqrt( stopper - ( ( nx / a ) ** 2 ) - ( ( ny / b ) ** 2 ) ) * c ) + 1) , 1):
                allval.append( ( ( ( np.pi * a * b * c ) / 3 ) ** ( 2 / 3 ) ) * ( ( ( nx / a ) ** 2 ) + ( ( ny / b ) ** 2 ) + ( ( nz / c ) ** 2 ) ) )
    allval = np.asarray(allval)   # List to array
    allval, val = np.unique(allval,return_counts=True) # Unique Fermi energies and number of unique states

    newval = []
    for i in range(len(val)):   # Loop over the number of unique Fermi energies to get the number of filled states under each Fermi energy 
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

# Logarithmic plots
def scatterLog(x1,y1,x2,y2):
    plot2 = plt.figure(1,figsize=(6,6))
    plt.xlabel("N")
    plt.ylabel("Dimensionless Fermi Energy")
    plt.plot(x1,y1)
    plt.plot(x2,y2,'--')
    plt.xscale('log')
    plt.yscale('log')
    return plot2   # Returns the info to make a logarithmic plot of numerical results

main()
