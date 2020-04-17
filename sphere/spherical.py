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
    t0 = process_time()
    kmax = float(input("Enter number to be multiplied by Pi: ")) * np.pi #Value of kappa_max
    file1 = open('nlm.txt','w')     # Creates a file for the nx, ny, and nz values
    file2 = open('kNk.txt','w')     # Creates a file for the kappa and N(kappa) values
    list1 = loop(kmax)     # List of lists containing 1) n, l, and m values and 2) kappa^2 values 
    kNk = pairGen(list1[1],False)     # List of lists containing 1) kappa and N(kappa) and 2) kappa and improved prescription N(kappa) [only calculated if set to True]
    numer = arrayNumer(kNk[0])     # List of lists containing numerical approximations containing 1) all terms, 2) volume term, and 3) volume and area terms
    regPlot(kNk[0],numer[0],numer[1],numer[2])     # Function call to make 2d plot of N(kappa) vs kappa including plots of numerical approximations
    scatterPlot(list1[0])     # Function call to make 3d plot of states contributing to N(kappa)
    t1 = process_time()
    plt.show()
    file1.write(" n  l  m \n")
    for i in range(len(list1[0])):
        file1.write(str(list1[0][i]) + "\n")
    file2.write(" kappa      Nk \n")
    for i in range(len(kNk[0])):
        file2.write(str(kNk[0][i]) + "\n")
    file1.close()
    file2.close()
    print("Number of states (N(k)):", int(kNk[0][-1][1]))
    print("Time elapsed:", t1 - t0, "seconds") 
    print("Job has completed successfully... Fare thee well")

# Numerically solves for kappa^2 and records the n, l, m values used to determine it
def loop(a):
    triplets = []
    k2 = []
    n = 1
    nok = True
    while nok != False:
        l = 0
        lok = True
        while lok != False:
            Bln = zeros.main(l,n)
            k2var = ( (((4 * np.pi) / 3) ** (2/3)) * (Bln[-1] ** 2) )
            if k2var > (a ** 2):     # Condition assuring that kappa^2 is smaller than kappa_max^2
                lok = False
                if l == 0:
                    nok = False
                else:
                    n += 1
            else:
                for i in range(-l,l+1):
                    temp = [n,l,i]
                    triplets.append(temp)
                    k2.append(round(k2var,5))
                l += 1
    triplets = np.array(triplets)
    k2 = np.array(k2)
    return [triplets,k2]     # List of lists containing 1) nx, ny, and nz values and 2) kappa^2 values 

# Creates arrays of ordered pairs
def pairGen(l,rx):
    k2, Nk = np.unique(l,return_counts=True)
    pairOG = [] 
    pairRx = []
    for i in range(len(k2)):    # Loop making the ordered pair, (kappa, N(kappa))
        k = math.sqrt(k2[i])
        if i > 0:
            N = temp[1] + Nk[i]
        else:
            N = Nk[i]
        temp = [k,N]
        pairOG.append(temp)
    if rx == True:     # Condition allowing the bypass of the improved prescription for N(kappa)
        for i in range(len(k2)):
            if i > 0:
                Rx1 = pairOG[i-1][1] + pairOG[i][1]
            else:
                Rx1 = pairOG[i][1]
            Rx2 = float(Rx1 / 2)
            temp = [pairOG[i][0],Rx2]
            pairRx.append(temp)
    pairOG = np.array(pairOG)
    pairRx = np.array(pairRx)
    return [pairOG,pairRx]     # List of lists containing 1) kappa and N(kappa) and 2) kappa and improved prescription N(kappa) [only calculated if set to True]
                            
# Solves for the numerical approximation to N(kappa) using the averaged asymptotic expression
def arrayNumer(l):
    allTerms = []
    vTerm = []
    vsTerms = []
    for i in range(len(l)):
        K = l[i][0]
        NK3 = ( (K ** 3) / (6 * (np.pi ** 2)) ) 
        NK2 = ( ((3 / (32 * np.pi)) ** (2/3)) * (K ** 2) )
        NK = ( ((2 / (9 * (np.pi ** 4))) ** (1/3)) * K )
        N = NK3 - NK2 + NK
        allTerms.append(N)
        vTerm.append(NK3)
        vsTerms.append(NK3 - NK2)
    allTerms = np.array(allTerms)
    vTerm = np.array(vTerm)
    vsTerms = np.array(vsTerms)
    return [allTerms,vTerm,vsTerms]     # List of lists containing numerical approximations containing 1) all terms, 2) volume term, and 3) volume and area terms

# Creates a 3d plot of the states contained within the geometry
def scatterPlot(l):
    plot1 = plt.figure(1,figsize=(12,4))
    ax = plot1.add_subplot(121, projection='3d')
    for i in range(len(l)):
        xs = l[i][2]
        ys = l[i][0]
        zs = l[i][1]
        ax.scatter(xs, ys, zs)
    ax.set_xlabel('m')
    ax.set_ylabel('n')
    ax.set_zlabel('l')
    return plot1

# Creates a 2d plot of kappa versus N(kappa), including plots of numerical approximations
def regPlot(l1,l2,l3,l4):
    xs = [0]
    y1 = [0]
    y2 = [0]
    y3 = [0]
    y4 = [0]
    plot2 = plt.figure(1,figsize=(12,4))
    plt.subplot(1,2,2)
    for i in range(len(l1)):
        xs.append(l1[i][0])
        y1.append(l1[i][1])
        y2.append(l2[i])
        y3.append(l3[i])
        y4.append(l4[i])
    plt.step(xs, y1, where='post', color='blue', label='Exact numerical')
    plt.plot(xs, y2, '--', color='grey', label='Approx. using all four terms')
    plt.plot(xs, y3, ':', color='orange', label='Approx. using volume')
    plt.plot(xs, y4, '-.', color='green', label='Approx. using volume and area')
    plt.legend()
    plt.xlabel("K")
    plt.ylabel("N(K)")
    return plot2

main()
