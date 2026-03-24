# Cosine calculator modeled after newtonDD.py

import numpy as np
import matplotlib.pyplot as pyp
import argparse

# this just makes it so that if you do "python3 cos1.py --test", it will run test environment
# where numbers are input through the command line.
parser = argparse.ArgumentParser(description="Provide cosine interpolation with Newton's Divided Differences")
parser.add_argument("-t","--test",action="store_true")
TEST = parser.parse_args().test

############################## FUNCTIONS #############################

# Evaluate divided difference interpolant
def newtonEval(t,coefs,x):
    n = len(coefs)
    value = coefs[n-1]
    for i in range(n-2,-1,-1): # same as n-2, n-3, n-4, ..., 0
        value = value*(t-x[i]) + coefs[i]
    return value

# Set up divided difference coefficients
def newtonDDsetup(x,y):
    n = len(x)
    if (len(y) != n): 
        print("ERROR CODE 1: x and y are different sizes")
        exit(1)

    # DD level 0
    # coefs[i] = y[i] for i=0,1,2,...,n-1
    coefs = [y[i] for i in range(n)] 
    
    # DD higher levels (bottom to top, overwrite lower entries as they are finished)
    for level in range(1,n): # 1,2,3,4, ... n-1
        for i in range(n-1,level-1,-1): #n-1, n-2, ..., level
            dx = x[i] - x[i-level] 
            if (dx==0): exit(2)
            coefs[i] = (coefs[i]-coefs[i-1])/dx
    return coefs

# x,y are 1d- arrays
def newtonDD(x,y):
    n = len(x)
    if (len(y) != n): exit(1)
    if TEST: print("x =",x)
    if TEST: print("y =",y)
    coefs = newtonDDsetup(x,y)
    if TEST: 
        print("The coefs are: ", end=" ")
        for i in range(n):
            print(" %g"%(coefs[i]),end=" ")
        print()

    m = 10*n
    minx = min(x)
    maxx = max(x)
    t = np.arange(minx,maxx+1,(maxx-minx)/m)
    val = newtonEval(t,coefs,x)

    # plot
    pyp.plot(t,val)
    for i in range(n):
        pyp.plot(x[i],y[i],'k*')
    pyp.xlabel("x")
    pyp.ylabel("y")
    pyp.xlim(minx,maxx)
    pyp.ylim(top=max(y))
    pyp.title("Newton DD interpolation")
    pyp.legend(["Newton DD interpolant","Data points"],loc="best")
    pyp.show()

########################################################
# cosine-specific functions
########################################################

# fundamental interpolation nodes on [0, pi/2]
xnodes = [0.0, np.pi/4, np.pi/2]
ynodes = [1.0, 1/np.sqrt(2), 0.0]
coefs_cos = newtonDDsetup(xnodes,ynodes)

def reduceCos(x):
    """
    Reduce x to xr in [0, pi/2] and return sign s
    so that cos(x) = s*cos(xr).
    """
    # reduce to [0,2pi)
    t = x % (2*np.pi)

    # reduce to [0,pi]
    if (t > np.pi):
        t = 2*np.pi - t

    # reduce to [0,pi/2]
    if (t <= np.pi/2):
        xr = t
        s = 1.0
    else:
        xr = np.pi - t
        s = -1.0

    return xr,s

def cos1(x):
    xr,s = reduceCos(x)
    return s*newtonEval(xr,coefs_cos,xnodes)

def printTable():
    xvals = [-1,1,-2,2,-3,3,-4,4,-14,14,-1000,1000]

    print("              x              cos(x)             cos1(x)            error(x)")
    print("--------------------------------------------------------------------------")
    for x in xvals:
        trueval = np.cos(x)
        approx = cos1(x)
        err = abs(trueval - approx)
        print("%15g %18.10f %18.10f %18.10f"%(x,trueval,approx,err))

def bonusPlot():
    t = np.arange(-np.pi,3*np.pi+0.01,0.01)
    trueval = np.cos(t)

    approx = np.zeros(len(t))
    for i in range(len(t)):
        approx[i] = cos1(t[i])

    pyp.plot(t,trueval)
    pyp.plot(t,approx,'--')

    # interpolation nodes
    for i in range(len(xnodes)):
        pyp.plot(xnodes[i],ynodes[i],'ro')

    pyp.annotate("(0,1)",(xnodes[0],ynodes[0]),textcoords="offset points",xytext=(5,8))
    pyp.annotate("(pi/4,1/sqrt(2))",(xnodes[1],ynodes[1]),textcoords="offset points",xytext=(5,8))
    pyp.annotate("(pi/2,0)",(xnodes[2],ynodes[2]),textcoords="offset points",xytext=(5,8))

    pyp.xlabel("x")
    pyp.ylabel("y")
    pyp.xlim(-np.pi,3*np.pi)
    pyp.title("cos(x) and cos1(x)")
    pyp.legend(["cos(x)","cos1(x)","Interpolation nodes"],loc="best")
    pyp.grid()
    pyp.show()

def testCos1():
    print("Interpolation nodes:")
    print("x =",xnodes)
    print("y =",ynodes)
    print()

    print("The coefs are: ", end=" ")
    for i in range(len(coefs_cos)):
        print(" %g"%(coefs_cos[i]),end=" ")
    print()
    print()

    printTable()
    print()

    ans = input("Make bonus plot? [y/n] ")
    if (ans[0] == 'y' or ans[0] == 'Y'):
        bonusPlot()

############################## MAIN #############################
if TEST:
    testCos1()


#I modeled my cos1 program after the structure of newtonDD.py. I kept the same Newton divided-difference helper functions, namely newtonEval and newtonDDsetup, and then built a cosine calculator on top of them.For cosine, the fundamental domain is [0, pi/2]. 
#This is because cosine is 2pi-periodic and even, so every x can first be reduced to [0, pi]. Then for x in [pi/2, pi], we use cos(x) = -cos(pi - x), which reduces the problem to [0, pi/2].

#On [0, pi/2], I interpolated cosine at the three nodes 0, pi/4, and pi/2,using the valuescos(0) = 1, cos(pi/4) = 1/sqrt(2), and cos(pi/2) = 0.

#The program computes cos(x), cos1(x), and the error |cos(x) - cos1(x)| for the required test values x = ±1, ±2, ±3, ±4, ±14, ±1000, and also includes the bonus plot.##