import numpy as np
import matplotlib.pyplot as pyp
import argparse

# this just makes it so that if you do "python3 cos12.py --test", it will run test environment
parser = argparse.ArgumentParser(description="Provide cosine interpolation with Newton's Divided Differences")
parser.add_argument("-t","--test",action="store_true")
TEST = parser.parse_args().test

############################## FUNCTIONS #############################

# Evaluate divided difference interpolant
def newtonEval(t,coefs,x):
    n = len(coefs)
    value = coefs[n-1]
    for i in range(n-2,-1,-1):
        value = value*(t-x[i]) + coefs[i]
    return value

# Set up divided difference coefficients
def newtonDDsetup(x,y):
    n = len(x)
    if (len(y) != n):
        print("ERROR CODE 1: x and y are different sizes")
        exit(1)

    coefs = [y[i] for i in range(n)]

    for level in range(1,n):
        for i in range(n-1,level-1,-1):
            dx = x[i] - x[i-level]
            if (dx == 0):
                exit(2)
            coefs[i] = (coefs[i]-coefs[i-1])/dx
    return coefs

########################################################
# cosine-specific functions
########################################################

# ------------------------------------------
# cos1(x): interpolation at 0, pi/4, pi/2
# ------------------------------------------
xnodes1 = [0.0, np.pi/4, np.pi/2]
ynodes1 = [1.0, 1/np.sqrt(2), 0.0]
coefs_cos1 = newtonDDsetup(xnodes1, ynodes1)

# ------------------------------------------
# cos2(x): 3 Chebyshev nodes on [0, pi/2]
# nodes on [a,b]:
# x_k = (a+b)/2 + (b-a)/2 * cos((2k+1)pi/(2n)), k=0,1,...,n-1
# here n=3, a=0, b=pi/2
# ------------------------------------------
a = 0.0
b = np.pi/2
n = 3

xnodes2 = []
for k in range(n):
    xk = (a+b)/2 + (b-a)/2 * np.cos((2*k+1)*np.pi/(2*n))
    xnodes2.append(xk)

# sort nodes from left to right
xnodes2 = sorted(xnodes2)
ynodes2 = [np.cos(xk) for xk in xnodes2]
coefs_cos2 = newtonDDsetup(xnodes2, ynodes2)

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
    return s*newtonEval(xr, coefs_cos1, xnodes1)

def cos2(x):
    xr,s = reduceCos(x)
    return s*newtonEval(xr, coefs_cos2, xnodes2)

def printTable():
    xvals = [1,-1,2,-2,3,-3,4,-4,14,-14,1000,-1000]

    print("                   x              cos(x)            cos1(x)        err cos1            cos2(x)        err cos2")
    print("----------------------------------------------------------------------------------------------------------------")
    for x in xvals:
        trueval = np.cos(x)
        approx1 = cos1(x)
        err1 = abs(trueval - approx1)
        approx2 = cos2(x)
        err2 = abs(trueval - approx2)

        print("%20g %18.10f %18.10f %15.10f %18.10f %15.10f"
              %(x,trueval,approx1,err1,approx2,err2))

def bonusPlots():
    t = np.arange(-np.pi,3*np.pi+0.01,0.01)
    trueval = np.cos(t)

    approx1 = np.zeros(len(t))
    approx2 = np.zeros(len(t))
    err1 = np.zeros(len(t))
    err2 = np.zeros(len(t))

    for i in range(len(t)):
        approx1[i] = cos1(t[i])
        approx2[i] = cos2(t[i])
        err1[i] = abs(trueval[i] - approx1[i])
        err2[i] = abs(trueval[i] - approx2[i])

    # Figure 1: cos(x) and cos1(x)
    pyp.figure()
    pyp.plot(t,trueval)
    pyp.plot(t,approx1,'--')

    for i in range(len(xnodes1)):
        pyp.plot(xnodes1[i],ynodes1[i],'ro')

    pyp.xlabel("x")
    pyp.ylabel("y")
    pyp.xlim(-np.pi,3*np.pi)
    pyp.title("cos(x) and cos1(x)")
    pyp.legend(["cos(x)","cos1(x)","Interpolation nodes"],loc="best")
    pyp.grid()
    pyp.show()

    # Figure 2: cos(x) and cos2(x)
    pyp.figure()
    pyp.plot(t,trueval)
    pyp.plot(t,approx2,'--')

    for i in range(len(xnodes2)):
        pyp.plot(xnodes2[i],ynodes2[i],'ro')

    pyp.xlabel("x")
    pyp.ylabel("y")
    pyp.xlim(-np.pi,3*np.pi)
    pyp.title("cos(x) and cos2(x)")
    pyp.legend(["cos(x)","cos2(x)","Chebyshev nodes"],loc="best")
    pyp.grid()
    pyp.show()

    # Figure 3: errors
    pyp.figure()
    pyp.plot(t,err1)
    pyp.plot(t,err2,'--')
    pyp.xlabel("x")
    pyp.ylabel("error")
    pyp.xlim(-np.pi,3*np.pi)
    pyp.title("Errors in cos1(x) and cos2(x)")
    pyp.legend(["|cos(x)-cos1(x)|","|cos(x)-cos2(x)|"],loc="best")
    pyp.grid()
    pyp.show()

def testCos():
    print("Interpolation nodes for cos1:")
    print("xnodes1 =", xnodes1)
    print("ynodes1 =", ynodes1)
    print()

    print("The coefficients for cos1 are: ", end="")
    for c in coefs_cos1:
        print(" %g" % c, end="")
    print("\n")

    print("Chebyshev interpolation nodes for cos2:")
    print("xnodes2 =", xnodes2)
    print("ynodes2 =", ynodes2)
    print()

    print("The coefficients for cos2 are: ", end="")
    for c in coefs_cos2:
        print(" %g" % c, end="")
    print("\n")

    printTable()
    print()

    ans = input("Make the three bonus plots? [y/n] ")
    if (ans[0] == 'y' or ans[0] == 'Y'):
        bonusPlots()

############################## MAIN #############################
if TEST:
    testCos()