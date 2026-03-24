#!/usr/bin/env python3
'''
 BISECTION METHOD

 Solves the problem
   f(x) = 0
 using the bisection algorithm.

 The main function is bisection:

 [state,root] = bisection(a, b, tolerance, maxIteration, debug)

 Inputs:
   a,b           The initial bounding interval, with a root between.
   tolerance     The convergence tolerance (must be > 0).
   maxIteration  The maximum number of iterations that can be taken.
   debug         Boolean for printing out information on every iteration.
 Outputs:
   root             The solution.
   state         An error status code.
     SUCCESS     Sucessful termination.
     WONT_STOP   Error: Exceeded maximum number of iterations.
     BAD_DATA    Error: The interval may not bracket a root.

 We assume the function f is given:
   double f(double);
'''

from math import exp,cos
from enum import Enum 
# Note: if don't have enum package, please use bisection2.py
############################## VARIABLES #############################
class STATE(Enum):
    SUCCESS = 0
    WONT_STOP = 1
    BAD_DATA = 2
############################## FUNCTIONS #############################

def f(x):
    return 3*x**5 - 5*x**3 - 2*x + 1# Example 1.2 in Sauer
    # return x*x 

def sgn(x):
    if (x>0):
        return 1
    elif (x<0):
        return -1
    elif (x==0):
        return 0
    else:
        return x #x is not a number

# Returns two values, the first is the error state and the second is the found root
def bisection(a,b,tolerance, maxIteration,debug):
    x = None
    # format string
    prec = 8
    fmt = f"Iter %d: x= %.{prec}g, dx = %.{prec}g, a = %.{prec}g, b = %.{prec}g, f(x) = %.{prec+4}g"
    
    # if necessary, swap a and b
    if ( a > b):
        c = a
        a = b
        b = c
    
    fa = f(a)
    fb = f(b)

    # make sure there is a root between a and b
    if(sgn(fa)*sgn(fb)>0.0):
        return STATE.BAD_DATA,x
    
    # iteration loop
    dx = b-a
    for iteration in range(maxIteration):
        dx/=2
        x = a+dx
        fx = f(x)
        if debug:
            print(fmt % (iteration, x, dx,a,b,fx))
        
        # Check error tolerance
        if (dx <= tolerance):
            return STATE.SUCCESS,x
        
        if(sgn(fa)*sgn(fx)>0.0):
            a = x
            fa = fx
        else:
            b = x
            fb = fx
    return STATE.WONT_STOP,x
################################ MAIN ###############################

### input
print("Solves the problem f(x) = 0 on interval [a,b] using the bisection algorithm")
a = float(input("Enter a: "))
b = float(input("Enter b: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = int(input("Monitor iterations? (1/0): ")) == 1

### Solve for a root
[s,root] = bisection(a,b,tol,maxIter,debug)
if s is STATE.SUCCESS:
    print("The root is %.12g"%(root))
    print("f(%.12g) = %.12g"%(root,f(root)))
    exit()
elif s is STATE.WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!"%(maxIter))
elif s is STATE.BAD_DATA:
    print("ERROR: Unsuitable interval!")
else:
    print("ERROR: Coding error!")
exit(1)
