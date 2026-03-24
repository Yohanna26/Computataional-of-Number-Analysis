#!/usr/bin/env python3
'''
 FIXED POINT (PICARD) ITERATION METHOD

 Solves the problem
   g(x) =     5x + 1 = e^(2x) + 6x^3
 using fixed point iteration. For a known true solution calculates errors

 The main function is fpi:
 
  [state, x, errors, iter] = fpi(@g, x0, tolerance, maxIteration, debug);

  Inputs:
    @g            Handle to function g
    x0            The initial guess at the fixed point
    tolerance     The convergence tolerance (must be > 0).
    maxIteration  The maximum number of iterations that can be taken.
    debug         Boolean for printing out information on every iteration.
  Outputs:
    x             The solution
    errors        Array with errors at each iteration
    iter          number of iterations to convergence
  Return:
    state         An error status code.
      SUCCESS     Sucessful termination.
      WONT_STOP   Error: Exceeded maximum number of iterations.
'''
from math import sqrt, exp
from numpy import zeros
############################## VARIABLES #############################
SUCCESS = 0
WONT_STOP = 1
x_true = 0 #negative root
############################## FUNCTIONS #############################
#f = 0.5*log(5*x + 1 - 6*x**3)

def f(x):
    return (exp(2*x) + 6*x**3 - 1)/5

def fpi(func,x0,TOL,MAX_ITERS,debug):
    global x_true, SUCCESS, WONT_STOP
    prec = 14
    fmt = f"Iter %3d: x= %.{prec}g, error = %.{prec}g, ratio = %.{prec}g"

    errors = zeros(MAX_ITERS+1)
    x = x0
    err = abs(x - x_true)
    errors[0] = err
    
    if debug:
      print(fmt % (0, x, err, float('nan')))
    
    ## FPI Loop
    for itn in range(1,MAX_ITERS+1):
        gx = func(x)
        dx = abs(x-gx)
        x = gx
        err = abs(x - x_true)
        errors[itn] = err
        ratio = err / errors[itn - 1] if errors[itn - 1] != 0 else float('nan')
        if debug:
            print(fmt % (itn, x, err, ratio))
        
        # Check error tolerance
        if dx <= TOL * (abs(x) + 1):
            return SUCCESS, x, errors, itn

    return WONT_STOP, x, errors, MAX_ITERS


################################ MAIN ###############################

###input
print("Solve the problem g(x)=x using fixed point iteration")
x0 = float(input("Enter guess at root: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = int(input("Monitor iterations? (1/0): ")) == 1

### Solve for the fixed point
[s,x,errors,iters] = fpi(f,x0,tol,maxIter,debug)
if s == SUCCESS:
    print("The root is %.16g"%(x))
    print("f(%.16g) = %.16g"%(x,f(x)))
    print("The number of iterations is %d"%(iters))

    errors = errors[:iters + 1]

    target = 5e-7
    first_ok = None
    for i, e in enumerate(errors):
        if e <= target:
            first_ok = i
            break

    if first_ok is not None:
        print("First iter with 6 correct decimals (error <= 5e-7): iter =", first_ok)
        print("Value at that iter rounded to 6 decimals:", format(x, ".6f"))
    else:
        print("Did not reach 6 correct decimals within the iterations run.")

    print("\nerrors =", errors)
else:
    print("ERROR: Failed to converge in %d iterations!" % maxIter)
