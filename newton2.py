#!/usr/bin/env python3
'''
 NEWTON'S METHOD

 Solves the problem f(x)=54x^6+45x^5-102x^4-69x^3+35x^2+16x-4 on the interval [-2,2]using Newton's method. For a known tru solution
 calculates errors.

 The main function is newton:
 
  [state,x,errors,iters] = newton(x0, tolerance, maxIteration, debug)

  Inputs:
    x0              The initial guess at the solution
    tolerance       The convergence tolerance (must be > 0).
    maxIteration    The maximum number of iterations that can be taken.
    debug           Boolean to set debugging output
  Outputs:
    x               The solution
    errors          Array with errors at each iteration
    iter            number of iterations to convergence
  Return:
    state           An error status code.
      SUCCESS       Sucessful termination.
      WONT_STOP     Error: Exceeded maximum number of iterations.
      BAD_ITERATE   Error: The function had a vanishing derivative

  Remark: We assume that we known the two functions
    f               The name of the function for which a root is sought
    df              The name of th derivative of the function. The derivative
                    of f must be computed by hand and coded correctly as df, or
                    newton will not work!

'''
from numpy import zeros

############################## VARIABLES #############################
SUCCESS = 0
WONT_STOP = 1
BAD_ITERATE = 2

############################## FUNCTIONS #############################
# The function for which a root is sought
def f(x):
    return 54*x**6 + 45*x**5 - 102*x**4 - 69*x**3 + 35*x**2 + 16*x - 4

# Derivative of f(x)
def df(x):
    return 324*x**5 + 225*x**4 - 408*x**3 - 207*x**2 + 70*x + 16

def newton(x0, TOL, MAX_ITERS, debug, x_true=None, multiplicity=1):
    global SUCCESS, WONT_STOP, BAD_ITERATE
    prec = 12
    eps = 1e-20

    # Print format (now includes rq, rl)
    fmt = f"Iter %d: x= %.{prec}g, dx= %.{prec}g, error = %.{prec}g, rq = %.{prec}g, rl = %.{prec}g"

    errors = zeros(MAX_ITERS+1)
    x = x0

    if x_true is not None:
        err = abs(x - x_true)
        errors[0] = err
    else:
        err = float('nan')
        errors[0] = float('nan')

    if debug:
        if x_true is None:
            print(f"Guess: x={x:.8g} (no true root provided, ratios disabled)")
        else:
            print("Guess: x=%.8g, error=%.8g" % (x, err))

    ## Newton Loop
    for itn in range(1, MAX_ITERS+1):
        dfx = df(x)
        if abs(dfx) < eps:
            state = BAD_ITERATE
            iters = itn
            return state, x, errors, iters

        dx = -f(x)/dfx
        dx = multiplicity * dx     # multiplicity=2 gives modified Newton
        x += dx

        # error and ratios
        rq = float('nan')
        rl = float('nan')

        if x_true is not None:
            err = abs(x - x_true)
            errors[itn] = err

            if itn >= 1 and errors[itn-1] != 0:
                rl = errors[itn] / errors[itn-1]
                rq = errors[itn] / (errors[itn-1]**2)
        else:
            errors[itn] = float('nan')

        if debug:
            # If x_true not provided, print NaNs for ratios
            print(fmt % (itn, x, dx, err, rq, rl))

        # Check tolerance (same rule as your original: |dx| <= TOL)
        if abs(dx) <= TOL:
            iter = itn
            state = SUCCESS
            return state, x, errors, iter

    state = WONT_STOP
    iter = itn
    return state, x, errors, iter


################################ MAIN ################################

print("Solve the problem f(x)=0 using Newton's method for Additional Problem 2")

x0 = float(input("Enter guess at root x0: "))
tol = float(input("Enter tolerance (e.g. 5e-7): "))
maxIter = int(input("Enter maxIteration (e.g. 50): "))
debug = bool(int(input("Monitor iterations? (1/0): ")))

# For error ratios you MUST provide a true root (use the known value for the root you're targeting)
x_true_in = input("Enter true root x_true for error calculation (or press Enter to skip): ").strip()
x_true = float(x_true_in) if x_true_in != "" else None

m_in = input("Enter multiplicity m for modified Newton (default 1, use 2 for double root): ").strip()
m = int(m_in) if m_in != "" else 1

[s, x, errors, iters] = newton(x0, tol, maxIter, debug, x_true=x_true, multiplicity=m)

if s == SUCCESS:
    print("The root is %.16g" % (x))
    print("The root to 6 decimals is %.6f" % (x))
    print("The number of iterations is %d" % (iters))
    errors = errors[:iters+1]
    print("errors =", errors)
elif s == WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!" % (maxIter))
elif s == BAD_ITERATE:
    print("ERROR: Obtained a vanishing derivative!")
else:
    print("ERROR: Coding error!")
