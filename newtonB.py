#!/usr/bin/env python3
'''
 NEWTON-BISECTION METHOD

 Solves the problem f(x)=0 using Newton-Bisection method. For a known true solution
 calculates errors.

 The main function is newton:
 
  [state,x,errors,iters] = newtonBIsection(x0, tolerance, maxIteration, debug)

  Inputs:
    a,b             The initial bounding interval, with a root between.
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
      BAD_DATA      Error: The interval may not bracket a root

  Remark: We assume that we known the two functions
    f               The name of the function for which a root is sought
    df              The name of the derivative of the function.
'''
#!/usr/bin/env python3
'''
 NEWTON-BISECTION METHOD (AP3)

 Solves f(x)=0 using Newton-Bisection method. If x_true is provided,
 prints error e_i and ratios:
   rq  = e_i / e_{i-1}^2
   rsl = e_i / (e_{i-1} e_{i-2})
   rl  = e_i / e_{i-1}

 NOTE about double root (-2/3): f does NOT change sign across an even-multiplicity root.
 Choose an interval endpoint equal (or extremely close) to the root so f(endpoint)=0,
 otherwise BAD_DATA may occur.
'''
from numpy import zeros, sign

############################## VARIABLES #############################
SUCCESS = 0
WONT_STOP = 1
BAD_DATA = 2

############################## FUNCTIONS #############################
def f(x):
    return 54*x**6 + 45*x**5 - 102*x**4 - 69*x**3 + 35*x**2 + 16*x - 4

def df(x):
    return 324*x**5 + 225*x**4 - 408*x**3 - 207*x**2 + 70*x + 16

def newtonBisection(a, b, TOL, MAX_ITERS, debug, x_true=None):
    prec = 12
    eps = 1e-20

    fmt = f"Iter %d: x= %.{prec}g, dx= %.{prec}g, error = %.{prec}g, "
    fmt += f"rq = %.{prec}g, rsl = %.{prec}g, rl = %.{prec}g, "
    fmt += f"interval = [%.{prec}g,%.{prec}g], Newton? %d"

    # ensure a < b
    if a > b:
        a, b = b, a

    fa = f(a)
    fb = f(b)

    # accept if it brackets OR includes a root at an endpoint (f(a)=0 or f(b)=0)
    endpoint_tol = TOL
    if sign(fa) * sign(fb) > 0.0 and abs(fa) > endpoint_tol and abs(fb) > endpoint_tol:
        return BAD_DATA, None, None, 0


    errors = zeros(MAX_ITERS + 1)

    # initial guess: midpoint
    x = a + (b - a) / 2.0
    err = abs(x - x_true) if x_true is not None else float("nan")
    errors[0] = err

    if debug:
        if x_true is None:
            print(f"Interval = [{a:.12g},{b:.12g}], initial x = {x:.12g} (no x_true; ratios will be nan)")
        else:
            print(f"Interval = [{a:.12g},{b:.12g}], initial x = {x:.12g}, error = {err:.12g}")

    fx = f(x)

    # shrink bracket using midpoint
    if sign(fa) * sign(fx) > 0.0:
        a = x
        fa = fx
    else:
        b = x
        fb = fx

    # Newton-Bisection loop
    for itn in range(1, MAX_ITERS + 1):
        dfx = df(x)
        usedNewton = True

        # Newton step if possible
        if abs(dfx) > eps:
            xNew = x - fx / dfx
            # if it leaves interval, revert to bisection
            if xNew < a or xNew > b:
                xNew = a + (b - a) / 2.0
                usedNewton = False
        else:
            xNew = a + (b - a) / 2.0
            usedNewton = False

        fxNew = f(xNew)

        # update bracket (or keep endpoint root)
        if sign(fa) * sign(fxNew) > 0.0:
            a = xNew
            fa = fxNew
        else:
            b = xNew
            fb = fxNew

        dx = xNew - x
        x = xNew
        fx = fxNew

        err = abs(x - x_true) if x_true is not None else float("nan")
        errors[itn] = err

        # ratios
        rq = float("nan")
        rsl = float("nan")
        rl = float("nan")
        if x_true is not None:
            if itn >= 1 and errors[itn-1] != 0:
                rl = errors[itn] / errors[itn-1]
                rq = errors[itn] / (errors[itn-1] ** 2)
            if itn >= 2 and errors[itn-1] != 0 and errors[itn-2] != 0:
                rsl = errors[itn] / (errors[itn-1] * errors[itn-2])

        if debug:
            print(fmt % (itn, x, dx, err, rq, rsl, rl, a, b, int(usedNewton)))

        # stopping rule
        if abs(dx) <= TOL:
            return SUCCESS, x, errors[:itn+1], itn

    return WONT_STOP, x, errors, MAX_ITERS


################################ MAIN ###############################
if __name__ == "__main__":
    print("Solve f(x)=0 on interval [a,b] using Newton-Bisection method (AP3)")

    a = float(input("Enter a: "))
    b = float(input("Enter b: "))
    tol = float(input("Enter tolerance (e.g. 5e-7): "))
    maxIter = int(input("Enter maxIteration (e.g. 50): "))
    debug = bool(int(input("Monitor iterations? (1/0): ")))

    x_true_in = input("Enter true root x_true for error calculation (or press Enter to skip): ").strip()
    x_true = float(x_true_in) if x_true_in != "" else None

    s, x, errors, iters = newtonBisection(a, b, tol, maxIter, debug, x_true)

    if s == SUCCESS:
        print("The root is %.16g" % (x))
        print("The root to 6 decimals is %.6f" % (x))
        print("The number of iterations is %d" % (iters))
        print("errors =", errors)
    elif s == WONT_STOP:
        print("ERROR: Failed to converge in %d iterations!" % (maxIter))
    elif s == BAD_DATA:
        print("ERROR: Unsuitable interval! (Must bracket a root or include it at an endpoint.)")
    else:
        print("ERROR: Coding error!")