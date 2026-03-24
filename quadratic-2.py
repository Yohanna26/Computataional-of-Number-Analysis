#!/usr/bin/env python
'''
    Solve ax^2 + bx + c = 0 for real roots

    return 0 in no error, 1 otherwise
'''

import argparse
from math import fabs,sqrt

parser = argparse.ArgumentParser(description="Solve ax^2 + bx + c =0 for real roots")
parser.add_argument("-d","--debug",action="store_true")
DEBUG = parser.parse_args().debug
r1 = None
r2 = None
############################## FUNCTIONS #############################

def quadraticFormula(a,b,c):
    global r1, r2
    discriminant = b*b - 4*a*c
    eps = 1e-20

    if DEBUG:
        print("a = %.20f" % a)
        print("b = %.20f" % b)
        print("c = %.20f" % c)
        print("D = %.20f" % discriminant)
    
    if(discriminant < 0 and fabs(discriminant) < eps):
        print("|abs(D)| =  %.6e; Setting D to 0" % fabs(discriminant))
        discriminant = 0
    
    if(discriminant < 0 or (fabs(a) < eps)):
        return 1
    
    if (b<0):
         # multiply through by -1 so b>=0
        a=-a
        b=-b
        c=-c
    
    r1 = -(b+sqrt(discriminant))/(2*a)
    r2 = -2*c / (b+sqrt(discriminant))
    return 0

############################## MAIN #############################

print("Solve ax^2 + bx + c = 0 for real roots.")
a = float(input("Enter a: "))
b = float(input("Enter b: "))
c = float(input("Enter c: "))

error = quadraticFormula(a,b,c)

if error:
    print("ERROR: a=0 or roots not real")
    exit(1)
print("Roots are %f and %f" % (r1,r2))