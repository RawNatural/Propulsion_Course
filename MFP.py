#! /usr/bin/env python
"""
Mass Flow Parameter calculation function.

The mass flow parameter is defined so that MFP(M,k,R)*Area is
constant for isentropic flow.  It can be used for solving nozzle
flow problems as follows:

If the desired outflow Mach number M2 and the inflow state 1 are
known, the required nozzle area ratio can be computed using

A2/A1 = MFP(M1,k,R)/MFP(M2,k,R)

If the area ratio and inflow state are known, the outflow Mach
number M2 can be found by iteratively solving

MFP(M2,k,R) =  MFP(M1,k,R)*A1/A2

In Python this iterative solution can be found using:

from scipy.optimize import fsolve
M2 = fsolve(lambda M: MFP(M,k,R) - MFP(M1,k,R)*A1/A2, M1 + 2.0)[0]

Refer to the scipy documentation for more information.
http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html

The initial guess of (M1 + 2.0) is appropriate when seeking a supersonic
value of M2.

Vince Wheatley's Matlab code ported to Python by:

Peter Blyton
Semester 1, 2011
MECH4450 Aerospace Propulsion
The University of Queensland
"""

import numpy as np

def MFP(M,kc,R):
    """
    Return the mass flow parameter.
    
    Arguments:
    M: (float) Flight Mach number.
    k: (float) Ratio of specific heats.
    R: (float) J/kg.K, gas constant.
    """
    return np.sqrt(kc/R)*M*(1.0 + 0.5*(kc - 1.0)*M**2)**(-0.5*(kc + 1.0)/(kc - 1.0))

if __name__ == "__main__":
    from scipy.optimize import fsolve
    kc = 1.4
    c_pc = 1004
    R = c_pc*(1 - 1/kc)
    M1 = 6.0
    A1 = 1.0
    A2 = 2.5
    M2 = fsolve(lambda M: MFP(M,kc,R) - MFP(M1,kc,R)*A1/A2, M1 + 2.0)[0]
    print("For a gas with k =", kc)
    print("R = ", R)
    print("inlet Mach number of ", M1)
    print("expanding through an area ratio of ", A2/A1)
    print("the exit Mach number is ", M2)
    print("done.")

