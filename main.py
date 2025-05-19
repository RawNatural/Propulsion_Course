import numpy as np
#m0 = something; alpha = 1; mC = m0/(1+alpha); mB = m0*alpha/(1+alpha)

""" General Assumptions """
c_pc = 1004 # [J/kg/K]
c_pt = 1239 # [J/kg/K]
c_ratio = c_pt/c_pc 
gamma_c = 1.4
gamma_t = 1.3
q_0 = 50 # [kPa]
T_0 = 250 # [K]
A_0 = 14.4 # [m2]
C_D = 0.03 # Drag Coefficient
A_ref = 382 # [m2]

""" Fuel Properties """
H = 120 # [MJ/kg]
f_st = 0.0291
phi = 1 # Maximum Equivalence Ratio

""" Inlet Modelling """
pi_dmax = 0.96
A_ratio =  2.5 # Diffuser Area Ratio # A_2/A_1

""" Core Modelling Parameters """
pi_c = 30 # Compressor Stagnation Pressure Ratio 
pi_f = 2.0 # Fan Stagnation Pressure Ratio
e_c = 0.91 # Compressor Polytropic Efficiency
e_f = 0.93 # Fan Polytropic Efficiency
pi_b = 0.99 # Burner Stagnation Pressure Ratio
pi_AB = 0.99 # Afterburner Stagnation Pressure Ratio
eff_b = 0.99 # Burner Combustion Efficiency
eff_AB = 0.99 # Afterburner Combustion Efficiency
e_t = 0.93 # Turbine Polytropic Efficiency
eff_m = 0.98 # Mechanical Transmission Efficiency
pi_n = 0.95 # Nozzle Stagnation Pressure Ratio
tau_n = 1 # Nozzle Stagnation Temperature Ratio
T_t4max = 2000 # [K] - Turbine Entry Stagnation Temperature
tau_lambda = (c_pt/c_pc)*(T_t4max/T_0)
#f = mf/mc # Burner Fuel Air Ratio
#f_AB = mf_AB/mc # Afterburner Fuel Air Ratio
tau_AB = 1.2 #T_t7/T_t5 #Also must initialise as something - design parameter. Can we figure this out based on combustion or something?


""" Bypass Modelling Parameters """
pi_f = 2.0 # Fan Stagnation Pressure Ratio
pi_fn = 0.95 # Nozzle Stagnation Pressure Ratio
tau_fn = 1 # Nozzle Stagnation Temperature Ratio
e_f = 2.0 # Fan Polytropic Efficiency
#f_B = mf_B/mB # Bypass Fuel Air Ratio
#pi_B = P_t14/P_t13
eff_B = eff_AB # Otherwise tau_B = 1

tau_B = 1.2 #T_t14/T_t13 # must initialise because it uses itself when finding itself.

Rc = c_pc*(1-1/gamma_c)
Rt = c_pt*(1-1/gamma_t)
a0 = np.sqrt(gamma_c*Rc*T_0)

alpha = 1

mC = None; mB = None;


tau_c = (pi_c)**e_t*((gamma_c - 1)/gamma_c)
tau_f = (pi_f)**e_t*((gamma_c - 1)/gamma_c)
tau_fn = 1


def tau_0(M0):
    return 1 + (gamma_c -1)/2 * M0**2

def pi_0(M0):
    return (1 + (gamma_c -1)/2 * M0**2)**(gamma_c/(gamma_c -1))

def tau_t(M0):
    return 1-1/(eff_m*(1+f_()))*tau['0']/tau['lambda']*(tau['c']-1)

"""Set all tau's and pi's in a dictionary"""
tau = {
    "c": tau_c,
    "lambda": tau_lambda,
    "n": tau_n,
    "fn": tau_fn,
    "f": tau_f,
    "0": None, #functions of M0
    "t": None,
    "AB": None,
    "B": None,
}

pi = {
    "c": pi_c,
    "d": pi_dmax,
    "f": pi_f,
    "b": pi_b,
    "AB": pi_AB,
    "n": pi_n,
    "fn": pi_fn,
    "f": pi_f,
    "t": None,
    "B": None,
    "0": None
}

f = {

}

"""Define Functions"""

def TB_ratio(M0):
    """ Bypass Temperature Ratio T19/T0 """ 
    TB_ratio = tau['fn']*tau_B_()*tau['f']*tau['0'] / (1*pi['0']*pi['B']*pi['f']*pi['fn'])
    return TB_ratio

def TC_ratio(M0):
    """ Core Temperature Ratio T9/T0 """
    TC_ratio = tau['lambda']*tau_t(M0)*(1/c_ratio) / ((1*pi['0']*pi['d']*pi['c']*pi['b']*pi['t']*pi['n']*pi['f'])**((gamma_t-1)/gamma_t))
    return TC_ratio

def PB_ratio_():
    """ Pressure stagnation ratio across Bypass Pt19/P19 """    
    return pi['fn']*pi['B']*pi['f']*pi['0']

def PC_ratio():
    """ Pressure stagnation ratio across Core Pt9/P9 """ 
    return pi['n']*pi['AB']*pi['t']*pi['b']*pi['c']*pi['f']*pi['0']

def M9():
    """ Mach number at point 9 """
    return (2/(gamma_t-1))*((1*pi['0']*pi['d']*pi['c']*pi['b']*pi['t']*pi['n'])**((gamma_t-1)/gamma_t)-1)

def M19(M0):
    """ Mach number at point 19 """
    return M0*1.5 #change for something later

def m0_(M0):
    rho0 = 100000/((M0*a0)**2)
    m_0 = rho0*a0*A_0*M0
    global mC, mB
    mC = m_0/(1+alpha); mB = m_0 * alpha / (1+alpha);
    return m_0

def f_():
    """ Calculates fuel air ratio for burner """
    numerator = tau['lambda'] - tau['c']*tau['0']
    denominator = (eff_b*H)/(c_pc*T_0)
    return numerator/denominator

def f_B_():
    """ Caclulates fuel air ratio for bypass """
    numerator = tau['B']*tau['f']*tau['0'] - tau['f']*tau['0']
    denominator = ((eff_B*H)/(c_pc*T_0)) - tau['B']*tau['f']*tau['0']
    return numerator/denominator

def tau_B_():
    """ Calculates tau['B'] """
    """
    tauB = 0; tB = 1;
    while abs(tB - tauB) < 0.1:
        numerator = f_B_(tB)*(eff_B*H/(c_pc*T_0)) + tau['f']*tau_0(M0)
        denominator = tau['f']*tau_0(M0) + f_B_()*tau['f']*tau_0(M0)
        tauB = numerator/denominator
    return tauB"""
    return tau_B

def pi_B():
    """ Assuming Isentropic across bypass burner """ #Task sheet says to use Rayleigh flow as in a ramjet (lecture 8 maybe)
    return 0.99 #use this for now because pi_AB and pi_b are 0.99

def f_AB_():
    """ Calculates fuel air ratio for after burner """
    numerator = (1+f_())*(tau['AB']*tau['t']*tau['lambda'] - tau['t']*tau['lambda']) 
    denominator = ((eff_AB*H)/(c_pt*T_0)) - tau['AB']*tau['t']*tau['lambda']
    return numerator/denominator

def tau_AB_():
    """tauAB = 0; tAB = 1;
    while abs(tAB - tauAB) < 0.1:
        numerator = f_AB_(tauAB) * (eff_AB*H/(c_pt*T_0))+tau['t']*tau['lambda']*(1 + f_())
        denominator = tau['t']*tau['lambda']*(1 + f_()) + f_AB_(tauAB)*tau['t']*tau['lambda']
        tauAB = numerator/denominator
    return tauAB"""
    return tau_AB #Actually supposed to calculate using f_st = f combustion things.

def pi_t(M0):
    return tau_t(M0) ** (gamma_t/((gamma_t-1)*e_t))

""" Key Notes
    From Assignment Task Sheet: P0 = P9 = P19, therefore P0/P9 = P0/P19 = 1 """
def F(M0):
    """ Calculates overall thrust """
    Fcore = m0_(M0)*(a0/(1+alpha))*((1+f_()+f_AB_())*(   np.sqrt(gamma_t*Rt*TC_ratio(M0)*gamma_c*Rc) * M9())-M0)
    Fbypass = m0_(M0)*(alpha*a0/(1+alpha))*((1+f_B_())*(   np.sqrt(gamma_t*Rt*TB_ratio(M0)*gamma_c*Rc) * M19(M0))-M0) # P0=P9, so last term all goes to 0.
    return Fcore+Fbypass

def ST():
    " Calculates specific thrust "
    return F(M0)/m0_(M0)

def SFC(M0):
    """ Calculates specifc fuel consumption """
    return ((m0_(M0)/(1 + alpha)) * (f_() + f_B_()*alpha + f_AB_())) / F(M0)


def setMode(mode, M0):
    if mode == 1:
        tau.update({"c": tau_c,
        "f": tau_f,
        "t": tau_t(M0),
        "0": tau_0(M0),
        "AB": 1,
        "B": 1,});
        pi.update({"c": pi_c,
        "AB": 1,
        "f": pi_f,
        "t": pi_t(M0),
        "B": 1,
        "0": pi_0(M0)})
        f_AB = 0
        f_B = 0

    elif mode == 2:
        tau.update({
        "c": tau_c,
        "f": tau_f,
        "t": tau_t(M0),
        "0": tau_0(M0),
        "AB": tau_AB,
        "B": 1})
        pi.update({"c": pi_c,
        "AB": pi_AB,
        "f": pi_f,
        "t": pi_t(M0),
        "B": 1,
        "0": pi_0(M0) } )
        f_AB = f_AB_()
        f_B = 0

    elif mode == 3:
        tau.update({"c": tau_c,
        "f": tau_f,
        "t": tau_t(M0),
        "0": tau_0(M0),
        "AB": tau_AB,
        "B": tau_B})
        pi.update({"c": pi_c,
        "AB": pi_AB,
        "f": pi_f,
        "t": pi_t(M0),
        "B": pi_B(),
        "0": pi_0(M0) })
        f = f_()
        f_B = f_B_()
        f_AB = f_AB_()

    elif mode == 4:
        # - Also in ramjet mode, there is shock losses, so Pt2/Pt0 is something else. See lecture 8.
        tau.update({"c": 1,
        "f": 1,
        "t": 1,
        "0": tau_0(M0),
        "AB": tau_AB,
        "B": tau_B})
        pi.update({"c": 1,
        "AB": pi_AB,
        "f": 1,
        "t": 1,
        "B": pi_B(),
        "0": pi_0(M0) })
        f = f_()
        f_B = f_B_()
        f_AB = f_AB_()

""" Main Code """

#M0 = np.linspace(0.1, 10, 10)
M0 = 3
# Initiate
tau['0'] = tau_0(M0)
tau['t'] = tau_t(M0)
tau['B'] = tau_B
tau['AB'] = tau_AB_()
pi["t"] = pi_t(M0)
pi["B"] = pi_B()
pi['0'] = pi_0(M0)
"""Note to selfs: When changing mode, use setMode(mode, M0) """

"""Make the plot for performance. - thrust? SFC? """

print(tau)

'''Tests'''
setMode(1, 3)
print("Mode 1")
print(tau)
print(pi)

print(SFC(M0))
print(ST())

"""
setMode(2, 3)
print("Mode 2")
print(tau)
print(pi)

setMode(4, 3)
print("Mode 4")
print(tau)
print(pi)
"""











