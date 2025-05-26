import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve
from scipy.optimize import root_scalar

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

"""Changes:
tauAB in TC ratio

"""

""" General Assumptions """
c_pc = 1004 # [J/kg/K]
c_pt = 1239 # [J/kg/K]
c_ratio = c_pt/c_pc 
gamma_c = 1.4
gamma_t = 1.3
q_0 = 50000 # [Pa]
T_0 = 250 # [K]
A_0 = 14.4 # [m2]
C_D = 0.03 # Drag Coefficient
A_ref = 382 # [m2]

""" Fuel Properties """
H = 120e6 # [J/kg]
#f_st = 2.38 
f_st = 0.0291 # Stoichiometric fuel ratio
phi = 1 # Maximum Equivalence Ratio

""" Inlet Modelling """
pi_dmax = 0.96
A_ratio =  2.5 # Diffuser Area Ratio # A_2/A_1

""" Core Modelling Parameters """
pi_c = 30 # Compressor Stagnation Pressure Ratio 
e_c = 0.91 # Compressor Polytropic Efficiency
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
#tau_AB = 1.2 #T_t7/T_t5
f_Bmax = f_st

""" Bypass Modelling Parameters """
pi_f = 2.0 # Fan Stagnation Pressure Ratio
pi_fn = 0.95 # Nozzle Stagnation Pressure Ratio
tau_fn = 1 # Nozzle Stagnation Temperature Ratio
e_f = 2.0 # Fan Polytropic Efficiency
#f_B = mf_B/mB # Bypass Fuel Air Ratio
#pi_B = P_t14/P_t13
eff_B = eff_AB # Otherwise tau_B = 1
#tau_B = 1.2 #T_t14/T_t13

""" Guesses """
#M19 = 5
#M_14 = 4
#M_13 = 3

Rc = c_pc*(1-1/gamma_c)
Rt = c_pt*(1-1/gamma_t)
a0 = np.sqrt(gamma_c*Rc*T_0)
alpha = 1

def m0_(M0):
    """ Calculates the mass flow rate m0  """
    #rho0 = 100000/((M0*a0)**2)
    #m_0 = rho0*a0*A_0*M0
    m_0 = 2* q_0 /(M0*a0) * A_0 * 2
    mC = m_0/(1+alpha); mB = m_0 * alpha / (1+alpha)
    return m_0

def tau_c_():
    """ Calculates temperature stagnation ratio Tt3/Tt2 across compressor """
    tau_c = pi['c']**((gamma_c - 1)/(gamma_c*e_c))
    #tau_c = 1.5
    return tau_c

def tau_f_():
    """ Calculates temperature stagnation ratio Tt13/Tt2 across fan """
    tau_f = pi['f']**((gamma_c - 1)/(gamma_c*e_f))
    return tau_f

def tau_0_(M0):
    """ Calculates stagnation temperature ratio across inlet """
    tau_0 = 1 + (((gamma_c - 1)/2) * M0**2)
    return tau_0

def pi_t_(): # Turbine Stagnation Pressure Ratio
    #tau_t = pi_t**(e_t*(gamma_t - 1)/gamma_t) #rearrange to:
    pi_t = tau['t']**((gamma_t)/(e_t*(gamma_t-1)))
    return pi_t

def pi_0_(M0):
    """ Calculates stagnation to normal pressure ratio at point 0: Pt0/P0 """
    pi_0 = tau_0_(M0)**(gamma_c/(gamma_c - 1))
    return pi_0

def pi_i_():
    #return (((gamma_c+1)*M0**2)/(2+(gamma_c-1)*M0**2))**(gamma_c/(gamma_c-1)) / (((2*gamma_c/(gamma_c+1))*M0**2 - (gamma_c-1)/(gamma_c+1))**(1/(gamma_c-1)))
    return pi['d']/pi_dmax

def pi_d_(M0):
    """ Calculates stagnation pressure ratio across the inlet Pt2/Pt0 """
    pi_ds = []
    if not isinstance(M0, float):
        for M in M0:
            eff_r = 1 if M <= 1 else 1-0.075*(M-1)**1.35 if M < 5 else 800/(M**4+985)
            pi_d = eff_r * pi_dmax
            pi_ds.append(pi_d)
        return np.array(pi_ds)
    else:
        eff_r = 1 if M0 <= 1 else 1-0.075*(M0-1)**1.35 if M0 < 5 else 800/(M0**4+985)
        pi_d = eff_r * pi_dmax
        return pi_d

"""Set all tau's and pi's in a dictionary"""
tau = {
    "c": None, #initialis later
    "lambda": tau_lambda,
    "n": tau_n,
    "fn": tau_fn,
    "f": None, #initialise later
    "0": None, #functions of M0 - iniitiate upon M0 creation
    "t": None, #initialise later
    "AB": None, #initialise later
    "B": None #Initialise later 
}


pi = {
    "c": pi_c,
    "f": pi_f,
    "b": pi_b,
    "AB": pi_AB,
    "n": pi_n,
    "fn": pi_fn,
    "f": pi_f,
    "i": 1,
    "t": None, #Initialise later
    "B": None, # Set later pi_B_(M_14, M_13),
    "d": None, #initialise as function of M0 later
    "0": None #functions of M0  - iniitiate upon M0 creation
}

f = {
    "f": None,
    "AB": None,
    "B": None,
    "Bmax": f_st
} #Initialise all later
current_mode = False

def TC_ratio(): #verified
    """ Core Temperature Ratio T9/T0 """ 
    TC_ratio = tau['lambda']*tau['t']*tau['AB']*(1/c_ratio) / ((1*pi['0']*pi['i']*pi['d']*pi['c']*pi['b']*pi['t']*pi['n'])**((gamma_t-1)/gamma_t))
    return TC_ratio

def PB_ratio_():
    """ Pressure stagnation ratio across Bypass Pt19/P19 """
    PB_ratio = pi['fn']*pi['B']*pi['f']*pi['0']*pi['d']*pi['i']
    return PB_ratio

def PC_ratio():
    """ Pressure stagnation ratio across Core Pt9/P9 """ 
    return pi['n']*pi['AB']*pi['t']*pi['b']*pi['c']*pi['f']*pi['0']*pi['i']

def M9():
    """ Mach number at point 9 """
    return np.sqrt(abs((2/(gamma_t-1))*((1*pi['0']*pi['d']*pi['c']*pi['b']*pi['t']*pi['n']*pi['i'])**((gamma_t-1)/gamma_t)-1)))

def MFP(M0,gamma_c,Rc):
    """
    Return the mass flow parameter.
    
    Arguments:
    M0: (float) Flight Mach number. #should be M0
    k: (float) Ratio of specific heats. #should be gamma_c
    R: (float) J/kg.K, gas constant. # Rc
    """
    return np.sqrt(gamma_c/Rc)*M0*(1.0 + 0.5*(gamma_c - 1.0)*M0**2)**(-0.5*(gamma_c + 1.0)/(gamma_c - 1.0))

def M2(M0):
    """ Calculates M2 """
    M1 = abs((1 + (gamma_c - 1)/2 * M0**2)/(gamma_c*M0**2 - (gamma_c - 1)/2))**0.5
    M2s = []
    if not isinstance(M1, float):
        for M_1 in M1:
            M_2 = fsolve(lambda M: MFP(M,gamma_c,Rc) - MFP(M_1,gamma_c,Rc)*(1/A_ratio), 0.2)[0]
            M2s.append(M_2)
        M2_ = np.array(M2s)
    else:
        M2_ = fsolve(lambda M: MFP(M,gamma_c,Rc) - MFP(M1,gamma_c,Rc)*(1/A_ratio), 0.2)[0]
    return M2_


def M19_(M0):
    """ Mach number at point 19 """
    M19 = np.sqrt(abs((((PB_ratio_())**((gamma_t - 1)/gamma_t) - 1)*2/(gamma_t - 1))))
    return M19 

def f_():
    """ Calculates fuel air ratio for burner """
    f_numerator = tau_lambda - tau['c']*tau['0']
    f_denominator = (eff_b*H)/(c_pc*T_0) - tau_lambda
    return f_numerator/f_denominator

def f_B_():
    """ Caclulates fuel air ratio for bypass """
    #numerator = tau['B']*tau['f']*tau['0'] - tau['f']*tau['0']
    #denominator = ((eff_B*H)/(c_pc*T_0)) - tau['B']*tau['f']*tau['0']
    #just set to f_st for now # limit to choking.
    #f_B = f_st if f_st <= f['Bmax'] else f['Bmax']
    f_B = f_st
    return f_B

def tau_t_(): #uses f()
    tau_t = 1 - tau['0']*((tau['c']-1)+alpha*(tau['f']-1)) / (tau_lambda*(1+f['f'])*eff_m)
    return tau_t

def tau_B_():
    """ Calculates tau['B'] """ 
    tau_B_num = f["B"]*eff_B*H/(c_pc*T_0) + tau['f']*tau['0']
    tau_B_den = tau['f']*tau['0'] + f["B"]*tau['f']*tau['0']
    tau_B = tau_B_num/tau_B_den
    return tau_B

def M13_(M0):
    #return np.sqrt(abs((((pi_f*pi['0'])**((gamma_t - 1)/gamma_t) - 1)*2/(gamma_t - 1))))
    return M2(M0)

'''def solve_M14(M14, M13, tau_B_val):
    if M14 <= 0:
        return 1e6  # Large penalty to keep M14 positive
    #M14 = M14[0] if isinstance(M14, (np.ndarray, list)) else M14
    M_14 = ((((1 + ((gamma_c - 1) / 2) * M14**2) / (1 + ((gamma_c - 1) / 2) * M13**2)) 
             * (M14**2 / M13**2) ) * ((1 + gamma_c * M13**2) / (1 + gamma_c * M14**2))**2) - tau_B_val#1.2#tau['B']
    if M_14 <= 0:
        return 1e6  # Large penalty to keep M14 positive
    return M_14

def M14_(M0):
    M13_vals = M13_(M0)
    tau_B_vals = tau['B']
    M14_guess = 0.2
    M14s = []

    # Convert scalars to 1-element lists for uniform handling
    if np.isscalar(M13_vals):
        M13_vals = [M13_vals]
    if np.isscalar(tau_B_vals):
        tau_B_vals = [tau_B_vals] * len(M13_vals)  # broadcast scalar to match M13

    for M_13, tauB in zip(M13_vals, tau_B_vals):
        #sol = root_scalar(solve_M14, args=(M_13, tauB), bracket=[0.05, 0.4])
        #M_14 = sol.root if sol.converged else np.nan
        M_14 = fsolve(solve_M14, M14_guess, args=(M_13, tauB))[0]
        if M_14 < 0.99:
            f['B'] = f_st
        while M_14 > 0.99:
            f['B'] -= 0.001
            #sol = root_scalar(solve_M14, args=(M_13, tauB), bracket=[0.05, 0.4])
            #M_14 = sol.root if sol.converged else np.nan
            M_14 = fsolve(solve_M14, M14_guess, args=(M_13, tauB))[0]
        M14s.append(M_14) 
        f['Bmax'] = f['B']
    return M14s[0] if len(M14s) == 1 else np.array(M14s)'''

from scipy.optimize import fsolve
import numpy as np

def solve_M14(M14, M13, tau_B_val):
    if M14 <= 0:
        return 1e6  # Penalize invalid guesses
    val = ((((1 + ((gamma_c - 1) / 2) * M14**2) / (1 + ((gamma_c - 1) / 2) * M13**2)) 
           * (M14**2 / M13**2)) * ((1 + gamma_c * M13**2) / (1 + gamma_c * M14**2))**2) - tau_B_val
    return val

def M14_(M0):
    M13_vals = M13_(M0)
    tau_B_vals = tau['B']
    M14_guess = 0.2
    M14s = []

    if np.isscalar(M13_vals):
        M13_vals = [M13_vals]
    if np.isscalar(tau_B_vals):
        tau_B_vals = [tau_B_vals] * len(M13_vals)

    for M_13, tauB in zip(M13_vals, tau_B_vals):
        def wrapped_f(M14): return solve_M14(M14, M_13, tauB)

        try:
            M_14 = fsolve(wrapped_f, M14_guess)[0]
        except Exception as e:
            print(f"fsolve failed for M13={M_13}, tauB={tauB}: {e}")
            M_14 = np.nan

        if not np.isfinite(M_14) or M_14 <= 0:
            print(f"Warning: Invalid M14 value returned: {M_14}")
            M_14 = np.nan


        """max_iters = 20
        iters = 0
        while M_14 > 0.99 and iters < max_iters:
            f['B'] -= 0.0000001
            M_14 = fsolve(wrapped_f, M14_guess)[0]
            iters += 1
"""
        M14s.append(M_14)
        f['Bmax'] = f['B']
        #print(f['Bmax'])

    return M14s[0] if len(M14s) == 1 else np.array(M14s)


def getCombustion_T_and_Ps(M13, M14):
    Tt13 = tau['0']*tau['f']*T_0
    T13 = Tt13/(1+(gamma_c-1)/2*M13**2)
    Tt14 = Tt13*tau['B']
    T14 = Tt14/(1+(gamma_c-1)/2*M14**2)
    Pt13 = pi['0']*pi['f']*pi['d']*T_0
    P13 = Pt13/((Tt13/T13)**(gamma_c/(gamma_c-1)))
    Pt14 = Pt13*pi['B']
    P14 = Pt14/((Tt14/T14)**(gamma_c/(gamma_c-1)))
    return [[T13, T14], [P13, P14]]

def pi_B_(M_13, M_14):
    """ Calculat2es stagnation pressure ratio Pt14/Pt13 across Bypass """
    pi_B = ((1 + gamma_c*M_13**2)/(1 + gamma_c*M_14**2))*((1 + (gamma_c - 1)/2 *M_14**2)/(1 + (gamma_c - 1)/2 * M_13**2))**(gamma_c/(gamma_c - 1))
    return pi_B

def TB_ratio():
    """ Bypass Temperature Ratio T19/T0 """ 
    TB_ratio = tau['fn']*tau['B']*tau['f']*tau['0'] / (1*pi['0']*pi['B']*pi['f']*pi['fn']*pi['d']*pi['i'])
    return TB_ratio

def f_AB_():
    """ Calculates fuel air ratio for after burner """
    f_AB = f_st - f_()
    return f_AB

def tau_AB_():
    tau_AB_num = f["AB"]*(eff_AB*H/(c_pt*T_0)) + tau['t']*tau_lambda*(1/c_ratio)*(1 + f_())
    tau_AB_den = tau['t']*tau_lambda*(1/c_ratio)*(1 + f['f'] + f["AB"])
    tau_AB = tau_AB_num/tau_AB_den
    #tau_AB = 1.2
    return tau_AB

def F(M0):
    """ Calculates overall thrust """
    Fcore = m0_(M0)*((a0/(1+alpha))*((1+f['f']+f['AB'])*(np.sqrt(gamma_t*Rt*TC_ratio()/(gamma_c*Rc)) * M9())-M0)) # P0=P9, so last term all goes to 0.
    Fbypass = m0_(M0)*((alpha*a0/(1+alpha))*((1+f["B"])*(np.sqrt(gamma_t*Rt*TB_ratio()/(gamma_c*Rc)) * M19_(M0))-M0)) # P0=P19, so last term all goes to 0.
    F_total = Fcore + Fbypass
    return F_total

def ST(M0):
    " Calculates specific thrust "
    ST = F(M0)/m0_(M0)
    return ST 

def SFC(M0):
    """ Calculates specifc fuel consumption """
    SFC = ((m0_(M0)/(1 + alpha)) * (f_() + f['B']*alpha + f["AB"])) / F(M0)
    return SFC

def setMode(mode, M0):
    global current_mode
    current_mode = mode
    if mode == 1:
        tau.update({"c": tau_c_(),
        "f": tau_f_(),
        "t": tau_t_(),
        "0": tau_0_(M0),
        "AB": 1,
        "B": 1,})
        pi.update({"c": pi_c,
        "AB": 1,
        "f": pi_f,
        "t": pi_t_(),
        "B": 1,
        "d": pi_d_(M0),
        "i": 1,
        "0": pi_0_(M0)})
        f.update({
            "f": f_(),
            "B": 0,
            "AB": 0
        })
    elif mode == 2:
        tau.update({
        "c": tau_c_(),
        "f": tau_f_(),
        "t": tau_t_(),
        "0": tau_0_(M0),
        "AB": tau_AB_(),
        "B": 1})
        pi.update({"c": pi_c,
        "AB": pi_AB,
        "f": pi_f,
        "t": pi_t_(),
        "B": 1,
        "d": pi_d_(M0),
        "i": 1,
        "0": pi_0_(M0)})
        f.update({
            "f": f_(),
            "B": 0,
            "AB": f_AB_()
        })
    elif mode == 3:
        tau.update({"c": tau_c_(),
        "f": tau_f_(),
        "t": tau_t_(),
        "0": tau_0_(M0),
        "AB": tau_AB_(),
        "B": tau_B_()})
        pi.update({"c": pi_c,
        "AB": pi_AB,
        "f": pi_f,
        "d": pi_d_(M0),
        "i": 1,
        "t": pi_t_(),
        "B": pi_B_(M14_(M0), M13_(M0)),
        "0": pi_0_(M0)})
        f.update({
            "f": f_(),
            "B": f_B_(),
            "AB": f_AB_()
        })
    elif mode == 4:
        # - Also in ramjet mode, there is shock losses, so Pt2/Pt0 is something else. See lecture 8.
        tau.update({"c": 1,
        "f": 1,
        "t": 1,
        "0": tau_0_(M0),
        "AB": tau_AB_(),
        "B": tau_B_()})
        pi.update({"c": 1,
        "AB": pi_AB,
        "f": 1,
        "t": 1,
        "d": 1,
        "i": pi_i_(),
        "B": pi_B_(M14_(M0), M13_(M0)),
        "0": pi_0_(M0)})
        f.update({
            "f": f_(),
            "B": f_B_(),
            "AB": f_AB_()
        })
"""Initialise pi_B and tau_B, and tau_t pi_t"""

"""Part 3"""
def M_max_():
    M_range = np.linspace(1, 10, 200)
    Tt0 = tau_0_(M_range)*T_0
    for i in range(len(M_range)):
        if Tt0[i] >= T_t4max:
            return M_range[i-1]
        
def M_turb_limit_(M0):
    f = f_()
    for i in range(len(M0)):
        if f[i] < 0.00005:
            return M0[i-1]
 
''' def initialise(M0):
    """Initilise other tau and pi which are dependant on M0 (or pi0)"""
    pi['0'] = pi_0_(M0)
    tau['0'] = tau_0_(M0)
    f.update({
            "f": f_(),
            "B": f_B_(),
            "AB": f_AB_()
        })
    tau['t'] = tau_t_()
    pi['t'] = pi_t_(tau['t'])
    tau['B'] = tau_B_()
    pi['B'] = pi_B_(M13_(M0), M14_(M0))
    pi['d'] = pi_d_(M0)
    tau['AB'] = tau_AB_()'''

if __name__ == "__main__":

    modes = [1, 2, 3, 4]
    mode_labels = {
        1: "Pure Turbofan",
        2: "Afterburning Turbofan",
        3: "Turboramjet",
        4: "Ramjet"
    }

    """Find M_max"""
    M_max = M_max_()

    """Set mach range"""
    M0 = np.linspace(1, M_max, 200)  # finer resolution

    """Initilise other tau and pi which are dependant on M0 (or pi0)"""
    pi['0'] = pi_0_(M0)
    tau['0'] = tau_0_(M0)
    tau['c'] = tau_c_()
    tau['f'] = tau_f_()
    f.update({
            "f": f_(),
            "B": f_B_(),
            "AB": f_AB_()
        })
    tau['t'] = tau_t_()
    pi['t'] = pi_t_()
    tau['B'] = tau_B_()
    pi['B'] = pi_B_(M13_(M0), M14_(M0))
    pi['d'] = pi_d_(M0)
    tau['AB'] = tau_AB_()
    setMode(3, M0)
    print(M14_(M0))
    print(f['Bmax'])

    """Find M_turb_limit"""
    M_turb_limit = M_turb_limit_(M0)

    print("Successfully Inititialised!")

    # Plot Specific Thrust (ST) vs Mach Number
    plt.figure(figsize=(10, 6))
    m_14s = []
    m0s = []
    for mode in modes:
        new = True
        setMode(mode, M0)
        ST_values = np.full(len(M0), np.nan)
        i = 0
        for M in M0:
            setMode(mode, M)
            if new:
                setMode(mode,M)
                f['Bmax'] = f_st
                new = False
            val = ST(M) if M < M_turb_limit or mode == 4 and M < M_max  else np.nan
            if hasattr(val, "__len__") and not isinstance(val, str):
                val = np.mean(val)
            if val < 0.01:
                break
            ST_values[i] = val
            i += 1
            if mode == 3:
                #print(M14_(M))
                m0s.append(M0); m_14s.append(np.float64(M14_(M)));
        plt.plot(M0, ST_values, label=mode_labels[mode])



    plt.grid(True)
    plt.xlabel("Flight Mach Number (M0)")
    plt.ylabel("Specific Thrust [Ns/kg]")
    plt.title("Specific Thrust vs Flight Mach Number for Different Modes")
    #plt.axhline(0, color='black', linewidth=2.5)  # bold y=0 line
    plt.axvline(x=M_turb_limit, color='gray', linestyle='--')
    plt.text(M_turb_limit, 700, 'Mach turbine limit ', rotation=0, va='bottom', ha='right')
    plt.axvline(x=M_max, color='gray', linestyle='--')
    plt.text(M_max, 700, 'Mach max ', rotation=0, va='bottom', ha='right')
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot Specific Fuel Consumption (SFC) vs Mach Number
    plt.figure(figsize=(10, 6))
    for mode in modes:
        new = True
        SFC_values = []
        for M in M0:
            setMode(mode, M)
            if new:
                setMode(mode,M)
                new = False
            val = SFC(M) if M < M_turb_limit or mode == 4 and M < M_max else np.nan
            val = np.nan if SFC(M) < 0 else val
            if hasattr(val, "__len__") and not isinstance(val, str):
                val = np.mean(val)
            SFC_values.append(val)

        plt.plot(M0, SFC_values, label=mode_labels[mode])

    plt.grid(True)
    plt.xlabel("Flight Mach Number (M0)")
    plt.ylabel("Specific Fuel Consumption [kg/(Ns)]")  # or your appropriate units
    plt.title("Specific Fuel Consumption vs Flight Mach Number for Different Modes")
    plt.axvline(x=M_turb_limit, color='gray', linestyle='--')
    plt.text(M_turb_limit, 0.0007, 'Mach turbine limit ', rotation=0, va='bottom', ha='right')
    plt.axvline(x=M_max, color='gray', linestyle='--')
    plt.text(M_max, 0.00025, 'Mach max ', rotation=0, va='bottom', ha='right')
    plt.legend()
    plt.tight_layout()
    plt.show()


    """ Task 2a """
    
    setMode(3, M0)
    setMode(3, M0)

    L_Bb = 4 # Bypass Burner Length

    def V_13_():
        """ Calculates """
        V_13 = np.sqrt(gamma_c*Rc*T13)
        return V_13

    [T_vals, P_vals] = getCombustion_T_and_Ps(M13_(M0), M14_(M0))
    T13, T14 = T_vals
    P13, P14 = P_vals
    L_Bb = 4 # Bypass Burner Length

    def V_13_():
        V_13 = np.sqrt(gamma_c*Rc*T13)*M13_(M0)
        return V_13

    def condition():
        """ L_Bb/V_13 Time for mass flow through bypass"""
        conditions = []
        for i in range(0, len(V_13_())):
            cond = L_Bb/V_13_()[i]
            conditions.append(cond)
        return conditions
    
    A = 2.4e19
    B = 0
    Ta = 30000
    P0 = 1e5
    P = P13
    #R = 8.314
    R = Rc#8.314
    n=1

    def getkf(T13):
        return (A*T13**B) * np.exp(-Ta/(T13))

    def getc(T13):
        return P/(Rc*T13)
    
    kf = getkf(T13) # Forward reaction coefficient

    c = getc(T13) #K = getK(T);
    X_H2 = 0.296; X_O2 = 0.148; X_N2 = 1.88
    X_H20p = 0.347; X_N2p = 0.653

    X_vals = [X_H2, X_O2, X_N2, X_H20p, X_N2p]

    r0 = 2*kf*c**1.5 * X_vals[0] * X_vals[1] ** 0.5 #- 2*kr * c**2 * X_vals[3]; #Reaction rate at time t0
    # subtraction of products negligible

    tau_comb = c * X_vals[3] / r0 # Time taken to reach equilibrium

    time_flow_bypass = condition()  # Call once to avoid redundant computation


    M_bypass_burn = False

    for i in range(len(tau_comb)):
        if time_flow_bypass[i] > tau_comb[i]:
            if not isinstance(M_bypass_burn, float):
                M_bypass_burn = M0[i - 1]
    print(f"M_bypass_burn = {M_bypass_burn:.2f}")
    
    # Plot Winning Graph - Task 4
    ST_values = []
    SFC_values = []
    T_margins = []
    fs = []
    fABs = []
    D = q_0*C_D*A_ref
    for M in M0:
        if M < M_bypass_burn:
            setMode(2, M)
            setMode(2, M) #set twice to confirm initialisation correct
        elif M < M_turb_limit:
            setMode(3, M)
            setMode(3, M)
        else:
            setMode(4, M)
            setMode(4, M)

        #setMode(1, M)
        #initialise()
        fs.append(f_())
        fABs.append(f_AB_())
        F_ = F(M)
        val = ST(M)
        SFC_ = SFC(M)
        T_margin = F_/D-1
        #if hasattr(val, "__len__"):
        #    val = np.mean(val)
        ST_values.append(np.atleast_1d(val)[0])
        SFC_values.append(np.atleast_1d(SFC_)[0])
        T_margins.append(np.atleast_1d(T_margin)[0])


    plt.figure(figsize=(10, 6))
    plt.plot(M0, ST_values, label="Specific Thrust")
    plt.grid(True)
    plt.xlabel("Flight Mach Number (M0)")
    plt.ylabel("Specific Thrust [Ns/kg]")
    plt.axvline(x=M_turb_limit, color='gray', linestyle='--')
    plt.text(M_turb_limit, 500, 'Mach turbine limit ', rotation=0, va='bottom', ha='left')
    plt.axvline(x=M_max, color='gray', linestyle='--')
    plt.text(M_max, 500, 'Mach max ', rotation=0, va='bottom', ha='right')
    plt.axvline(x=M_bypass_burn, color='gray', linestyle='--')
    plt.text(M_bypass_burn, 500, 'Mach Bypass Burn ', rotation=0, va='bottom', ha='right')
    plt.title("Specific Thrust vs Flight Mach Number for Transitioning Modes")
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(M0, SFC_values, label="Specific Fuel Consumption")
    plt.grid(True)
    plt.xlabel("Flight Mach Number (M0)")
    plt.ylabel("Specific Fuel Consumption [N/(kg/s)]") 
    plt.title("Specific Fuel Consumption vs Flight Mach Number for Transitioning Modes")
    #plt.ylim(0, 0.0012)
    plt.axvline(x=M_turb_limit, color='gray', linestyle='--')
    plt.text(M_turb_limit, 0.0007, 'Mach turbine limit ', rotation=0, va='bottom', ha='left')
    plt.axvline(x=M_max, color='gray', linestyle='--')
    plt.text(M_max, 0.00025, 'Mach max ', rotation=0, va='bottom', ha='right')
    plt.axvline(x=M_bypass_burn, color='gray', linestyle='--')
    plt.text(M_bypass_burn, 0.00025, 'Mach Bypass Burn ', rotation=0, va='bottom', ha='right')
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(M0, T_margins, label="Thrust Margin")
    plt.grid(True)
    plt.xlabel("Flight Mach Number (M0)")
    plt.ylabel("Thrust Margin") 
    plt.title("Thrust Margin vs Flight Mach Number for Transitioning Modes")
    plt.legend()
    plt.tight_layout()
    plt.show()

    
    """plt.figure()
    plt.plot(T13, tau_comb, "-")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Reaction time (s)")
    plt.show()"""

    plt.figure(figsize=(10,6))
    plt.plot(M0, np.log(tau_comb), "-")
    plt.xlabel("Mach Number")
    plt.ylabel("Log of reaction time")
    plt.show()

    """Part 2b"""

    print(f"Temp14 range = {T14[-1]-T14[0]}")
    print(np.min(T14))
    print(np.max(T14))


    def getK(T_): # between mach 3-5
        expv = -27.5+(-8.145+27.5)*(T_-T14[0])/(T14[-1]-T14[0])     
        # Linear interpolation from -27.5 to -8.145 w.r.t. T14 between Mach 3- Mach 5. Used to find progress variable
        return np.exp(expv)

    def progress_variable():
        return (getK(T14)*np.sqrt(2)*(P0/P14))**(2/3)

    X_H2_dissociated = progress_variable()
    percentage_hydrogen_left = X_H2_dissociated / X_vals[0] * 100

    plt.figure(figsize=(10,6))
    plt.plot(M0, percentage_hydrogen_left, "-")
    plt.xlabel("Mach number")
    plt.ylabel("Percentage Hydrogen not burned")
    plt.show()

   
def margin_lim(M0):
    for i in range(len(M0)):
        if T_margins[i] < 1:
            return M0[i-1]

""" Task 5: Getting the margin limit value (which is above M_Max """
def margin_lim():
    M0 = np.linspace(M_max, 10, 100)
    setMode(4, M0)
    thrust_margin = (F(M0)/D) - 1
    for i, tm in enumerate(thrust_margin):
        if tm < 1:
            return M0[i-1]

print(f" Thrust Margin Limit = {margin_lim():.2f}")
