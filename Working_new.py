import numpy as np
import matplotlib.pyplot as plt

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
f_st = 2.38 # Stoichiometric fuel ratio
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
tau_AB = 1.2 #T_t7/T_t5
pi_t = 0.99 # Turbine Stagnation Pressure Ratio Assumed

""" Bypass Modelling Parameters """
pi_f = 2.0 # Fan Stagnation Pressure Ratio
pi_fn = 0.95 # Nozzle Stagnation Pressure Ratio
tau_fn = 1 # Nozzle Stagnation Temperature Ratio
e_f = 2.0 # Fan Polytropic Efficiency
#f_B = mf_B/mB # Bypass Fuel Air Ratio
#pi_B = P_t14/P_t13
eff_B = eff_AB # Otherwise tau_B = 1
tau_B = 1.2 #T_t14/T_t13

""" Guesses """
M19 = 5
M_14 = 4
M_13 = 3

Rc = c_pc*(1-1/gamma_c)
Rt = c_pt*(1-1/gamma_t)
a0 = np.sqrt(gamma_c*Rc*T_0)
alpha = 1
M0 = np.linspace(0.1, 10, 10)

def m0_(M0):
    rho0 = 100000/((M0*a0)**2)
    m_0 = rho0*a0*A_0*M0
    global mC, mB
    mC = m_0/(1+alpha); mB = m_0 * alpha / (1+alpha)
    return m_0

def tau_c_():
    """ Calculates temperature stagnation ratio Tt3/Tt2 across compressor """
    tau_c = pi_c**((gamma_c - 1)/(gamma_c*e_c))
    
    return tau_c
print(f"tau_c = {tau_c_()}")

def tau_f_():
    """ Calculates temperature stagnation ratio Tt13/Tt2 across fan """
    tau_f = pi_f**((gamma_c - 1)/(gamma_c*e_f))
    return tau_f
print(f"tau_f = {tau_f_()}")

def tau_t_():
    tau_t = pi_t**(e_t*(gamma_t - 1)/gamma_t)
    return tau_t
print(f"tau_t = {tau_t_()}")

def tau_0_(M0):
    tau_0 = 1 + (gamma_c - 1)/2 * M0**2
    return tau_0
print(f"tau_0 = {tau_0_(M0)}")

def pi_B_(M_13, M_14):
    """ Calculates stagnation pressure ratio Pt14/Pt13 across Bypass """
    pi_B = ((1 + gamma_c*M_13**2)/(1 + gamma_c*M_14**2))*((1 + (gamma_c - 1)/2 *M_14**2)/(1 + (gamma_c - 1)/2 * M_13**2))**(gamma_c/(gamma_c - 1))
    return pi_B
print(f"pi_B = {pi_B_(M_13, M_14)}")

def pi_0_(M0):
    """ Calculates stagnation pressure ratio across the inlet Pt0/P0 """
    pi_0 = tau_0_(M0)**(gamma_c/(gamma_c - 1))
    return pi_0
print(f"pi_0 = {pi_0_(M0)}")

"""Set all tau's and pi's in a dictionary"""
tau = {
    "c": tau_c_(),
    "lambda": tau_lambda,
    "n": tau_n,
    "fn": tau_fn,
    "f": tau_f_(),
    "0": tau_0_(M0), #functions of M0
    "t": tau_t_(),
    "AB": tau_AB,
    "B": None,
}
print(f"Tau values: {tau}")

pi = {
    "c": pi_c,
    "d": pi_dmax,
    "f": pi_f,
    "b": pi_b,
    "AB": pi_AB,
    "n": pi_n,
    "fn": pi_fn,
    "f": pi_f,
    "t": pi_t,
    "B": pi_B_(M_14, M_13),
    "0": pi_0_(M0)
}
print(f"Pi values = {pi}")

def TC_ratio():
    """ Core Temperature Ratio T9/T0 """
    TC_ratio = tau['lambda']*tau_t_()*(1/c_ratio) / ((1*pi['0']*pi['d']*pi['c']*pi['b']*pi['t']*pi['n']*pi['f'])**((gamma_t-1)/gamma_t))
    return TC_ratio
print(f"TC_ratio = {TC_ratio()}")

def PB_ratio_():
    """ Pressure stagnation ratio across Bypass Pt19/P19 """
    PB_ratio = pi['fn']*pi['B']*pi['f']*pi['0']
    return PB_ratio
print(f"PB_ratio = {PB_ratio_()}")    

def PC_ratio():
    """ Pressure stagnation ratio across Core Pt9/P9 """ 
    print(f"PC_ratio = {PC_ratio}")
    return pi['n']*pi['AB']*pi['t']*pi['b']*pi['c']*pi['f']*pi['0']

def M9():
    """ Mach number at point 9 """
    return (2/(gamma_t-1))*((1*pi['0']*pi['d']*pi['c']*pi['b']*pi['t']*pi['n'])**((gamma_t-1)/gamma_t)-1)

def M19(M0):
    """ Mach number at point 19 """
    return M0*1.5 #change for something later

def f_():
    """ Calculates fuel air ratio for burner """
    f_numerator = tau['lambda'] - tau['c']*tau['0']
    f_denominator = (eff_b*H)/(c_pc*T_0)
    global F
    f = f_numerator/f_denominator
    return f
print(f" f = {f_()}")

def f_B_():
    """ Caclulates fuel air ratio for bypass """
    # numerator = tau['B']*tau['f']*tau['0'] - tau['f']*tau['0']
    # denominator = ((eff_B*H)/(c_pc*T_0)) - tau['B']*tau['f']*tau['0']
    global f_B
    f_B = np.linspace(0.1, 1, 10)
    return f_B
print(f"f_B = {f_B_()}")

def tau_B_():
    """ Calculates tau['B'] """
    tau_B_num = f_B_()*eff_B*H/(c_pt*T_0) + tau_f_()*tau_0_(M0)
    tau_B_den = tau_f_()*tau_0_(M0) + f_B_()*tau_f_()*tau_0_(M0)
    tau_B = tau_B_num/tau_B_den
    return tau_B
print(f"tau_b = {tau_B_()}")

def TB_ratio():
    """ Bypass Temperature Ratio T19/T0 """ 
    TB_ratio = tau['fn']*tau_B_()*tau['f']*tau['0'] / (1*pi['0']*pi['B']*pi['f']*pi['fn'])
    return TB_ratio
print(f"TB_ratio = {TB_ratio()}")

def f_AB_():
    """ Calculates fuel air ratio for after burner """
    global f_AB
    f_AB = f_st - f_()
    return f_AB
print(f"f_AB = {f_AB_()}")

def tau_AB_():
    tau_AB_num = f_AB_()*(eff_AB*H/(c_pt*T_0)) + tau_t_()*tau_lambda*(1 + f_())
    tau_AB_den = tau_t_()*tau_lambda*(1 + f_()) + f_AB_()*tau_t_()*tau_lambda
    tau_AB = tau_AB_num/tau_AB_den
    return tau_AB
print(f"tau_AB = {tau_AB_()}")

def F(M0):
    """ Calculates overall thrust """
    Fcore = m0_(M0)*(a0/(1+alpha))*((1+f_()+f_AB_())*(np.sqrt(gamma_t*Rt*TC_ratio()*gamma_c*Rc) * M9())-M0) # P0=P9, so last term all goes to 0.
    Fbypass = m0_(M0)*(alpha*a0/(1+alpha))*((1+f_B_())*(np.sqrt(gamma_t*Rt*TB_ratio()*gamma_c*Rc) * M19(M0))-M0) # P0=P19, so last term all goes to 0.
    F_total = Fcore + Fbypass
    return F_total
print(f"F_total = {F(M0)}")

def ST(M0):
    " Calculates specific thrust "
    ST = F(M0)/m0_(M0)
    return ST 
print(f"Specific Thrust = {ST(M0)}")

def SFC(M0):
    """ Calculates specifc fuel consumption """
    SFC = ((m0_(M0)/(1 + alpha)) * (f_() + f_B_()*alpha + f_AB_())) / F(M0)
    return SFC
print(f"Specific Fuel Consumption = {SFC(M0)}")

def setMode(mode, M0):
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
        "t": pi_t,
        "B": 1,
        "0": pi_0_(M0)})
        f_AB = 0
        f_B = 0
    elif mode == 2:
        tau.update({
        "c": tau_c_(),
        "f": tau_f_(),
        "t": tau_t_(),
        "0": tau_0_(M0),
        "AB": tau_AB,
        "B": 1})
        pi.update({"c": pi_c,
        "AB": pi_AB,
        "f": pi_f,
        "t": pi_t,
        "B": 1,
        "0": pi_0_(M0)})
        f_AB = f_AB_()
        f_B = 0
    elif mode == 3:
        tau.update({"c": tau_c_(),
        "f": tau_f_(),
        "t": tau_t_(),
        "0": tau_0_(M0),
        "AB": tau_AB,
        "B": tau_B})
        pi.update({"c": pi_c,
        "AB": pi_AB,
        "f": pi_f,
        "t": pi_t,
        "B": pi_B_(M_14, M_13),
        "0": pi_0_(M0)})
        f = f_()
        f_B = f_B_()
        f_AB = f_AB_()
    elif mode == 4:
        # - Also in ramjet mode, there is shock losses, so Pt2/Pt0 is something else. See lecture 8.
        tau.update({"c": 1,
        "f": 1,
        "t": 1,
        "0": tau_0_(M0),
        "AB": tau_AB,
        "B": tau_B})
        pi.update({"c": 1,
        "AB": pi_AB,
        "f": 1,
        "t": 1,
        "B": pi_B_(M_14, M_13),
        "0": pi_0_(M0)})
        f = f_()
        f_B = f_B_()
        f_AB = f_AB_()

if __name__ == "__main__":
    M0_range = np.linspace(0.1, 10, 50)  # finer resolution
    modes = [1, 2, 3, 4]
    mode_labels = {
        1: "Pure Turbofan",
        2: "Afterburning Turbofan",
        3: "Turboramjet",
        4: "Ramjet"
    }

    # Plot Specific Thrust (ST) vs Mach Number
    plt.figure(figsize=(10, 6))
    for mode in modes:
        ST_values = []
        for M in M0_range:
            setMode(mode, M)
            val = ST(M)
            if hasattr(val, "__len__") and not isinstance(val, str):
                val = np.mean(val)
            ST_values.append(val)

        plt.plot(M0_range, ST_values, label=mode_labels[mode])

    plt.grid(True)
    plt.xlabel("Flight Mach Number (M0)")
    plt.ylabel("Specific Thrust [Ns/kg]")
    plt.title("Specific Thrust vs Flight Mach Number for Different Modes")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Plot Specific Fuel Consumption (SFC) vs Mach Number
    plt.figure(figsize=(10, 6))
    for mode in modes:
        SFC_values = []
        for M in M0_range:
            setMode(mode, M)
            val = SFC(M)
            if hasattr(val, "__len__") and not isinstance(val, str):
                val = np.mean(val)
            SFC_values.append(val)

        plt.plot(M0_range, SFC_values, label=mode_labels[mode])

    plt.grid(True)
    plt.xlabel("Flight Mach Number (M0)")
    plt.ylabel("Specific Fuel Consumption [kg/(Ns)]")  # or your appropriate units
    plt.title("Specific Fuel Consumption vs Flight Mach Number for Different Modes")
    plt.legend()
    plt.tight_layout()
    plt.show()
