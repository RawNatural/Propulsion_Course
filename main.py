import numpy as np
#m0 = something; alpha = 1; mC = m0/(1+alpha); mB = m0*alpha/(1+alpha)

""" General Assumptions """
c_pc = 1004 # [J/kg/K]
c_pt = 1239 # [J/kg/K]
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
A_2/A_1 = 2.5 # Diffuser Area Ratio

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
f = mf/mc # Burner Fuel Air Ratio
tau_AB = T_t7/T_t5
f_AB = mf_AB/mc # Afterburner Fuel Air Ratio

""" Bypass Modelling Parameters """
pi_f = 2.0 # Fan Stagnation Pressure Ratio
e_f = 2.0 # Fan Polytropic Efficiency
f_B = mf_B/mB # Bypass Fuel Air Ratio
tau_B = T_t14/T_t13
eff_B = eff_AB # Otherwise tau_B = 1

R = c_pc*(1-1/gamma_c)
a0 = np.sqrt(gamma_c*R*T_0)
print(a0)

# For now
alpha = 1

m0 = something; alpha = 1; mC = m0/(1+alpha); mB = m0*alpha/(1+alpha)

#x = False

#tau_0 = x; tau_b = x; tau_lambda = x; tau_lambAB = x; tau_B = x; tau_AB = x; tau_f = x;

#tau_c = pi_C*(gamma-1)/(gamma*e_c); tau_t = pi_T*(e_t*(gamma-1)/gamma);  #which gamma?
#tau_c = pi_C*(gamma-1)/(gamma*e_c); 

""" Key Notes
    From Assignment Task Sheet: P0 = P9 = P19, therefore P0/P9 = P0/P19 = 1 """
M0 = varies;

def TB_ratio():
    """ Bypass Temperature Ratio T19/T0 """ 
    TB_ratio = tau_fn*tau_B*tau_f*tau_0
    return TB_ratio

def TC_ratio_():
def TC_ratio():
    """ Core Temperature Ratio T9/T0 """
    TC_ratio = tau_lambda*tau_t*tau_AB*(1/c_ratio)
    TC_ratio = tau_lambda*tau_t*(1/c_ratio) / ((1*pi0*pi_dmax*pi_C*pi_b*pi_t*pi_n)**((gamma_t-1)/gamma_t))
    return TC_ratio

def PB_ratio_():
    """ Pressure stagnation ratio across Bypass Pt19/P19 from """
    """ Pressure stagnation ratio across Bypass Pt19/P19 """
    return 

def PC_ratio():
    """ Pressure stagnation ratio across Core Pt9/P9 """ 
    return 

def M9():
    """ Mach number at point 9 """
    return (2/(gamma_t-1))*((1*pi0*pi_dmax*pi_C*pi_b*pi_t*pi_n)**((gamma_t-1)/gamma_t)-1)

def M19():
    """ Mach number at point 19 """
    return 

def m0(M0):
    rho0 = 100000/((M0*a0)**2)
    return rho0*a0*A_0*M0

def f():
    """ Calculates fuel air ratio for burner """
    f = (((eff_b*H)/(c_pc*T_0))**-1)*(tau_lambda - tau_c*tau_0)
    return f

def f_B():
    """ Caclulates fuel air ratio for bypass """
    f_B = ((((eff_B*H)/(c_pc*T_0)) - tau_B*tau_f*tau_0)**-1)*(tau_B*tau_f*tau_0 - tau_f*tau_0)
    return f_B

def f_AB():
    """ Calculates fuel air ratio for after burner """
    f_AB = ((((eff_AB*H)/(c_pt*T_0)) - tau_AB*tau_t*tau_lambda)**-1)*(1+f_())*(tau_AB*tau_t*tau_lambda - tau_t*tau_lambda)
    f_AB = ((((eff_AB*H)/(c_pt*T_0)) - tau_AB*tau_t*tau_lambda)**-1)*(1+f())*(tau_AB*tau_t*tau_lambda - tau_t*tau_lambda)
    return f_AB

def tau_t():
    return 1-1/(eff_m*(1+f()))*tau_0/tau_lambda*(tau_c-1);

def pi_t():
    return tau_t() ** (gamma_t/((gamma_t-1)*e_t))

def F():
    """ Calculates overall thrust """
    Fcore = m0*(a0/(1+alpha)*((1+f_())*(   np.sqrt(gamma_t*Rt*TC_ratio_()*gamma_c*Rc) * M9 + np.sqrt(gamma_c*Rt*TC_ratio_()*gamma_t*Rc) * (1-PC_ratio_())/(M9*gamma_c)  )-M0))
    Fbypass = m0*(alpha*a0/(1+alpha)*((1+f_B_())*(   np.sqrt(gamma_t*Rt*TB_ratio_()*gamma_c*Rc) * M19 + np.sqrt(gamma_c*Rt*TB_ratio_()*gamma_t*Rc) * (1-P0/P19)/(M19*gamma_c)  )-M0))
    Fcore = m0*(a0/(1+alpha))*((1+f()+f_AB())*(   np.sqrt(gamma_t*Rt*TC_ratio()*gamma_c*Rc) * M9)-M0)
    Fbypass = m0*(alpha*a0/(1+alpha))*((1+f_B())*(   np.sqrt(gamma_t*Rt*TB_ratio()*gamma_c*Rc) * M19)-M0) # P0=P9, so last term all goes to 0.
    return Fcore+Fbypass

def ST():
    " Calculates specific thrust "
    return F_()/m0

def SFC():
    """ Calculates specifc fuel consumption """
    return ((m0/(1 + alpha)) * (f() + f_B()*alpha + f_AB())) / F()
