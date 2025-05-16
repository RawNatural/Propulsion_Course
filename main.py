import numpy as np
import matplotlib.pyplot as plt


pi_b = 0.99; pi_AB = 0.9; pi_n = 0.95; pi_B = 1; pi_dmax = 0.96; pi_C = 30;
eff_b = 0.99; eff_AB = 0.99; eff_B = eff_AB; eff_m = 0.98
e_t = 0.93; e_f = 0.93; e_c = 0.91
c_pc = 1004; c_pt = 12239; gamma_c = 1.4; gamma_t = 1.3;
c_ratio = c_pt/c_pc

m0 = 500; #change this value

alpha = 1; mC = m0/(1+alpha); mB = m0*alpha/(1+alpha)

#x = False

#tau_0 = x; tau_b = x; tau_lambda = x; tau_lambAB = x; tau_B = x; tau_AB = x; tau_f = x;

#tau_c = pi_C*(gamma-1)/(gamma*e_c); 

#M0 = varies;

def TB_ratio():
    """ Bypass Temperature Ratio T19/T0 """ 
    TB_ratio = tau_fn*tau_B*tau_f*tau_0
    return TB_ratio

def TC_ratio():
    """ Core Temperature Ratio T9/T0 """
    TC_ratio = tau_lambda*tau_t*(1/c_ratio) / ((1*pi0*pi_dmax*pi_C*pi_b*pi_t*pi_n)**((gamma_t-1)/gamma_t))
    return TC_ratio

def PB_ratio():
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
    f_AB = ((((eff_AB*H)/(c_pt*T_0)) - tau_AB*tau_t*tau_lambda)**-1)*(1+f())*(tau_AB*tau_t*tau_lambda - tau_t*tau_lambda)
    return f_AB

def tau_t():
    return 1-1/(eff_m*(1+f()))*tau_0/tau_lambda*(tau_c-1);

def pi_t():
    return tau_t() ** (gamma_t/((gamma_t-1)*e_t))

def F():
    """ Calculates overall thrust """
    Fcore = m0*(a0/(1+alpha))*((1+f()+f_AB())*(   np.sqrt(gamma_t*Rt*TC_ratio()*gamma_c*Rc) * M9)-M0)
    Fbypass = m0*(alpha*a0/(1+alpha))*((1+f_B())*(   np.sqrt(gamma_t*Rt*TB_ratio()*gamma_c*Rc) * M19)-M0) # P0=P9, so last term all goes to 0.
    return Fcore+Fbypass

def SFC():
    """ Calculates specifc fuel consumption """
    return ((m0/(1 + alpha)) * (f() + f_B()*alpha + f_AB())) / F()
