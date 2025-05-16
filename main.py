import numpy as np
import matplotlib.pyplot as plt


pi_b = 0.99; pi_AB = 0.9; pi_n = 0.95; pi_B = 1; pi_dmax = 0.96; pi_C = 30;
eff_b = 0.99; eff_AB = 0.99; eff_B = eff_AB; eff_m = 0.98
e_t = 0.93; e_f = 0.93; e_c = 0.91
c_pc = 1004; c_pt = 12239; gamma_c = 1.4; gamma_t = 1.3;
c_ratio = c_pt/c_pc

m0 = something; alpha = 1; mC = m0/(1+alpha); mB = m0*alpha/(1+alpha)

#x = False

#tau_0 = x; tau_b = x; tau_lambda = x; tau_lambAB = x; tau_B = x; tau_AB = x; tau_f = x;

#tau_c = pi_C*(gamma-1)/(gamma*e_c); tau_t = pi_T*(e_t*(gamma-1)/gamma);  #which gamma?

M0 = varies;

def TB_ratio_():
    """ Bypass Temperature Ratio T19/T0 """ 
    TB_ratio = tau_fn*tau_B*tau_f*tau_0
    return TB_ratio

def TC_ratio_():
    """ Core Temperature Ratio T9/T0 """
    TC_ratio = tau_lambda*tau_t*tau_AB*(1/c_ratio)
    return TC_ratio

def PB_ratio_():
    """ Pressure stagnation ratio across Bypass Pt19/P19 """
    return 

def PC_ratio_():
    """ Pressure stagnation ratio across Core Pt9/P9 """ 
    return 

def M9_():
    """ Need to define Pt19/P0 """
    return (2/(yt-1))*(())

def M19_():
    """ Need to define Pt9/P0 """
    return 

def f_():
    """ Calculates fuel air ratio for burner """
    f = (((eff_b*H)/(c_pc*T_0))**-1)*(tau_lambda - tau_c*tau_0)
    return f

def f_B_():
    """ Caclulates fuel air ratio for bypass """
    f_B = ((((eff_B*H)/(c_pc*T_0)) - tau_B*tau_f*tau_0)**-1)*(tau_B*tau_f*tau_0 - tau_f*tau_0)
    return f_B

def f_AB_():
    """ Calculates fuel air ratio for after burner """
    f_AB = ((((eff_AB*H)/(c_pt*T_0)) - tau_AB*tau_t*tau_lambda)**-1)*(1+f_())*(tau_AB*tau_t*tau_lambda - tau_t*tau_lambda)
    return f_AB

def F_():
    """ Calculates overall thrust """
    Fcore = m0*(a0/(1+alpha)*((1+f_())*(   np.sqrt(gamma_t*Rt*TC_ratio_()*gamma_c*Rc) * M9 + np.sqrt(gamma_c*Rt*TC_ratio_()*gamma_t*Rc) * (1-PC_ratio_())/(M9*gamma_c)  )-M0))
    Fbypass = m0*(alpha*a0/(1+alpha)*((1+f_B_())*(   np.sqrt(gamma_t*Rt*TB_ratio_()*gamma_c*Rc) * M19 + np.sqrt(gamma_c*Rt*TB_ratio_()*gamma_t*Rc) * (1-P0/P19)/(M19*gamma_c)  )-M0))
    return Fcore+Fbypass

def SFC_():
    """ Calculates specifc fuel consumption """
    return ((m0/(1 + alpha)) * (f_() + f_B_()*alpha + f_AB_())) / F_()



