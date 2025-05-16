import numpy as np



pi_b = 0.99; pi_AB = 0.9; pi_n = 0.95; pi_B = 1; pi_dmax = 0.96; pi_C = 30;
eff_b = 0.99; eff_AB = 0.99; eff_B = eff_AB; eff_m = 0.98
e_t = 0.93; e_f = 0.93; e_c = 0.91
c_pc = 1004; c_pt = 12239; yc = 1.4; yt = 1.3;

m0 = something; alpha = 1; mC = m0/(1+alpha); mB = m0*alpha/(1+alpha)

x = False

tau_0 = x; tau_b = x; tau_lambda = x; tau_lambAB = x; tau_B = x; tau_AB = x; tau_f = x;

tau_c = pi_C*(gamma-1)/(gamma*e_c); tau_t = pi_T*(e_t*(gamma-1)/gamma);  #which gamma?

M0 = varies;



def f_():
    f = (((eff_b*H)/(c_pc*T_0))**-1)*(tau_lambda - tau_c*tau_0)
    return f

def f_B_():
    f_B = ((((eff_B*H)/(c_pc*T_0)) - tau_B*tau_f*tau_0)**-1)*(tau_B*tau_f*tau_0 - tau_f*tau_0)
    return f_B

def f_AB_():
    f_AB = ((((eff_AB*H)/(c_pt*T_0)) - tau_AB*tau_t*tau_lambda)**-1)*(1+f_())*(tau_AB*tau_t*tau_lambda - tau_t*tau_lambda)
    return f_AB

def SFC_():
    return m0/2 * f_() + f_B_() + f_AB_()



def F_():
    return m0*(a0/(1+alpha))*(((2/(gamma-1))*tau_lambAB*(1-1/(tau_0*tau_c*tau_t)))**0.5-M0) \
    + alpha/(1+alpha)* a0 * ((tau_lambAB/(tau_0*tau_f))**0.5 * (2/(gamma-1)*(tau_0*tau_f-1)**0.5 - M0))