import numpy as np
import matplotlib.pyplot as plt

""" Flight Conditions """
M0 = 1.5
T0 = 217  # [K]
# P0 = P9 (assume nozzle is perfectly expanded)

""" Engine Limits """
T_t4 = 1900  # [K]

""" Gas Properties """
H = 42.8e6  # [J/kg], converted from MJ/kg
c_pc = 1004  # [J/kg·K]
gamma_c = 1.4
c_pt = 1239  # [J/kg·K]
gamma_t = 1.3
R_c = c_pc * (1 - 1 / gamma_c)
R_t = c_pt * (1 - 1 / gamma_t)

print(f"R_c = {R_c:.2f} J/kg·K, R_t = {R_t:.2f} J/kg·K")
f_st = 0.0291
""" Technology Level 1 """
pi_d = 0.85
e_c = 0.8
pi_b = 0.9
eff_b = 0.85
e_t = 0.8
pi_n = 0.9
eff_m = 0.98

# Define compression ratio range
pi_c = np.linspace(0.01, 40, 100)  # Start from 1.01 to avoid division by zero
a_0 = np.sqrt(gamma_c * R_c * T0)

def tau_0(M0):
    return 1 + ((gamma_c - 1) / 2) * M0**2

def pi_0(M0):
    return tau_0(M0)**(gamma_c / (gamma_c - 1))

def tau_lambda():
    return (c_pt/c_pc)*(T_t4 / T0)

def tau_c(pi_c_val):
    return pi_c_val**(e_t*((gamma_c - 1) / (gamma_c * e_c)))

def f(M0, pi_c_val):
    return ((tau_lambda() - tau_0(M0) * tau_c(pi_c_val))) / ((eff_b * H) / (c_pc * T0) - tau_lambda())

def tau_t(M0, pi_c_val):
    return 1 - 1 / (eff_m * (1 + f(M0, pi_c_val))) * (tau_0(M0) / tau_lambda()) * (tau_c(pi_c_val) - 1)

def pi_t(M0, pi_c_val):
    return tau_t(M0, pi_c_val)**(gamma_t / ((gamma_t - 1) * e_t))

def T9_T0(M0, pi_c_val):
    #total_pi = (pi_0(M0) * pi_d * pi_c_val * pi_b * pi_t(M0, pi_c_val) * pi_n)**((gamma_t-1)/gamma_t)
    #return tau_lambda() * tau_t(M0, pi_c_val) * (c_pc / c_pt) * (1 / total_pi)
    print(f"TC: {tau_lambda()*tau_t(M0, pi_c_val)*(c_pc / c_pt) / ((1*pi_0(M0)*pi_d*pi_c_val*pi_b*pi_t(M0, pi_c_val)*pi_n)**((gamma_t-1)/gamma_t))}")
    return tau_lambda()*tau_t(M0, pi_c_val)*(c_pc / c_pt) / ((1*pi_0(M0)*pi_d*pi_c_val*pi_b*pi_t(M0, pi_c_val)*pi_n)**((gamma_t-1)/gamma_t))


M9s = []
def M9(M0, pi_c_val):
    total_pi = pi_0(M0) * pi_d * pi_c_val * pi_b * pi_t(M0, pi_c_val) * pi_n
    M_9 = np.sqrt((2 / (gamma_t - 1)) * (total_pi**((gamma_t - 1) / gamma_t) - 1))
    M9s.append(M_9)
    return M_9
TCratios = []
def ST(M0, pi_c_val):
    TCratios.append(T9_T0(M0, pi_c_val))
    return a_0 * ((1 + f(M0,pi_c_val)) * np.sqrt((gamma_t * R_t / (gamma_c * R_c)) * T9_T0(M0, pi_c_val)) * M9(M0, pi_c_val) - M0)

# Compute ST for all pi_c
ST_vals = ST(M0, pi_c)
print(ST_vals)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(pi_c, ST_vals, label='Specific Thrust')
plt.xlabel('Compressor Pressure Ratio (π_c)')
plt.ylabel('Specific Thrust [m/s]')
plt.title('Specific Thrust vs Compressor Pressure Ratio (Tech Level 1)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

print(len(M9s))
print(len(pi_c))

plt.figure(figsize=(8, 5))
plt.plot(pi_c, M9s[0], label='Specific Thrust')
plt.xlabel('Compressor Pressure Ratio (π_c)')
plt.ylabel('M9')
plt.title('Specific Thrust vs Compressor Pressure Ratio (Tech Level 1)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()