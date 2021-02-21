#### Functions to work out size and temperature dependency

import numpy as np

# B0 for uptake
def size_growth(B_U, Ma):
    B0_U = B_U * Ma**(-0.25)
    return B0_U

# B0 for maintenance respiration
def size_resp(B_R, Ma):
    B0_R = B_R * Ma**(-0.25)
    return B0_R

# Arrhenius/Sharpe-Schoolfield for maintenance growth
def temp_growth(k, T, Tref, T_pk,N, B_U, Ma, Ea_U, Ea_D):
    Sharpe = (size_growth(B_U, Ma) * np.exp((-Ea_U/k) * ((1/T)-(1/Tref))))#/(1 + (Ea_U/(Ea_D - Ea_U)) * np.exp(Ea_D/k * (1/T_pk - 1/T)))
    return Sharpe

# Arrhenius/Sharpe-Schoolfield for maintenance respiration
def temp_resp(k, T, Tref, T_pk,N, B_R, Ma,  Ea_R, Ea_D):
    Sharpe = (size_resp(B_R, Ma) * np.exp((-Ea_R/k) * ((1/T)-(1/Tref))))#/(1 + (Ea_R/(Ea_D - Ea_R)) * np.exp(Ea_D/k * (1/T_pk - 1/T)))
    return Sharpe
