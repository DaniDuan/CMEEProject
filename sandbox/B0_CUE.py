### getting Ea0_CUE
import numpy as np
import parameters as par

M = N = 1000
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_diff = 0.6
lf = 0.4 # Leakage
k = 0.0000862 # Boltzman constant

Ea_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:N] # Ea for uptake
Ea_R = Ea_U - Ea_diff # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
B_U = (10**(2.84 + (-4.96 * Ea_U))) + 30 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration
T_pk_U = Tref + np.random.normal(32, 3, size = N)
T_pk_R = T_pk_U + 3
U0 = par.params(N, M, 273.15, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[0]
R0 = par.params(N, M, 273.15, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[1]
l = par.params(N, M, 273.15, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[2]
l_sum = np.sum(l, axis=1)
B0_CUE = np.mean((U0 @ (1 - l_sum)*0.1 - R0)/(np.sum(U0, axis = 1)*0.1))
np.round(B0_CUE, 3) ## = 0.016

