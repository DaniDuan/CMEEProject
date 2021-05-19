### getting B0_CUE
import numpy as np
import parameters as par

M = N = 1000
Tref = 283.15 # Reference temperature Kelvin, 0 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_CUE = 0.3
lf = 0.4 # Leakage
k = 0.0000862 # Boltzman constant7

B_U = np.repeat(2, N)
B_R = B_U * (1 - lf) - 0.2 * B_U # B0 for respiration
Ea_U = np.random.beta(25, ((25 - 1/3) / 0.8) + 2/3 - 25, N)
Ea_R = Ea_U - Ea_CUE * (B_U * (1 - lf) - B_R) / B_R
# B_U = (10**(2.84 + (-4.96 * Ea_U))) + 3 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
# B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration
T_pk_U = Tref + np.random.normal(35, 3, size = N)
T_pk_R = T_pk_U + 3

U0 = par.params(N, M, Tref, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[0]
R0 = par.params(N, M, Tref, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[1]
l = par.params(N, M, Tref, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[2]
l_sum = np.sum(l, axis=1)
B0_CUE = np.mean((U0 @ (1 - l_sum) - R0)/(np.sum(U0, axis = 1)))
np.round(B0_CUE, 3) ## = 0.016

### Testing
T = 273.15
U = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[0] # Uptake
R = par.params(N, M, T, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[1] # Respiration
CUE = (U @ (1 - l_sum) - R)/(np.sum(U, axis = 1))
np.round(np.mean(CUE), 3)


B0_CUE = 0.1 * np.exp((-Ea_CUE/k) * ((1/Tref)-(1/273.15)))
