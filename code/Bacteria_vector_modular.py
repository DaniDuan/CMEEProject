import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import sys
import importlib
import math
from sklearn.utils import shuffle
from random import randint 
import size_temp_funcs as st
import parameters as par
import model_func as mod
import random

######### Main Code ###########

######## Set up parameters ###########

N = 50 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 25 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
# pk = 20 # Peak above Tref, degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
t_n = 30 # Number of temperatures to run the model at, model starts at 20

# Assembly
ass = 1 # Assembly number, i.e. how many times the system can assemble
tv = 10 # immigration times inside one assembly
t_fin = 100 # Number of time steps
x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (0.1)))) # Starting concentration for resources and consumers
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant


##### Intergrate system forward #####

def ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, x0, Ea_D, typ, K):
    # pars_out = np.empty((t_n-20, 19)

    # Setted Parameters
    k = 0.0000862 # Boltzman constant
    # T_pk = Tref + pk # Peak above Tref, Kelvin

    result_array = np.empty((0,N+M)) # Array to store data in for plotting
    CUE_out = np.empty((0,N))
    rich_seires = np.empty((0,tv))
    t_point = 0
    


    for i in range(ass):

        # Set up Ea (activation energy) and B0 (normalisation constant)
        # Based on Tom Smith's observations
        Ea_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:N] # Ea for uptake
        Ea_R = Ea_U - 0.8 # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
        B_U = (10**(2.84 + (-4.96 * Ea_U))) + 4 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
        B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration

        pk_U = np.random.normal(25, 3, size = N)
        pk_R = pk_U + 2
        T_pk_R = Tref + pk_R 
        T_pk_U = Tref + pk_U

        p = np.concatenate((np.array([1]), np.repeat(1, M-1)))  # Resource input

        # Set up model
        U = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D)[0] # Uptake
        R = par.params(N, M, T, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D)[1] # Respiration
        l = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D)[2] # Leakage
        l_sum = np.sum(l, axis=1)


        rich = np.empty((0, tv))


        for j in range(tv):

            t = np.linspace(0,t_fin-1,t_fin) # resetting 't' if steady state not reached (see below)
            pars = (U, R, l, p, l_sum, N, M, typ, K) # Parameters to pass onto model

            # Run model
            pops = odeint(mod.metabolic_model, y0=x0, t=t, args = pars) # Integrate

            # # Steady state test
            # ss_test = np.round(abs((pops[t_fin-1,0:N]) - (pops[t_fin-50,0:N])),3) # Find difference between last timestep and timestep at t_fin - 50 (i.e. 150 if t_fin = 200)
            # while True:
            #     if  np.any(ss_test > 0):    # If there is a difference then that consumer not yet reached steady state
            #         t = np.linspace(0,99,100) # reset t so shorter, only 100 timesteps
            #         pops2 = odeint(mod.metabolic_model, y0=pops[t_fin-1,:], t=t, args=pars) # re-run model using last timestep concentrations of consumers and resources
            #         pops = np.append(pops, pops2, axis=0) # append results of additional run to orginial run
            #         t_fin = t_fin + 100 # adjust timesteps number
            #         ss_test = np.round(abs((pops[t_fin-1,0:N]) - (pops[t_fin-50,0:N])),3) # Find again if consumers reached steady state now
            #     elif np.all(ss_test == 0):
            #         break # Once no difference in consumer concentration then stop performing additional model runs 
            #     else:
            #         pops=pops # If at steady state then nothing happens
            
            # t_fin = 100
            # t_point = t_point + t_fin
            pops = np.round(pops, 7)


            # G_i = [i for i in range(len(xc)-1) if np.any((xc[i+1,]-xc[i,]) > 0)] # before reaching steady state


            ###Assembly###

            # Find which consumers have gone extinct
            rem_find = pops[t_fin-1,0:N] # Get the consumer concentrations from last timestep
            ext = np.where(rem_find<0.01) # Find which consumers have a concentration < 0.01 g/mL, i.e. are extinct
            rem_find = np.where(rem_find<0.01,0.1,rem_find) # Replace extinct consumers with a new concentration of 0.1

            # Richness
            rich = np.append(rich, N - len(rem_find[ext])) # Richness
            # t_rich = np.append(t_rich, t_point)

            
            x0 = np.concatenate((rem_find, pops[t_fin-1,N:N+M])) # Join new concentrations for consumers with those of resources
            

            # Invasion

            # New Ea_ and Ea_R
            Ea_tmp_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:len(ext[0])] # Ea for uptake cut to length(ext),i.e. the number of extinct consumers
            Ea_U[ext] = Ea_tmp_U # Replace removed Eas with new Eas
            Ea_R = Ea_U - 0.8

            # New B0
            B_U = 10**(2.84 + (-4.96 * Ea_U)) + 4
            B_R = 10**(1.29 + (-1.25 * Ea_R))

            # New Tpeak
            pk_U = np.random.normal(25, 3, size = len(rem_find[ext]))
            pk_R = pk_U + 2
            T_pk_R = Tref + pk_R
            T_pk_U = Tref + pk_U

            # New U&R
            U_new = par.params(len(rem_find[ext]), M, T, k, Tref, T_pk_U, B_U[ext], B_R[ext],Ma, Ea_U[ext], Ea_R[ext], Ea_D[ext])[0]
            U[ext] = U_new
            R_new = par.params(len(rem_find[ext]), M, T, k, Tref, T_pk_U, B_U[ext], B_R[ext],Ma, Ea_U[ext], Ea_R[ext], Ea_D[ext])[1]
            R[ext] = R_new


            result_array = np.append(result_array, pops, axis=0)


            # # CUE
            # xc =  pops[:,0:N] # consumer
            # r =  pops[:,N:N+M] # resources
            # if typ == 2:
            #     xr = r /(K + r) # type 2, monod function
            # else:
            #     xr = r #type 1 
            # SL = (1 - l_sum) * xr
            # C = np.einsum('ij,kj->ik', SL, U) - R
            # dCdt = xc * C
            # CUE = dCdt / (xc*np.einsum('ij,kj->ik', xr, U)) # CUE of single species
            # # CUE = C / np.einsum('ij,kj->ik', xr, U)
            # CUE_out = np.append(CUE_out,np.round(CUE, 5), axis = 0)
            # # CUE_out = np.nan_to_num(CUE, nan=0)

            # # Richness & Community CUE
            # for a in range(len(pops)):
            #     dCdt_com = np.sum(np.nan_to_num(dCdt[a,:]/xc[a,:], nan=0))
            #     U_com = np.sum(xr[a,:]*U[np.where(pops[a,0:N] >= 0.01)]) # Community level uptake
            #     R_com = np.sum(R[np.where(pops[a,0:N] >= 0.01)]) # Community level respiration
            #     CUE_com = dCdt_com/U_com # Community level CUE
            #     CUE_com_out = np.append(CUE_com_out, np.round(CUE_com, 5))

        # x0 = pops[len(pops)-1,:]
        # x0[N:N+M] = x0[N:N+M] + 0.1
        # p = p + 1
        rich_seires = np.append(rich_seires, [rich], axis = 0)

    return result_array, rich_seires, # CUE_out # CUE_com_out

ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, x0, Ea_D, typ, K)
