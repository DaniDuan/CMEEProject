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

N = 25 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 30 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
# pk = 20 # Peak above Tref, degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
t_n = 25 # Number of temperatures to run the model at, model starts at 20

# Assembly
ass = 1 # Assembly number, i.e. how many times the system can assemble
tv = 10 # immigration times inside one assembly
t_fin = 50 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant


##### Intergrate system forward #####

def ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, typ, K):
    # pars_out = np.empty((t_n-20, 19)

    # Setted Parameters
    k = 0.0000862 # Boltzman constant
    # T_pk = Tref + pk # Peak above Tref, Kelvin

    result_array = np.empty((0,N+M)) # Array to store data in for plotting
    CUE_out = np.empty((0,N))
    rich_seires = np.empty((0,tv))
    # U_out = np.empty((0,M))
    U_out_total = np.empty((0,M))
    sur_rate = np.empty((0,4))


    for i in range(ass):

        x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (0.1)))) # Starting concentration for resources and consumers

        # Set up Ea (activation energy) and B0 (normalisation constant)
        # Based on Tom Smith's observations
        Ea_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:N] # Ea for uptake
        Ea_R = Ea_U - 0.6 # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
        B_U = (10**(2.84 + (-4.96 * Ea_U))) + 4 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
        B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration

        T_pk_U = Tref + np.random.normal(32, 5, size = N)
        T_pk_R = T_pk_U + 3

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
            pops = np.round(pops, 7)


            # G_i = [i for i in range(len(xc)-1) if np.any((xc[i+1,]-xc[i,]) > 0)] # before reaching steady state


            ###Assembly###

            # Find which consumers have gone extinct
            rem_find = pops[t_fin-1,0:N] # Get the consumer concentrations from last timestep
            ext = np.where(rem_find<0.01) # Find which consumers have a concentration < 0.01 g/mL, i.e. are extinct
            sur = np.where(rem_find>0.01)
            # U_out = np.append(U_out, U[sur[0]], axis = 0)
            U_out_total = np.append(U_out_total, U, axis = 0)


            jaccard = np.zeros((N,N)) # Competition
            np.fill_diagonal(jaccard,1)
            for i in range(N):
                for j in range(N):
                   jaccard[i,j] = np.sum(np.minimum(U[i,],U[j,]))/np.sum(np.maximum(U[i,],U[j,])) 
            comp = np.mean(jaccard, axis = 0)
            U_ac = comp*np.sum(U,axis=1)

            U_range = np.arange(0,400,100)
            # U_range = np.arange(np.floor(np.min(U_ac)/100)*100, np.ceil(np.max(U_ac)/100)*100, 100)
            s_total = np.empty((0))
            s_sur = np.empty((0))
            for i in range(len(U_range)):
                 s_total = np.append(s_total, ((U_ac >= U_range[i]) & (U_ac < U_range[i]+100)).sum())
                 s_sur = np.append(s_sur, ((U_ac[sur[0]] >= U_range[i]) & (U_ac[sur[0]] < U_range[i]+100)).sum())
            sur_rate = np.append(sur_rate, [np.array(np.nan_to_num(s_sur/s_total, nan=0))], axis = 0)


            # Richness
            rich = np.append(rich, N - len(rem_find[ext])) # Richness
            # t_rich = np.append(t_rich, t_point)

            rem_find = np.where(rem_find<0.01,0.1,rem_find) # Replace extinct consumers with a new concentration of 0.1
            x0 = np.concatenate((rem_find, pops[t_fin-1,N:N+M])) # Join new concentrations for consumers with those of resources
            

            # Invasion

            # New Ea_ and Ea_R
            Ea_tmp_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:len(ext[0])] # Ea for uptake cut to length(ext),i.e. the number of extinct consumers
            Ea_U[ext] = Ea_tmp_U # Replace removed Eas with new Eas
            Ea_R = Ea_U - 0.6

            # New B0
            B_U = 10**(2.84 + (-4.96 * Ea_U)) + 4
            B_R = 10**(1.29 + (-1.25 * Ea_R))

            # New Tpeak
            pk_U = np.random.normal(32, 5, size = len(rem_find[ext]))
            pk_R = pk_U + 3
            T_pk_R = Tref + pk_R
            T_pk_U = Tref + pk_U

            # New U&R
            U_new = par.params(len(rem_find[ext]), M, T, k, Tref, T_pk_U, B_U[ext], B_R[ext],Ma, Ea_U[ext], Ea_R[ext], Ea_D[ext])[0]
            U[ext] = U_new
            R_new = par.params(len(rem_find[ext]), M, T, k, Tref, T_pk_U, B_U[ext], B_R[ext],Ma, Ea_U[ext], Ea_R[ext], Ea_D[ext])[1]
            R[ext] = R_new
            
            result_array = np.append(result_array, pops, axis=0)

            # # CUE
            # dCdt = pops[:,0:N] * ((1 - l_sum) * pops[:,N:N+M] @ U.transpose() - R)
            # CUE = dCdt / (pops[:,0:N] * (pops[:,N:N+M] @ U.transpose())) # CUE of single species
            # # CUE = C / np.einsum('ij,kj->ik', xr, U)
            # CUE_out = np.append(CUE_out,np.round(CUE, 5), axis = 0)


            # # Community CUE
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

    return result_array, rich_seires, l, U_out_total, sur_rate # , CUE_out

B = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, typ, K)
