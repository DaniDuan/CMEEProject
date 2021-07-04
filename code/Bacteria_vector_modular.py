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

N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 0 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin, 10 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_CUE = 0.3
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly number, i.e. how many times the system can assemble
t_fin = 3000 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant


##### Intergrate system forward #####

def ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K):
    '''
    Main function for the simulation of resource uptake and growth of microbial communities.
    '''

    # Setted Parameters
    k = 0.0000862 # Boltzman constant
    a = 15 # The alpha value for beta distribution in Ea
    B_U = 2
    B_R = 1
    # B_U = 0.48 # mean = 0.4810934, N = 377
    # B_R = 0.18 # mean_CUE0 = 0.2203849, SE = 0.04119449, N = 26


    ### Creating empty array for storing data ###
    result_array = np.empty((0,N+M)) # Array to store data in for plotting
    rich_series = np.empty((0))
    U_out_total = np.empty((0,M))
    # U_ac_total = np.empty((0,N))
    R_out = np.empty((0, N))
    CUE_out = np.empty((0,N))
    Ea_CUE_out = np.empty((0,N))
    overlap = np.empty((0,N))
    crossf = np.empty((0,N))


    for i in range(ass):

        ### Resetting values for every assembly ###

        x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (1)))) # Starting concentration for resources and consumers

        # Set up Ea (activation energy) and B0 (normalisation constant) based on Tom Smith's observations
        T_pk_U = 273.15 + np.random.normal(35, 5, size = N)
        T_pk_R = T_pk_U + 3
        # B_U = np.random.exponential(0.48, N) # mean = 0.4810934, N = 377
        # CUE_0 = np.random.uniform(0.18, 0.26, N) # mean = 0.2203849, SE = 0.04119449, N = 26
        # B_R = B_U * (1 - lf) - CUE_0 * B_U
        Ea_U = np.random.beta(a, ((a - 1/3) / (0.82/4)) + 2/3 - a, N)*4
        Ea_R = np.random.beta(a, ((a - 1/3) / (0.67/4)) + 2/3 - a, N)*4
        # Ea_R = Ea_U - Ea_CUE * (B_U * (1 - lf) - B_R) / B_R
        
        p = np.repeat(p_value, M)  # Resource input

        # Set up model
        U = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[0] # Uptake
        R = par.params(N, M, T, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[1] # Respiration
        l = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[2] # Leakage
        l_sum = np.sum(l, axis=1)

        # Integration
        t = np.linspace(0,t_fin-1,t_fin) # resetting 't' if steady state not reached (see below)
        pars = (U, R, l, p, l_sum, N, M, typ, K) # Parameters to pass onto model

        pops = odeint(mod.metabolic_model, y0=x0, t=t, args = pars) # Integrate
        pops = np.round(pops, 7)


        ### Storing simulation results ###
        result_array = np.append(result_array, pops, axis=0)
        U_out_total = np.append(U_out_total, U, axis = 0)
        R_out = np.append(R_out, [R], axis = 0)


        ### Analysis ###

        # Competition for resources
        jaccard = np.array([[np.sum(np.minimum(U[i,],U[j,]))/np.sum(np.maximum(U[i,],U[j,])) for j in range(N)] for i in range(N)])
        comp = np.mean(jaccard, axis = 0)
        overlap = np.append(overlap, [comp], axis = 0)
        # U_ac_total = np.append(U_ac_total, [comp*np.sum(U,axis=1)], axis = 0)


        # Cross-feeding
        leak = U@l
        cf = np.array([[np.sum(np.minimum(leak[i], U[j]))/np.sum(np.maximum(leak[i], U[j])) for j in range(N)] for i in range(N)])
        np.fill_diagonal(cf, 0)
        crossf = np.append(crossf, [np.mean(cf, axis = 1)], axis = 0)


        # Richness
        rich_series = np.append(rich_series, [len(np.where(pops[t_fin-1,0:N])[0])]) # Get the consumer concentrations from last timestep and find out survivors


        # CUE
        CUE = (U @ (1 - l_sum) - R)/(np.sum(U, axis = 1))
        CUE_out = np.append(CUE_out, [CUE], axis = 0)
        Ea_CUE_out = np.append(Ea_CUE_out, [B_R*(Ea_U - Ea_R)/(B_U*(1 - lf) - B_R)], axis = 0)


        ### Previous code for calculating simulated CUE ###

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

    return result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf
