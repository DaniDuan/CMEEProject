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
from scipy.optimize import fsolve

######### Main Code ###########

######## Set up parameters ###########

N = 15 # Number of consumers
M = 15 # Number of resources

# Temperature params
T = 273.15 + 35 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_diff = 0.2
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 2 # Assembly number, i.e. how many times the system can assemble
tv = 100 # immigration times inside one assembly
t_fin = 50 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant


##### Function for calculating CUE #####
def CUE_calc(U, R, T_pk_U, l_sum): 
    k = 0.0000862 # Boltzman constant
    B0_CUE = 0.175 
    T_pk_CUE = T_pk_U #CUE paper
    Ea_D_CUE = np.repeat(3.5,N)
    CUE = (U @ (1 - l_sum)*0.1 - R)/(np.sum(U, axis = 1)*0.1)
    EaCUE_values = np.empty(0)
    for i in range(len(CUE)):
        Ea_CUE = 0
        def f(Ea_CUE):
            return(CUE[i] - B0_CUE * 2.71828**((-Ea_CUE/k) * ((1/T)-(1/Tref)))/(1 + (Ea_CUE/(Ea_D[i] - Ea_CUE)) * 2.71828**(Ea_D[i]/k * (1/T_pk_CUE[i] - 1/T))))
        answer = fsolve(f, 0.5)
        EaCUE_values = np.append(EaCUE_values, answer)

    return CUE, EaCUE_values


##### Intergrate system forward #####

def ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, Ea_diff, lf, p_value, typ, K):
    '''
    Main function for the simulation of resource uptake and growth of microbial communities.
    '''

    # Setted Parameters
    k = 0.0000862 # Boltzman constant

    ### Creating empty array for storing data ###
    result_array = np.empty((0,N+M)) # Array to store data in for plotting
    rich_series = np.empty((0,tv))
    U_out_total = np.empty((0,M))
    U_ac_total = np.empty((0,N))
    R_out = np.empty((0))
    CUE_out = np.empty((0,N))
    Ea_CUE_out = np.empty((0,N))


    for i in range(ass):

        ### Resetting values for every assembly ###

        x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (0.1)))) # Starting concentration for resources and consumers

        # Set up Ea (activation energy) and B0 (normalisation constant) based on Tom Smith's observations
        Ea_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:N] # Ea for uptake
        Ea_R = Ea_U - Ea_diff # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
        # B_U = (10**(2.84 + (-4.96 * Ea_U))) + 4
        # B_R = (10**(1.29 + (-1.25 * Ea_R))) + 1.2
        B_U = (10**(2.84 + (-4.96 * Ea_U))) + 30 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
        B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration
        T_pk_U = Tref + np.random.normal(32, 5, size = N)
        T_pk_R = T_pk_U + 3

        p = np.repeat(p_value, M)  # Resource input
        # p = np.arange(M,0,-1)
        # p = np.concatenate((np.array([1,1,1]), np.repeat(0, M-3)))

        # Set up model
        U = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[0] # Uptake
        R = par.params(N, M, T, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[1] # Respiration
        l = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[2] # Leakage
        l_sum = np.sum(l, axis=1)

        rich = np.empty((0, tv))


        for j in range(tv):

            ### Integration ###
            
            t = np.linspace(0,t_fin-1,t_fin) # resetting 't' if steady state not reached (see below)
            pars = (U, R, l, p, l_sum, N, M, typ, K) # Parameters to pass onto model

            pops = odeint(mod.metabolic_model, y0=x0, t=t, args = pars) # Integrate
            pops = np.round(pops, 7)

            # G_i = [i for i in range(len(xc)-1) if np.any((xc[i+1,]-xc[i,]) > 0)] # before reaching steady state

            ### Analysis ###

            # Find which consumers have gone extinct
            rem_find = pops[t_fin-1,0:N] # Get the consumer concentrations from last timestep
            ext = np.where(rem_find<0.01) # Find which consumers have a concentration < 0.01 g/mL, i.e. are extinct
            sur = np.where(rem_find>0.01)
            U_out_total = np.append(U_out_total, U, axis = 0)
            R_out = np.append(R_out, np.mean(R))
            
            # Competition for resources
            jaccard = np.zeros((N,N)) # Competition
            np.fill_diagonal(jaccard,1)
            jaccard = np.array([[np.sum(np.minimum(U[i,],U[j,]))/np.sum(np.maximum(U[i,],U[j,])) for j in range(N)] for i in range(N)])
            comp = np.mean(jaccard, axis = 0)
            U_ac_total = np.append(U_ac_total, [comp*np.sum(U,axis=1)], axis = 0)

            # Richness
            rich = np.append(rich, N - len(rem_find[ext])) # Richness

            # CUE and Ea_CUE
            CUE_result, Ea_CUE_result =  CUE_calc(U, R, T_pk_U, l_sum)
            CUE_out = np.append(CUE_out, [CUE_result], axis = 0)
            Ea_CUE_out = np.append(Ea_CUE_out, [Ea_CUE_result], axis = 0)

            ### Invasion ###

            rem_find = np.where(rem_find<0.01,0.1,rem_find) # Replace extinct consumers with a new concentration of 0.1
            x0 = np.concatenate((rem_find, pops[t_fin-1,N:N+M])) # Join new concentrations for consumers with those of resources
            
            # New Ea_ and Ea_R
            Ea_tmp_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:len(ext[0])] # Ea for uptake cut to length(ext),i.e. the number of extinct consumers
            Ea_U[ext] = Ea_tmp_U # Replace removed Eas with new Eas
            Ea_R = Ea_U - Ea_diff

            # New B0
            B_U = 10**(2.84 + (-4.96 * Ea_U)) + 30
            B_R = 10**(1.29 + (-1.25 * Ea_R)) # + 1.2

            # New Tpeak
            pk_U = np.random.normal(32, 5, size = len(rem_find[ext]))
            pk_R = pk_U + 3
            T_pk_R[ext] = Tref + pk_R
            T_pk_U[ext] = Tref + pk_U

            # New U&R
            U_new = par.params(len(rem_find[ext]), M, T, k, Tref, T_pk_R[ext], B_U[ext], B_R[ext],Ma, Ea_U[ext], Ea_R[ext], Ea_D[ext], lf)[0]
            U[ext] = U_new
            R_new = par.params(len(rem_find[ext]), M, T, k, Tref, T_pk_U[ext], B_U[ext], B_R[ext],Ma, Ea_U[ext], Ea_R[ext], Ea_D[ext], lf)[1]
            R[ext] = R_new
            
            ### Storing simulation results ###
            result_array = np.append(result_array, pops, axis=0)

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

        # x0 = pops[len(pops)-1,:]
        # x0[N:N+M] = x0[N:N+M] + 0.1
        # p = p + 1
        
        rich_series = np.append(rich_series, [rich], axis = 0)

    return result_array, rich_series, l, U_out_total, U_ac_total, R_out, CUE_out, Ea_CUE_out


# result_array, rich_series, l, U_out_total, U_ac_total, R_out, CUE_out, Ea_CUE_out = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, Ea_diff, lf, p_value, typ, K)
# U_out = np.array([np.mean(U_out_total[N*i:N*(i+1)-1,:]) for i in range(tv*ass)])
# fig, ax1 = plt.subplots()
# ax2 = ax1.twiny()
# ax1.scatter(U_out, np.ndarray.flatten(rich_series), c = 'r')
# ax2.scatter(R_out, np.ndarray.flatten(rich_series), c = "b")
# plt.show()
