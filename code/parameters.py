### Define params

import numpy as np
from numpy import random
import size_temp_funcs as st
import scipy as sc

# Parameters
def params(N, M, T, k, Tref, T_pk, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D):
    # Uptake
    # U = np.random.sample([N,M])+300

    # u1 = np.zeros([M,M])
    # u_temp = st.temp_growth(k, T, Tref, T_pk, N, B_U, Ma, Ea_U, Ea_D) # uptake rates and maintenance are temp-dependant
    # np.fill_diagonal(u1,u_temp[0:M]) # fill in temp dependant uptake on diagonals
    # u2 = np.zeros([M,M])
    # u3 = np.zeros([M,M])
    # if N == M:
    #     U = u1
    # elif M == N/2:
    #     np.fill_diagonal(u2,u_temp[M:M*2]) # fill in temp dependant uptake on diagonals
    #     U = np.concatenate((u1, u2), axis=0)
    # else:
    #     np.fill_diagonal(u2,u_temp[M:M*2]) # fill in temp dependant uptake on diagonals
    #     np.fill_diagonal(u3,u_temp[M+1:M*3])
    #     U = np.concatenate((u1, u2, u3), axis=0)
    

    # Give random uptake for each resources, sum up to the total uptake of the bacteria
    U_sum = st.temp_growth(k, T, Tref, T_pk, N, B_U, Ma, Ea_U, Ea_D)
    u_1 = np.empty((0,N)) 
    for i in range(M-1):
        mean = U_sum[i]/N
        random.seed(i)
        a = np.array([np.random.uniform(0, mean, size = N)])
        u_1 = np.append(u_1,a,axis = 0)
    u_2 = np.array([U_sum[0:N] - np.sum(u_1, axis = 0)])
    # b = np.where(u_2<0)
    # u_1[:,b] = u_1[:,b] + u_2[b]/N-1
    # u2 = np.where(u2<0, 0, u2)
    U_raw = np.append(u_1, u_2, axis = 0)
    for i in range(U_raw.shape[1]): 
        random.seed(i)
        np.random.shuffle(U_raw[:,i])
    U = np.transpose(U_raw)


    # Respiration
    ar_rm = st.temp_resp(k, T, Tref,T_pk, N, B_R, Ma, Ea_R, Ea_D) # find how varies with temperature (ar = arrhenius)
    # #Rm = sc.full([N], (0.3))
    R = ar_rm

    
    ### Calc  R as function  of U/Rm rule for competitors. Then need to back-calc Ea.

    # R_2 = np.array([])
    # if M == N:
    #     R = st.temp_resp(k, T, Tref,T_pk, N, B_r, Ma, Ea_R, Ea_D) # find how varies with temperature (ar = arrhenius)
    # elif M == N/2:
    #     R_1 = st.temp_resp(k, T, Tref,T_pk, N, B_r, Ma, Ea_R, Ea_D)[0:M] # find how varies with temperature (ar = arrhenius)
    #     for g in range(0, M):
    #         #print(g)
    #         R_1x = range(0,M) # upper diag row
    #         R_1y = range(0,M) #np.where(np.diag(U)) # upper diag column
    #         R_2x = np.array(np.where(np.diag(U[M:N]))) + M # lower diag row
    #         R_2y = range(0,M) # lower diag column
    #         R_i = U[R_2x[0,g], R_2y[g]] * R_1[g] * Meta_ratio[g]/U[R_1x[g], R_1y[g]]
    #         #print(R_i)
    #         R_2 = np.append(R_2, R_i)
    #         #print(R_2)
    #     R = np.concatenate((R_1, R_2))
 



    # Excretion
    l = np.zeros([M,M])
    for i in range(M-1): 
        l[i,i+1] =  0.4
    l[M-1,0] = 0.4
    

    # External resource input
    p = np.concatenate((np.array([1]), np.repeat(1, M-1)))  #np.repeat(1, M) #np.ones(M) #
    return U, R, l, p
