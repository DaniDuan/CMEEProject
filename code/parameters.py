### Define params

import numpy as np
import size_temp_funcs as st

# Parameters
def params(N, M, T, k, Tref, T_pk, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf):
    '''
    Returning matrices and vectors of uptake, respiration and excretion.
    '''
    # Uptake
    # Give random uptake for each resources, sum up to the total uptake of the bacteria
    U_sum = st.temp_growth(k, T, Tref, T_pk, N, B_U, Ma, Ea_U, Ea_D)
    diri = np.transpose(np.random.dirichlet(np.ones(M),N))
    U = np.transpose(diri*U_sum)

    # Respiration
    R = st.temp_resp(k, T, Tref,T_pk, N, B_R, Ma, Ea_R, Ea_D) # find how varies with temperature (ar = arrhenius)

    # Excretion
    # SUMMING UP TO 0.4  
    # l_raw = np.array([[np.random.normal(1/(i-1),0.005)* lf if i-1>0 else np.random.normal(lf, 0.005) for i in range(M,0,-1)] for i in range(1,M+1)])
    # fix = [[1 if j>i else 0 for j in range(M)] for i in range(M)]
    # fix[M-1][0] = 1
    # l = np.transpose(l_raw) * fix

    l_raw = np.array([[np.random.normal(1/(i),0.005)* lf for i in range(M,0,-1)] for i in range(1,M+1)])
    l = [[l_raw[j,i] if j>=i else 0 for j in range(M)] for i in range(M)]

    # l = np.zeros([M,M])
    # for i in range(M-1): 
    #     l[i,i+1] = lf
    # # l[M-1,0] = lf

    return U, R, l
