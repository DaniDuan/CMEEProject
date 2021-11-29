### Define params

import numpy as np
import size_temp_funcs as st

# Parameters
def params(N, M, T, k, Tref, T_pk_U, T_pk_R, B_U, B_R, Ma, Ea_U, Ea_R, Ea_D, lf):
    '''
    Returning matrices and vectors of uptake, respiration and excretion.
    '''
    # Uptake
    # Give random uptake for each resources, sum up to the total uptake of the bacteria
    U_sum = st.temp_growth(k, T, Tref, T_pk_U, N, B_U, Ma, Ea_U, Ea_D)
    diri = np.transpose(np.random.dirichlet(np.full(M,1),N))
    U = np.transpose(diri*U_sum)

    # jaccard = np.array([[np.sum(np.minimum(U[i,],U[j,]))/np.sum(np.maximum(U[i,],U[j,])) for j in range(N)] for i in range(N)])
    # np.mean(np.mean(jaccard, axis = 0))

    # Respiration
    R = st.temp_resp(k, T, Tref,T_pk_R, N, B_R, Ma, Ea_R, Ea_D) # find how varies with temperature (ar = arrhenius)

    # Excretion
    # SUMMING UP TO 0.4 
    # l_raw = np.array([[np.random.normal(1/(i),0.005)* lf for i in range(M,0,-1)] for i in range(1,M+1)])
    # l = [[l_raw[j,i] if j>=i else 0 for j in range(M)] for i in range(M)]

    # Allowing maximum of 2 metabolic products
    l_raw = np.array([[np.random.normal(lf/3,0.005) if i >= 3 and i > M-3 else np.random.normal(1/(i),0.005)* lf for i in range(M,0,-1)] for i in range(1,M+1)])
    l = [[l_raw[j,i] if j>=i and j-i <3 else 0 for j in range(M)] for i in range(M)]
    # np.round(l,3) 

    return U, R, l


