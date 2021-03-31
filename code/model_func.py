#### Model equations
import numpy as np

def metabolic_model(x,t, U, R, l, p, l_sum, N, M, K, typ):
    
    xc =  x[0:N] # consumer
    r =  x[N:N+M] # resources

    # Functional response
    if typ == 2:
        xr = r /(K + r) # type 2, monod function
    else:
        xr = r #type 1 

    ## Consumers
    # calculate leakeage
    SL = (1 - l_sum) * xr
    #uptake rate and maintenance
    C = np.sum(SL * U, axis=1) - R
    #dCdt
    dCdt = xc * C
    
    ## Resources
    dSdt = p - np.multiply((xc @ U).transpose(), xr) + np.einsum('i,k,ik,kj->j', xc, xr, U, l)

    return np.array(np.concatenate((dCdt, dSdt)))