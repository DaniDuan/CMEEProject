#### Model equations
import numpy as np

def metabolic_model(x,t, U, R, l, p, l_sum, N, M, K, typ):
    A = np.empty((N+M))
    A[0:N] = x[0:N] * (np.sum((1 - l_sum) * x[N:N+M] * U, axis=1) - R)
    # A[N:N+M] = p - np.multiply((x[0:N] @ U).transpose(), x[N:N+M]) + np.einsum('i,k,ik,kj->j', x[0:N], x[N:N+M], U, l, optimize = True)
    A[N:N+M] = p - x[0:N] @ (U * x[N:N+M]) + x[0:N] @ ((U * x[N:N+M]) @ l)
    
    # xc =  x[0:N] # consumer
    # xr =  x[N:N+M] # resources

    # Functional response
    # if typ == 2:
    #     xr = r /(K + r) # type 2, monod function
    # else:
        # xr = r #type 1 

    # ## Consumers
    # # calculate leakeage
    # SL = (1 - l_sum) * xr
    # #uptake rate and maintenance
    # C = np.sum(SL * U, axis=1) - R
    # #dCdt
    # dCdt = xc * C
    
    # dCdt = x[0:N] * (np.sum((1 - l_sum) * x[N:N+M] * U, axis=1) - R)

    ## Resources
    # dSdt = p - np.multiply((xc @ U).transpose(), xr) + np.einsum('i,k,ik,kj->j', xc, xr, U, l)

    return A #np.array(np.concatenate((dCdt, dSdt)))