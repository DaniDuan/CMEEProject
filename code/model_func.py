#### Model equations
import numpy as np

def metabolic_model(x,t, U, R, l, p, l_sum, N, M, typ, K):
    '''
    ODEs for integration.
    '''
    
    if typ == 1:
        Mo = x[N:N+M]
    if typ == 2:
        Mo = x[N:N+M]/(K + x[N:N+M]) # Monod function for type 2

    A = np.empty((N+M))
    A[0:N] = x[0:N] * (np.sum((1 - l_sum) * Mo * U, axis=1) - R)
    # x[0:N] * (((1-l_sum)*U) @ x[N:N+M] - R)
    A[N:N+M] = p - x[0:N] @ (U * Mo * x[N:N+M]) + x[0:N] @ ((U * Mo * x[N:N+M]) @ l)
    
    return A