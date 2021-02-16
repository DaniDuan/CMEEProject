# x = x0

# x = result_array[5]
# xc =  x[0:N] # consumer
# r =  x[N:N+M] # resources
# xr = r/(K + r) # type 2
# SL = (1 - l_sum) * xr
# C = np.sum(SL * U, axis=1) - R
# dCdt = xc * C
# CUE = dCdt / np.sum(xr * U, axis=1)


x = result_array
xc =  x[:,0:N] # consumer
r =  x[:,N:N+M] # resources

# Functional response
xr = r/(K + r) # type 2

SL = (1 - l_sum) * xr
#uptake rate and maintenance
C = np.einsum('ij,kj->ik', SL, U) - R
#dCdt
dCdt = xc * C

CUE = dCdt / np.einsum('ij,kj->ik', xr, U)


t_plot = sc.linspace(0,len(result_array),len(result_array))
plt.plot(t_plot, CUE, 'r-', label = 'Consumers', linewidth=0.7)
plt.show()

