import Bacteria_vector_modular as bvm
import parameters as par
import numpy as np 
import scipy as sc
import matplotlib.pylab as plt

############# Testing ###########
# x = x0

# x = result_array[50]
# xc =  x[0:N] # consumer
# r =  x[N:N+M] # resources
# # xr = r
# xr = r/(K + r) # type 2
# SL = (1 - l_sum) * xr
# C = np.sum(SL * U, axis=1) - R
# dCdt = xc * C
# CUE = dCdt / np.sum(xr * U, axis=1)


########## Giving Parameters ###########
N = 10 # Number of consumers
M = 5 # Number of resources

# Temperature params
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
pk = 20 # Peak above Tref, degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
t_n = 21 # Number of temperatures to run the model at, model starts at 20

# Assembly
ass = 1 # Assembly number, i.e. how many times the system can assemble
t_fin = 100 # Number of time steps
x0 = np.concatenate((sc.full([N], (0.1)),sc.full([M], (0.1)))) # Starting concentration for resources and consumers
typ = 2 # Functional response, Type I or II
K = 0.5 # Half saturation constant


def CUE_cal(t_fin, N, M, t_n,  Tref, Ma, ass, x0, pk, Ea_D, typ):

    # Setted Parameters
    k = 0.0000862 # Boltzman constant
    T_pk = Tref + pk # Peak above Tref, Kelvin
    Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
    t = sc.linspace(0,t_fin-1,t_fin) # Time steps

    Ea_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:N] # ?Ea for uptake
    Ea_R = Ea_U - 0.8 # ?Ea for respiration, which should always be lower than Ea_U so 'peaks' later
    B_U = (10**(2.84 + (-4.96 * Ea_U))) + 4 # ?B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
    B_R = (10**(1.29 + (-1.25 * Ea_R))) # ?B0 for respiration

    l_sum = np.full((1,M),0.4) # Total leakage for each resource

    x = bvm.ass_temp_run(t_fin, N, M, t_n,  Tref, Ma, ass, x0, pk, Ea_D, typ)[0] # Result array
    xc =  x[:,0:N] # consumer
    r =  x[:,N:N+M] # resources

    # Functional response
    if typ == 2:
        xr = r /(K + r) # type 2, monod function
    else:
        xr = r #type 1 


    #uptake rate and maintenance
    CUE_out = np.empty((0,N))

    for i in range(20, t_n):
        T = 273.15 + i
        U = par.params(N, M, T, k, Tref, T_pk, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D)[0] # Uptake
        R = par.params(N, M, T, k, Tref, T_pk, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D)[1] # Respiration
        SL = (1 - l_sum) * xr[t_fin*(i-20)*ass:t_fin*(i-19)*ass,:]
        C = np.einsum('ij,kj->ik', SL, U) - R
        dCdt = xc[t_fin*(i-20)*ass:t_fin*(i-19)*ass,:] * C
        
        CUE = dCdt / np.einsum('ij,kj->ik', xr[t_fin*(i-20)*ass:t_fin*(i-19)*ass,:], U)
        CUE_out = np.append(CUE_out,CUE, axis = 0)

    
    return CUE_out


def plot_CUE(CUE_out):
    t_plot = sc.linspace(0,len(CUE_out),len(CUE_out))
    plt.plot(t_plot, CUE_out, 'r-', label = 'CUE', linewidth=0.7)
    plt.ylabel('CUE')
    plt.xlabel('Time')
    plt.title('Carbon Use Efficiency dynamics')
    plt.show()
    return


CUE_out = CUE_cal(t_fin, N, M, t_n,  Tref, Ma, ass, x0, pk, Ea_D, typ)
plot_CUE(CUE_out)