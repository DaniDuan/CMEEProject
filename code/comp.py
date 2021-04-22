from Bacteria_vector_modular import ass_temp_run
import numpy as np
import matplotlib.pylab as plt

######## Set up parameters ###########

N = 25 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 35 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
# pk = 20 # Peak above Tref, degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance

# Assembly
tv = 20 # immigration times inside one assembly
t_fin = 50 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant


def plot_comp(t_fin, N, M, T, Tref, Ma, tv, Ea_D, typ, K):

    ass = 1
    result_array, rich_seires, l, U_out_total, U_ac_total = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, typ, K)
    
    # uptake = np.array([result_array[i*t_fin:(i+1)*t_fin,0:N] @ U_out_total[i*N:(i+1)*N,:] * result_array[i*t_fin:(i+1)*t_fin,N:N+M] for i in range(tv)])
    # uptake_end = np.array([uptake[i][t_fin-1] for i in range(tv)])
    
    sur = [np.where(result_array[(i+1)*t_fin-1, 0:N] > 0.01) for i in range(tv)]# sur 

    U_range = [np.arange(0, np.ceil(np.max(U_ac_total)/100)*100, 100) for i in range(tv)]
    sur_rate = np.empty((0,int(np.ceil(np.max(U_ac_total)/100))))
    s_total = np.empty((0))
    s_sur = np.empty((0))
    for j in range(tv):
        for i in range(len(U_range[j])):
            s_total = np.append(s_total, ((U_ac_total[j,:] >= U_range[j][i]) & (U_ac_total[j,:] < U_range[j][i]+100)).sum())
            s_sur = np.append(s_sur, ((U_ac_total[j,:][sur[j][0]] >= U_range[j][i]) & (U_ac_total[j,:][sur[j][0]] < U_range[j][i]+100)).sum())
    
    sur_rate = np.array(s_sur/s_total).reshape(tv,int(np.ceil(np.max(U_ac_total)/100)))
    sur_rate_filtered = [i[j] for i, j in zip(sur_rate.T, (~np.isnan(sur_rate)).T)]

    plt.boxplot(sur_rate_filtered)
    plt.title(T-273.15)
    plt.show()

    return sur_rate
