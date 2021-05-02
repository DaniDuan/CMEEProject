from Bacteria_vector_modular import ass_temp_run
import numpy as np
import matplotlib.pylab as plt

######## Set up parameters ###########

N = 25 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 25 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
# pk = 20 # Peak above Tref, degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_diff = 0.6
lf = 0.4

# Assembly
tv = 20 # immigration times inside one assembly
t_fin = 50 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant


def plot_comp(t_fin, N, M, T, Tref, Ma, tv, Ea_D, lf, typ, K):
    '''
    Calculating the actual species resource uptake accounting competition, 
    and returning the survival rate of species with the actual uptake values at different temperatures.
    '''

    ass = 1
    result_array, rich_seires, l, U_out_total, U_ac_total = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, Ea_diff, lf, typ, K)
    
    # uptake = np.array([result_array[i*t_fin:(i+1)*t_fin,0:N] @ U_out_total[i*N:(i+1)*N,:] * result_array[i*t_fin:(i+1)*t_fin,N:N+M] for i in range(tv)])
    # uptake_end = np.array([uptake[i][t_fin-1] for i in range(tv)])
    
    sur = [np.where(result_array[(i+1)*t_fin-1, 0:N] > 0.01) for i in range(tv)]# sur 

    U_range = [np.arange(0, np.ceil(np.max(U_ac_total)/10)*10, 10) for i in range(tv)]
    sur_rate = np.empty((0,int(np.ceil(np.max(U_ac_total)/10))))
    s_total = np.empty((0))
    s_sur = np.empty((0))
    for j in range(tv):
        for i in range(len(U_range[j])):
            s_total = np.append(s_total, ((U_ac_total[j,:] >= U_range[j][i]) & (U_ac_total[j,:] < U_range[j][i]+10)).sum())
            s_sur = np.append(s_sur, ((U_ac_total[j,:][sur[j][0]] >= U_range[j][i]) & (U_ac_total[j,:][sur[j][0]] < U_range[j][i]+10)).sum())
    
    sur_rate = np.array(s_sur/s_total).reshape(tv,int(np.ceil(np.max(U_ac_total)/10)))
    sur_rate_filtered = [i[j] for i, j in zip(sur_rate.T, (~np.isnan(sur_rate)).T)]
    
    return sur_rate_filtered



# sur_rate_filtered = plot_comp(t_fin, N, M, T, Tref, Ma, tv, Ea_D, typ, K)

# sur_rate_mean = []
# sur_rate_ci = []
# # using np.mean()
# for i in range(len(sur_rate_filtered)):
#     sur_rate_mean.append(np.mean(sur_rate_filtered[i]))
#     sur_rate_ci.append(1.96 * np.std(sur_rate_filtered[i],axis = 0)/(len(sur_rate_filtered[i])**0.5))
    
# plt.scatter(range(0,len(sur_rate_mean)*10,10),sur_rate_mean, c = "r")
# plt.errorbar(range(0,len(sur_rate_mean)*10,10),sur_rate_mean, yerr=sur_rate_ci, fmt="o", color='g')
# plt.legend(range(20,40,5))
# plt.title(T-273.15)
# plt.show()

