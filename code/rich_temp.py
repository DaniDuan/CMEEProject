from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import numpy as np

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
Tref = 273.15 + 0 # Reference temperature Kelvin, 10 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_CUE = 0.3
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly number, i.e. how many times the system can assemble
tv = 1 # immigration times inside one assembly
t_fin = 3000 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant

T_c = 6 # How many temperatures to cover (how many cycles to run)

############# Defining a Function for Running Model ##########
def funcs_with_temp(T_c, t_fin, N, M, Tref, Ma, ass, tv, Ea_D, lf, p_value, typ, K): 
    '''
    Returning community richness at different temperatures.
    '''
    
    # rich_temp_mean = np.empty((0))
    # rich_temp_ci = np.empty((0))
    CUE_sur = np.empty((0, ass*tv))
    U_sur = np.empty((0, ass*tv))
    R_sur = np.empty((0, ass*tv))
    rich_out = np.empty((0, ass*tv))


    for i in range(T_c):
        T = 273.15 + 5 * i # Temperature
        result_array, rich_series, l, U_out_total, U_ac_total, R_out, CUE_out, Ea_CUE_out, overlap = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, lf, p_value, typ, K)
        sur = [np.where(result_array[(i+1)*t_fin-1, 0:N]) for i in range(tv*ass)]
        CUE_sur = np.append(CUE_sur, [np.array([np.mean(CUE_out[i][sur[i]]) for i in range(len(sur))])], axis = 0)
        U = np.sum(U_out_total, axis = 1).reshape(tv*ass, N)
        U_sur = np.append(U_sur, [np.array([np.mean(U[i][sur[i]]) for i in range(len(sur))])], axis = 0)
        R_sur = np.append(R_sur, [np.array([np.mean(R_out[i][sur[i]]) for i in range(len(sur))])], axis = 0)
        rich_out = np.append(rich_out, [rich_series.flatten()], axis = 0)
        # rich_mean = np.mean(rich_series, axis = 0)
        # rich_ci = 1.96 * np.std(rich_series,axis = 0)/(ass**0.5)
        # rich_temp_mean = np.append(rich_temp_mean, rich_mean[tv-1])
        # rich_temp_ci = np.append(rich_temp_ci, rich_ci[tv-1])
        
    return CUE_sur, U_sur, R_sur, rich_out

CUE_sur, U_sur, R_sur, rich_out = funcs_with_temp(T_c, t_fin, N, M, Tref, Ma, ass, tv, Ea_D, lf, p_value, typ, K)

# CUE_mean = np.nanmean(CUE_sur_out, axis = 1)
# CUE_CI = 1.96 * np.nanstd(CUE_sur_out,axis = 1)/((CUE_sur_out.shape[1])**0.5)

def mean_ci(output_matrix, T_c, ass, tv):
    mean = []
    ci = []
    for i in range(T_c):
        om_mean = np.mean(output_matrix[i].reshape(ass,tv), axis = 0)  
        om_ci = 1.96 * np.std(output_matrix[i].reshape(ass,tv),axis = 0)/(ass**0.5)
        mean = np.append(mean, om_mean[tv-1])
        ci = np.append(ci, om_ci[tv-1])
    return mean, ci

# rich_mean, rich_ci = mean_ci(rich_out, T_c, ass, tv)
# CUE_mean, CUE_ci = mean_ci(CUE_sur, T_c, ass, tv)
# U_mean, U_ci = mean_ci(U_sur, T_c, ass, tv)
# R_mean, R_ci = mean_ci(R_sur, T_c, ass, tv)

# plt.scatter(np.mean(rich_out, axis = 0), np.mean(CUE_sur, axis = 0))
# plt.show()

# T_plot = range(0, 5*T_c, 5)
# plt.plot(T_plot, CUE_mean, 'r-', linewidth=0.7, label = 'CUE')
# plt.plot(T_plot, Ea_CUE_mean, 'maroon', linewidth=0.7, label = 'Ea')
# plt.fill_between(T_plot, CUE_mean - CUE_CI, CUE_mean + CUE_CI, color='r', alpha=.1)
# plt.fill_between(T_plot, Ea_CUE_mean - Ea_CUE_CI, Ea_CUE_mean + Ea_CUE_CI, color='maroon', alpha=.1)
# plt.ylabel('CUE & Ea')
# plt.xlabel('Temp')
# plt.legend()
# plt.show()

# T_plot = range(0, 5*T_c, 5)
# plt.plot(T_plot, rich_mean, linewidth=0.7)
# plt.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, alpha=.1)
# plt.ylabel('Richness')
# plt.xlabel('Temp')
# plt.show()


# plt.plot(rich_temp_mean, CUE_mean)
# plt.show()