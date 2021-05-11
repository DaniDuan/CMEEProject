from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import numpy as np

########## Setting Parameters ###########
N = 15 # Number of consumers
M = 15 # Number of resources

# Temperature params
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
# pk = 20 # Peak above Tref, degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_diff = 0.6
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration


# Assembly
ass = 6 # Assembly number, i.e. how many times the system can assemble
tv = 300 # immigration times inside one assembly
t_fin = 50 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant

T_c = 9 # How many temperatures to cover (how many cycles to run)

############# Defining a Function for Running Model ##########
def funcs_with_temp(T_c, t_fin, N, M, Tref, Ma, ass, tv, Ea_D, Ea_diff, lf, p_value, typ, K): 
    '''
    Returning community richness at different temperatures.
    '''
    
    rich_temp_mean = np.empty((0))
    rich_temp_ci = np.empty((0))
    simp = np.empty((0))
    CUE_sur_out = np.empty((0, ass*tv))
    Ea_CUE_sur_out = np.empty((0,ass*tv))

    for i in range(T_c):
        T = 273.15 + 5 * i # Temperature
        result_array, rich_series, l, U_out_total, U_ac_total, R_out, CUE_out, Ea_CUE_out = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, Ea_diff, lf, p_value, typ, K)
        sur = [np.where(result_array[(i+1)*t_fin-1, 0:N] > 0.01) for i in range(tv*ass)]
        CUE_sur = np.array([np.mean(CUE_out[i][sur[i]]) for i in range(len(sur))])
        CUE_sur_out = np.append(CUE_sur_out, [CUE_sur], axis = 0)
        Ea_CUE_sur = np.array([np.mean(Ea_CUE_out[i][sur[i]]) for i in range(len(sur))])
        Ea_CUE_sur_out = np.append(Ea_CUE_sur_out, [Ea_CUE_sur], axis = 0)
        rich_mean = np.mean(rich_series, axis = 0)
        rich_ci = 1.96 * np.std(rich_series,axis = 0)/(ass**0.5)
        rich_temp_mean = np.append(rich_temp_mean, rich_mean[tv-1])
        rich_temp_ci = np.append(rich_temp_ci, rich_ci[tv-1])
        
    return rich_temp_mean, rich_temp_ci, CUE_sur_out, Ea_CUE_sur_out

rich_temp_mean, rich_temp_ci, CUE_sur_out, Ea_CUE_sur_out = funcs_with_temp(T_c, t_fin, N, M, Tref, Ma, ass, tv, Ea_D, Ea_diff, lf, p_value, typ, K)


# T_plot = range(0, 5*T_c, 5)
# plt.plot(T_plot, rich_temp_mean, 'b-', linewidth=0.7)
# plt.fill_between(T_plot, rich_temp_mean - rich_temp_ci, rich_temp_mean + rich_temp_ci, color='b', alpha=.1)
# plt.ylabel('Richness')
# plt.xlabel('Temp')
# plt.show()

# plt.legend(title = "Leakage")
# plt.ylabel('Richness')
# plt.xlabel('Temp')
# plt.show()