from Bacteria_vector_modular import ass_temp_run
import numpy as np
import matplotlib.pylab as plt
from matplotlib.lines import Line2D

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 0 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin, 10 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
# Ea_CUE = 0.3
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly number, i.e. how many times the system can assemble
tv = 1 # immigration times inside one assembly
t_fin = 3000 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant


def plot_con_res_CUE(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, lf, p_value, typ, K):
    '''
    Plotting time series of resource, species biomass and community richness.
    '''
    result_array, rich_series, l, U_out_total, U_ac_total, R_out, CUE_out, Ea_CUE_out = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, lf, p_value, typ, K)

    t_plot = np.linspace(0,len(result_array),len(result_array))
    
    plt.plot(t_plot, result_array[:,N:N+M], 'b-', linewidth=0.7)
    plt.ylabel('Resources')
    plt.xlabel('Time')
    # plt.title('Consumer-Resource population dynamics')
    plt.show()

    plt.plot(t_plot, result_array[:,0:N], 'g-', linewidth=0.7)
    plt.ylabel('Consumers')
    plt.xlabel('Time')
    plt.show()

    # plt.plot(t_plot, CUE_out, 'r-', label = 'Species level', linewidth=0.7)
    # # plt.plot(t_plot, CUE_com_out, 'k-', label = 'Community level', linewidth = 0.5)
    # # plt.ylim(bottom = -1)
    # plt.ylabel('CUE')
    # plt.xlabel('Time')
    # plt.title('Carbon Use Efficiency dynamics')
    # # plt.legend([Line2D([0], [0], color='red', lw=2), Line2D([0], [0], color='black', lw=2)], ['Species level', 'Community level'])
    # plt.show()

    # t_rich = np.linspace(t_fin, tv*t_fin, tv)
    # rich_mean = np.mean(rich_series, axis = 0)
    # rich_ci = 1.96 * np.std(rich_series,axis = 0)/(ass**0.5)
    # plt.plot(t_rich, rich_mean, 'c-', linewidth=0.7)
    # plt.fill_between(t_rich, rich_mean - rich_ci, rich_mean + rich_ci, color='b', alpha=.1)
    # plt.ylabel('Richness')
    # plt.xlabel('Time')
    # plt.show()


    # # CUE_mean = np.mean(CUE_series[[np.arange(tv*ass, step = tv)+i for i in range(tv)], :], axis = 1)
    # t_CUE = np.linspace(t_fin, ass*tv*t_fin, ass*tv)
    # plt.plot(t_CUE, CUE_series, 'r-', linewidth=0.7)
    # plt.ylabel('CUE')
    # plt.xlabel('Time')
    # plt.show()

    # CUE_c_mean = np.mean(CUE_c_series, axis = 0)
    # # CUE_ci = 1.96 * np.std(CUE_c_series,axis = 0)/CUE_c_mean
    # plt.plot(t_rich, CUE_c_mean, 'c-', linewidth=0.7)
    # # plt.fill_between(t_rich, CUE_c_mean - CUE_ci, CUE_c_mean + CUE_ci, color='b', alpha=.1)
    # plt.ylabel('CUE(Community level)')
    # plt.xlabel('Time')
    # plt.show()

    return  #CUE_out

# import time
# start = time.time()

# plot_con_res_CUE(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, Ea_diff, lf, p_value, typ, K)

# print((time.time() - start)/60)