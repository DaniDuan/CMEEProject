from Bacteria_vector_modular import ass_temp_run
import numpy as np
import matplotlib.pylab as plt
from matplotlib.lines import Line2D

########## Setting Parameters ###########
N = 50 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 25 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy

# Assembly
ass = 2 # Assembly times at each temperature
tv = 5 # Turnovers inside one assembly
t_fin = 100 # Number of time steps for each temperature
x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (0.1)))) # Starting concentration for resources and consumers
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)


def plot_con_res_CUE(t_fin, N, M, T, Tref, Ma, ass, tv, x0, Ea_D, typ, K):
    
    rich, t_rich, result_array, CUE_t = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, x0, Ea_D, typ, K)

    t_plot = np.linspace(0,len(result_array),len(result_array))
    
    plt.plot(t_plot, result_array[:,N:N+M], 'b-', linewidth=0.7)
    plt.ylabel('Resources')
    plt.xlabel('Time')
    # plt.title('Consumer-Resource population dynamics')
    # plt.legend([Line2D([0], [0], color='green', lw=2), Line2D([0], [0], color='blue', lw=2)], ['Consumer', 'Resources'])
    plt.show()

    plt.plot(t_plot, result_array[:,0:N], 'g-', linewidth=0.7)
    plt.ylabel('Consumers')
    plt.xlabel('Time')
    plt.show()

    # plt.plot(t_plot, CUE_out, 'r-', label = 'Species level', linewidth=0.7)
    # plt.plot(t_plot, CUE_com_out, 'k-', label = 'Community level', linewidth = 0.5)
    # # plt.ylim(bottom = -1)
    # plt.ylabel('CUE')
    # plt.xlabel('Time')
    # plt.title('Carbon Use Efficiency dynamics')
    # plt.legend([Line2D([0], [0], color='red', lw=2), Line2D([0], [0], color='black', lw=2)], ['Species level', 'Community level'])
    # plt.show()

    # plt.plot(t_plot, CUE_out_U, 'r-', label = 'Species level', linewidth=0.7)
    # plt.plot(t_plot, CUE_com_U_out, 'k-', label = 'Community level', linewidth = 0.5)
    # # plt.ylim(bottom = -1)
    # plt.ylabel('CUE (U & R)')
    # plt.xlabel('Time')
    # plt.title('Carbon Use Efficiency dynamics(U & R)')
    # plt.legend([Line2D([0], [0], color='red', lw=2), Line2D([0], [0], color='black', lw=2)], ['Species level', 'Community level'])
    # plt.show()

    plt.figure()
    plt.plot(t_rich, rich, 'c-', linewidth=0.7)
    plt.ylabel('Richness')
    plt.xlabel('Time')
    plt.show()

    # plt.figure()
    # plt.plot(t_rich, CUE_t, 'r-', linewidth=0.7)
    # plt.ylabel('CUE')
    # plt.xlabel('Time')
    # plt.show()


    return

plot_con_res_CUE(t_fin, N, M, T, Tref, Ma, ass, tv, x0, Ea_D, typ, K)