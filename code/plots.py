from Bacteria_vector_modular import ass_temp_run
import numpy as np
import matplotlib.pylab as plt
from matplotlib.lines import Line2D

########## Setting Parameters ###########
N = 1 # Number of consumers
M = 1 # Number of resources

# Temperature params
T = 273.15 + 25 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
pk_U = np.random.normal(25, 3, size = N)
pk_R = pk_U + 2
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy

# Assembly
ass = 5 # Assembly times at each temperature
t_fin = 100 # Number of time steps for each temperature
x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (0.1)))) # Starting concentration for resources and consumers
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)

cc = 2 # Times for increasing resource concentration for the next assembly

def plot_con_res_CUE(t_fin, N, M, T,  Tref, Ma, ass, x0, pk_R, pk_U, Ea_D, typ, K, cc):
    
    result_array, CUE_out, CUE_com_out, CUE_out_U, CUE_com_U_out, rich_d = ass_temp_run(t_fin, N, M, T,  Tref, Ma, ass, x0, pk_R, pk_U, Ea_D, typ, K, cc)[3:9]

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

    plt.plot(t_plot, CUE_out, 'r-', label = 'Species level', linewidth=0.7)
    plt.plot(t_plot, CUE_com_out, 'k-', label = 'Community level', linewidth = 0.5)
    # plt.ylim(bottom = -1)
    plt.ylabel('CUE')
    plt.xlabel('Time')
    plt.title('Carbon Use Efficiency dynamics')
    plt.legend([Line2D([0], [0], color='red', lw=2), Line2D([0], [0], color='black', lw=2)], ['Species level', 'Community level'])
    plt.show()

    plt.plot(t_plot, CUE_out_U, 'r-', label = 'Species level', linewidth=0.7)
    plt.plot(t_plot, CUE_com_U_out, 'k-', label = 'Community level', linewidth = 0.5)
    # plt.ylim(bottom = -1)
    plt.ylabel('CUE (U & R)')
    plt.xlabel('Time')
    plt.title('Carbon Use Efficiency dynamics(U & R)')
    plt.legend([Line2D([0], [0], color='red', lw=2), Line2D([0], [0], color='black', lw=2)], ['Species level', 'Community level'])
    plt.show()

    plt.plot(t_plot, rich_d, 'c-', linewidth=0.7)
    plt.ylabel('Richness')
    plt.xlabel('Time')
    plt.show()

    return 

plot_con_res_CUE(t_fin, N, M, T,  Tref, Ma, ass, x0, pk_R, pk_U, Ea_D, typ, K, cc)