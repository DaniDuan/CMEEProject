from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import numpy as np

########## Setting Parameters ###########
N = 10 # Number of consumers
M = 5 # Number of resources

# Temperature params
T = 273.15 + 25 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
# pk = 20 # Peak above Tref, degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
t_n = 30 # Number of temperatures to run the model at, model starts at 20

# Assembly
ass = 6 # Assembly number, i.e. how many times the system can assemble
tv = 10 # immigration times inside one assembly
t_fin = 50 # Number of time steps
x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (0.1)))) # Starting concentration for resources and consumers
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant

T_c = 6 # How many temperatures to cover (how many cycles to run)

############# Defining a Function for Running Model ##########
def funcs_with_temp(T_c, t_fin, N, M, Tref, Ma, ass, tv, x0, Ea_D, typ, K): 
    rich_temp_mean = np.empty((0))
    rich_temp_ci = np.empty((0))

    for i in range(T_c):
        T = 273.15 + 10 + 5 * i # Temperature
        rich_seires = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, x0, Ea_D, typ, K)[1]
        rich_mean = np.mean(rich_seires, axis = 0)
        rich_ci = 1.96 * np.std(rich_seires,axis = 0)/(ass**0.5)
        rich_temp_mean = np.append(rich_temp_mean, rich_mean[tv-1])
        rich_temp_ci = np.append(rich_temp_ci, rich_ci[tv-1])
    
    T_plot = range(10, 10 + 5*T_c, 5)
    plt.plot(T_plot, rich_temp_mean, 'b-', linewidth=0.7)
    plt.fill_between(T_plot, rich_temp_mean - rich_temp_ci, rich_temp_mean + rich_temp_ci, color='b', alpha=.1)
    plt.ylabel('Richness')
    plt.xlabel('Temp')
    plt.show()
    
    return rich_temp_mean

    
funcs_with_temp(T_c, t_fin, N, M, Tref, Ma, ass, tv, x0, Ea_D, typ, K)