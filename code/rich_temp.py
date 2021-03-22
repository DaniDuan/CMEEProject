from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import numpy as np

########## Setting Parameters ###########
N = 10 # Number of consumers
M = 5 # Number of resources

# Temperature params
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
pk_U = np.random.normal(25, 3, size = N)
pk_R = pk_U + 2
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy

# Assembly
ass = 8 # Assembly times at each temperature
t_fin = 100 # Number of time steps for each temperature
x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (0.1)))) # Starting concentration for resources and consumers
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)

cc = 0 # Times for increasing resource concentration for the next assembly
T_c = 6 # How many temperatures to cover (how many cycles to run)

############# Defining a Function for Running Model ##########
def funcs_with_temp(T_c, t_fin, N, M, Tref, Ma, ass, x0, pk_R, pk_U, Ea_D, typ, K, cc): 
    rich_st = np.empty((0))

    for i in range(T_c):
        T = 273.15 + 20 + 2 * i # Temperature
        results = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, x0, pk_R, pk_U, Ea_D, typ, K, cc)
        rich_temp = results[0]
        rich_st = np.append(rich_st, rich_temp[ass - 1])
    
    T_plot = range(20, 20 + 2*T_c, 2)
    plt.plot(T_plot, rich_st, 'b-', linewidth=0.7)
    plt.ylabel('Richness')
    plt.xlabel('Temp')
    plt.show()
    
    return rich_st

    
funcs_with_temp(T_c, t_fin, N, M, Tref, Ma, ass, x0, pk_R, pk_U, Ea_D, typ, K, cc)