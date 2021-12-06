import numpy as np
from scipy.integrate import solve_ivp
# from scipy.integrate import odeint
# import matplotlib.pylab as plt
# from matplotlib.lines import Line2D
import size_temp_funcs as st
import parameters as par
import model_func as mod


######### Main Code ###########

######## Set up parameters ###########

N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 28 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_CUE = 0.3
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly number, i.e. how many times the system can assemble
t_fin = 2000 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 5 # Half saturation constant
# strc = 1

##### Intergrate system forward #####

def ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K):
    '''
    Main function for the simulation of resource uptake and growth of microbial communities.
    '''

    # Setted Parameters
    k = 0.0000862 # Boltzman constant
    a = 15 # The alpha value for beta distribution in Ea
    B_R0 = 1.70 * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref))) # Using CUE0 = 0.22, mean growth rate = 0.48
    B_U0 = (1.70/(1 - lf - 0.22)) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref)))

    ### Creating empty array for storing data ###
    result_array = np.empty((0,N+M)) # Array to store data in for plotting
    rich_series = np.empty((0))
    U_out_total = np.empty((0,M))
    R_out = np.empty((0, N))
    CUE_out = np.empty((0,N))
    Sr = np.empty((0,N))
    Ea_CUE_out = np.empty((0,N))
    overlap = np.empty((0,N))
    crossf = np.empty((0,N))

    for i in range(ass):

        ### Resetting values for every assembly ###

        x0 = np.concatenate((np.full([N], 0.1), np.full([M], 1))) # Starting concentration for resources and consumers

        # Set up Ea (activation energy) and B0 (normalisation constant) based on Tom Smith's observations
        B_U = np.random.normal(B_U0, 0.1*B_U0, N) # Adding variation into B0
        B_R = np.random.normal(B_R0, 0.1*B_R0, N)
        T_pk_U = 273.15 + np.random.normal(35, 5, size = N)
        T_pk_R = T_pk_U + 3
        Ea_U = np.random.beta(a, ((a - 1/3) / (0.82/4)) + 2/3 - a, N)*4
        # Ea_R = Ea_U
        Ea_R = np.random.beta(a, ((a - 1/3) / (0.67/4)) + 2/3 - a, N)*4

        p = np.repeat(p_value, M)  # Resource input

        # Set up model
        U, R, l = par.params(N, M, T, k, Tref, T_pk_U, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, np.repeat(Ea_D,N), lf) # Uptake
        l_sum = np.sum(l, axis=1)
        
        # With cost
        # R_cost = R*(1-lf)*np.sum(U/np.sum(U), axis = 1)*150

        # Integration
        t = np.linspace(0,t_fin-1,t_fin) 
        pars = (U, R, l, p, l_sum, N, M, typ, K) # Parameters to pass onto model
        # pars = (U, R_cost, l, p, l_sum, N, M, typ, K) # Parameters to pass onto model

        # pops, infodict = odeint(mod.metabolic_model, y0=x0, t=t, args = pars, full_output=1) # Integrate
        pops = solve_ivp(mod.metabolic_model, t_span= [0,t_fin], y0=x0, t_eval = t, args = pars, method = 'BDF') # Integrate
        pops.y = np.transpose(np.round(pops.y, 7))

        ### Storing simulation results ###
        result_array = np.append(result_array, pops.y, axis=0)
        U_out_total = np.append(U_out_total, U, axis = 0)
        R_out = np.append(R_out, [R], axis = 0)


        ### Analysis ###

        # Competition for resources
        jaccard = np.array([[np.sum(np.minimum(U[i,],U[j,]))/np.sum(np.maximum(U[i,],U[j,])) for j in range(N)] for i in range(N)])
        comp = np.mean(jaccard, axis = 0)
        overlap = np.append(overlap, [comp], axis = 0)

        # Cross-feeding
        leak = U@l
        cf = np.array([[np.sum(np.minimum(leak[i], U[j]))/np.sum(np.maximum(leak[i], U[j])) for j in range(N)] for i in range(N)])
        crossf = np.append(crossf, [np.nanmean(cf, axis = 1)], axis = 0)

        # Richness
        rich_series = np.append(rich_series, [len(np.where(pops.y[t_fin-1,0:N])[0])]) # Get the consumer concentrations from last timestep and find out survivors

        # CUE
        CUE = (U @ (1 - l_sum) - R)/(np.sum(U, axis = 1))
        CUE_out = np.append(CUE_out, [CUE], axis = 0)
        Ea_CUE_out = np.append(Ea_CUE_out, [B_R*(Ea_U - Ea_R)/(B_U*(1 - lf) - B_R)], axis = 0)

        # S*
        Sr = np.append(Sr, [R/(np.sum(U, axis = 1)*(1-lf))], axis = 0)


    return result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr

result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)



# rich = np.array([len(np.where(result_array[i,0:N])[0]) for i in range(len(result_array))]).reshape(ass,t_fin)
# rich_mean = np.mean(rich, axis = 0)
# rich_ci = 1.96 * np.std(rich,axis = 0)/(ass**0.5)

# t_plot = np.linspace(0,t_fin,t_fin)

# plt.plot(t_plot, rich_mean)
# plt.fill_between(t_plot, rich_mean - rich_ci, rich_mean + rich_ci, color = 'b', alpha=.1)
# plt.ylabel('Richness')
# plt.xlabel('Time')
# plt.show()


# t_plot = np.linspace(0,t_fin,t_fin)

# plt.plot(t_plot, result_array[:,N:N+M], 'b-', linewidth=0.7)
# plt.ylabel('Resources')
# plt.xlabel('Time')
# # plt.title('Consumer-Resource population dynamics')
# plt.show()

# plt.plot(t_plot, result_array[:,0:N], 'g-', linewidth=0.7)
# plt.ylabel('Consumers')
# plt.xlabel('Time')
# plt.show()