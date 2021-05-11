from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 15 # Number of consumers
M = 15 # Number of resources

# Temperature params
T = 273.15 + 25 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance


# Assembly
ass = 6 # Assembly number, i.e. how many times the system can assemble
tv = 300 # immigration times inside one assembly
t_fin = 50 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant


p_value = 1 # External input resource concentration
rich_mean = np.empty((0,9))
for i in range(10):
    Ea_diff = 0.1*i
    rich_clct = []
    for j in range(9):
        lf = (j + 1)*0.1 # Leakage
        rich_series = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, tv, Ea_D, Ea_diff, lf, p_value, typ, K)[1]
        rich_clct = np.append(rich_clct, np.mean(rich_series[:,tv-1]))
    rich_mean = np.append(rich_mean, [rich_clct], axis = 0)


Xlabel = np.round(np.arange(0.1,1,0.1), 1)
Ylabel = np.round(1.5 - np.arange(0,1,0.1), 1)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yticks(range(len(Ylabel)))
ax.set_xticks(range(len(Xlabel)))
ax.set_yticklabels(Ylabel)
ax.set_xticklabels(Xlabel)
plt.xlabel("leakage")
plt.ylabel("Respiration Ea")
im = plt.imshow(rich_mean, cmap = 'OrRd')
plt.colorbar(im)
plt.title( "2-D Heat Map" )
plt.show()