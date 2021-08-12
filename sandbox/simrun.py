import time
start = time.time()

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
Tref = 273.15 + 15 # Reference temperature Kelvin, 10 degrees C !!!
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly times at each temperature
t_fin = 4000 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)
T_c = 31 # How many temperatures to cover (how many cycles to run)


rich = np.empty((0, ass))
sU = []
sR = []
eU = []
eR = []
U_var = []
extU_var = []
all_CUE = np.empty((0))
all_CUE_ci = np.empty((0))
sur_var = np.empty((0))
sur_var_ci = np.empty((0))
ext_var = np.empty((0))
ext_var_ci = np.empty((0))
all_Ea = np.empty((0))
all_Ea_ci = np.empty((0))
sur_CUE = []
sur_Ea = []
ext_CUE = []
ext_Ea = []
sur_overlap = []
ext_overlap = []
sur_crossf = []
ext_crossf = []


for i in range(T_c):
    T = 273.15 + i # Temperature
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
    rich = np.append(rich, [rich_series.flatten()], axis = 0)
    
    sur = [np.where(result_array[(i+1)*t_fin-1, 0:N]) for i in range(ass)]
    ext = [np.where(result_array[(i+1)*t_fin-1, 0:N] == 0) for i in range(ass)]
    
    sU.append(np.concatenate([np.sum(U_out_total, axis = 1).reshape(ass, N)[i][sur[i]] for i in range(len(sur))]).ravel())
    sR.append(np.concatenate([R_out[i][sur[i]] for i in range(len(sur))]).ravel())
    eU.append(np.concatenate([np.sum(U_out_total, axis = 1).reshape(ass, N)[i][ext[i]] for i in range(len(ext))]).ravel())
    eR.append(np.concatenate([R_out[i][ext[i]] for i in range(len(ext))]).ravel())
    U_var.append(np.concatenate([np.var(U_out_total,axis = 1).reshape(ass, N)[i][sur[i]] for i in range(len(sur))]).ravel())
    extU_var.append(np.concatenate([np.var(U_out_total,axis = 1).reshape(ass, N)[i][ext[i]] for i in range(len(ext))]).ravel())
    
    all_CUE = np.append(all_CUE, [np.mean(CUE_out)])
    all_CUE_ci = np.append(all_CUE_ci, [1.96 * np.std(CUE_out)/((ass*N)**0.5)])
    
    sur_CUE.append(np.concatenate([CUE_out[i][sur[i]] for i in range(len(sur))]))
    sur_var = np.append(sur_var, [np.nanmean([np.var(CUE_out[i][sur[i]]) for i in range(len(sur))])])
    sur_var_ci = np.append(sur_var_ci, [1.96 * np.nanstd([np.var(CUE_out[i][sur[i]]) for i in range(len(sur))])/(np.sum(rich_series)**0.5)])
    sur_Ea.append(np.concatenate([Ea_CUE_out[i][sur[i]] for i in range(len(sur))]))
    ext_CUE.append(np.concatenate([CUE_out[i][ext[i]] for i in range(len(ext))]))
    ext_Ea.append(np.concatenate([Ea_CUE_out[i][ext[i]] for i in range(len(ext))]))
    ext_var = np.append(ext_var, [np.mean([np.var(CUE_out[i][ext[i]]) for i in range(len(ext))])])
    ext_var_ci = np.append(ext_var_ci, [1.96 * np.std([np.var(CUE_out[i][ext[i]]) for i in range(len(ext))])/(np.sum(N-rich_series)**0.5)])
    
    all_Ea = np.append(all_Ea, [np.mean(Ea_CUE_out)])
    all_Ea_ci = np.append(all_Ea_ci, [1.96 * np.std(Ea_CUE_out)/((ass*N)**0.5)])
    
    sur_overlap.append(np.concatenate([overlap[i][sur[i]] for i in range(len(sur))]))
    ext_overlap.append(np.concatenate([overlap[i][ext[i]] for i in range(len(ext))]))
    
    sur_crossf.append(np.concatenate([crossf[i][sur[i]] for i in range(len(sur))]))
    ext_crossf.append(np.concatenate([crossf[i][ext[i]] for i in range(len(ext))]))
    
print((time.time() - start)/60)


rich_mean = np.nanmean(rich, axis = 1)
rich_ci =  1.96 * np.nanstd(rich,axis = 1)/(ass**0.5)

CUE_mean = np.array([np.mean(sur_CUE[i]) for i in range(T_c)])
CUE_ci = np.array([1.96 * np.std(sur_CUE[i])/(len(ext_CUE[i])**0.5) for i in range(T_c)])
CUE_ext_mean = np.array([np.mean(ext_CUE[i]) for i in range(T_c)])
CUE_ext_ci = np.array([1.96 * np.std(ext_CUE[i])/(len(ext_CUE[i])**0.5) for i in range(T_c)])

Ea_mean = np.array([np.mean(sur_Ea[i]) for i in range(T_c)])
Ea_ci = np.array([1.96 * np.std(sur_Ea[i])/(len(sur_Ea[i])**0.5) for i in range(T_c)])

sU_mean = np.array([np.mean(sU[i]) for i in range(T_c)])
sU_ci = np.array([1.96 * np.std(sU[i])/(len(sU[i])**0.5) for i in range(T_c)])
sR_mean = np.array([np.mean(sR[i]) for i in range(T_c)])
sR_ci = np.array([1.96 * np.std(sR[i])/(len(sR[i])**0.5) for i in range(T_c)])

U_var_mean = np.array([np.mean(U_var[i]) for i in range(T_c)])
U_var_ci = np.array([1.96 * np.std(U_var[i])/(len(U_var[i])**0.5) for i in range(T_c)])
extU_var_mean = np.array([np.mean(extU_var[i]) for i in range(T_c)])
extU_var_ci = np.array([1.96 * np.std(extU_var[i])/(len(extU_var[i])**0.5) for i in range(T_c)])

eU_mean = np.array([np.mean(eU[i]) for i in range(T_c)])
eU_ci = np.array([1.96 * np.std(eU[i])/(len(eU[i])**0.5) for i in range(T_c)])
eR_mean = np.array([np.mean(eR[i]) for i in range(T_c)])
eR_ci = np.array([1.96 * np.std(eR[i])/(len(eR[i])**0.5) for i in range(T_c)])

rich_sur = [[np.repeat(rich[i][j], rich[i][j]) for j in range(ass)] for i in range(T_c)]
rich_ext = [[np.repeat(rich[i][j], N - rich[i][j]) for j in range(ass)] for i in range(T_c)]

overlap_sur_mean = np.array([np.mean(sur_overlap[i]) for i in range(T_c)])
overlap_sur_ci = np.array([1.96 * np.std(sur_overlap[i])/(len(sur_overlap[i])**0.5) for i in range(T_c)])
overlap_ext_mean = np.array([np.mean(ext_overlap[i]) for i in range(T_c)])
overlap_ext_ci = np.array([1.96 * np.std(ext_overlap[i])/(len(ext_overlap[i])**0.5) for i in range(T_c)])

crossf_sur_mean = np.array([np.mean(sur_crossf[i]) for i in range(T_c)])
crossf_sur_ci = np.array([1.96 * np.std(sur_overlap[i])/(len(sur_overlap[i])**0.5) for i in range(T_c)])
crossf_ext_mean = np.array([np.mean(ext_crossf[i]) for i in range(T_c)])
crossf_ext_ci = np.array([1.96 * np.std(ext_overlap[i])/(len(ext_overlap[i])**0.5) for i in range(T_c)])

Ea = []
for i in range(T_c): 
    n = 0
    for j in rich[i]:
        j = int(j)
        Ea.append(sur_Ea[i][n:n+j])
        n = n + j

A = list(zip(rich.flatten(),Ea))
A.sort(key = lambda x: x[0])

rich_v = []
for x in rich.flatten():
    if x not in rich_v:
        rich_v.append(x)
rich_v = np.sort(rich_v)


Ea_sorted = []
for i in rich_v:
    sorting = np.empty((0))
    for j in range(len(A)):
        if [x[0] for x in A][j] == i:
            sorting = np.append(sorting, A[j][1])
    Ea_sorted.append(sorting)

    
meanEa = []
for i in range(len(Ea_sorted)):
    meanEa.append(np.mean(Ea_sorted[i]))
    


T_plot = range(0, T_c, 1)
T_sur = [[np.repeat(T_plot[i], rich[i][j]) for j in range(ass)] for i in range(T_c)]

plt.plot(T_plot, rich_mean)
plt.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, color = 'b', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Richness')
# plt.text(-5,28,'A',fontsize= 'xx-large')
plt. show()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ln1 = ax1.plot(T_plot, sU_mean,'darkorange', label = "Survivor Uptake Rate")
ax1.fill_between(T_plot, sU_mean - sU_ci, sU_mean + sU_ci, color='darkorange', alpha=.1)
ln3 = ax1.plot(T_plot, eU_mean,'tan', label = "Extinct Uptake Rate")
ax1.fill_between(T_plot, eU_mean - eU_ci, eU_mean + eU_ci, color='tan', alpha=.1)

ln2 = ax2.plot(T_plot, sR_mean, 'g', label = "Survivor Respiration Rate")
ax2.fill_between(T_plot, sR_mean - sR_ci, sR_mean + sR_ci, color='g', alpha=.1)
ln4 = ax2.plot(T_plot, eR_mean, 'darkseagreen', label = "Extinct Respiration Rate")
ax2.fill_between(T_plot, eR_mean - eR_ci, eR_mean + eR_ci, color='darkseagreen', alpha=.1)

ax1.set_xlabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Uptake rate')
ax2.set_ylabel('Repiration Rate')

lns = ln1+ln2+ln3+ln4
ax1.legend(lns, [i.get_label() for i in lns], loc = 2)
# ax1.text(-5,400,'B',fontsize= 'x-large')
plt. show()

plt.plot(T_plot, CUE_mean, 'r', label = "Survivor")
plt.fill_between(T_plot, CUE_mean - CUE_ci, CUE_mean + CUE_ci, color='r', alpha=.1)
plt.plot(T_plot, CUE_ext_mean, 'dimgrey', label = "Extinct")
plt.fill_between(T_plot,  CUE_ext_mean - CUE_ext_ci,  CUE_ext_mean + CUE_ext_ci, color='dimgrey', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('CUE')
plt.legend()
# plt.text(-5,0.57,'C',fontsize= 'x-large')
plt. show()


n = 0
for i in range(T_c):
    plt.scatter(ext_Ea[i], np.concatenate(rich_ext[i]), color = 'lightgrey', alpha = 0.5, s = 10)
    if n == 0: 
        plt.scatter(sur_Ea[i], np.concatenate(rich_sur[i]), color = 'b', alpha = 0.7, s = 10, label = 'Survivors')
        n = 1
    else:
        plt.scatter(sur_Ea[i], np.concatenate(rich_sur[i]), color = 'b', alpha = 0.7, s = 10)

m, b = np.polyfit(np.delete(meanEa,np.where(np.isnan(meanEa))), np.delete(rich_v,np.where(np.isnan(meanEa))), 1)
x = np.arange(0.07, 1.45, 0.01)
# plt.plot(x, m*x + b, color = 'maroon',linewidth=2)
plt.xlabel('Thermal Sensitivity(Ea) of CUE')
plt.ylabel('Richness')
plt.legend(loc = 1)
# plt.text(-2.2,27,'D',fontsize= 'x-large')
plt.show()


plt.plot(T_plot, overlap_sur_mean, color = 'b', alpha = 0.7, label = 'Resource Overlap - Survivor')
plt.plot(T_plot, overlap_ext_mean, color = 'lightsteelblue', label = 'Resource Overlap - Extinct')
plt.fill_between(T_plot, overlap_sur_mean - overlap_sur_ci, overlap_sur_mean + overlap_sur_ci, color = 'b', alpha=.1)
plt.fill_between(T_plot, overlap_ext_mean - overlap_ext_ci, overlap_ext_mean + overlap_ext_ci, color = 'lightsteelblue', alpha=.1)
plt.plot(T_plot, crossf_sur_mean, color = 'r', alpha = 0.7, label = 'Facilitation - Survivor')
plt.plot(T_plot, crossf_ext_mean, color = 'rosybrown', label = 'Facilitation - Extinct')
plt.fill_between(T_plot, crossf_sur_mean - crossf_sur_ci, crossf_sur_mean + crossf_sur_ci, color = 'r', alpha=.1)
plt.fill_between(T_plot, crossf_ext_mean - crossf_ext_ci, crossf_ext_mean + crossf_ext_ci, color = 'rosybrown', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Pair-wise Interactions')
plt.legend()
# plt.text(-4.5,0.35,'E',fontsize= 'x-large')
plt.show()


plt.plot(T_plot, sur_var,'orangered', label = "Survivor")
plt.fill_between(T_plot, sur_var - sur_var_ci, sur_var + sur_var_ci, color='orangered', alpha=.1)
plt.plot(T_plot, ext_var,'grey', label = "Extinct")
plt.fill_between(T_plot, ext_var - ext_var_ci, ext_var + ext_var_ci, color='grey', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Variance of CUE')
plt.legend()
# plt.text(-4.5,0.31,'F',fontsize= 'x-large')
plt.show()
print(sur_var)

plt.plot(T_plot, U_var_mean,'darkorange', label = "Survivor")
plt.fill_between(T_plot, U_var_mean - U_var_ci, U_var_mean + U_var_ci, color='darkorange', alpha=.1)
plt.plot(T_plot, extU_var_mean,'grey', label = "Extinct")
plt.fill_between(T_plot, extU_var_mean - extU_var_ci, extU_var_mean + extU_var_ci, color='grey', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Variance of Uptake')
plt.legend(loc = 2)
# plt.text(-4.5,0.31,'F',fontsize= 'x-large')
plt.show()



########################################################################################################

import time
start = time.time()

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly times at each temperature
t_fin = 4000 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)
T_c = 7 # How many temperatures to cover (how many cycles to run)


rich_mean = np.empty((0,T_c))
var_Ea = np.empty((0,T_c))
for i in range(T_c):
    Tref = 273.15 + 5 * i
    all_var_Ea = []
    rich = np.empty((0,ass))
    for j in range(T_c):
        T = 273.15 + 5 * j # Temperature
        result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
        rich = np.append(rich, [rich_series.flatten()], axis = 0)
        sur = [np.where(result_array[(i+1)*t_fin-1, 0:N]) for i in range(ass)]
        all_var_Ea = np.append(all_var_Ea, [np.mean([np.var(Ea_CUE_out[i][sur[i]]) for i in range(len(sur))])])
        # all_var = np.append(all_var, [np.mean(np.var(CUE_out, axis = 1))])
    var_Ea = np.append(var_Ea, [all_var_Ea], axis = 0)
    rich_mean = np.append(rich_mean, [np.mean(rich, axis = 1)], axis = 0)
    
print((time.time() - start)/60)

Xlabel = np.round(np.arange(0, 5 * T_c, 5), 1)
Ylabel = np.round(np.arange(30, -5, -5),1)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yticks(range(len(Ylabel)))
ax.set_xticks(range(len(Xlabel)))
ax.set_yticklabels(Ylabel)
ax.set_xticklabels(Xlabel)
plt.xlabel("Temperature")
plt.ylabel("Reference Temperature")
im = plt.imshow(np.flip(rich_mean,axis = 0), cmap = 'YlOrBr')
plt.colorbar(im)
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yticks(range(len(Ylabel)))
ax.set_xticks(range(len(Xlabel)))
ax.set_yticklabels(Ylabel)
ax.set_xticklabels(Xlabel)
plt.xlabel("Temperature")
plt.ylabel("Reference Temperature")
im = plt.imshow(np.flip(var_Ea,axis = 0), cmap = 'YlOrBr')
plt.colorbar(im)
plt.show()


###########################################################################################################
import time
start = time.time()

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
Tref = 273.15
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy
p_value = 1 # External input resource concentration

# Assembly
ass = 50 # Assembly times at each temperature
t_fin = 4000 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)
T_c = 7 # How many temperatures to cover (how many cycles to run)


rich_mean = np.empty((0,T_c-1))
rich_ci = np.empty((0,T_c-1))
cf = np.empty((0,T_c-1))

for i in range(5):
    if i == 0:
        lf = 0
    else:
        lf = i*0.2-0.1
    all_cf = []
    rich = np.empty((0,ass))
    for j in range(T_c-1):
        T = 273.15 + 5 * (j+1) # Temperature
        result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
        rich = np.append(rich, [rich_series.flatten()], axis = 0)
        sur = [np.where(result_array[(i+1)*t_fin-1, 0:N]) for i in range(ass)]
        all_cf = np.append(all_cf, [np.mean(np.concatenate([crossf[i][sur[i]] for i in range(len(sur))]))])
    cf = np.append(cf, [all_cf], axis = 0)
    rich_mean = np.append(rich_mean, [np.mean(rich, axis = 1)], axis = 0)
    rich_ci = np.append(rich_ci, [1.96 * np.std(rich, axis = 1)/(ass**0.5)])
    
print((time.time() - start)/60)

# Xlabel = np.round(np.arange(5, 5 * T_c, 5), 1)
# Ylabel = np.array((0.7, 0.5, 0.3, 0.1, 0))
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_yticks(range(len(Ylabel)))
# ax.set_xticks(range(len(Xlabel)))
# ax.set_yticklabels(Ylabel)
# ax.set_xticklabels(Xlabel)
# plt.xlabel("Temperature")
# plt.ylabel("Leakage")
# im = plt.imshow(np.flip(rich_mean,axis = 0), cmap = 'YlOrBr')
# plt.colorbar(im)
# plt.show()
rich_ci = rich_ci.reshape(5,6)
for i in range(5):
    if i == 0:
        lf = 0
    else:
        lf = i*0.2-0.1
    plt.plot(np.arange(5, 5 * T_c, 5), rich_mean[i], label = 'l = %s' %np.round(lf,1), alpha = 0.8)
    plt.fill_between(np.arange(5, 5 * T_c, 5), rich_mean[i] - rich_ci[i], rich_mean[i] + rich_ci[i], alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Richness')
plt.legend()
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yticks(range(len(Ylabel)))
ax.set_xticks(range(len(Xlabel)))
ax.set_yticklabels(Ylabel)
ax.set_xticklabels(Xlabel)
plt.xlabel("Temperature")
plt.ylabel("Leakage")
im = plt.imshow(np.flip(cf,axis = 0), cmap = 'YlOrBr')
plt.colorbar(im)
plt.show()


######################################################################################################################

import time
start = time.time()

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
Tref = 273.15 + 0 # Reference temperature Kelvin, 10 degrees C !!!
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly times at each temperature
t_fin = 4000 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)
T_c = 31 # How many temperatures to cover (how many cycles to run)

rich = np.empty((0, ass))
sur_Sr = []
sur_Ea = []

for i in range(T_c):
    T = 273.15 + i # Temperature
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
    rich = np.append(rich, [rich_series.flatten()], axis = 0)
    sur = [np.where(result_array[(i+1)*t_fin-1, 0:N]) for i in range(ass)]
    ext = [np.where(result_array[(i+1)*t_fin-1, 0:N] == 0) for i in range(ass)]
    
    sur_Sr.append(np.concatenate([Sr[i][sur[i]] for i in range(len(sur))]))
    sur_Ea.append(np.concatenate([Ea_CUE_out[i][sur[i]] for i in range(len(sur))]))

    
print((time.time() - start)/60)


rich_mean = np.nanmean(rich, axis = 1)
rich_ci =  1.96 * np.nanstd(rich,axis = 1)/(ass**0.5)

Sr_mean = np.array([np.mean(sur_Sr[i]) for i in range(T_c)])
Sr_ci = np.array([1.96 * np.std(sur_Sr[i])/(len(sur_Sr[i])**0.5) for i in range(T_c)])

Ea_mean = np.array([np.mean(sur_Ea[i]) for i in range(T_c)])
Ea_ci = np.array([1.96 * (np.std(sur_Ea[i])* ((0.6*B_U) - B_R)/B_R)/(len(sur_Ea[i])**0.5) for i in range(T_c)])
DEa = Ea_mean * ((0.6*B_U) - B_R)/B_R

T_plot = range(0, T_c, 1)

plt.plot(T_plot, rich_mean)
plt.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, color = 'b', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Richness')
# plt.text(-5,28,'A',fontsize= 'xx-large')
plt. show()


plt.scatter(Sr_mean, rich_mean, s=20, color = 'k', alpha = 0.7)
plt.errorbar(Sr_mean, rich_mean, xerr=Sr_ci, fmt=',', color = 'k', alpha = 0.7)
plt.xlabel('Community average S*')
plt.ylabel('Richness')

a,b = np.polyfit(Sr_mean, np.log(rich_mean), 1)
x = np.arange(np.min(Sr_mean), p.max(Sr_mean), 0.001)
y = np.exp(b) * np.exp(a*x)
plt.plot(x, y, color = 'k', alpha = 0.7)

# model = np.poly1d(np.polyfit(Sr_mean, rich_mean, 5))
# plt.plot(x, model(x), color = 'k', alpha = 0.7)

plt.show()