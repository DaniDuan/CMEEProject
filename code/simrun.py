import time
start = time.time()

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np
import scipy as sc
from scipy import stats

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
Tref = 273.15 + 0 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 50 # Assembly times at each temperature
t_fin = 4000 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)
T_c = 31 # How many temperatures to cover (how many cycles to run)


rich = np.empty((0, ass))
eq = np.empty((0, ass))
eq_sur = np.empty((0, ass))
eq_sur_ci = np.empty((0, ass))
dEa = np.empty((0, ass))
dEa_sur = np.empty((0, ass))
dEa_sur_ci = np.empty((0, ass))
sU = []
sR = []
eU = []
eR = []
U_var = []
extU_var = []
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
sur_Sr = []
all_Sr = np.empty((0))
sur_overlap = []
ext_overlap = []
sur_crossf = []
ext_crossf = []


for i in range(T_c):
    T = 273.15 + i # Temperature
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
    rich = np.append(rich, [rich_series.flatten()], axis = 0)
    
    sur = [np.where(result_array[(i+1)*t_fin-1, 0:N]) for i in range(ass)]
    ext = [np.where(result_array[(i+1)*t_fin-1, 0:N] == 0) for i in range(ass)]
    
    sU.append(np.concatenate([np.sum(U_out_total, axis = 1).reshape(ass, N)[i][sur[i]] for i in range(len(sur))]).ravel())
    sR.append(np.concatenate([R_out[i][sur[i]] for i in range(len(sur))]).ravel())
    eU.append(np.concatenate([np.sum(U_out_total, axis = 1).reshape(ass, N)[i][ext[i]] for i in range(len(ext))]).ravel())
    eR.append(np.concatenate([R_out[i][ext[i]] for i in range(len(ext))]).ravel())
    U_var.append(np.concatenate([np.var(U_out_total,axis = 1).reshape(ass, N)[i][sur[i]] for i in range(len(sur))]).ravel())
    extU_var.append(np.concatenate([np.var(U_out_total,axis = 1).reshape(ass, N)[i][ext[i]] for i in range(len(ext))]).ravel())
        
    sur_CUE.append(np.concatenate([CUE_out[i][sur[i]] for i in range(len(sur))]))
    sur_var = np.append(sur_var, [np.nanmean([np.var(CUE_out[i][sur[i]]) for i in range(len(sur))])])
    sur_var_ci = np.append(sur_var_ci, [1.96 * np.nanstd([np.var(CUE_out[i][sur[i]]) for i in range(len(sur))])/(np.sum(rich_series)**0.5)])
    sur_Ea.append(np.concatenate([Ea_CUE_out[i][sur[i]] for i in range(len(sur))]))
    ext_CUE.append(np.concatenate([CUE_out[i][ext[i]] for i in range(len(ext))]))
    ext_Ea.append(np.concatenate([Ea_CUE_out[i][ext[i]] for i in range(len(ext))]))
    ext_var = np.append(ext_var, [np.mean([np.var(CUE_out[i][ext[i]]) for i in range(len(ext))])])
    ext_var_ci = np.append(ext_var_ci, [1.96 * np.std([np.var(CUE_out[i][ext[i]]) for i in range(len(ext))])/(np.sum(N-rich_series)**0.5)])
    
    sur_Sr.append(np.concatenate([Sr[i][sur[i]] for i in range(len(sur))]))
    all_Sr = np.append(all_Sr, [np.mean(Sr)])
    
    eq_sur = np.append(eq_sur, [np.array([np.mean(np.abs(Sr[i][sur[i]] - np.mean(Sr,axis=1)[i])) for i in range (ass)])], axis = 0)
    eq_sur_ci = np.append(eq_sur_ci, [np.array([1.96*np.std(np.abs(Sr[i][sur[i]] - np.mean(Sr,axis=1)[i]))/(len(sur[i])**0.5) for i in range (ass)])],axis =0)
    eq = np.append(eq, [np.mean([np.abs(Sr[i,:] - np.mean(Sr,axis=1)[i]) for i in range (ass)], axis = 1)], axis= 0)

    all_Ea = np.append(all_Ea, [np.mean(Ea_CUE_out)])
    all_Ea_ci = np.append(all_Ea_ci, [1.96 * np.std(Ea_CUE_out)/((ass*N)**0.5)])
    dEa_sur = np.append(dEa_sur, [np.array([np.mean(np.abs(Ea_CUE_out[i][sur[i]] - np.mean(Ea_CUE_out,axis=1)[i])) for i in range (ass)])], axis = 0)
    dEa_sur_ci = np.append(dEa_sur_ci, [np.array([1.96*np.std(np.abs(Ea_CUE_out[i][sur[i]] - np.mean(Ea_CUE_out,axis=1)[i]))/(len(sur[i])**0.5) for i in range (ass)])],axis =0)
    dEa = np.append(dEa, [np.mean([np.abs(Ea_CUE_out[i,:] - np.mean(Ea_CUE_out,axis=1)[i]) for i in range (ass)], axis = 1)], axis= 0)
    
    sur_overlap.append(np.concatenate([overlap[i][sur[i]] for i in range(len(sur))]))
    ext_overlap.append(np.concatenate([overlap[i][ext[i]] for i in range(len(ext))]))
    
    sur_crossf.append(np.concatenate([crossf[i][sur[i]] for i in range(len(sur))]))
    ext_crossf.append(np.concatenate([crossf[i][ext[i]] for i in range(len(ext))]))
    
print((time.time() - start)/60)


# temp_rich = {'rich':rich, 'eq':eq, 'eq_sur':eq_sur, 'eq_sur_ci': eq_sur_ci, 'dEa': dEa, 'dEa_sur': dEa_sur, 'dEa_sur_ci': dEa_sur_ci, \
#      'sU': sU, 'sR': sR, 'eU': eU, 'eR':eR, 'U_var':U_var, 'extU_var':extU_var, 'sur_var':sur_var, 'sur_var_ci':sur_var_ci,\
#      'ext_var': ext_var, 'ext_var_ci':ext_var_ci, 'all_Ea': all_Ea, 'all_Ea_ci': all_Ea_ci, 'sur_CUE':sur_CUE, \
#      'sur_Ea':sur_Ea, 'ext_CUE':ext_CUE, 'ext_Ea':ext_Ea, 'sur_Sr': sur_Sr, 'all_Sr':all_Sr, 'sur_overlap':sur_overlap, \
#      'ext_overlap':ext_overlap, 'sur_crossf':sur_crossf, 'ext_crossf':ext_crossf}
# np.save('../data/temp_rich.npy', temp_rich) 

# temp_rich = np.load('../data/temp_rich.npy',allow_pickle='TRUE').item()

rich_mean = np.nanmean(rich, axis = 1)
rich_ci =  1.96 * np.nanstd(rich,axis = 1)/(ass**0.5)

CUE_mean = np.array([np.mean(sur_CUE[i]) for i in range(T_c)])
CUE_ci = np.array([1.96 * np.std(sur_CUE[i])/(len(ext_CUE[i])**0.5) for i in range(T_c)])
CUE_ext_mean = np.array([np.mean(ext_CUE[i]) for i in range(T_c)])
CUE_ext_ci = np.array([1.96 * np.std(ext_CUE[i])/(len(ext_CUE[i])**0.5) for i in range(T_c)])

Ea_mean = np.array([np.mean(sur_Ea[i]) for i in range(T_c)])
Ea_ci = np.array([1.96 * np.std(sur_Ea[i])/(len(sur_Ea[i])**0.5) for i in range(T_c)])

Sr_mean = np.array([np.mean(sur_Sr[i]) for i in range(T_c)])
Sr_ci = np.array([1.96 * np.std(sur_Sr[i])/(len(sur_Sr[i])**0.5) for i in range(T_c)])

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

surEa_mean = np.array([np.mean(sur_Ea[i]) for i in range(T_c)])
surEa_ci = np.array([1.96 * np.std(sur_Ea[i])/(len(sur_Ea[i])**0.5) for i in range(T_c)])

Ea = []
Ea_e = []
for i in range(T_c): 
    n = 0
    for j in rich[i]:
        j = int(j)
        Ea.append(sur_Ea[i][n:n+j])
        Ea_e.append(ext_Ea[i][n:n+j])
        n = n + j

A = list(zip(rich.flatten(),Ea))
B = list(zip(rich.flatten(),Ea_e))
A.sort(key = lambda x: x[0])
B.sort(key = lambda x: x[0])

rich_v = []
for x in rich.flatten():
    if x not in rich_v:
        rich_v.append(x)
rich_v = np.sort(rich_v)


Ea_sorted = []
Ea_e_sorted = []
for i in rich_v:
    sorting = np.empty((0))
    sorting_e = np.empty((0))
    for j in range(len(A)):
        if [x[0] for x in A][j] == i:
            sorting = np.append(sorting, A[j][1])
            sorting_e = np.append(sorting_e, B[j][1])
    Ea_sorted.append(sorting)
    Ea_e_sorted.append(sorting_e)

    
meanEa = []
ciEa = []
meanEa_e = []
ciEa_e = []
for i in range(len(Ea_sorted)):
    meanEa.append(np.mean(Ea_sorted[i]))
    ciEa.append(1.96 * np.std(Ea_sorted[i])/(len(Ea_sorted[i])**0.5))
    meanEa_e.append(np.mean(Ea_e_sorted[i][np.where(Ea_e_sorted[i]>-500)]))
    ciEa_e.append(1.96 * np.std(Ea_e_sorted[i][np.where(Ea_e_sorted[i]>-500)])/(len(Ea_e_sorted[i][np.where(Ea_e_sorted[i]>-500)])**0.5))

T_plot = range(0, T_c, 1)
T_sur = [[np.repeat(T_plot[i], rich[i][j]) for j in range(ass)] for i in range(T_c)]


# plt.rcParams["figure.figsize"] = (15,9)
# # plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
# plt.rcParams.update({'font.size': 25})
# plt.rc('xtick', labelsize=25) 
# plt.rc('ytick', labelsize=25) 
# # plt.rcParams.update(mpl.rcParamsDefault)

plt.rcParams["figure.figsize"] = (12,9)
# plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
plt.rcParams.update({'font.size': 30})
plt.rc('xtick', labelsize=30) 
plt.rc('ytick', labelsize=30) 
# plt.rcParams.update(mpl.rcParamsDefault)

plt.plot(T_plot, rich_mean)
plt.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, alpha=0.1,linewidth=2.5)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Richness')
plt.text(-5,14,'A',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/selectingEaCUE.png')
plt. show()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ln1 = ax1.plot(T_plot, sU_mean,'darkorange', label = "Survivor Uptake Rate",linewidth=2.5)
ax1.fill_between(T_plot, sU_mean - sU_ci, sU_mean + sU_ci, color='darkorange', alpha=.1)
ln3 = ax1.plot(T_plot, eU_mean,'tan', label = "Extinct Uptake Rate", alpha=.7,linewidth=2.5)
ax1.fill_between(T_plot, eU_mean - eU_ci, eU_mean + eU_ci, color='tan', alpha=.3)
ln2 = ax2.plot(T_plot, sR_mean, 'g', label = "Survivor Respiration Rate",linewidth=2.5)
ax2.fill_between(T_plot, sR_mean - sR_ci, sR_mean + sR_ci, color='g', alpha=.1)
ln4 = ax2.plot(T_plot, eR_mean, 'darkseagreen', label = "Extinct Respiration Rate", alpha=.5,linewidth=2.5)
ax2.fill_between(T_plot, eR_mean - eR_ci, eR_mean + eR_ci, color='darkseagreen', alpha=.1)
ax1.set_xlabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Uptake Rate (1/Time)')
ax2.set_ylabel('Repiration Rate (Mass/Volume*Time)')
lns = ln1+ln2+ln3+ln4
ax1.legend(lns, [i.get_label() for i in lns], loc = 2)
# text(-5,400,'B',fontsize= 'x-large')
ax1.text(-5,900,'A',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/selectingEaCUE_1.png')
plt. show()

plt.plot(T_plot, CUE_mean, 'r', label = "Survivor",linewidth=2.5)
plt.fill_between(T_plot, CUE_mean - CUE_ci, CUE_mean + CUE_ci, color='r', alpha=.1)
plt.plot(T_plot, CUE_ext_mean, 'dimgrey', label = "Extinct",linewidth=2.5)
plt.fill_between(T_plot,  CUE_ext_mean - CUE_ext_ci,  CUE_ext_mean + CUE_ext_ci, color='dimgrey', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('CUE')
plt.legend()
plt.text(-5,0.6,'B',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/selectingEaCUE_2.png')
plt. show()


# n = 0
# for i in range(T_c):
#     plt.scatter(ext_Ea[i], np.concatenate(rich_ext[i]), color = 'lightgrey', alpha = 0.5, s = 50)
#     if n == 0: 
#         plt.scatter(sur_Ea[i], np.concatenate(rich_sur[i]), color = 'b', alpha = 0.5, s =50, label = 'Survivors')
#         n = 1
#     else:
#         plt.scatter(sur_Ea[i], np.concatenate(rich_sur[i]), color = 'b', alpha = 0.5, s = 50)
# m, b = np.polyfit(np.delete(meanEa,np.where(np.isnan(meanEa))), np.delete(rich_v,np.where(np.isnan(meanEa))), 1)
# # x = np.arange(0.13, 1.5, 0.01)

plt.scatter(meanEa, rich_v, color = 'b', alpha = 0.7, s = 50, label = 'Survivors')
plt.errorbar(meanEa, rich_v, xerr=ciEa, fmt=',', color = 'b', alpha = 0.7)
plt.scatter(meanEa_e, rich_v, color = 'grey', alpha = 0.7, s = 50, label = 'Extinct')
plt.errorbar(meanEa_e, rich_v, xerr=ciEa_e, fmt=',', color = 'grey', alpha = 0.7)

m, b, r_value, p_value, std_err = sc.stats.linregress(meanEa, rich_v)
print(r_value**2)
x = np.arange((np.max(rich_v)-b)/m, (np.min(rich_v)-b)/m, 0.01)
plt.plot(x, m*x + b, color = 'k',linewidth=3)
plt.text(0.7, 12.5, '$R^2 = $%s' %np.round(r_value**2, 3), fontsize = 22)
plt.xlabel('Thermal Sensitivity of CUE')
plt.ylabel('Richness')
plt.legend()
plt.text(-0.35,20,'B',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/EaCUE_richness.png')
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
plt.legend(fontsize = 'x-small', framealpha = 0.4)
# plt.tight_layout()
plt.text(-6,0.4,'A',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/Resource_overlap.png')
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

plt.plot(T_plot, np.log(U_var_mean),'darkorange', label = "Survivor",linewidth=2.5)
plt.fill_between(T_plot, np.log(U_var_mean-U_var_ci), np.log(U_var_mean+U_var_ci), color='darkorange', alpha=.2)
plt.plot(T_plot, np.log(extU_var_mean),'grey', label = "Extinct")
plt.fill_between(T_plot, np.log(extU_var_mean-extU_var_ci), np.log(extU_var_mean+extU_var_ci), color='grey', alpha=.2)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Variance of Uptake (log)')
plt.legend(loc = 2,fontsize = 'x-small')
# plt.tight_layout()
plt.text(-6,7.2,'B',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/VarU.png')
plt.show()

plt.scatter(Sr_mean, rich_mean, s=20, color = 'k', alpha = 0.7)
plt.errorbar(Sr_mean, rich_mean, xerr=Sr_ci, fmt=',', color = 'k', alpha = 0.7)
plt.xlabel('Community average S*')
plt.ylabel('Richness')
plt.show()


eq_mean = np.nanmean(eq_sur, axis = 1)
eq_ci =  1.96 * np.nanstd(eq_sur,axis = 1)/(ass**0.5)


plt.plot(T_plot, eq_mean, 'r', label = 'Competitive exclusion', alpha = 0.7)
plt.fill_between(T_plot, eq_mean - eq_ci, eq_mean + eq_ci, color='r', alpha=.1)
# plt.plot(T_plot, (1 - overlap_sur_mean)/99, 'r', label = 'Resource partitioning', alpha = 0.7)
# plt.fill_between(T_plot, (1 - overlap_sur_mean)/99 - overlap_sur_ci, (1 - overlap_sur_mean)/99 + overlap_sur_ci, color = 'r', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('|Survivor $S_i^* - \overline{S^*}$|')
# plt.legend()
# plt.savefig('../thesis/Figures/eq_st.png')
plt.show()

dEa_mean = np.nanmean(dEa_sur, axis = 1)
dEa_ci = 1.96 * np.nanstd(dEa_sur,axis = 1)/(ass**0.5)

m, b, r_value, p_value, std_err = sc.stats.linregress(dEa_mean, rich_mean)
print(r_value**2)
x = np.arange(np.min(dEa_mean), np.max(dEa_mean), 0.001)
plt.plot(x, m*x + b, color = 'k',linewidth=1, alpha = 0.7)
plt.scatter(dEa_mean, rich_mean, s=20, color = 'k', alpha = 0.7)
plt.errorbar(dEa_mean, rich_mean, xerr=dEa_ci, fmt=',', color = 'k', alpha = 0.3)
plt.errorbar(dEa_mean, rich_mean, yerr=rich_ci, fmt=',', color = 'k', alpha = 0.3)
plt.xlabel('Ea Difference')
plt.ylabel('Richness')
plt.show()

m, b, r_value, p_value, std_err = sc.stats.linregress(surEa_mean, rich_mean)
print(r_value**2)
x = np.arange(np.min(surEa_mean), np.max(surEa_mean), 0.001)
plt.plot(x, m*x + b, color = 'k',linewidth=1, alpha = 0.7)
plt.scatter(surEa_mean, rich_mean, s=20, color = 'k', alpha = 0.7)
plt.errorbar(surEa_mean, rich_mean, xerr=surEa_ci, fmt=',', color = 'k', alpha = 0.3)
plt.errorbar(surEa_mean, rich_mean, yerr=rich_ci, fmt=',', color = 'k', alpha = 0.3)
plt.xlabel('Average survivor $E_{a_{CUE}}$')
plt.ylabel('Richness')
plt.text(1, 13, '$R^2 = $%s' %np.round(r_value**2, 3))
plt.show()


plt.scatter(eq_mean, rich_mean, s=20, color = 'k', alpha = 0.7)
# plt.errorbar(dEa_sur, rich_mean, xerr=dEa_sur_ci, fmt=',', color = 'k', alpha = 0.7)
plt.xlabel('Fitness Difference')
plt.ylabel('Richness')
plt.show()

# fig, ax1 = plt.subplots()
# ax2 = ax1.twiny()
# ax1.scatter(dEa_sur, rich_mean, s=20, color = 'b', alpha = 0.7, label = 'Ea difference')
# ax1.set_xlabel('Ea Difference')
# ax2.scatter(eq_sur, rich_mean, s=20, color = 'k', alpha = 0.7, label = 'Fitness difference')
# ax2.set_xlabel('Fitness Difference')
# ax1.set_ylabel('Richness')
# plt.show()

plt.scatter(surEa_mean, Sr_mean, s=100, color = 'k', alpha = 0.7)
plt.errorbar(surEa_mean, Sr_mean, xerr=surEa_ci, fmt=',', color = 'k', alpha = 0.3)
plt.errorbar(surEa_mean, Sr_mean, yerr=Sr_ci, fmt=',', color = 'k', alpha = 0.3)
plt.xlabel('Average survivor $E_{a_{CUE}}$')
plt.ylabel('Average survivor $S^*$')
plt.text(-0.1,0.57,'A',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/Eadiff_fitdiff.png')
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
        result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
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
plt.text(1,11.5,'A',fontsize= 'x-large')
plt.savefig('../data/leakage_rich.png')
plt.show()


Xlabel = np.round(np.arange(5, 5 * T_c, 5), 1)
Ylabel = np.array((0.7, 0.5, 0.3, 0.1, 0))
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
plt.text(-1.3,-0.9,'B',fontsize= 'x-large')
plt.savefig('../data/cf_lf.png')
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
plt.show()

a,b = np.polyfit(Sr_mean, np.log(rich_mean), 1)
x = np.arange(np.min(Sr_mean), p.max(Sr_mean), 0.001)
y = np.exp(b) * np.exp(a*x)
plt.plot(x, y, color = 'k', alpha = 0.7)

# model = np.poly1d(np.polyfit(Sr_mean, rich_mean, 5))
# plt.plot(x, model(x), color = 'k', alpha = 0.7)

plt.show()


###################################################################################################################
import numpy as np
from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt

######### Main Code ###########

######## Set up parameters ###########

N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 0 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
# lf = 0 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly number, i.e. how many times the system can assemble
t_fin = 4000 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 5 # Half saturation constant


rich_mean = np.empty((0,5))
rich_ci = np.empty((0,5))

for i in range(5):
    if i == 0:
        lf = 0
    else:
        lf = i*0.2-0.1
    all_cf = []
    rich = np.empty((0,ass))
    for i in range(5):
        strc = 10**(i-2)
        result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K, strc)
        rich = np.append(rich, [rich_series.flatten()], axis = 0)
    rich_mean = np.append(rich_mean, [np.mean(rich, axis = 1)], axis = 0)
    rich_ci = np.append(rich_ci, [1.96 * np.std(rich, axis = 1)/(ass**0.5)])

Xlabel = np.array((0.01, 0.1, 1, 10, 100))
Ylabel = np.array((0.7, 0.5, 0.3, 0.1, 0))
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yticks(range(len(Ylabel)))
ax.set_xticks(range(len(Xlabel)))
ax.set_yticklabels(Ylabel)
ax.set_xticklabels(Xlabel)
plt.xlabel("Alpha value for Dirichlet distribution")
plt.ylabel("Leakage")
im = plt.imshow(np.flip(rich_mean,axis = 0), cmap = 'binary')
plt.colorbar(im)
plt.show()

#####################################################################################################
import time
start = time.time()

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 0 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin, 10 degrees C !!!
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy
Ea_CUE = 0.3
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 30 # Assembly times at each temperature
tv = 1 # immigration times inside one assembly
t_fin = 4000 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)

########## Running Model ###########
result_array_1, rich_series_1, l_1, U_out_total_1, R_out_1, CUE_out_1, Ea_CUE_out_1, overlap_1, crossf_1, Sr_1 = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
print(np.mean(rich_series_1))

rich_1 = np.array([len(np.where(result_array_1[i,0:N])[0]) for i in range(len(result_array_1))]).reshape(ass,t_fin)
rich_mean_1 = np.mean(rich_1, axis = 0)
rich_ci_1 = 1.96 * np.std(rich_1,axis = 0)/(ass**0.5)
t_plot = np.linspace(0,t_fin,t_fin)

plt.plot(t_plot[0:2000], rich_mean_1[0:2000])
plt.fill_between(t_plot[0:2000], rich_mean_1[0:2000] - rich_ci_1[0:2000], rich_mean_1[0:2000] + rich_ci_1[0:2000], alpha=.2)
plt.ylabel('Richness')
plt.xlabel('Time')
# plt.xlim([0,2000])
plt.text(-400,110,'B',fontsize= 'x-large')
plt.savefig('../thesis/Figures/rich_decay.png')
plt.show()


print((time.time() - start)/60)

##########################################################################################
x = np.arange(0, 2, 0.1)
plt.plot(x, 4.47*x, linewidth = 1.2, label = 'Linear function')
for i in range(6):
    K = 2**(i - 3) # Half saturation constant for Monod equation(Type II)
    plt.plot(x, 4.47*x/(K+x), linewidth = 1.2, label = 'K = %s' %K)
plt.legend(framealpha = 0.3)
plt.ylabel('Uptake rate (Mass/Volume*Time)')
plt.xlabel('Resource concentration (Mass/Volume)')
plt.text(-0.3, 9.5,'A',fontsize= 'x-large')
plt.savefig('../data/K.png')
plt.show()

N = 1 # Number of consumers
M = 1 # Number of resources

# Temperature params
T = 273.15 + 0 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin, 10 degrees C !!!
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy
# Ea_CUE = 0.3
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 1 # Assembly times at each temperature
tv = 1 # immigration times inside one assembly
t_fin = 200 # Number of time steps for each temperature


########## Running Model ###########

t_plot = np.linspace(0,t_fin,t_fin)

typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant for Monod equation(Type II)
np.random.seed(0)
result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
plt.plot(t_plot, result_array[:,N:N+M], linewidth= 1.2, label = 'Linear function')
typ = 2 # Functional response, Type I or II
for i in range(6):
    K = 2**(i - 3) # Half saturation constant for Monod equation(Type II)
    np.random.seed(0)
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
    plt.plot(t_plot, result_array[:,N:N+M], linewidth= 1.2, label = 'K = %s' %K)
plt.legend(framealpha = 0.3, loc = 1)
plt.ylabel('Resources concentration (Mass/Volume)')
plt.xlabel('Time')
plt.text(-35, 20,'C',fontsize= 'x-large')
plt.savefig('../data/K_res.png')
plt.show()

typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant for Monod equation(Type II)
np.random.seed(0)
result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
plt.plot(t_plot, result_array[:,0:N], linewidth= 1.2, label = 'Linear function')
typ = 2 # Functional response, Type I or II
for i in range(6):
    K = 2**(i - 3) # Half saturation constant for Monod equation(Type II)
    np.random.seed(0)
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
    plt.plot(t_plot, result_array[:,0:N], linewidth= 1.2, label = 'K = %s' %K)
plt.legend(framealpha = 0.3, loc = 1)
plt.ylabel('Consumer biomass (Mass/Volume)')
plt.xlabel('Time')
plt.text(-35, 3.2,'D',fontsize= 'x-large')
plt.savefig('../data/K_con.png')
plt.show()

#############################################################################################################

import time
start = time.time()

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
T = 273.15 + 0 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin, 10 degrees C !!!
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy
# Ea_CUE = 0.3
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 1 # Assembly times at each temperature
tv = 1 # immigration times inside one assembly
t_fin = 1000 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)

########## Running Model ###########
np.random.seed(1)
result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
print(np.mean(rich_series))

print((time.time() - start)/60)

t_plot = np.linspace(0,t_fin,t_fin)

plt.plot(t_plot, result_array[:,N:N+M], 'b-', linewidth=1)
plt.ylabel('Resources concentration')
plt.xlabel('Time')
# plt.title('Consumer-Resource population dynamics')
# plt.text(-150,3.2,'A',fontsize= 'x-large')
plt.savefig('../../pre/resource_con_example.png')
plt.show()

plt.plot(t_plot, result_array[:,0:N], 'g-', linewidth=1)
plt.ylabel('Consumer biomass')
plt.xlabel('Time')
# plt.text(-150, 4.2,'B',fontsize= 'x-large')
plt.savefig('../../pre/consumer_con_example.png')
plt.show()

###########################################################################################
import time
start = time.time()

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 2 # Number of consumers
M = 1 # Number of resources

# Temperature params
T = 273.15 + 25 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin, 10 degrees C !!!
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
ass = 1 # Assembly times at each temperature
tv = 1 # immigration times inside one assembly
t_fin = 10 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)

########## Running Model ###########
result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
print(np.mean(rich_series))

from matplotlib.lines import Line2D

rich = np.array([len(np.where(result_array[i,0:N])[0]) for i in range(len(result_array))]).reshape(ass,t_fin)
rich_mean = np.mean(rich, axis = 0)
rich_ci = 1.96 * np.std(rich,axis = 0)/(ass**0.5)

t_plot = np.linspace(0,len(result_array),len(result_array))

plt.plot(t_plot, result_array[:,N:N+M], 'b-', linewidth=2.5, label = "Resources")
plt.plot(t_plot, result_array[:,0:N], 'g-', linewidth=2.5, label = "Consumers")
plt.ylabel('Consumer & Resource Concentration')
plt.xlabel('Time')
# plt.title('Consumer-Resource Concentration Dynamics')
plt.legend([Line2D([0], [0], color='green', lw=2), Line2D([0], [0], color='blue', lw=2)], ['Consumers', 'Resources'])
plt.text(-1.75,1.13,'B',fontsize= 'x-large')
plt.savefig('../data/2_consumer_25C.png')
plt.show()

######################################################################################################
import numpy as np
import matplotlib.pyplot as plt

N = 100
k = 0.0000862 # Boltzman constant
Tref = 273.15 + 0 # Reference temperature Kelvin, 0 degrees C
T = 273.15 + np.linspace(0,60,61) # Temperatures
Ea_D = 3.5
lf = 0.4

B_R = 1.70 * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref))) # Using CUE0 = 0.22, mean growth rate = 0.48
B_U = (1.70/(1 - lf - 0.22)) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref)))


np.random.seed(11) # 11 19 20
T_pk_U = 273.15 + np.random.normal(35, 5, size = N)
T_pk_R = T_pk_U + 3
a = 15
Ea_U = np.random.beta(a, ((a - 1/3) / (0.82/4)) + 2/3 - a, N)*4
Ea_R = np.random.beta(a, ((a - 1/3) / (0.67/4)) + 2/3 - a, N)*4

for i in range(N):
    U_Sharpe = B_U * np.exp((-Ea_U[i]/k) * ((1/T)-(1/Tref)))/(1 + (Ea_U[i]/(Ea_D - Ea_U[i])) * np.exp(Ea_D/k * (1/T_pk_U[i] - 1/T))) 
    plt.plot(T - 273.15, U_Sharpe, color = 'darkorange')

plt.annotate(s='', xy=(30,1500), xytext=(0,1500), arrowprops=dict(arrowstyle='<->',facecolor='black', lw=2))
plt.axvline(x = 30, c = 'k', ls = '--')
plt.text(10,1350,'OTR')
plt.xlabel('Temperature ($^\circ$C)') 
plt.ylabel('Uptake Rate')
# plt.text(-12,1650,'B',fontsize= 'x-large')
plt.savefig('../../pre/simulated_uptake_rate.png')
plt.show()

for i in range(N):
    R_Sharpe = B_R * np.exp((-Ea_R[i]/k) * ((1/T)-(1/Tref)))/(1 + (Ea_R[i]/(Ea_D - Ea_R[i])) * np.exp(Ea_D/k * (1/T_pk_R[i] - 1/T))) 
    plt.plot(T - 273.15, R_Sharpe, 'darkgreen')

plt.annotate(s='', xy=(30,400), xytext=(0,400), arrowprops=dict(arrowstyle='<->',facecolor='black', lw=2))
plt.axvline(x = 30, c = 'k', ls = '--')
plt.text(10,360,'OTR')
plt.xlabel('Temperature ($^\circ$C)') 
plt.ylabel('Repiration Rate')
# plt.text(-12,460,'C',fontsize= 'x-large')
plt.savefig('../../pre/simulated_res_rate.png')
plt.show()


##################################################################################################
import numpy as np
import matplotlib.pyplot as plt
N = 50000
T_pk_U = 273.15 + np.random.normal(35, 5, size = N)
T_pk_R = T_pk_U + 3
plt.hist(T_pk_U, 30, color = "orangered", density = True, alpha = 0.5, label = "Uptake rate")
plt.hist(T_pk_R, 30, color = "g", density = True, alpha = 0.5, label = "Respiration rate")
plt.axvline(273.15+35, color="orangered", linestyle='dashed', linewidth = 3) # Mean
plt.axvline(273.15+38, color='darkgreen', linestyle='dashed', linewidth = 3) # Mean
plt.xlabel("Peak Temperature")
plt.ylabel("Density")
plt.legend(framealpha = 0.3)
plt.text(275,0.09,'A',fontsize= 'x-large')
plt.savefig('../data/SI_UR.png')
plt.show()

a = 20 # The alpha value for beta distribution in Ea
Ea_U = np.random.beta(a, ((a - 1/3) / (0.82/4)) + 2/3 - a, N)*4
Ea_R = np.random.beta(a, ((a - 1/3) / (0.67/4)) + 2/3 - a, N)*4

plt.hist(Ea_U, 30, color = "orangered", density = True, alpha = 0.5, label = "Uptake rate")
plt.hist(Ea_R, 30, color = "g", density = True, alpha = 0.5, label = "Respiration rate")
plt.axvline(0.82, color="orangered", linestyle='dashed', linewidth = 3) # Median
plt.axvline(0.67, color='darkgreen', linestyle='dashed', linewidth = 3) # Median
plt.xlabel("Activation Energy (Ea)")
plt.ylabel("Density")
plt.legend(framealpha = 0.3)
plt.text(-0.1,3.4,'B',fontsize= 'x-large')
plt.savefig('../data/SI_UR_1.png')
plt.show()

####################################################################################################
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
tv = 1 # immigration times inside one assembly
t_fin = 4000 # Number of time steps for each temperature
T_c = 26 # How many temperatures to cover (how many cycles to run)


T_plot = range(0, T_c, 1)

typ = 1 # Functional response, Type I or II
K = 0
rich = np.empty((0, ass))
for i in range(T_c):
    T = 273.15 + i # Temperature
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
    rich = np.append(rich, [rich_series.flatten()], axis = 0)

rich_mean = np.mean(rich, axis = 1)
rich_ci =  1.96 * np.std(rich,axis = 1)/(ass**0.5)
plt.plot(T_plot, rich_mean, linewidth= 1, label = 'Linear function')
plt.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, alpha=.1)


typ = 2 # Functional response, Type I or II
for j in range(6):
    K = 2**(j - 3) # Half saturation constant for Monod equation(Type II)
    rich = np.empty((0, ass))

    for i in range(T_c):
        T = 273.15 + i # Temperature
        result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K)
        rich = np.append(rich, [rich_series.flatten()], axis = 0)

    rich_mean = np.mean(rich, axis = 1)
    rich_ci =  1.96 * np.std(rich,axis = 1)/(ass**0.5)
    plt.plot(T_plot, rich_mean, label = 'K = %s' %K)
    plt.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, alpha=.1)

plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Richness')
plt.legend()
plt.show()
    
print((time.time() - start)/60)


##################################################################################################################
import numpy as np
import matplotlib.pyplot as plt

k = 0.0000862 # Boltzman constant
Tref = 273.15 + 0 # Reference temperature Kelvin, 0 degrees C
lf = 0.4
T = 273.15 + np.linspace(0,30,100) # Temperatures
T_pk_U = 273.15 + 35
T_pk_R = T_pk_U + 3

Ea_D = 3.5 # Deactivation energy
Ea_CUE = 0.3
B0_CUE = 0.1 * np.exp((-Ea_CUE/k) * ((1/Tref)-(1/273.15)))
B_U = 4.47
B_R = 1.70
Ea_U = 0.82 # Ea for uptake
Ea_R = 0.67
B_R*(Ea_U-Ea_R)/(B_U*(1-lf)-B_R)
U_Sharpe = B_U * np.exp((-Ea_U/k) * ((1/T)-(1/Tref)))/(1 + (Ea_U/(Ea_D - Ea_U)) * np.exp(Ea_D/k * (1/T_pk_U - 1/T))) 
R_Sharpe = B_R * np.exp((-Ea_R/k) * ((1/T)-(1/Tref)))/(1 + (Ea_R/(Ea_D - Ea_R)) * np.exp(Ea_D/k * (1/T_pk_R - 1/T)))
CUE_Sharpe = 0.22* np.exp((-B_R*(Ea_U-Ea_R)/(B_U*(1-lf)-B_R)/k) * ((1/T)-(1/Tref)))

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ln1 = ax1.plot(T - 273.15, U_Sharpe, 'darkorange', linewidth=2, label = "Uptake Rate")
ln2 = ax1.plot(T - 273.15, R_Sharpe, 'darkgreen', linewidth=2, label = "Respiration Rate")
ln3 = ax2.plot(T - 273.15, CUE_Sharpe, 'r', linewidth=2, label = 'CUE')
ax1.set_xlabel('Temperature ($^\circ$C)', fontsize = 30)
ax1.set_ylabel('Uptake & Respiration Rates', fontsize = 30)
ax2.set_ylabel('CUE', fontsize = 30)
# plt.title('Modified Sharpe-Schoolfield Temperature Performance Curve')
lns = ln1+ln2+ln3
ax1.legend(lns, [i.get_label() for i in lns], loc = 2, fontsize = 25)
ax1.annotate(s='', xy=(9,35), xytext=(14,56), arrowprops=dict(arrowstyle='<-', lw=2))
ax1.text(11,62,'$E_{a_{CUE}}$')
plt.tight_layout()
plt.savefig('../../pre/Figures/URCUE.png')
plt.show()


###############################################################################################################################
import numpy as np

M = 5
lf = 0.4
alpha = 100
# l_raw = np.array([[np.random.normal(lf/3,0.005) if i >= 3 and i > M-3 else np.random.normal(1/(i),0.005)* lf for i in range(M,0,-1)] for i in range(1,M+1)])
# l = [[l_raw[j,i] if j>=i and j-i <3 else 0 for j in range(M)] for i in range(M)]
l = np.stack([np.concatenate((np.zeros(i),(np.random.dirichlet(np.full(3,alpha),1)*lf).flatten(), np.zeros(M-3-i))) if i <= M-3 else np.concatenate((np.zeros(i),(np.random.dirichlet(np.full(M-i,alpha),1)*lf).flatten())) for i in range(M)])

im = plt.imshow(l, cmap = 'binary')
ax = plt.gca()
ax.set_xticks(np.arange(-.5,5,1), minor = True)
ax.set_yticks(np.arange(-.5,5,1), minor = True)
plt.xticks([])
plt.yticks([])
ax.grid(which = 'minor', color = 'k', linestyle='-', linewidth=1)
plt.title( "" )
plt.grid()
plt.show()

np.round(l,3)

import parameters as par
import size_temp_funcs as st
import numpy as np
import matplotlib.pyplot as plt


N = 5 # Number of consumers
M = 5 # Number of resources

# Temperature params
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
pk_U =  273.15 + np.random.normal(35, 5, size = N)
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
k = 0.0000862 # Boltzman constant
T_pk = Tref + pk_U # Peak above Tref, Kelvin
T = 273.15 + 25 # Temperature
Ea_U = np.random.beta(4, ((4 - 1/3) / 0.82) + 2/3 - 4, N) # Ea for uptake
B_U0 = (1.70/(1 - lf - 0.22)) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref)))
B_U = np.random.normal(B_U0, 0.1*B_U0, N) # Adding variation into B0

U_sum = st.temp_growth(k, T, Tref, T_pk, N, B_U, Ma, Ea_U, Ea_D)
diri = np.transpose(np.random.dirichlet(np.ones(M),N))
U = np.transpose(diri*U_sum)
print(U_sum)
print(U)

Xlabel = ''
Ylabel = ''
fig = plt.figure()
ax = fig.add_subplot(111)
im = plt.imshow(U, cmap = 'binary')
ax = plt.gca()
ax.set_xticks(np.arange(-.5,5,1), minor = True)
ax.set_yticks(np.arange(-.5,5,1), minor = True)
plt.xticks([])
plt.yticks([])
ax.grid(which = 'minor', color = 'k', linestyle='-', linewidth=1)
plt.title( "" )
plt.show()

plt.rcParams["figure.figsize"] = (15,9)

T = 273.15 + np.linspace(0,50,51) # Temperatures
for i in range(N):
    U_Sharpe = B_U[i] * np.exp((-Ea_U[i]/k) * ((1/T)-(1/Tref)))/(1 + (Ea_U[i]/(Ea_D - Ea_U[i])) * np.exp(Ea_D/k * (1/pk_U[i] - 1/T))) 
    plt.plot(T - 273.15, U_Sharpe, 'darkgreen', linewidth=4)
    plt.axvline(25, linewidth = 2.5, color = 'darkblue')
    plt.yticks([])
    plt.ylim(0,165)
    plt.xticks([])
    plt.show()




plt.hist(Ea_CUE_out.flatten()[(Ea_CUE_out.flatten()>-10) & (Ea_CUE_out.flatten()<10)],100)
plt.xlabel("$E_{a_\u03B5}$")
plt.ylabel("Density")
plt.savefig('../thesis/Figures/d_EaCUE.png')
plt.show()
