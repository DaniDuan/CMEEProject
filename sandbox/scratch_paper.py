import numpy as np
import matplotlib.pylab as plt


s = np.random.dirichlet((10, 5, 3, 1), 10).transpose()
s.shape
np.mean(s, axis = 0)
np.sum(s, axis = 0)
np.mean(s, axis = 1)
np.sum(s, axis = 1)


w = np.random.dirichlet()


plt.barh(range(20), s[0])
plt.barh(range(20), s[1], left=s[0], color='g')
plt.barh(range(20), s[2], left=s[0]+s[1], color='r')

f1 = np.random.dirichlet(np.full(10,1), 100000) 
plt.hist(plt.hist(f1, 30, density = True)[0])
plt.show()

plt.hist(np.random.beta(1, 1, 50000)*4, 30, density = True)
plt.show()

w = np.random.dirichlet(np.ones(M),N)

A = np.transpose(np.array([[1,2,3],[3,4,5],[5,6,7],[10,20,30]]))
B = np.transpose(np.array([1,2,3]))
A*B


    # U_sum = st.temp_growth(k, T, Tref, T_pk, N, B_U, Ma, Ea_U, Ea_D)
    # u_1 = np.empty((0,N)) 
    # for i in range(M-1):
    #     mean = U_sum[i]/N
    #     random.seed(i)
    #     a = np.array([np.random.uniform(0, mean, size = N)])
    #     u_1 = np.append(u_1,a,axis = 0)
    # u_2 = np.array([U_sum[0:N] - np.sum(u_1, axis = 0)])
    # # b = np.where(u_2<0)
    # # u_1[:,b] = u_1[:,b] + u_2[b]/N-1
    # # u2 = np.where(u2<0, 0, u2)
    # U_raw = np.append(u_1, u_2, axis = 0)
    # for i in range(U_raw.shape[1]): 
    #     random.seed(i)
    #     np.random.shuffle(U_raw[:,i])
    # U = np.transpose(U_raw)

random.seed(0)
l_raw = np.array([[np.random.normal(1/(i-1),0.005)* 0.4 if i-1>0 else np.random.normal(0.4, 0.005) for i in range(M,0,-1)] for i in range(1,M+1)])
l = np.transpose(l_raw) * fix


A = np.array([1,2,2,3,3,4,5,6,7,8])
B = [A[i*2] for i in range(5)]
B = np.array([[5,10],[10,20],[30,40]])
A*B
C = np.empty((0))
for i in range(5):
    C = np.append(C,i)] for i in range(5)

a = np.array(np.where(B>10))
len(B[np.where(B==10)])

rich = np.append(rich, N - len(rem_find[ext]))


import parameters as par
import size_temp_funcs as st
import numpy as np

N = 10 # Number of consumers
M = 5 # Number of resources

# Temperature params
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
np.random.seed(0)
pk_U = np.random.normal(25, 3, size = N)
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
k = 0.0000862 # Boltzman constant
T_pk = Tref + pk_U # Peak above Tref, Kelvin
T = 273.15 + 21 # Model temperature
Ea_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:N] # Ea for uptake
Ea_R = Ea_U - 0.8 # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
B_U = (10**(2.84 + (-4.96 * Ea_U))) + 4 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration

U_sum = st.temp_growth(k, T, Tref, T_pk, N, B_U, Ma, Ea_U, Ea_D)
# np.random.seed(0)
diri = np.transpose(np.random.dirichlet(np.full(M,1/M),N))
U = np.transpose(diri*U_sum)
np.round(U)
print(U_sum)
print(U)


for i in range(5):
    print(i)


A = np.array([[1,2,2,np.nan,np.nan,np.nan]])
A = np.nan_to_num(A, nan=0)
np.mean(A)


np.random.seed(0)
A = np.random.normal(1,1,10)
B = np.random.normal(6,6,6)
A
B

C = np.einsum('ij,kj->ik', SL, U) - R
dCdt = xc * C
np.sum(dCdt,axis = 0)/np.sum(xc*np.einsum('ij,kj->ik', xr, U),axis = 0)



from matplotlib import pyplot as plt
import numpy as np

#some example data
x= np.linspace(0.1, 9.9, 20)
y = 3.0 * x
#some confidence interval
ci = 1.96 * np.std(A,axis = 1)/np.mean(A,axis = 1)

fig, ax = plt.subplots()
ax.plot(t_rich,A)
ax.fill_between(t_rich, (A-ci), (A+ci), color='b', alpha=.1)
plt.show()



a = 10 #per group
b = 5 #group
c = 8
A = np.arange(a*b*c).reshape(a*b,c)

np.arange(a*b, step = a)+i
sub = np.array([np.sum(A[np.arange(a*b, step = a)+i, :], axis =0) for i in range(b)])/b
np.mean(A[[np.arange(a*b, step = a)+i for i in range(a)], :], axis = 1)
A = CUE_series
a = tv
b = ass

A = np.array([[1,2,3],[2,3,1],[-1,5,0]])
B = A < 0
A[np.where(~B.any(axis = 1))]
[i for i in range(len(A)-1) if any(A[i+1,]-A[i,] > 0)]

A[np.where(A[2,]>0)]
A[2,]>0

[i for i in range(len(xc)-1) if np.sum(xc[i+1,])-np.sum(xc[i,])!=0]

B = np.random.normal(6,6,6)



(np.sum(np.einsum('ij,kj->ik', SL[2:5], U), axis = 0) - R * len(range(2,5)))*xc[2:5]
pops[5,0:N]-pops[2,0:N]

            CUE = 1 - xc*R/(dCdt + xc*R)


A = np.empty((0,N,tv*ass))
A = np.append(A, [[1,2],[1,2],[1,2]])


np.random.randint(tv, size = 1)


a = np.array([[1,2],[2,2],[3,2]])
b = np.array([[2],[3],[4]])
np.dstack((a,b))

U_out = np.empty((0))
U_out = np.array([U_out,U[np.where(rem_find > 0.01)[0]]])
U_out = np.append(U_out, U[np.where(rem_find > 0.01)[0]], axis = 2)

a = np.ones((3,))
b = np.ones((2,))
c = np.array([a, b])


func = [np.where(U_out[i,:] == np.max(U_out, axis = 1)[i]) for i in range(U_out.shape[0])]
func_ab

k = 0.0000862 # Boltzman constant
T = np.linspace(0,50,51) # Temperatures
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
T_pk_U = Tref + np.random.normal(32, 5, size = 1) 
T_pk_R = T_pk_U + 3
Ea_U = np.round(np.random.normal(1.5, 0.01, 1),3) # Ea for uptake
Ea_R = Ea_U - 0.6 # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
B_U = (10**(2.84 + (-4.96 * Ea_U))) + 4 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration
Ea_D = 3.5 # Deactivation energy


U_Sharpe = B_U * np.exp((-Ea_U/k) * ((1/T)-(1/Tref)))/(1 + (Ea_U/(Ea_D - Ea_U)) * np.exp(Ea_D/k * (1/T_pk_U - 1/T))) 
R_Sharpe = B_R * np.exp((-Ea_R/k) * ((1/T)-(1/Tref)))/(1 + (Ea_R/(Ea_D - Ea_R)) * np.exp(Ea_D/k * (1/T_pk_R - 1/T)))


B = np.array([[0,1,2,4],[0,2,2,4]])
np.sum(np.min(B,axis = 0))/np.sum(np.max(B,axis = 0))

np.ceil(130/100)*100


a = np.array([1, 2,3])
b = np.array([10, 20, 30])

x = np.array(list(zip(a, b)))





#######################################################
import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import sys
import importlib
import math
from sklearn.utils import shuffle
from random import randint 
import size_temp_funcs as st
import parameters as par
import model_func as mod
import random
N = 25 # Number of consumers
M = 50 # Number of resources
T = 273.15 + 30 # Temperature
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
t_n = 25 # Number of temperatures to run the model at, model starts at 20
ass = 1 # Assembly number, i.e. how many times the system can assemble
tv = 10 # immigration times inside one assembly
t_fin = 50 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 1 # Half saturation constant
k = 0.0000862 # Boltzman constant
result_array = np.empty((0,N+M)) # Array to store data in for plotting
CUE_out = np.empty((0,N))
rich_seires = np.empty((0,tv))
U_out_total = np.empty((0,M))
U_ac_total = np.empty((0,N))
x0 = np.concatenate((np.full([N], (0.1)),np.full([M], (0.1)))) # Starting concentration for resources and consumers
Ea_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:N] # Ea for uptake
Ea_R = Ea_U - 0.4 # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
B_U = (10**(2.84 + (-4.96 * Ea_U))) + 4 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration
T_pk_U = Tref + np.random.normal(32, 5, size = N)
T_pk_R = T_pk_U + 3
p = np.concatenate((np.array([1]), np.repeat(1, M-1)))  # Resource input
U = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D)[0] # Uptake
R = par.params(N, M, T, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D)[1] # Respiration
l = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D)[2] # Leakage
l_sum = np.sum(l, axis=1)
rich = np.empty((0, tv))
t = np.linspace(0,t_fin-1,t_fin) # resetting 't' if steady state not reached (see below)
pars = (U, R, l, p, l_sum, N, M, typ, K) # Parameters to pass onto model
pops = odeint(mod.metabolic_model, y0=x0, t=t, args = pars) # Integrate
pops = np.round(pops, 7)
rem_find = pops[t_fin-1,0:N] # Get the consumer concentrations from last timestep
ext = np.where(rem_find<0.01) # Find which consumers have a concentration < 0.01 g/mL, i.e. are extinct
sur = np.where(rem_find>0.01)
# U_out = np.append(U_out, U[sur[0]], axis = 0)
U_out_total = np.append(U_out_total, U, axis = 0)
jaccard = np.zeros((N,N)) # Competition
np.fill_diagonal(jaccard,1)
for i in range(N):
    for j in range(N):
        jaccard[i,j] = np.sum(np.minimum(U[i,],U[j,]))/np.sum(np.maximum(U[i,],U[j,])) 
comp = np.mean(jaccard, axis = 0)
np.mean(comp)
np.mean(comp[sur[0]])
np.mean(comp[ext[0]])

import numpy as np
import matplotlib.pyplot as plt

k = 0.0000862 # Boltzman constant
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
T = Tref + np.linspace(0,50,51) # Temperatures
T_pk_U = Tref + 32
T_pk_R = T_pk_U + 3
Ea_U = 1.5 # Ea for uptake
Ea_R = Ea_U - 0.2 # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
B_U = (10**(2.84 + (-4.96 * Ea_U))) + 4 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration
Ea_D = 3.5 # Deactivation energy


U_Sharpe = B_U * np.exp((-Ea_U/k) * ((1/T)-(1/Tref)))/(1 + (Ea_U/(Ea_D - Ea_U)) * np.exp(Ea_D/k * (1/T_pk_U - 1/T))) 
R_Sharpe = B_R * np.exp((-Ea_R/k) * ((1/T)-(1/Tref)))/(1 + (Ea_R/(Ea_D - Ea_R)) * np.exp(Ea_D/k * (1/T_pk_R - 1/T)))

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(T - Tref, U_Sharpe, 'r-', linewidth=0.7)
ax2.plot(T - Tref, R_Sharpe, 'b-', linewidth=0.7)
ax1.set_xlabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Uptake Rate', color = 'r')
ax2.set_ylabel('Respiration Rate', color = 'b')
plt.title('Modified Sharpe-Schoolfield Temperature Dependence')
plt.show()


    colors = ['r','b']


################################## PCA code from somewhere #################################################
import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.preprocessing import StandardScaler

iris = datasets.load_iris()
X = iris.data
y = iris.target

# In general, it's a good idea to scale the data prior to PCA.
scaler = StandardScaler()
scaler.fit(X)
X=scaler.transform(X)    
pca = PCA()
x_new = pca.fit_transform(X)

def myplot(score,coeff,labels=None):
    xs = score[:,0]
    ys = score[:,1]
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    plt.scatter(xs * scalex,ys * scaley, c = y)
    for i in range(n):
        plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
        if labels is None:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+str(i+1), color = 'g', ha = 'center', va = 'center')
        else:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.xlabel("PC{}".format(1))
    plt.ylabel("PC{}".format(2))
    plt.grid()

#Call the function. Use only the 2 PCs.
myplot(x_new[:,0:2],np.transpose(pca.components_[0:2, :]))
plt.show()



############################# how to make a heatmap #######################################
data = np.normal.random((15,15))
plt.imshow( data , cmap = 'RdBu')
plt.title( "2-D Heat Map" )
plt.show()



#################################################################################################
import numpy as np
import parameters as par

N = 10 # Number of consumers
M = 5 # Number of resources
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
Ea_diff = 0.6
lf = 0.4 # Leakage
k = 0.0000862 # Boltzman constant
T = 273.15 + 35 # Temperature
Ea_U = np.round(np.random.normal(1.5, 0.01, N),3)[0:N] # Ea for uptake
Ea_R = Ea_U - Ea_diff # Ea for respiration, which should always be lower than Ea_U so 'peaks' later
B_U = (10**(2.84 + (-4.96 * Ea_U))) + 30 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
B_R = (10**(1.29 + (-1.25 * Ea_R))) # B0 for respiration
T_pk_U = Tref + np.random.normal(32, 3, size = N)
T_pk_R = T_pk_U + 3

U = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[0] # Uptake
R = par.params(N, M, T, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[1] # Respiration
l = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[2] # Leakage
l_sum = np.sum(l, axis=1)


B0_CUE = 0.175 
T_pk_CUE = T_pk_U #CUE paper
Ea_D_CUE = np.repeat(3.5,N)

CUE = (U @ (1 - l_sum)*0.1 - R)/(np.sum(U, axis = 1)*0.1)
CUE

from scipy.optimize import fsolve
EaCUE_values = []
for i in range(len(CUE)):
    def f(Ea_CUE):
        return(CUE[i] - B0_CUE * 2.71828**((-Ea_CUE/k) * ((1/T)-(1/Tref)))/(1 + (Ea_CUE/(Ea_D[i] - Ea_CUE)) * 2.71828**(Ea_D[i]/k * (1/T_pk_CUE[i] - 1/T))))
    answer = fsolve(f, 0.5)
    EaCUE_values = np.append(EaCUE_values, answer)

EaCUE_values


##############################################################################################
import numpy as np
from scipy.optimize import fsolve

Tref = 273.15
T_pk_CUE = Tref + 32
Ea_D = 3.5
B0_CUE = 0.114
k = 0.0000862
CUE = np.arange(0.2,0.65, 0.05)
T = Tref + np.arange(0, 45, 5)
EaCUE_values = np.empty(0)
for i in range(len(CUE)):
    Ea_CUE = 0
    def f(Ea_CUE):
        return(CUE[i] - B0_CUE * 2.71828**((-Ea_CUE/k) * ((1/T[i])-(1/Tref)))/(1 + (Ea_CUE/(Ea_D - Ea_CUE)) * 2.71828**(Ea_D/k * (1/T_pk_CUE - 1/T[i]))))
    answer = fsolve(f, 0.5)
    EaCUE_values = np.append(EaCUE_values, answer)
np.round(EaCUE_values, 3)


np.arange(0.1, 0.46, 0.07)




import numpy as np

M = 5
lf = 0.4
l_raw = np.array([[np.random.normal(1/(i),0.005)* lf for i in range(M,0,-1)] for i in range(1,M+1)])
fix = [[1 if j>=i else 0 for j in range(M)] for i in range(M)]
# fix[M-1][M-1] = 1
l = np.transpose(l_raw) * fix
print(l)
print(np.sum(l, axis = 1))



from scipy.integrate import odeint
from scipy.special import gamma, airy
y1_0 = 1.0 / 3**(2.0/3.0) / gamma(2.0/3.0)
y0_0 = -1.0 / 3**(1.0/3.0) / gamma(1.0/3.0)
y0 = [y0_0, y1_0]
t = np.arange(0, 4.0, 0.01)

def RHS(y, t):
    return [t*y[1],y[0]]

%timeit odeint(RHS, y0, t)

def RHS(y, t):
    y[0],y[1] = t*y[1],y[0]
    return y

%timeit odeint(RHS, y0, t)

from numba import jit
@jit
def RHS(y, t):
    y[0],y[1] = t*y[1],y[0]
    return y


a, b > 1
a > 3
(a - 1/3) / (a + b - 2/3) = 0.8

a = 15
b = ((a - 1/3) / (0.67/4)) + 2/3 - a
# np.mean(np.random.beta(a, ((a - 1/3) / 0.67) + 2/3 - a, 5000))
b
plt.hist(np.random.beta(a, b, 50000)*4)
plt.show()


U0 = 2
R0 = 1
T = 273.15 + 0 # Temperature
Tref = 273.15 + 10 # Reference temperature Kelvin, 10 degrees C
B_U = np.repeat(U0, N) # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
B_R = np.repeat(R0, N) # B0 for respiration
B0_CUE = (U0*(1 - lf) - R0)/U0 # B0 value for CUE
Ea_U = np.random.beta(35, ((35 - 1/3) / 0.8) + 2/3 - 35, N)
Ea_R = np.random.beta(35, ((35 - 1/3) / 0.8) + 2/3 - 35, N)
Ea_CUE = (R0 * (Ea_U - Ea_R))/(U0 * (1 - lf) - R0)
T_pk_U = Tref + np.random.normal(32, 5, size = N)
T_pk_R = T_pk_U + 3
U = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[0] # Uptake
R = par.params(N, M, T, k, Tref, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[1] # Respiration
l = par.params(N, M, T, k, Tref, T_pk_U, B_U, B_R,Ma, Ea_U, Ea_R, Ea_D, lf)[2] # Leakage
l_sum = np.sum(l, axis=1)

Ea_CUE = (R0 * (Ea_U - Ea_R))/(U0 * (1 - lf) - R0)
np.mean((U @ (1 - l_sum) - R)/(np.sum(U, axis = 1)))
np.mean(Ea_CUE)

result_array[:, 0:N]

rich = np.array([len(np.where(result_array[i,0:N])[0]) for i in range(len(result_array))])
len(np.where(result_array[i,0:N]>0.01)[0])


np.concatenate([np.sum(U_out_total, axis = 1).reshape(ass, N)[i][sur[i]] for i in range(len(sur))]).ravel()


all_CUE.append(np.concatenate([CUE_out[i][sur[i]] for i in range(len(sur))]).ravel())


import numpy as np
import matplotlib.pyplot as plt

k = 0.0000862 # Boltzman constant
Tref = 273.15 + 10 # Reference temperature Kelvin, 0 degrees C
T = 273.15 + np.linspace(0,30,31) # Temperatures
T_pk_U = 273.15 + 35
T_pk_R = T_pk_U + 3
lf = 0.4
Ea_D = 3.5 # Deactivation energy
# B_U = 2 # B0 for uptake - ** NB The '+4' term is added so B_U>> B_R, otherwise often the both consumers can die and the resources are consumed
B_U = 2 * np.exp((-0.8/k) * ((1/Tref)-(1/273.15)))/(1 + (0.8/(Ea_D - 0.8)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref))) # U is always 2 at 0 degree
B_R = 1 * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref)))
Ea_U = 0.8
Ea_R = 0.67
Ea_CUE = B_R*(Ea_U - Ea_R)/(B_U*(1 - lf) - B_R)

B0_CUE = 0.1 * np.exp((-Ea_CUE/k) * ((1/Tref)-(1/273.15)))
U_Sharpe = B_U * np.exp((-Ea_U/k) * ((1/T)-(1/Tref)))/(1 + (Ea_U/(Ea_D - Ea_U)) * np.exp(Ea_D/k * (1/T_pk_U - 1/T))) 
R_Sharpe = B_R * np.exp((-Ea_R/k) * ((1/T)-(1/Tref)))/(1 + (Ea_R/(Ea_D - Ea_R)) * np.exp(Ea_D/k * (1/T_pk_R - 1/T)))

CUE = (U_Sharpe*0.6-R_Sharpe)/U_Sharpe
CUE_Sharpe = B0_CUE * np.exp((-Ea_CUE/k) * ((1/T)-(1/Tref)))/(1 + (Ea_CUE/(Ea_D - Ea_CUE)) * np.exp(Ea_D/k * (1/(273.15+27) - 1/T)))
np.round(CUE,3)
np.round(CUE_Sharpe,3)

plt.plot(T - 273.15, CUE, 'g-', linewidth=0.7, label = 'CUE')
plt.plot(T - 273.15, CUE_Sharpe, 'b-', linewidth=0.7, label = 'CUE_Sharpe')
plt.legend()
plt.show()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(T - 273.15, U_Sharpe, 'r-', linewidth=0.7)
ax2.plot(T - 273.15, R_Sharpe, 'b-', linewidth=0.7)
ax1.set_xlabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Uptake Rate', color = 'r')
ax2.set_ylabel('Respiration Rate', color = 'b')
plt.title('Modified Sharpe-Schoolfield Temperature Dependence')
plt.show()

# 10 degrees, U & R
# 6.63991734506992  3.3199649753526126

################################################################
N = 100
Tref = 273.15
Ea_CUE = -0.3
k = 0.0000862
Ea_D = 3.5
lf = 0.4
B0_CUE = 0.1 * np.exp((-Ea_CUE/k) * ((1/Tref)-(1/273.15)))    
B_U = 2 * np.exp((-0.8/k) * ((1/Tref)-(1/273.15)))/(1 + (0.8/(Ea_D - 0.8)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref))) # U is always 2 at 0 degree
B_R = B_U * (1 - lf) - B0_CUE * B_U
Ea_U = np.random.beta(25, ((25 - 1/3) / 0.8) + 2/3 - 25, N)
Ea_R = Ea_U - Ea_CUE * (B_U * (1 - lf) - B_R) / B_R
# 0.675 0.990 0.849
# 0.437 0.876 0.726



from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits import axisartist
import matplotlib.pyplot as plt

host = host_subplot(111, axes_class=axisartist.Axes)
plt.subplots_adjust(right=0.75)

par1 = host.twinx()
par2 = host.twinx()

par2.axis["right"] = par2.new_fixed_axis(loc="right", offset=(60, 0))

par1.axis["right"].toggle(all=True)
par2.axis["right"].toggle(all=True)

p1, = host.plot([0, 1, 2], [0, 1, 2], label="Density")
p2, = par1.plot([0, 1, 2], [0, 3, 2], label="Temperature")
p3, = par2.plot([0, 1, 2], [50, 30, 15], label="Velocity")

host.set_xlim(0, 2)
host.set_ylim(0, 2)
par1.set_ylim(0, 4)
par2.set_ylim(1, 65)

host.set_xlabel("Distance")
host.set_ylabel("Density")
par1.set_ylabel("Temperature")
par2.set_ylabel("Velocity")

host.legend()

host.axis["left"].label.set_color(p1.get_color())
par1.axis["right"].label.set_color(p2.get_color())
par2.axis["right"].label.set_color(p3.get_color())

plt.show()


import parameters as par
import size_temp_funcs as st
import numpy as np

N = 10 # Number of consumers
M = 5 # Number of resources

# Temperature params
Tref = 273.15 # Reference temperature Kelvin, 0 degrees C
pk_U =  273.15 + np.random.normal(35, 5, size = N)
Ma = 1 # Mass
Ea_D = np.repeat(3.5,N) # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
k = 0.0000862 # Boltzman constant
T_pk = Tref + pk_U # Peak above Tref, Kelvin
T = 273.15 + 30 # Temperature
Ea_U = np.random.beta(4, ((4 - 1/3) / 0.82) + 2/3 - 4, N) # Ea for uptake
B_U = 2 * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref))) # B0 for uptake

U_sum = st.temp_growth(k, T, Tref, T_pk, N, B_U, Ma, Ea_U, Ea_D)
np.random.seed(0)
diri = np.transpose(np.random.dirichlet(np.ones(M),N))
U = np.transpose(diri*U_sum)
print(U_sum/np.sum(U_sum))
print(np.std(U_sum/np.sum(U_sum)))
U/np.repeat(U_sum, M).reshape(N, M)



import numpy as np

M = 5
lf = 0.4

l_raw = np.array([[np.random.normal(1/(i),0.005)* lf for i in range(M,0,-1)] for i in range(1,M+1)])
l = [[l_raw[j,i] if j>=i else 0 for j in range(M)] for i in range(M)]
l = np.round(l, 3)
print(np.round(l, 3))
print(np.sum(l, axis = 1))

A = np.array((1,2,3,0,5))
plt.scatter(np.repeat(1,5), A)
plt.scatter(np.repeat(2,5), A)
plt.show()


Ea_D = 3.5
Tref = 273.15 + np.linspace(0,50,51) # Temperatures
B_U = np.random.uniform(1, 3, N) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref))) # U is always 2 at 0 degree
B_R = np.random.uniform(0, 2, N) * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref))) 



a = 1.4
b = ((a - 1/3) / (0.5326/3)) + 2/3 - a
# np.mean(np.random.beta(a, ((a - 1/3) / 0.67) + 2/3 - a, 5000))
b
plt.hist(np.random.beta(a, b, 50000)*3)
plt.show()

lam = np.log(2)/0.089 # median value 0.089
plt.hist(np.random.exponential(1/lam, 10000))
plt.show()

plt.hist(np.random.exponential(0.4810934, 10000))
plt.show()

lam = np.log(2)/0.1790591 # median value 0.089
plt.hist(np.random.exponential(1/lam, 10000))
plt.show()

b = 1
np.median(np.random.beta(b, ((b - 1/3) / 0.175) + 2/3 - b, 50000))
plt.hist(np.random.beta(b, ((b - 1/3) / 0.175) + 2/3 - b, 50000))
plt.show()

def selection_sort(x):
    for i in range(len(x)):
        swap = i + np.argmin(x[i:])
        (x[i], x[swap]) = (x[swap], x[i])
    return x

array = (B_U*0.6 - B_R)/B_U
selection_sort(array)


plt.hist(np.random.uniform(0.05230875, 0.3884611, 10000))
plt.show()


K = 1
x = np.arange(0, 2, 0.1)
plt.plot(x, 4.47*x/(K+x))
plt.plot(x, 4.47*x)
plt.show()


for i in range(6):
    K = 2**(i - 3)
    print(K)



B = np.diag(x[N:N+M])

1.96 * np.std(np.sum(U_out_total, axis = 1))/((ass*N)**0.5)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
for i in range(6):
    T = 273.15 + 5*i
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, invasion_result, inv_U, inv_R = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K, inv, sp)
    y = 1 - np.mean(R_out)/np.mean(np.sum(U_out_total, axis = 1))
    ax1.scatter(i*5,y)
    ax2.scatter(i*5, np.mean(CUE_out),color = "red")
plt.show()


U = R/(1 - lf - CUE_0)

np.var(np.array((0,101,1)))

im = plt.imshow(cf, cmap = 'binary')
ax = plt.gca()
ax.set_xticks(np.arange(-.5,5,1), minor = True)
ax.set_yticks(np.arange(-.5,5,1), minor = True)
plt.xticks([])
plt.yticks([])
ax.grid(which = 'minor', color = 'k', linestyle='-', linewidth=1)
plt.title( "" )
plt.grid()
plt.show()

im = plt.imshow(jaccard, cmap = 'binary')
ax = plt.gca()
ax.set_xticks(np.arange(-.5,5,1), minor = True)
ax.set_yticks(np.arange(-.5,5,1), minor = True)
plt.xticks([])
plt.yticks([])
ax.grid(which = 'minor', color = 'k', linestyle='-', linewidth=1)
plt.title( "" )
plt.grid()
plt.show()

while np.any(infodict.get('nfe')< 0):
    print(1)

n = 0
while n<2:
    n = n+1
print(n)

np.any(infodict.get('tolsf')<0)

for i in range(5):
    if i == 0:
        l = 0
    else:
        l = i*0.2-0.1
    print(l)


U_s = U_out_total[(ass-1)*N-1:ass*N-1,:][sur[ass-1]]
R_s = R_out[ass-1][sur[ass-1]]
CUE_s = CUE_out[ass-1][sur[ass-1]]
l
l_sum = np.sum(l, axis=1)
t = np.linspace(0,t_fin-1,t_fin) 
for i in range(rich_series[ass-1]):
    pars = (U_s[i,:], R_s[i], l, p, l_sum, 1, M, typ, K) # Parameters to pass onto model
    pops, infodict = odeint(mod.metabolic_model, y0=x0, t=t, args = pars, full_output=1) # Integrate

