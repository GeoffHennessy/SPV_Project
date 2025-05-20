import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

# Constants
e = 1.602 * 10 ** - 19  # Fundamental unit of charge
k = 8.617 * 10 ** (-5)
eps = 11.5 * e / 100  # Electrical permittivity of silicon in correct units


# These are the default parameters that will be used in each function. Change them here or when you call the function
Nt = 1.0 * 10 ** -2  # Available surface states per unit area
E0 = 0.2
T0 = 293.15
pb0 = 1.0 * 10 ** 15
nb0= 1.0 * 10 ** 10
V = np.arange(-0.11, 1, 0.0001)


def ohm(voltage,T):
    return voltage / (k * T)

def LD (T,eps,pb):
    return (np.sqrt(k * T * eps / (pb * e **2)))

def Qss(voltage,E,T,Nt):
    return Nt * (1 - 1 / (1 + np.exp((E - voltage) / (k * T))))

def Qsc(voltage,T,eps,nb,pb):
    # function 1 is ohm(), function2 is LD()
    return (eps * np.sqrt(2) * k * T / (e * LD(T,eps,pb)) * ((np.exp(-ohm(voltage,T)) + ohm(voltage,T) - 1) +
    (nb / pb) * (np.exp(ohm(voltage,T)) - ohm(voltage,T) - 1)) ** (1 / 2))

def intersection(voltage, *args):
    E,T,Nt,nb,pb,eps = args
    return Qsc(voltage,T,eps,nb,pb) - Qss(voltage,E,T,Nt)

# For Graph functions, args will always be default values
def VvsE(E):
    return opt.bisect(intersection,0,5, args = (E,T0,Nt,nb0,pb0,eps), xtol = 1e-5, maxiter = 100)
    
def VvsT(T):
    val = opt.bisect(intersection, 0,5, args = (E0,T,Nt,nb0,pb0,eps), xtol = 1e-5, maxiter = 100)
    return val

def VvsNb(nb,ratio): # These have a ratio because the ratio of Nb, Pb are often fixed for a material
    return opt.bisect(intersection,0,5, args = (E0,T0,Nt,nb,nb*ratio,eps), xtol = 1e-5, maxiter = 100)

def VvsPb(pb,ratio): # These have a ratio because the ratio of Nb, Pb are often fixed for a material
    val = opt.bisect(intersection,0,5, args = (E0,T0,Nt,pb/ratio,pb,eps), xtol = 1e-5, maxiter = 100)
    return val


# Generating Data
Evals = np.linspace(-0.4,0.4,200)
Pbvals = np.logspace(9,21,200)
Nbvals = np.logspace(6,18,200)
Tvals = np.linspace(100,400,200)

#%%%%%%%%%%%%%%%%%%%%%
# Graph 1: E vs T heat map
grid1 = np.zeros((len(Evals),len(Tvals)))

for i in range(len(Evals)):
    E0 = Evals[i]
    for j in range(len(Tvals)):
        T = Tvals[j]
        Vval = VvsT(T)
        grid1[i,j] = Vval

E0, T0 = 0.2, 293.15

# Graph 2: Pb and Nb heat map
grid2 = np.zeros((len(Pbvals),len(Nbvals)))

for i in range(len(Pbvals)):
    pb0 = Pbvals[i]
    for j in range(len(Nbvals)):
        Nb = Nbvals[j]
        Vval = VvsNb(Nb,pb0/Nb)
        grid2[i,j] = Vval

pb0, nb0 = 1e15, 1e10

# Graph 3: Pb and T heat map
grid3 = np.zeros((len(Tvals),len(Pbvals)))

for i in range(len(Tvals)):
    T0 = Tvals[i]
    for j in range(len(Pbvals)):
        pb = Pbvals[j]
        Vval = VvsPb(pb,1000)
        grid3[i,j] = Vval

pb0,T0 = 1e15, 293.15

# Graph 4: Pb and E heat map
grid4 = np.zeros((len(Evals),len(Pbvals)))

for i in range(len(Evals)):
    E0 = Evals[i]
    for j in range(len(Pbvals)):
        Pb = Pbvals[j]
        Vval = VvsPb(Pb,1000)
        grid4[i,j] = Vval

pb0, E0 = 1e15, .2


# Graph 5: Nb and T heat map
grid5 = np.zeros((len(Tvals),len(Nbvals)))

for i in range(len(Tvals)):
    T0 = Tvals[i]
    for j in range(len(Nbvals)):
        Nb = Nbvals[j]
        Vval = VvsNb(Nb,1000)
        grid5[i,j] = Vval

nb0, T0 = 1e10, 293.15

# Graph 6: Nb and E heat map
grid6 = np.zeros((len(Evals),len(Nbvals)))

for i in range(len(Evals)):
    E0 = Evals[i]
    for j in range(len(Nbvals)):
        Nb = Nbvals[j]
        Vval = VvsNb(Nb,1000)
        grid6[i,j] = Vval

pb0, E0 = 1e15, 0.2

# Graph 7: Nb and T heat map not constant ratio
grid7 = np.zeros((len(Tvals),len(Nbvals)))

for i in range(len(Tvals)):
    T0 = Tvals[i]
    for j in range(len(Nbvals)):
        Nb = Nbvals[j]
        Vval = VvsNb(Nb,pb0/Nb)
        grid7[i,j] = Vval

nb0, T0 = 1e10, 293.15

# Graph 8: Nb and E heat map not constant ratio
grid8 = np.zeros((len(Evals),len(Nbvals)))

for i in range(len(Evals)):
    E0 = Evals[i]
    for j in range(len(Nbvals)):
        Nb = Nbvals[j]
        Vval = VvsNb(Nb,pb0/Nb)
        grid8[i,j] = Vval

pb0, E0 = 1e15, 0.2

# Graph 9: Pb and T heat map not constant ratio
grid9 = np.zeros((len(Tvals),len(Pbvals)))

for i in range(len(Tvals)):
    T0 = Tvals[i]
    for j in range(len(Pbvals)):
        pb = Pbvals[j]
        Vval = VvsPb(pb,pb/nb0)
        grid9[i,j] = Vval

pb0,T0 = 1e15, 293.15

# Graph 10: Pb and E heat map not constant ratio
grid10 = np.zeros((len(Evals),len(Pbvals)))

for i in range(len(Evals)):
    E0 = Evals[i]
    for j in range(len(Pbvals)):
        Pb = Pbvals[j]
        Vval = VvsPb(Pb,pb/nb0)
        grid10[i,j] = Vval

pb0, E0 = 1e15, .2

#%%%%%%%%%%%
# Graphing results

# #Graph 1:
# fig1 = plt.figure(figsize = (10,10))
# ax1 = plt.subplot()
# ax1.set_title("E Vs T")
# xticks = np.arange(-0.4-0.4,0.2)
# yticks = np.arange(100,400,50)
# graph1 = ax1.pcolormesh(Evals,Tvals,grid1)
# bar = fig1.colorbar(graph1)

# #Graph 2:
# fig2 = plt.figure(figsize = (10,10))
# ax2 = plt.subplot()
# ax2.set_title("Pb vs Nb")
# ax2.set_xscale("log")
# ax2.set_yscale("log")
# graph2 = ax2.pcolormesh(Pbvals,Nbvals,grid2)
# bar = fig2.colorbar(graph2)


# #Graph 3:
# fig3 = plt.figure(figsize = (10,10))
# ax3= plt.subplot()
# ax3.set_title("T Vs Pb")
# ax3.set_yscale("log")
# graph3 = ax3.pcolormesh(Tvals,Pbvals,grid3)
# bar = fig3.colorbar(graph3)


# #Graph 4:
# fig4 = plt.figure(figsize = (10,10))
# ax4= plt.subplot()
# ax4.set_title("E Vs Pb")
# ax4.set_yscale("log")
# graph4 = ax4.pcolormesh(Evals,Pbvals,grid4)
# bar = fig4.colorbar(graph4)


# #Graph 5:
# fig5 = plt.figure(figsize = (10,10))
# ax5= plt.subplot()
# ax5.set_title("T vs Nb (Pb/Nb is constant)")
# ax5.set_yscale("log")
# graph5 = ax5.pcolormesh(Tvals,Nbvals,grid5)
# bar = fig5.colorbar(graph5)


# #Graph 6:
# fig6 = plt.figure(figsize = (10,10))
# ax6= plt.subplot()
# ax6.set_title("E vs Nb (Pb/Nb is constant)")
# ax6.set_yscale("log")
# graph6 = ax6.pcolormesh(Evals,Nbvals,grid6)
# bar = fig6.colorbar(graph6)


# #Graph 7:
# fig7 = plt.figure(figsize = (10,10))
# ax7= plt.subplot()
# ax7.set_title("T vs Nb (Pb/Nb is not constant)")
# ax7.set_yscale("log")
# graph7 = ax7.pcolormesh(Tvals,Nbvals,grid7)
# bar = fig7.colorbar(graph7)


# #Graph 8:
# fig8 = plt.figure(figsize = (10,10))
# ax8= plt.subplot()
# ax8.set_title("E vs Nb (Pb/Nb is not constant)")
# ax8.set_yscale("log")
# graph8 = ax8.pcolormesh(Evals,Nbvals,grid8)
# bar = fig8.colorbar(graph8)

#Graph 9:
fig9 = plt.figure(figsize = (10,10))
ax9= plt.subplot()
ax9.set_title("T vs Pb (Pb/Nb is not constant)")
ax9.set_yscale("log")
graph9 = ax9.pcolormesh(Tvals,Nbvals,grid9)
bar = fig9.colorbar(graph9)


#Graph 10:
fig10 = plt.figure(figsize = (10,10))
ax10= plt.subplot()
ax10.set_title("E vs Pb (Pb/Nb is not constant)")
ax10.set_yscale("log")
graph10 = ax10.pcolormesh(Evals,Nbvals,grid10)
bar = fig10.colorbar(graph10)
plt.show()