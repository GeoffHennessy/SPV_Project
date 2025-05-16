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

def intersection(voltage,E,T,Nt,nb,pb):
    Plot1 = Qss(voltage,E,T,Nt)
    Plot2 = Qsc(voltage,T,eps,nb,pb)
    return(Plot2-Plot1)

def intersection(voltage, *args):
    E,T,Nt,nb,pb,eps = args
    return Qsc(voltage,T,eps,nb,pb) - Qss(voltage,E,T,Nt)

# For Graph functions, args will always be default values
def VvsE(E):
    return opt.bisect(intersection,0,5, args = (E,T0,Nt,nb0,pb0,eps), xtol = 1e-5, maxiter = 100)
    
def VvsT(T):
    val = opt.bisect(intersection, 0,5, args = (E0,T,Nt,nb0,pb0,eps), xtol = 1e-5, maxiter = 100)
    return val

def VvsNb(nb):
    return opt.minimize(intersection,0,5, args = (E0,T0,Nt,nb,pb0,eps), xtol = 1e-5, maxiter = 100)

def VvsPb(pb):
    val = opt.minimize(intersection,0,5, args = (E0,T0,Nt,nb0,pb,eps), xtol = 1e-5, maxiter = 100)
    return val








#%%%%%%%%%%%%%%%%%%%%%%%%
# Generating the data
VVals = np.arange(-0.2,0.5,1e-5)
x = opt.minimize(intersection, 0.2, args = (E0,T0,Nt,nb0,pb0,eps))
#Graph 1: (Qss, Qsc0 vs V
QssVals = Qss(VVals,E0,T0,Nt)
QscVals = Qsc(VVals,T0,eps,nb0,pb0)

VT_Vals = []
TVals = np.arange(100,401,1)
for i in TVals:
    V_Val = VvsT(i)
    VT_Vals.append(V_Val)


#%%%%%%%%%%%%%%%%%%%%%%%%
# Generating Graph
fig1 = plt.figure(figsize = (10,10))
ax1 = plt.subplot()
ax1.plot(VVals,QssVals, c = "orange", label = "Qss")
ax1.plot(VVals,QscVals, c = "blue", label = "Qsc")
ax1.set_yscale("log")
ax1.legend()

fig2 = plt.figure(figsize = (10,10))
ax2 = plt.subplot()
ax2.plot(TVals,VT_Vals)
plt.show()
