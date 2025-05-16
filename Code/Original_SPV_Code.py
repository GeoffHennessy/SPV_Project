import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

# Constants
e = 1.602 * 10 ** - 19  # Fundamental unit of charge
k = 8.617 * 10 ** (-5)
eps = 11.5 * e / 100  # Electrical permittivity of silicon in correct units




def ohm(voltage,T):
    return voltage / (k * T)

def LD (T,eps,pb):
    return (np.sqrt(k * T * eps / (pb * e **2)))

def Qss(voltage,E,T,Nt):
    return Nt * (1 - 1 / (1 + np.exp((E - voltage) / (k * T))))

def Qsc(voltage,T,eps,nb,pb,function1,function2):
    # function 1 is ohm(), function2 is LD()
    return (eps * np.sqrt(2) * k * T / (e * function2(T,eps,pb)) * ((np.exp(-function1(voltage,T)) + function1(voltage,T) - 1) +
    (nb / pb) * (np.exp(function1(voltage,T)) - function1(voltage,T) - 1)) ** (1 / 2))

def intersection(voltage,E,T,Nt,nb,pb,function1,function2,function3,function4):
    # function1 is Qss, function2 is Qsc, function3 is ohm, function4 is LD
    Plot1 = function1(voltage,E,T,Nt)
    Plot2 = function2(voltage,T,eps,nb,pb,function3,function4)
    return(Plot2-Plot1)

def IntersectionGraph(voltage,E,T,Nt,nb,pb,eps,function1,function2,function3,function4):
    # function 1 is Qss(), function 2 is Qsc(), function3() is ohm, function 4 is LD
    Qss = function1(voltage,E,T,Nt)
    Qsc = function2(voltage,T,eps,nb,pb,function3,function4)

    plt.plot(voltage, Qss, label='Qss')
    plt.plot(voltage, Qsc, label='Qsc')
    plt.yscale("log")
    plt.ylim(1*10**-6,1*10**-1)
    plt.xlim(-0.12,0.6)
    plt.legend()
    plt.xlabel('Surface Potential (V)', fontsize = 22)
    plt.ylabel(r'Surface Charge C$\cdot$ cm$^{-2}$', fontsize = 22)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.title(r'Tracking $Q_{ss}$ amd $Q_{sc}$ as functions of voltage', fontsize = 26)
    plt.show()
    return



def VvsE(Resolution,EStart,T,Nt,nb,pb,function1,function2,function3,function4):
    # function1 is Qss, function2 is Qsc, function3 is ohm, function4 is LD
    VArray = np.zeros(Resolution)
    EArray = np.zeros(Resolution)

    for x in range(1, Resolution):
        E = EStart + 0.01 * x
        Intersect = scipy.optimize.bisect(intersection, 0, 2, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray[x] = Intersect
        EArray[x] = E

    plt.plot(EArray, VArray)
    plt.xlabel(r'Energy Difference (V)', fontsize=15)
    plt.ylabel('Surface Potential (V)', fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title('Relationship Between E and Surface Potential', fontsize=17)
    plt.show()
    return

def VvsPb(Resolution,pbStart,E,T,Nt,function1,function2,function3,function4):
    VArray = np.array([])
    pbArray = np.array([])

    T = 100
    for x in range(Resolution):
        nb = 1 * 10 ** (pbStart - 3)
        pb = 10 ** (pbStart + 0.1 * x) # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        pbArray = np.append(pbArray,pb)
    plt.plot(pbArray,VArray, label = 'T = 100')

    VArray = np.array([])
    pbArray = np.array([])
    T = 150
    for x in range(Resolution):
        nb = 1 * 10 ** (pbStart - 3)
        pb = 10 ** (pbStart + 0.1 * x) # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        pbArray = np.append(pbArray,pb)
    plt.plot(pbArray,VArray, label = 'T = 150')

    VArray = np.array([])
    pbArray = np.array([])
    T = 200
    for x in range(Resolution):
        nb = 1 * 10 ** (pbStart - 3)
        pb = 10 ** (pbStart + 0.1 * x) # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        pbArray = np.append(pbArray,pb)
    plt.plot(pbArray,VArray, label = 'T = 200')

    VArray = np.array([])
    pbArray = np.array([])
    T = 250
    for x in range(Resolution):
        nb = 1 * 10 ** (pbStart - 3)
        pb = 10 ** (pbStart + 0.1 * x) # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        pbArray = np.append(pbArray,pb)
    plt.plot(pbArray,VArray, label = 'T = 250')

    VArray = np.array([])
    pbArray = np.array([])
    T = 300
    for x in range(Resolution):
        nb = 1 * 10 ** (pbStart - 3)
        pb = 10 ** (pbStart + 0.1 * x) # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        pbArray = np.append(pbArray,pb)
    plt.plot(pbArray,VArray, label = 'T = 300')

    VArray = np.array([])
    pbArray = np.array([])
    T = 350
    for x in range(Resolution):
        nb = 1 * 10 ** (pbStart - 3)
        pb = 10 ** (pbStart + 0.1 * x) # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        pbArray = np.append(pbArray,pb)
    plt.plot(pbArray,VArray, label = 'T = 350')


    plt.legend(fontsize = 15)
    plt.xscale('log')
    plt.xlabel(r'$n_b$', fontsize = 30)
    plt.ylabel('Surface Potential (V)', fontsize = 30)
    plt.xticks(fontsize = 26)
    plt.yticks(fontsize = 26)
    plt.title(r'Relationship between Electron Density $(n_b)$ and Surface Potential', fontsize = 35)
    plt.show()
    return

def VvsNb(Resolution,nbStart,E,T,Nt,function1,function2,function3,function4):
    VArray = np.array([])
    nbArray = np.array([])

    for x in range(Resolution):
        nb = 1 * 10 ** (nbStart + 0.1 * x)
        pb = nb * 1000 # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        nbArray = np.append(nbArray,nb)

    plt.plot(nbArray,VArray)
    plt.xscale('log')
    plt.xlabel(r'$n_b$', fontsize = 30)
    plt.ylabel('Surface Potential (V)', fontsize = 30)
    plt.xticks(fontsize = 26)
    plt.yticks(fontsize = 26)
    plt.title(r'Relationship between Electron Density $(n_b)$ and Surface Potential', fontsize = 35)
    plt.show()
    return

def VvsT(Resolution,TStart,E,Nt,nb,pb,function1,function2,function3,function4):
    VArray = np.array([])
    TArray = np.array([])

    for y in range(5):
        E = y * 0.1
        for x in range(Resolution):
            T = TStart + x #This will be set to start at 100 and end at 600 to match the cross in the later EvsT heat map
            Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
            VArray = np.append(VArray, Intersect)
            TArray = np.append(TArray,T)


    plt.plot(TArray,VArray)
    plt.xlabel(r'$Temperature (k)$', fontsize = 30)
    plt.ylabel('Surface Potential (V)', fontsize = 30)
    plt.xticks(fontsize = 26)
    plt.yticks(fontsize = 26)
    plt.title(r'Relationship between Temperature and Surface Potential', fontsize = 35)
    plt.show()
    return

# Heat maps are here
def EAndTMap(Estart,Eend,Tstart,Tend,Nt,nb,pb,function1,function2,function3,function4):
    DE = 0.005
    DT = 1
    ELen = int((Eend-Estart) / DE)
    TLen = int((Tend - Tstart) / DT)
    EArray = np.linspace(Estart,Eend,ELen)
    TArray = np.linspace(Tstart,Tend,TLen)
    VArray = np.zeros((ELen,TLen))
    for x in range(ELen):
        E = EArray[x]

        for y in range(TLen):
            T = TArray[y]
            intersect = scipy.optimize.bisect(intersection, 0, 1, args=(E,T,Nt,nb,pb,function1,function2,function3,function4),  xtol=1.0e-5, maxiter=400)
            VArray[x,y] = intersect

    plt.pcolormesh(TArray, EArray, VArray, shading='auto')
    clb = plt.colorbar()
    clb.set_label('Surface Potential (V)', size = 40)
    clb.ax.tick_params(labelsize= 34)
    plt.xlabel('Temperature (k)', fontsize= 40)
    plt.xticks(fontsize= 34)
    # plt.xticks(np.arange(Estart,Eend,step=0.1))
    plt.ylabel('Energy Difference (V)', fontsize= 40)
    plt.yticks(fontsize= 34)
    # plt.title('Surface Potential as a Function of Temperature and Energy Difference', fontsize=25)
    plt.show()
    return

def NbvsPb(nbStart,nbEnd,pbRes,E,T,Nt,function1,function2,function3,function4):
    Dnb = 0.1
    Dpb = 0.01
    nbLen = int((nbEnd-nbStart) / Dnb)
    pbLen = int((pbRes) / Dpb)
    nbArray = np.logspace(nbStart, nbEnd, nbLen)
    pbArray = np.logspace(nbStart + 3, nbEnd + pbRes, pbLen)
    VArray = np.zeros((nbLen, pbLen))

    for x in range(nbLen):
        nb = nbArray[x]

        for y in range(pbLen):
            pb = pbArray[y]
            intersect = scipy.optimize.bisect(intersection, 0, 10, args=(E,T,Nt,nb,pb,function1,function2,function3,function4),  xtol=1.0e-5, maxiter=400)
            V = intersect
            VArray[x,y] = V

    print(VArray[25,25])
    plt.pcolormesh(pbArray,nbArray,VArray, shading = 'auto')
    clb = plt.colorbar()
    clb.ax.tick_params(labelsize=34)
    clb.set_label('Surface Potential (V)', fontsize = 40)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$p_b$',fontsize = 40)
    plt.xticks(fontsize = 34)
    plt.ylabel(r'$n_b$', fontsize = 40)
    plt.yticks(fontsize = 34)
    plt.show()
    return

def TvsPb(TStart,TEnd,pbStart,pbEnd,E,Nt,function1,function2,function3,function4):
    DT = 1
    Dpb = 0.05

    TLen = int((TEnd - TStart) / DT)
    pbLen = int((pbEnd - pbStart) / Dpb)
    TArray = np.linspace(TStart, TEnd, TLen)
    pbArray = np.logspace(pbStart, pbEnd, pbLen)
    VArray = np.zeros((TLen, pbLen))

    for x in range(TLen):
        T = TArray[x]

        for y in range(pbLen):
            pb = pbArray[y]
            nb = pb/1000 # Note that the ratio of nb and pb is kept constant here
            intersect = scipy.optimize.bisect(intersection, 0, 10, args=(E,T,Nt,nb,pb,function1,function2,function3,function4),  xtol=1.0e-5, maxiter=400)
            V = intersect
            VArray[x,y] = V
    yticks = np.arange(TStart,TEnd,50)
    plt.pcolormesh(pbArray,TArray,VArray, shading = 'auto')
    clb = plt.colorbar()
    clb.set_label('Surface Potential (V)', size = 40)
    clb.ax.tick_params(labelsize= 34)
    plt.xscale('log')
    plt.xlabel(r'$p_b$',fontsize = 40)
    plt.xticks(fontsize = 34)
    plt.ylabel('Temperature (K)', fontsize = 40)
    plt.yticks(yticks,fontsize = 34)
    plt.show()
    return

def EvsPb(EStart,EEnd,pbStart,pbEnd,T,Nt,function1,function2,function3,function4):
    DE = 0.01
    Dpb = 0.05

    ELen = int((EEnd - EStart) / DE)
    pbLen = int((pbEnd - pbStart) / Dpb)
    EArray = np.linspace(EStart, EEnd, ELen)
    pbArray = np.logspace(pbStart, pbEnd, pbLen)
    VArray = np.zeros((ELen, pbLen))

    for x in range(ELen):
        E = EArray[x]

        for y in range(pbLen):
            pb = pbArray[y]
            nb = pb/1000 # Note that the ratio of nb and pb is kept constant here
            intersect = scipy.optimize.bisect(intersection, 0, 10, args=(E,T,Nt,nb,pb,function1,function2,function3,function4),  xtol=1.0e-5, maxiter=400)
            V = intersect
            VArray[x,y] = V

    plt.pcolormesh(pbArray,EArray,VArray, shading = 'auto')
    clb = plt.colorbar()
    clb.set_label('Surface Potential (V)', size = 40)
    clb.ax.tick_params(labelsize=34)
    plt.xscale('log')
    plt.xlabel(r'$p_b$',fontsize = 40)
    plt.xticks(fontsize = 34)
    plt.ylabel('Energy Difference', fontsize = 40)
    plt.yticks(fontsize = 34)
    plt.show()
    return

# These are the default parameters that will be used in each function. Change them here or when you call the function
Nt = 1.0 * 10 ** -2  # Available surface states per unit area
E = 0.2
T = 293.15  # Temperature in kelvin (Used STP...need's updating)
pb = 1.0 * 10 ** 15
nb = 1.0 * 10 ** 10
V = np.arange(-0.11, 1, 0.0001)

# Call the functions here
IntersectionGraph(V,0.2,T,Nt,nb,pb,eps,Qss,Qsc,ohm,LD)
# VvsNb(100,5,E,T,Nt,Qss,Qsc,ohm,LD) # nb start is the power of nb. So nbStart = 7 would set nb to 10^7
# VvsE(150,-0.4,100,Nt,nb,pb,Qss,Qsc,ohm,LD)
VvsT(500,100,E,Nt,nb,pb,Qss,Qsc,ohm,LD) #VvsT(Resolution,TStart,E,Nt,nb,pb,function1,function2,function3,function4)
# VvsPb(150,8,E,T,Nt,Qss,Qsc,ohm,LD)
# EAndTMap(-0.5,1,100,600,Nt,nb,pb,Qss,Qsc,ohm,LD) # Paramerters (Estart,Eend,Tstart,Tend...)
# NbvsPb(5,18,3,E,T,Nt,Qss,Qsc,ohm,LD)
# TvsPb(100,351,8,18,E,Nt,Qss,Qsc,ohm,LD)
# EvsPb(-0.5,1.0,8,18,T,Nt,Qss,Qsc,ohm,LD)