import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
# Try ImShow for graphing
# Turn the heat maps into a function. Also finish and check the EvNb function.
# Constants
e = 1.602 * 10 ** - 19  # Fundamental unit of charge
k = 8.617 * 10 ** (-5)
# ohm = V / (k * T)  # This normally has e in the numerator. But when we express V as eV, the e cancels
eps = 11.5 * e / 100  # Electrical permittivity of silicon in correct units

# Formulas used to generate graph
def ohm(voltage,T):
    return voltage / (k * T)

def LD(T,eps,pb):
    return (np.sqrt(k * T * eps / (pb * e ** 2)))

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

# def IntersectionGraph(voltage,E,T,Nt,nb,pb,eps,function1,function2,function3,function4):
#     # function 1 is Qss(), function 2 is Qsc(), function3() is ohm, function 4 is LD
#     Qss = function1(voltage,E,T,Nt)
#     Qsc = function2(voltage,T,eps,nb,pb,function3,function4)
#
#     plt.plot(voltage, Qss, label='Qss')
#     plt.plot(voltage, Qsc, label='Qsc')
#     plt.yscale("log")
#     plt.ylim(1*10**-6,1*10**-1)
#     plt.legend()
#     plt.xlabel('Surface Potential (V)', fontsize = 16)
#     plt.ylabel(r'Surface Charge C$\cdot$ $cm^{-2}$', fontsize = 16)
#     plt.xticks(fontsize = 15)
#     plt.yticks(fontsize = 15)
#     plt.title(r'Tracking $Q_{ss}$ amd $Q_{sc}$ as functions of voltage', fontsize = 18)
#     plt.show()
#     return

# def VvsE(Resolution,EStart,T,Nt,nb,pb,function1,function2,function3,function4):
#     # function1 is Qss, function2 is Qsc, function3 is ohm, function4 is LD
#     VArray = []
#     EArray = []
#
#     for x in range(1, Resolution):
#         E = EStart + 0.01 * x
#         Intersect = scipy.optimize.bisect(intersection, 0, 2, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
#         VArray = np.append(VArray, Intersect)
#         EArray = np.append(EArray, E)
#
#     plt.plot(EArray, VArray)
#     plt.xlabel(r'Energy Difference (V)', fontsize=15)
#     plt.ylabel('Surface Potential (V)', fontsize=15)
#     plt.xticks(fontsize=12)
#     plt.yticks(fontsize=12)
#     plt.title('Relationship Between E and Surface Potential', fontsize=17)
#     plt.show()
#     return




# #This graph tracks V as n_b grows exponentially. Note the ratio of n_b/p_b remains the same

# def VvsNb(Resolution,nbStart,E,T,Nt,function1,function2,function3,function4):
#     VArray = np.array([])
#     nbArray = np.array([])
#
#     for x in range(Resolution):
#         nb = 1 * 10 ** (nbStart + 0.1 * x)
#         pb = nb * 1000 # This keeps the ratio between nb and pb the same
#
#         Intersect = scipy.optimize.bisect(intersection, 0, 2, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
#         VArray = np.append(VArray, Intersect)
#         nbArray = np.append(nbArray,nb)
#
#     plt.plot(nbArray,VArray)
#     plt.xscale('log')
#     plt.xlabel(r'$n_b$', fontsize = 15)
#     plt.ylabel('Surface Potential (V)', fontsize = 15)
#     plt.xticks(fontsize = 12)
#     plt.yticks(fontsize = 12)
#     plt.title(r'Relationship between Electron Density $(n_b)$ and Surface Potential', fontsize = 17)
#     plt.show()
#     return



# This section generates the E and T matrices that will be used in pcolormesh
# EArray,TArray = np.mgrid[-0.2:1.80:100j,100:501:50j]

# This section generates the values of V as functions of E and T and stores that in a matrix.


# def EAndTMap(Estart,Eend,Tstart,Tend,Nt,nb,pb,function1,function2,function3,function4):
#     DE = 0.005
#     DT = 1
#     ELen = int((Eend-Estart) / DE)
#     TLen = int((Tend - Tstart) / DT)
#     EArray = np.linspace(Estart,Eend,ELen)
#     TArray = np.linspace(Tstart,Tend,TLen)
#     VArray = np.zeros((ELen,TLen))
#     for x in range(ELen):
#         E = EArray[x]
#
#         for y in range(TLen):
#             T = TArray[y]
#             intersect = scipy.optimize.bisect(intersection, 0, 1, args=(E,T,Nt,nb,pb,function1,function2,function3,function4),  xtol=1.0e-5, maxiter=400)
#             VArray[x,y] = intersect
#
#     plt.pcolormesh(TArray, EArray, VArray, shading='auto')
#     # plt.pcolormesh(TArray,EArray, VArray, shading='auto')
#     clb = plt.colorbar()
#     clb.set_label('Surface Potential (V)', fontsize=15)
#     plt.xlabel('Temperature (k)', fontsize=15)
#     plt.xticks(fontsize=12)
#     plt.ylabel('Energy Difference (V)', fontsize=15)
#     plt.yticks(fontsize=12)
#     plt.title('Equilibrium Potential as a Function of Temperature and Energy Difference', fontsize=17)
#     plt.show()
#     return

# This is the older version
# def EAndTMap(ERes,TRes,EStart,TStart,Nt,nb,pb,function1,function2,function3,function4):
#     VArray = np.zeros((ERes, TRes))
#     EArray = np.arange(EStart,EStart + 0.02 * ERes, 0.02)
#     TArray = np.arange(TStart,TStart + 8 * TRes, 8)
#     for x in range(ERes):
#         E = EStart + 0.002 * x
#
#         for y in range(TRes):
#             T = 8 * y + TStart
#             intersect = scipy.optimize.bisect(intersection, 0, 10, args=(E,T,Nt,nb,pb,function1,function2,function3,function4),  xtol=1.0e-5, maxiter=400)
#             VArray[x,y] = intersect
#     plt.pcolormesh(TArray,EArray, VArray, shading='auto')
#     clb = plt.colorbar()
#     clb.set_label('Surface Potential (V)', fontsize=15)
#     plt.xlabel('Temperature (k)', fontsize=15)
#     plt.xticks(fontsize=12)
#     plt.ylabel('Energy Difference (V)', fontsize=15)
#     plt.yticks(fontsize=12)
#     plt.title('Equilibrium Potential as a Function of Temperature and Energy Difference', fontsize=17)
#     plt.show()
#
#     return



# def NbvsPb(nbStart,nbEnd,pbRes,E,T,Nt,function1,function2,function3,function4):
#     Dnb = 0.1
#     Dpb = 0.01
#     nbLen = int((nbEnd-nbStart) / Dnb)
#     pbLen = int((pbRes) / Dpb)
#     nbArray = np.logspace(nbStart, nbEnd, nbLen)
#     pbArray = np.logspace(nbStart + 3, nbEnd + pbRes, pbLen)
#     VArray = np.zeros((nbLen, pbLen))
#
#     for x in range(nbLen):
#         nb = nbArray[x]
#
#         for y in range(pbLen):
#             pb = pbArray[y]
#             intersect = scipy.optimize.bisect(intersection, 0, 10, args=(E,T,Nt,nb,pb,function1,function2,function3,function4),  xtol=1.0e-5, maxiter=400)
#             V = intersect
#             VArray[x,y] = V
#
#     print(VArray[25,25])
#     plt.pcolormesh(pbArray,nbArray,VArray, shading = 'auto')
#     clb2 = plt.colorbar()
#     clb2.set_label('Surface Potential (V)')
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.xlabel(r'$p_b$',fontsize = 16)
#     plt.xticks(fontsize = 12)
#     plt.ylabel(r'$n_b$', fontsize = 16)
#     plt.yticks(fontsize = 12)
#     plt.title(r'Surface Potential as a Function of $p_b$ and $n_b$', fontsize = 17)
#     plt.show()
#     return

def VvsPb(Resolution,pbStart,E,T,Nt,function1,function2,function3,function4):
    VArray = np.array([])
    pbArray = np.array([])

    for x in range(Resolution):
        E = -0.2
        pb = 10 ** (pbStart + 0.1 * x)
        nb = pb/1000 # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        pbArray = np.append(pbArray,pb)


    VArray2 = np.array([])
    pbArray2 = np.array([])

    E = 0.0
    for x in range(Resolution):
        pb = 10 ** (pbStart + 0.1 * x)
        nb = pb/1000 # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray2 = np.append(VArray2, Intersect)
        pbArray2 = np.append(pbArray2,pb)

    VArray3 = np.array([])
    pbArray3 = np.array([])

    E = 0.2
    for x in range(Resolution):
        pb = 10 ** (pbStart + 0.1 * x)
        nb = pb/1000 # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray3 = np.append(VArray3, Intersect)
        pbArray3 = np.append(pbArray3,pb)

    VArray4 = np.array([])
    pbArray4 = np.array([])

    E = 0.4
    for x in range(Resolution):
        pb = 10 ** (pbStart + 0.1 * x)
        nb = pb / 1000  # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20,args=(E, T, Nt, nb, pb, function1, function2, function3, function4),xtol=1.0e-5, maxiter=400)
        VArray4 = np.append(VArray4, Intersect)
        pbArray4 = np.append(pbArray4, pb)

    VArray5 = np.array([])
    pbArray5 = np.array([])

    E = 0.6
    for x in range(Resolution):
        pb = 10 ** (pbStart + 0.1 * x)
        nb = pb / 1000  # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20,args=(E, T, Nt, nb, pb, function1, function2, function3, function4),xtol=1.0e-5, maxiter=400)
        VArray5 = np.append(VArray5, Intersect)
        pbArray5 = np.append(pbArray5, pb)

    VArray0 = np.array([])
    pbArray0 = np.array([])

    E = -0.4
    for x in range(Resolution):
        pb = 10 ** (pbStart + 0.1 * x)
        nb = pb / 1000  # This keeps the ratio between nb and pb the same

        Intersect = scipy.optimize.bisect(intersection, 0, 20,args=(E, T, Nt, nb, pb, function1, function2, function3, function4),xtol=1.0e-5, maxiter=400)
        VArray0 = np.append(VArray0, Intersect)
        pbArray0 = np.append(pbArray0, pb)


    plt.plot(pbArray0, VArray0, label='E = -0.4')
    plt.plot(pbArray, VArray, label='E = -0.2')
    plt.plot(pbArray2, VArray2, label='E = 0.0')
    plt.plot(pbArray3, VArray3, label='E = 0.2')
    plt.plot(pbArray4, VArray4, label='E = 0.4')
    plt.plot(pbArray5, VArray5, label='E = 0.6')

    plt.legend(fontsize = 15)
    plt.xscale('log')
    plt.xlabel(r'$p_b$', fontsize = 15)
    plt.ylabel('Surface Potential (V)', fontsize = 15)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.title(r'Cross sections of Electron Density $(p_b)$ and Surface Potential $(E)$', fontsize = 17)
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
    clb.ax.tick_params(labelsize=34)
    plt.xscale('log')
    plt.xlabel(r'$p_b$',fontsize = 40)
    plt.xticks(fontsize = 24)
    plt.ylabel('Temperature (K)', fontsize = 28)
    plt.yticks(yticks,fontsize = 24)
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
    clb.set_label('Surface Potential (V)', size = 28)
    clb.ax.tick_params(labelsize=24)
    plt.xscale('log')
    plt.xlabel(r'$p_b$',fontsize = 28)
    plt.xticks(fontsize = 24)
    plt.ylabel('Energy Difference', fontsize = 28)
    plt.yticks(fontsize = 24)
    plt.show()
    return

def VvsT(Resolution,TStart,E,Nt,nb,pb,function1,function2,function3,function4):
    #manually adjusting E and graphing it. E starts at -0.2, graphed, then repeated in iterations of 0.2
    E = -0.2
    TArray = np.array([])
    VArray = np.array([])
    for x in range(Resolution):
        T = TStart + x  # This will be set to start at 100 and end at 600 to match the cross in the later EvsT heat map
        Intersect = scipy.optimize.bisect(intersection, 0, 20,
                                          args=(E, T, Nt, nb, pb, function1, function2, function3, function4),
                                          xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        TArray = np.append(TArray, T)
    plt.plot(TArray, VArray, label='E = -0.2')

    E = 0.0
    TArray = np.array([])
    VArray = np.array([])
    for x in range(Resolution):
        T = TStart + x  # This will be set to start at 100 and end at 600 to match the cross in the later EvsT heat map
        Intersect = scipy.optimize.bisect(intersection, 0, 20,
                                          args=(E, T, Nt, nb, pb, function1, function2, function3, function4),
                                          xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        TArray = np.append(TArray, T)
    plt.plot(TArray, VArray, label='E = 0.0')

    E = 0.2
    TArray = np.array([])
    VArray = np.array([])
    for x in range(Resolution):
        T = TStart + x #This will be set to start at 100 and end at 600 to match the cross in the later EvsT heat map
        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        TArray = np.append(TArray,T)
    plt.plot(TArray,VArray, label = 'E = 0.2')

    E = 0.4
    TArray = np.array([])
    VArray = np.array([])
    for x in range(Resolution):
        T = TStart + x #This will be set to start at 100 and end at 600 to match the cross in the later EvsT heat map
        Intersect = scipy.optimize.bisect(intersection, 0, 20, args=(E,T,Nt,nb,pb,function1,function2,function3,function4), xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        TArray = np.append(TArray,T)
    plt.plot(TArray,VArray, label = 'E = 0.4')

    E = 0.6
    TArray = np.array([])
    VArray = np.array([])
    for x in range(Resolution):
        T = TStart + x  # This will be set to start at 100 and end at 600 to match the cross in the later EvsT heat map
        Intersect = scipy.optimize.bisect(intersection, 0, 20,
                                          args=(E, T, Nt, nb, pb, function1, function2, function3, function4),
                                          xtol=1.0e-5, maxiter=400)
        VArray = np.append(VArray, Intersect)
        TArray = np.append(TArray, T)
    plt.plot(TArray, VArray, label='E = 0.6')


    plt.legend(fontsize = 15)
    plt.xlabel(r'$Temperature (k)$', fontsize = 30)
    plt.ylabel('Surface Potential (V)', fontsize = 30)
    plt.xticks(fontsize = 26)
    plt.yticks(fontsize = 26)
    plt.title(r'Relationship between Temperature and Surface Potential', fontsize = 35)
    plt.show()
    return

# These are the default parameters that will be used in each function
Nt = 1.0 * 10 ** -2  # Available surface states per unit area
E = 0.5
T = 293.15  # Temperature in kelvin (Used STP...need's updating)
pb = 1.0 * 10 ** 15
nb = 1.0 * 10 ** 10
V = np.arange(-0.11, 1, 0.0001)

# Execute functions here:

# VvsNb(T,5,E,T,Nt,Qss,Qsc,ohm,LD) # nb start is the power of nb. So nbStart = 7 would set nb to 10^7
# EAndTMap(-0.5,1,100,600,Nt,nb,pb,Qss,Qsc,ohm,LD)
# NbvsPb(5,18,3,E,T,Nt,Qss,Qsc,ohm,LD)
# TvsPb(100,351,8,18,E,Nt,Qss,Qsc,ohm,LD)
# EvsPb(-0.5,0.5,8,18,T,Nt,Qss,Qsc,ohm,LD)
VvsPb(150,6,E,T,Nt,Qss,Qsc,ohm,LD)
VvsT(500,100,0.5,Nt,nb,pb,Qss,Qsc,ohm,LD) #VvsT(Resolution,TStart,E,Nt,nb,pb,function1,function2,function3,function4)