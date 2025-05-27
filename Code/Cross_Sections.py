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








#%%%%%%%%%%%%%%%%%%%%%%%%
# Generating the data
VVals = np.arange(-0.2,0.5,1e-5)

#Graph 1: (Qss, Qsc0 vs V)
QssVals = Qss(VVals,E0,T0,Nt)
QscVals = Qsc(VVals,T0,eps,nb0,pb0)

# Graph 2: (V vs T)
VT_Vals = []
TVals = np.arange(100,401,1)
for i in TVals:
    V_Val = VvsT(i)
    VT_Vals.append(V_Val)

# Graph 3: (V vs E)
VE_Vals = []
EVals = np.linspace(-0.4,1.0,1000)
for i in EVals:
    V_Val = VvsE(i)
    VE_Vals.append(V_Val)

# Graph 4: (V vs Pb)
VPb_Vals = []
PbVals = np.logspace(6,23,1000)
for i in PbVals:
    V_Val = VvsPb(i,1e3)
    VPb_Vals.append(V_Val)

# Graph 5: (V vs Nb)
VNb_Vals = []
NbVals = np.logspace(6,23,10000)
for i in PbVals:
    V_Val = VvsPb(i,1e3)
    VNb_Vals.append(V_Val)

# Graph 6: (Cross secctions of E vs T)
VE_data1, VE_data2, VE_data3, VE_data4, VE_data5, VE_data6 = [], [], [], [], [], []
EVals = np.linspace(-0.4,1.0,1000)

T0 = 100
for i in EVals:
    V_Val = VvsE(i)
    VE_data1.append(V_Val)

T0 = 150
for i in EVals:
    V_Val = VvsE(i)
    VE_data2.append(V_Val)

T0 = 200
for i in EVals:
    V_Val = VvsE(i)
    VE_data3.append(V_Val)

T0 = 250
for i in EVals:
    V_Val = VvsE(i)
    VE_data4.append(V_Val)
    
T0 = 300
for i in EVals:
    V_Val = VvsE(i)
    VE_data5.append(V_Val)

T0 = 350
for i in EVals:
    V_Val = VvsE(i)
    VE_data6.append(V_Val)

T0 = 293.15 # Setting T0 back to it's original value

# Graph 7: (Cross secctions of Pb vs T)
VPBT_data1, VPBT_data2, VPBT_data3, VPBT_data4, VPBT_data5, VPBT_data6 = [], [], [], [], [], []
PbVals = np.logspace(6,23,1000)

T0 = 100
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBT_data1.append(V_Val)

T0 = 150
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBT_data2.append(V_Val)

T0 = 200
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBT_data3.append(V_Val)

T0 = 250
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBT_data4.append(V_Val)
    
T0 = 300
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBT_data5.append(V_Val)

T0 = 350
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBT_data6.append(V_Val)

T0 = 293.15 # Setting T0 back to it's original value

# Graph 8: (Cross sections of Pb vs E)

VPBE_data1, VPBE_data2, VPBE_data3, VPBE_data4, VPBE_data5, VPBE_data6 = [], [], [], [], [], []
PbVals = np.logspace(6,23,1000)

E0 = -0.4
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBE_data1.append(V_Val)

E0 = -0.2
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBE_data2.append(V_Val)

E0 = 0 
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBE_data3.append(V_Val)

E0 = 0.2
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBE_data4.append(V_Val)
    
E0= 0.4
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBE_data5.append(V_Val)

E0= 0.6
for i in PbVals:
    V_Val = VvsPb(i,1000)
    VPBE_data6.append(V_Val)

E0 = 0.2 # Resetting E0

# Graph 9: (Cross sections of T and Nb)
VNbT_data1, VNbT_data2, VNbT_data3, VNbT_data4, VNbT_data5, VNbT_data6 = [], [], [], [], [], []
NbVals = np.logspace(6,23,1000)

T0 = 100
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbT_data1.append(V_Val)

T0 = 150
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbT_data2.append(V_Val)

T0 = 200
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbT_data3.append(V_Val)

T0 = 250
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbT_data4.append(V_Val)
    
T0 = 300
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbT_data5.append(V_Val)

T0 = 350
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbT_data6.append(V_Val)

T0 = 293.15 # Setting T0 back to it's original value

# Graph 10: (Cross sections of E and Nb)
VNbE_data1, VNbE_data2, VNbE_data3, VNbE_data4, VNbE_data5, VNbE_data6 = [], [], [], [], [], []
NbVals = np.logspace(6,23,1000)

E0 = -0.4
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbE_data1.append(V_Val)

E0 = -0.2
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbE_data2.append(V_Val)

E0 = 0.0
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbE_data3.append(V_Val)

E0 = 0.2
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbE_data4.append(V_Val)
    
E0 = 0.4
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbE_data5.append(V_Val)

E0 = 0.6
for i in NbVals:
    V_Val = VvsNb(i,1000)
    VNbE_data6.append(V_Val)

E = 0.2

# Graph 11 (Nb vs Pb cross sections)
VNbPb_data1, VNbPb_data2, VNbPb_data3, VNbPb_data4, VNbPb_data5 = [], [], [], [], []
NbVals = np.logspace(6,18,1000)

pb0 = 1e9
for i in NbVals:
    V_Val = VvsNb(i,pb0/i)
    VNbPb_data1.append(V_Val)

pb0 = 1e12
for i in NbVals:
    V_Val = VvsNb(i,pb0/i)
    VNbPb_data2.append(V_Val)

pb0 = 1e15
for i in NbVals:
    V_Val = VvsNb(i,pb0/i)
    VNbPb_data3.append(V_Val)

pb0 = 1e18
for i in NbVals:
    V_Val = VvsNb(i,pb0/i)
    VNbPb_data4.append(V_Val)
    
pb0 = 1e21
for i in NbVals:
    V_Val = VvsNb(i,pb0/i)
    VNbPb_data5.append(V_Val)

#%%%%%%%%%%%%%%%%%%%%%%%%
# Generating Graph

# Graph 1
fig1 = plt.figure(figsize = (10,10))
ax1 = plt.subplot()
ax1.plot(VVals,QssVals, c = "orange", label = "Qss")
ax1.plot(VVals,QscVals, c = "blue", label = "Qsc")
ax1.set_yscale("log")
ax1.set_xlabel("Surface Potential (V)")
ax1.set_ylabel("Charge Densities (C/cm)")
fig1.savefig("Code\Figures\Intersection_Graph.png")
ax1.legend()

# Graph 2
fig2 = plt.figure(figsize = (10,10))
ax2 = plt.subplot()
ax2.plot(TVals,VT_Vals)
ax2.set_title("Surface Voltage vs Temperature")
ax2.set_xlabel("Temperature (k)")
ax2.set_ylabel("Surface Voltage (V)")
fig2.savefig("Code\Figures\Voltage_Temp.png")

# Graph 3
fig3 = plt.figure(figsize = (10,10))
ax3 = plt.subplot()
ax3.plot(EVals, VE_Vals)
ax3.set_title("Surface Voltage vs Energy Difference")
ax3.set_xlabel("Energy (J)")
ax3.set_ylabel("Surface Voltage (V)")
fig3.savefig("Code\Figures\Voltage_E.png")

# Graph 4
fig4 = plt.figure(figsize = (10,10))
ax4 = plt.subplot()
ax4.set_xscale("log")
ax4.plot(PbVals, VPb_Vals)
ax4.set_title("Surface Voltage vs Hole Density")
ax4.set_xlabel("Hole Density (C/cm)")
ax4.set_ylabel("Surface Voltage (V)")
fig4.savefig("Code\Figures\V_Pb.png")

# Graph 5
fig5 = plt.figure(figsize = (10,10))
ax5 = plt.subplot()
ax5.set_xscale("log")
ax5.plot(NbVals, VNb_Vals)
ax5.set_title("Surface Voltage vs Carrier Density")
ax5.set_xlabel("Carrier Density (C/cm)")
ax5.set_ylabel("Surface Voltage (V)")
fig5.savefig("Code\Figures\V_Nb.png")

# Graph 6

fig6 = plt.figure(figsize = (10,10))
ax6 = plt.subplot()
ax6.set_title("V vs E")
ax6.set_ylabel("Surface Voltage")
ax6.set_xlabel("Energy Difference")
ax6.plot(EVals,VE_data1, label = "T = 100")
ax6.plot(EVals,VE_data2, label = "T = 150")
ax6.plot(EVals,VE_data3, label = "T = 200")
ax6.plot(EVals,VE_data4, label = "T = 250")
ax6.plot(EVals,VE_data5, label = "T = 300")
ax6.plot(EVals,VE_data6, label = "T = 350")
fig6.savefig("Code\Figures\E_T_Cross_Sections.png")

# Graph 7

fig7 = plt.figure(figsize = (10,10))
ax7 = plt.subplot()
ax7.set_title("V vs Pb")
ax7.set_ylabel("Surface Voltage")
ax7.set_xlabel("Hole Density")
ax7.set_xscale("log")
ax7.plot(PbVals,VPBT_data1, label = "T = 100")
ax7.plot(PbVals,VPBT_data2, label = "T = 150")
ax7.plot(PbVals,VPBT_data3, label = "T = 200")
ax7.plot(PbVals,VPBT_data4, label = "T = 250")
ax7.plot(PbVals,VPBT_data5, label = "T = 300")
ax7.plot(PbVals,VPBT_data6, label = "T = 350")
ax7.legend()
fig7.savefig("Code\Figures\Pb_T_Cross_Sections.png")
# Graph 8

fig8 = plt.figure(figsize = (10,10))
ax8 = plt.subplot()
ax8.set_title("V vs Pb")
ax8.set_ylabel("Surface Voltage")
ax8.set_xlabel("Hole Density")
ax8.set_xscale("log")
ax8.plot(PbVals,VPBE_data1, label = "E = -0.4")
ax8.plot(PbVals,VPBE_data2, label = "E = -0.2")
ax8.plot(PbVals,VPBE_data3, label = "E = 0.0")
ax8.plot(PbVals,VPBE_data4, label = "E = 0.2")
ax8.plot(PbVals,VPBE_data5, label = "E = 0.4")
ax8.plot(PbVals,VPBE_data6, label = "E = 0.6")
ax8.legend()
fig8.savefig("Code\Figures\Pb_E_Cross_Sections.png")

# Graph 9
fig9 = plt.figure(figsize = (10,10))
ax9 = plt.subplot()
ax9.set_xscale("log")
ax9.set_title("V vs Nb")
ax9.set_ylabel("Surface Voltage")
ax9.set_xlabel("Carrier Denisty")
ax9.plot(NbVals,VNbT_data1, label = "T = 100")
ax9.plot(NbVals,VNbT_data2, label = "T = 150")
ax9.plot(NbVals,VNbT_data3, label = "T = 200")
ax9.plot(NbVals,VNbT_data4, label = "T = 250")
ax9.plot(NbVals,VNbT_data5, label = "T = 300")
ax9.plot(NbVals,VNbT_data6, label = "T = 350")
ax9.legend()
fig9.savefig("Code\Figures\Carrier_T_Cross_Sections.png")

# Graph 10

fig10 = plt.figure(figsize = (10,10))
ax10 = plt.subplot()
ax10.set_xscale("log")
ax10.set_title("V vs Nb")
ax10.set_ylabel("Surface Voltage")
ax10.set_xlabel("Carrier Density")
ax10.plot(NbVals,VNbE_data1, label = "E = -0.4")
ax10.plot(NbVals,VNbE_data2, label = "E = -0.2")
ax10. plot(NbVals,VNbE_data3, label = "E = 0")
ax10. plot(NbVals,VNbE_data4, label = "E = 0.2")
ax10. plot(NbVals,VNbE_data5, label = "E = 0.4")
ax10. plot(NbVals,VNbE_data6, label = "E = 0.6")
ax10.legend()
fig10.savefig("Code\Figures\Carrier_E_Cross_Sections.png")

# Graph 11

fig11 = plt.figure(figsize = (10,10))
ax11 = plt.subplot()
ax11.set_xscale("log")
ax11.set_title("V vs Nb")
ax11.set_ylabel("Surface Voltage")
ax11.set_xlabel("Carrier Density")
ax11. plot(NbVals,VNbPb_data1, label = "Pb = 1e9")
ax11. plot(NbVals,VNbPb_data2, label = "Pb = 1e12")
ax11. plot(NbVals,VNbPb_data3, label = "Pb = 1e15")
ax11. plot(NbVals,VNbPb_data4, label = "Pb = 1e18")
ax11. plot(NbVals,VNbPb_data5, label = "Pb = 1e21")
ax11.legend()
fig11.savefig("Code\Figures\Carrier_Pb_Cross_Sections.png")
# plt.show()