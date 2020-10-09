import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

cm = 0.01
mm = 0.001
ft = 0.3048
# file = 'methane_pox_on_pt_2.cti'
# composition = 'CH4:1, O2:1.5, AR:0.1'
# phase_cat = 'Pt_surf'
# T = 800 + 273.15
# P = ct.one_atm
# l = 0.3 * cm
# A = 1*cm**2
# Sb = 1000 / cm
# Q = 2/3 * cm**3
# e = 0.3

# file = 'khai_test_cata_2.cti'
# composition = 'CH4:2, O2:1'
# phase_cat = 'Ni_surf'
# T = 1000
# P = ct.one_atm
# l = 1
# d = 1 * ft
# A = np.pi * d**2 /4
# Sb = 0
# Q = 0.1 * A
# e = 1

# file = 'khai_test_cata_2.cti'
# composition = 'CH4:1.33, O2:0.81, N2:97.86'#'CH4:2, O2:1'
# phase_cat = 'Ni_surf'
# T = 800 + 273.15
# P = ct.one_atm
# l = 0.3 * cm
# A = 1*cm**2
# Sb = 1000 / cm
# Q = 2/3 * cm**3
# e = 0.3

file = 'Ni_surf_mech.cti'
composition = 'CH4:1.33, O2:0.81, N2:97.86'#'CH4:2, CO2:2, N2:96' #'CH4:1, O2:1.5, AR:0.1' #'CH4:1.33, O2:0.81, N2:97.86'
phase_cat = 'Ni_surf'   # Phase of catalyst
T = 973         # studying temperature
t0 = 298.15     # initial temperature of gas mixture
P = 0.986923267* ct.one_atm # studying pressure 1 bar
l = 27 * mm     # length of reactor
d = 10 * mm     # diameter of reactor
A = np.pi * d**2 / 4
Spv = 9.85e6     # surface area per unit volume
Q = 200/3 * cm**3   # flow rate of volume
e = 0.42         # porosity

A_eq = A*e      # equivalent area
N = 1000         # Number of step

gas1 = ct.Solution(file,'gas')
gas1.TPX = t0,ct.one_atm,composition # initial conditions for calculate the mass flow rate

gas1()
m_dot = Q*gas1.density_mass
print('mass flow rate',m_dot)
gas1.TPX = T,P,composition
sur1 = ct.Interface(file,phase_cat,[gas1])
sur1.TP = T,P

r1 = ct.IdealGasConstPressureReactor(gas1,energy = 'on')
#r1.volume = A_eq * l    # reacting volume
#print(r1.volume)
rsur1 = ct.ReactorSurface(sur1,r1, A = l * A_eq * Spv)
sim1 = ct.ReactorNet([r1])

TDY = gas1.TDY
cov = rsur1.coverages
gas1.TDY = TDY
gas1()
print(sur1.report())
#t_total1 = l * A_eq / Q    # Estimate the residence time
t_total1 = r1.mass / m_dot
print(t_total1)
dt = t_total1 / N
t1 = (np.arange(N) + 1) * dt
z1 = np.zeros_like(t1)
u1 = np.zeros_like(t1)
state1 = ct.SolutionArray(r1.thermo)
NI = np.zeros_like(z1)
O = np.zeros_like(z1)
CO = np.zeros_like(z1)
H = np.zeros_like(z1)
C = np.zeros_like(z1)
H2O = np.zeros_like(z1)
OH = np.zeros_like(z1)
CO2 = np.zeros_like(z1)
for n1, t_i in enumerate(t1):
    # perform time integration
    sim1.advance(t_i)
    # compute velocity and transform into space
    u1[n1] = m_dot / A_eq / r1.thermo.density
    z1[n1] = z1[n1 - 1] + u1[n1] * dt
    state1.append(r1.thermo.state)
    NI[n1] = sur1['NI(S)'].X*100
    O[n1] = sur1['O(S)'].X*100
    CO[n1] = sur1['CO(S)'].X*100
    H[n1] = sur1['H(S)'].X*100
    C[n1] = sur1['C(S)'].X*100
    H2O[n1] = sur1['H2O(S)'].X*100
    OH[n1] = sur1['OH(S)'].X*100
    CO2[n1] = sur1['CO2(S)'].X*100
print(u1)
plt.figure()
plt.subplot(121)
plt.plot(z1*100, state1.X[:, gas1.species_index('CH4')]*100, label='CH4')
plt.plot(z1*100, state1.X[:, gas1.species_index('O2')]*100, label=' O2')
plt.plot(z1*100, state1.X[:, gas1.species_index('CO2')]*100, label=' CO2')
plt.plot(z1*100, state1.X[:, gas1.species_index('CO')]*100, label='CO')
plt.plot(z1*100, state1.X[:, gas1.species_index('H2O')]*100, label='H2O')
plt.plot(z1*100, state1.X[:, gas1.species_index('H2')]*100, label='H2')
plt.title('Lagrangian Particle ')
plt.xlabel('$z$ [cm]')
plt.ylabel('Molar fraction %')
plt.legend(loc = 0)
#plt.show()

#plt.figure()
plt.subplot(122)
plt.plot(t1, state1.X[:, gas1.species_index('CH4')]*100, label='CH4')
plt.plot(t1, state1.X[:, gas1.species_index('O2')]*100, label=' O2')
plt.plot(t1, state1.X[:, gas1.species_index('CO2')]*100, label=' CO2')
plt.plot(t1, state1.X[:, gas1.species_index('CO')]*100, label='CO')
plt.plot(t1, state1.X[:, gas1.species_index('H2O')]*100, label='H2O')
plt.plot(t1, state1.X[:, gas1.species_index('H2')]*100, label='H2')
plt.title('Lagrangian Particle ')
plt.xlabel('$t$ [s]')
plt.ylabel('Molar fraction %')
plt.legend(loc = 0)
plt.show()
# plt.figure()
# plt.plot(z1*100,NI,label='NI(S)')
# plt.plot(z1*100,O,label='O(S)')
# plt.plot(z1*100,CO,label='CO(S)')
# plt.plot(z1*100,H,label='H(S)')
# plt.xlabel('axial position [cm]')
# plt.ylabel('Surface coverage')
# plt.legend(loc=0)
# plt.show()

# plt.figure()
# plt.plot(z1*100,C,label='C(S)')
# plt.plot(z1*100,H2O,label='H2O(S)')
# plt.plot(z1*100,CO2,label='CO2(S)')
# plt.plot(z1*100,OH,label='OH(S)')
# plt.xlabel('axial position [cm]')
# plt.ylabel('Surface coverage')
# plt.legend(loc=0)
# plt.show()
