import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

"""Calc analytical stuff"""
N = 4
T = np.logspace(-2, 2, 100) # in units of J/k_B
beta = 1/T

Z = 2*np.exp(-8*beta) + 2*np.exp(8*beta) + 12

E1 = 16/Z*(np.exp(-8*beta) - np.exp(8*beta)) # first moment
E2 = 128/Z*(np.exp(8*beta) + np.exp(-8*beta))
Cv = (E2 - E1**2)/(N*T**2)

M1 = (8*np.exp(8*beta) + 16)/Z # (4*np.exp(8*beta) + 4*np.exp(-8*beta) + 16)/Z
M2 = 32/Z*(np.exp(8*beta) + 1) #16/Z*(np.exp(8*beta) + 1)
chi = (M2 - M1**2)/(N*T)

"""Load numerical stuff"""
chi_n = pa.mat()
Cv_n = pa.mat()
Temp_n = pa.mat()

chi_n.load('chi.txt')
Cv_n.load('Cv.txt')
Temp_n.load('Temp.txt')

chi_n = np.array(chi_n)
Cv_n = np.array(Cv_n)
Temp_n = np.array(Temp_n)

#print(beta)
#print(1/Temp_n)

fig, (ax1, ax2) = plt.subplots(ncols=2)
ax1.plot(1/Temp_n, Cv_n, label='Numerical heat capacity')
ax1.plot(beta, Cv, label='Analytical heat capacity', ls='--')

ax2.plot(1/Temp_n, chi_n, label='Numerical suseptibility')
ax2.plot(beta, chi, label='Analytical suseptibility', ls='--')

[[axi.grid(), axi.legend(), axi.set_xlim(0,7), axi.set_xlabel(r'$\beta$ [$J/k_B$]')] for axi in (ax1, ax2)]
ax1.set_ylim(0,0.5)
ax1.set_ylabel(r'$C_v$ [$k_B/T$]')

ax2.set_ylim(0,0.18)
ax2.set_ylabel(r'\chi [J]')
plt.savefig('T_L2')
plt.show()

"""
fig, (ax1, ax2) = plt.subplots(ncols=2)
ax1.plot(T, (E2-E1**2))
ax2.plot(T, (M2-M1**2))

ax1.legend()
ax2.legend()
plt.show()
"""


""" find value of single temp"""
T = 1 # in units of J/k_B
beta = 1/T

Z = 2*np.exp(-8*beta) + 2*np.exp(8*beta) + 12

E1 = 16/Z*(np.exp(-8*beta) - np.exp(8*beta)) # first moment
E2 = 128/Z*(np.exp(8*beta) + np.exp(-8*beta))

print('Heat capacity')
print('First moment: ', E1)
print('Second moment: ', E2)
print('Std: ', (E2 - E1**2))

Cv = (E2 - E1**2)/(4*T**2)
print('Heat capacity: ', Cv)

M1 = (8*np.exp(8*beta) + 16)/Z # (4*np.exp(8*beta) + 4*np.exp(-8*beta) + 16)/Z
M2 = 32/Z*(np.exp(8*beta) + 1) #16/Z*(np.exp(8*beta) + 1)

chi = (M2 - M1**2)/(4*T)

print('')
print('Magnetisation')
print('First moment: ', M1)
print('Second moment: ', M2)
print('Std: ', (M2 - M1**2))

print('Magnetisation: ', chi)

