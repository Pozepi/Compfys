import numpy as np
import matplotlib.pyplot as plt

N = 4
T = np.linspace(1, 10, 100) # in units of J/k_B
beta = 1/T

Z = 2*np.exp(-4*beta) + 2*np.exp(4*beta) + 12

E1 = 8/Z*(np.exp(-4*beta) - np.exp(4*beta)) # first moment
E2 = 32/Z*(np.exp(4*beta) + np.exp(-4*beta))

Cv = (E2 - E1**2)/(N*T**2)

M1 = (4*np.exp(4*beta) + 4*np.exp(-4*beta) + 16)/Z
M2 = 16/Z*(np.exp(4*beta) + np.exp(-4*beta) + 2)

chi = (M2 - M1**2)/(4*T)

fig, (ax1, ax2) = plt.subplots(ncols=2)
ax1.plot(T, Cv, label='Heat capacity')
ax2.plot(T, chi, label='Magnetization thingie')

ax1.legend()
ax2.legend()
plt.show()

fig, (ax1, ax2) = plt.subplots(ncols=2)
ax1.plot(T, (E2-E1**2))
ax2.plot(T, (M2-M1**2))

ax1.legend()
ax2.legend()
plt.show()



""" find value of single temp"""
T = 1 # in units of J/k_B
beta = 1/T

Z = 2*np.exp(-4*beta) + 2*np.exp(4*beta) + 12

E1 = 8/Z*(np.exp(-4*beta) - np.exp(4*beta)) # first moment
E2 = 32/Z*(np.exp(4*beta) + np.exp(-4*beta))

print('Heat capacity')
print('First moment: ', E1)
print('Second moment: ', E2)
print('Std: ', (E2 - E1**2))

Cv = (E2 - E1**2)/(4*T**2)
print('Heat capacity: ', Cv)

M1 = (4*np.exp(4*beta) + 4*np.exp(-4*beta) + 16)/Z
M2 = 16/Z*(np.exp(4*beta) + np.exp(-4*beta) + 2)

chi = (M2 - M1**2)/(4*T)

print('')
print('Magnetisation')
print('First moment: ', M1)
print('Second moment: ', M2)
print('Std: ', (M2 - M1**2))

print('Magnetisation: ', chi)

