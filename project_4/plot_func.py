import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

eps = pa.mat()
eps.load("eps_burn_in_test_ordered_T1.txt")
eps = np.array(eps[:,0])
x = np.linspace(0,len(eps), len(eps))
print(np.std(eps))

m = pa.mat()
m.load("m_burn_in_test_ordered_T1.txt")
m = np.array(m[:,0])
print(eps)
print(m)
for i in range(len(eps)):
    eps[i] /= i+1
    m[i] /= i+1
fig, ax = plt.subplots(ncols=2)
ax[0].plot(x, eps)
ax[1].plot(x, m)
plt.show()

eps = pa.mat()
eps.load("eps_burn_in_test_unordered_T1.txt")
eps = np.array(eps[:,0])
x = np.linspace(0,len(eps), len(eps))
print(np.std(eps))

m = pa.mat()
m.load("m_burn_in_test_unordered_T1.txt")
m = np.array(m[:,0])
print(eps)
print(m)
for i in range(len(eps)):
    eps[i] /= i+1
    m[i] /= i+1
fig, ax = plt.subplots(ncols=2)
ax[0].plot(x, eps)
ax[1].plot(x, m)
plt.show()

eps = pa.mat()
eps.load("eps_burn_in_test_ordered_T24.txt")
eps = np.array(eps[:,0])
x = np.linspace(0,len(eps), len(eps))
print(np.std(eps))

m = pa.mat()
m.load("m_burn_in_test_ordered_T24.txt")
m = np.array(m[:,0])
print(eps)
print(m)
for i in range(len(eps)):
    eps[i] /= i+1
    m[i] /= i+1
fig, ax = plt.subplots(ncols=2)
ax[0].plot(x, eps)
ax[1].plot(x, m)
plt.show()

eps = pa.mat()
eps.load("eps_burn_in_test_unordered_T24.txt")
eps = np.array(eps[:,0])
x = np.linspace(0,len(eps), len(eps))
print(np.std(eps))

m = pa.mat()
m.load("m_burn_in_test_unordered_T24.txt")
m = np.array(m[:,0])
print(eps)
print(m)
for i in range(len(eps)):
    eps[i] /= i+1
    m[i] /= i+1
fig, ax = plt.subplots(ncols=2)
ax[0].plot(x, eps)
ax[1].plot(x, m)
plt.show()