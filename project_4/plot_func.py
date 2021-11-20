import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt


print("Welcome to the plotting program for MCMC, please input a choice:")
print("[1]: Burn in time plot")

choice = input("Pick an option: ")
if choice == "1":
    path_to_file = "Burn_in_files/"
    eps = pa.mat()
    eps.load(path_to_file+"eps_burn_in_test_ordered_T1.txt")
    eps = np.array(eps[:,0])
    x = np.linspace(0,len(eps), len(eps))

    m = pa.mat()
    m.load(path_to_file+"m_burn_in_test_ordered_T1.txt")
    m = np.array(m[:,0])

    for i in range(len(eps)):
        eps[i] /= i+1
        m[i] /= i+1
    fig, ax = plt.subplots(ncols=2)
    ax[0].plot(x, eps, label = 'ϵ for T=1, ordered')
    ax[1].plot(x, m, label = '|m| for T=1, ordered')

    eps2 = pa.mat()
    eps2.load(path_to_file+"eps_burn_in_test_unordered_T1.txt")
    eps2 = np.array(eps2[:,0])

    m2 = pa.mat()
    m2.load(path_to_file+"m_burn_in_test_unordered_T1.txt")
    m2 = np.array(m2[:,0])

    for i in range(len(eps)):
        eps2[i] /= i+1
        m2[i] /= i+1
    ax[0].plot(x, eps2, label = 'ϵ for T=1, unordered')
    ax[1].plot(x, m2, label = '|m| for T=1, unordered')

    eps3 = pa.mat()
    eps3.load(path_to_file+"eps_burn_in_test_ordered_T24.txt")
    eps3 = np.array(eps3[:,0])

    m3 = pa.mat()
    m3.load(path_to_file+"m_burn_in_test_ordered_T24.txt")
    m3 = np.array(m3[:,0])

    for i in range(len(eps)):
        eps3[i] /= i+1
        m3[i] /= i+1
    ax[0].plot(x, eps3, label = 'ϵ for T=2.4, ordered')
    ax[1].plot(x, m3, label = '|m| for T=2.4, ordered')

    eps4 = pa.mat()
    eps4.load(path_to_file+"eps_burn_in_test_unordered_T24.txt")
    eps4 = np.array(eps4[:,0])

    m4 = pa.mat()
    m4.load(path_to_file+"m_burn_in_test_unordered_T24.txt")
    m4 = np.array(m4[:,0])

    for i in range(len(eps)):
        eps4[i] /= i+1
        m4[i] /= i+1
    ax[0].plot(x, eps4, label = 'ϵ for T=2.4, unordered')
    ax[1].plot(x, m4, label = '|m| for T=2.4, unordered')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].legend()
    ax[1].legend()
    ax[0].set_title('Burn in time for ϵ measured in MC cycles')
    ax[1].set_title('Burn in time for |m| measured in MC cycles')
    ax[0].set_xlabel('Cycles')
    ax[1].set_xlabel('Cycles')
    ax[0].grid()
    ax[1].grid()
    ax[0].set_ylabel('ϵ [J]')
    ax[1].set_ylabel('|m|')
    plt.show()

elif choice == '2':
    eps = pa.mat()
    eps2 = pa.mat()
    eps.load("approximate_eps_T1.txt")
    eps = np.array(eps[:,0])
    eps2.load("approximate_eps_T24.txt")
    eps2 = np.array(eps2[:,0])
    for i in eps:
        print(i)
    fig, ax = plt.subplots(ncols=2)
    ax[0].hist(eps, bins=20)
    ax[1].hist(eps2, bins=20)
    plt.show()



"""
chi = pa.mat()
chi.load("chi_L20.txt")
chi = np.array(chi[:,0])
cv = pa.mat()
cv.load("Cv_L20.txt")
cv = np.array(cv[:,0])
eps = pa.mat()
eps.load("expect_eps_L20.txt")
eps = np.array(eps[:,0])
T = np.linspace(2.1, 2.4, len(eps))

plt.plot(T, chi)
plt.plot(T, cv)
plt.plot(T, eps)
plt.show()
"""