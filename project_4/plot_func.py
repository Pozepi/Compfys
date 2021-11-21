import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt


print("Welcome to the plotting program for MCMC, please input a choice:")
print("[1]: Burn in time plot")
print("[2]")

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
    """for i in eps:
        print(i)"""
    print(eps)
    fig, ax = plt.subplots(ncols=2)
    ax[0].hist(eps, bins=20)
    ax[1].hist(eps2, bins=20)
    plt.show()

elif choice == '3':
    eps = pa.mat()
    m = pa.mat()
    chi = pa.mat()
    Cv = pa.mat()
    temp = pa.mat()

    eps.load("expect_eps_L2.txt")
    m.load("expect_m_L2.txt")
    chi.load("chi_L2.txt")
    Cv.load("Cv_L2.txt")
    temp.load("Temp_L2.txt")

    eps = np.array(eps[:,0])
    m = np.array(m[:,0])
    chi = np.array(chi[:,0])
    Cv = np.array(Cv[:,0])
    temp = np.array(temp[:,0])

    N = 4
    T = np.linspace(2.1, 2.4, 100) # in units of J/k_B
    beta = 1/T

    Z = 2*np.exp(-8*beta) + 2*np.exp(8*beta) + 12

    E1 = 16/Z*(np.exp(-8*beta) - np.exp(8*beta)) # first moment
    E2 = 128/Z*(np.exp(8*beta) + np.exp(-8*beta))
    Cv_a = (E2 - E1**2)/(N*T**2)

    M1 = (8*np.exp(8*beta) + 16)/Z # (4*np.exp(8*beta) + 4*np.exp(-8*beta) + 16)/Z
    M2 = 32/Z*(np.exp(8*beta) + 1) #16/Z*(np.exp(8*beta) + 1)
    chi_a = (M2 - M1**2)/(N*T)

    eps_a = E1/N
    m_a = M1/N

    fig, ax = plt.subplots(2,2)
    ax[0,0].plot(1/temp, eps, label = 'Numerical <ϵ>')
    ax[0,0].plot(beta, eps_a, label = 'Analytical <ϵ>')
    ax[0,0].grid()
    ax[0,0].legend()
    ax[0,1].plot(1/temp, m, label = 'Numerical <|m|>')
    ax[0,1].plot(beta, m_a, label = 'Analytical <|m|>')
    ax[0,1].grid()
    ax[0,0].legend()
    ax[1,0].plot(1/temp, chi, label = 'Numerical χ')
    ax[1,0].plot(beta, chi_a, label = 'Analytical χ')
    ax[1,0].grid()
    ax[1,0].legend()
    ax[1,1].plot(1/temp, Cv, label = 'Numerical Cv')
    ax[1,1].plot(beta, Cv_a, label = 'Analytical Cv')
    ax[1,1].grid()
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