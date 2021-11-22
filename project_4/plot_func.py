import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt


print("Welcome to the plotting program for MCMC, please input a choice:")
print("[1]: Burn in time plot")
print("[2]: Oppgave 6")
print("[3]: Analytical vc numerical")

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
    ax[0].plot(x, eps, label = '<ϵ> for T=1, ordered')
    ax[1].plot(x, m, label = '<|m|> for T=1, ordered')

    eps2 = pa.mat()
    eps2.load(path_to_file+"eps_burn_in_test_unordered_T1.txt")
    eps2 = np.array(eps2[:,0])

    m2 = pa.mat()
    m2.load(path_to_file+"m_burn_in_test_unordered_T1.txt")
    m2 = np.array(m2[:,0])

    for i in range(len(eps)):
        eps2[i] /= i+1
        m2[i] /= i+1
    ax[0].plot(x, eps2, label = '<ϵ> for T=1, unordered')
    ax[1].plot(x, m2, label = '<|m|> for T=1, unordered')

    eps3 = pa.mat()
    eps3.load(path_to_file+"eps_burn_in_test_ordered_T24.txt")
    eps3 = np.array(eps3[:,0])

    m3 = pa.mat()
    m3.load(path_to_file+"m_burn_in_test_ordered_T24.txt")
    m3 = np.array(m3[:,0])

    for i in range(len(eps)):
        eps3[i] /= i+1
        m3[i] /= i+1
    ax[0].plot(x, eps3, label = '<ϵ> for T=2.4, ordered')
    ax[1].plot(x, m3, label = '<|m|> for T=2.4, ordered')

    eps4 = pa.mat()
    eps4.load(path_to_file+"eps_burn_in_test_unordered_T24.txt")
    eps4 = np.array(eps4[:,0])

    m4 = pa.mat()
    m4.load(path_to_file+"m_burn_in_test_unordered_T24.txt")
    m4 = np.array(m4[:,0])

    for i in range(len(eps)):
        eps4[i] /= i+1
        m4[i] /= i+1
    ax[0].plot(x, eps4, label = '<ϵ> for T=2.4, unordered')
    ax[1].plot(x, m4, label = '<|m|> for T=2.4, unordered')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].legend()
    ax[1].legend()
    ax[0].set_title('Burn in time for <ϵ> measured in MC cycles')
    ax[1].set_title('Burn in time for <|m|> measured in MC cycles')
    ax[0].set_xlabel('Cycles')
    ax[1].set_xlabel('Cycles')
    ax[0].grid()
    ax[1].grid()
    ax[0].set_ylabel('<ϵ> [J]')
    ax[1].set_ylabel('<|m|>')

    plt.savefig('burn_in.pdf')
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
    #print(eps)
    fig, ax = plt.subplots(figsize=(11,6), ncols=2)
    #ax[0].plot(eps)
    #ax[1].plot(eps2)
    ax[0].hist(eps, bins='auto')
    ax[0].set_title(r'$T = 1 \; J/k_B$')
    ax[1].hist(eps2, bins='auto')
    ax[1].set_title(r'$T = 2.4 \; J/k_B$')
    [[axi.set_ylabel('Count'), axi.set_xlabel(r"$\epsilon$ [J]"), axi.grid()] for axi in ax]
    plt.tight_layout()
    plt.savefig('eps_dist.pdf')
    plt.show()

elif choice == '3':
    eps = pa.mat()
    m = pa.mat()
    chi = pa.mat()
    Cv = pa.mat()
    temp = pa.mat()
    eps2 = pa.mat()
    m2 = pa.mat()
    chi2 = pa.mat()
    Cv2 = pa.mat()

    eps.load("eps.txt")
    m.load("m.txt")
    chi.load("chi.txt")
    Cv.load("Cv.txt")
    temp.load("Temp.txt")

    eps2.load("eps50.txt")
    m2.load("m50.txt")
    chi2.load("chi50.txt")
    Cv2.load("Cv50.txt")

    eps = np.array(eps[:,0])
    m = np.array(m[:,0])
    chi = np.array(chi[:,0])
    Cv = np.array(Cv[:,0])
    temp = np.array(temp[:,0])

    eps2 = np.array(eps2[:,0])
    m2 = np.array(m2[:,0])
    chi2 = np.array(chi2[:,0])
    Cv2 = np.array(Cv2[:,0])


    N = 4
    T = temp #np.logspace(-1, 1, 100) # in units of J/k_B
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

    fig, ax = plt.subplots(figsize=(11,6), 2,2)
    plt.suptitle('Analytical vs numerical, 10 000 and 250 000 MC cycles')
    ax[0,0].plot(beta, eps_a, label = 'Analytical <ϵ>')
    ax[0,0].plot(1/temp, eps, ls = '--',label = 'Numerical <ϵ> 250 000 cycles')
    ax[0,0].plot(1/temp, eps2, ls = '--',label = 'Numerical <ϵ> 10 000 cycles')
    ax[0,0].grid()
    ax[0,0].legend()
    ax[0,0].set_xlabel('β = $1/k_BT$  $[\quad J^{-1}]$')
    ax[0,0].set_ylabel('<ϵ> [J]')
    ax[0,0].set_title('<ϵ> analytical vs numerical')

    ax[0,1].plot(beta, m_a, label = 'Analytical <|m|>')
    ax[0,1].plot(1/temp, m, ls = '--',label = 'Numerical <|m|> 250 000 cycles')
    ax[0,1].plot(1/temp, m2, ls = '--',label = 'Numerical <|m|> 10 000 cycles')
    ax[0,1].grid()
    ax[0,1].legend()
    ax[0,1].set_xlabel('β = $1/k_BT$  $[\quad J^{-1}]$')
    ax[0,1].set_ylabel('<|m|>')
    ax[0,1].set_title('<|m|> analytical vs numerical')

    ax[1,0].plot(beta, chi_a, label = 'Analytical χ')
    ax[1,0].plot(1/temp, chi, ls = '--',label = 'Numerical χ 250 000 cycles')
    ax[1,0].plot(1/temp, chi2, ls = '--',label = 'Numerical χ 10 000 cycles')
    ax[1,0].grid()
    ax[1,0].legend()
    ax[1,0].set_xlabel('β = $1/k_BT$  $[\quad J^{-1}]$')
    ax[1,0].set_ylabel('χ')
    ax[1,0].set_title('χ analytical vs numerical')

    ax[1,1].plot(beta, Cv_a, label = 'Analytical Cv')
    ax[1,1].plot(1/temp, Cv, ls = '--',label = 'Numerical Cv 250 000 cycles')
    ax[1,1].plot(1/temp, Cv, ls = '--',label = 'Numerical Cv 10 000 cycles')
    ax[1,1].grid()
    ax[1,1].legend()
    ax[1,1].set_xlabel('β = $1/k_BT$  $[\quad J^{-1}]$')
    ax[1,1].set_ylabel('Cv [Kb]')
    ax[1,1].set_title('Cv analytical vs numerical')

    #[[[axi.set_xscale('symlog')] for axi in axij] for axij in ax]
    plt.savefig('anavsnum.pdf')
    plt.show()

elif choice == '4':
    from scipy.ndimage import gaussian_filter1d

    files = ['_L20.txt', '_L40.txt', '_L60.txt', '_L80.txt', '_L100.txt']
    color = ['r', 'g', 'b', 'k', 'y']

    for file, c in zip(files, color):
        cvn = pa.mat()
        cvn.load('Cv'+file)
        t = pa.mat()
        t.load('Temp'+file)
        cvn = np.array(cvn)
        smooth =  gaussian_filter1d(cvn, 5, axis=0)
        plt.plot(t, cvn, c+'x', label='Numerical data L = '+file[2:-4])
        plt.plot(t, smooth, c, label='Gaussian filter L = '+file[2:-4])

    plt.legend()
    plt.savefig('data.pdf')
    plt.show()


    L = np.array([float(file[2:-4]) for file in files])
    T_c = []
    for file in files:
        cv = pa.mat()
        cv.load('Cv'+file)
        cv = np.array(cv)[:,0]
        smooth =  gaussian_filter1d(cv, 5)
        i = np.where(smooth==np.max(smooth))[0][0]
        T_c.append(cv[i])
        #t = pa.mat()
        #t.load('Temp'+file)
        #t = np.array(t)[:,0]
        #T_c.append(np.trapz(cv*t, t))

    from scipy import stats
    slope, intercept, r, p, se = stats.linregress(L, T_c*L)
    print(slope, intercept)
    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot(L, T_c*L, label='$L\cdot T_c$')
    ax.plot(L, L*slope + intercept, label='Linear fit')
    ax.legend()
    ax.grid()
    plt.savefig('linear_fit.pdf')
