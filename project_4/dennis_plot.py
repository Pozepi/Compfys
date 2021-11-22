import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
from scipy.ndimage import gaussian_filter1d


"""cvn = pa.mat()
cvn.load('Cv_L100.txt')
t = pa.mat()
t.load('Temp_L100.txt')

plt.plot(t,cvn,'rx')
cvn = np.array(cvn)
smooth = gaussian_filter1d(cvn, 3, axis=0)
plt.plot(t,smooth)
plt.show()"""

files = ['_L100.txt', '_L20.txt', '_L40.txt', '_L80.txt', 'L_60.txt']
color = ['r', 'g', 'b', 'k']

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
plt.show()

"""
Cv = pa.mat()
Cv.load("Cv_L20.txt")

Cv2 = pa.mat()
Cv2.load("Cv_L40.txt")

T = pa.mat()
T.load("Temp_L20.txt")

eps = pa.mat()
eps.load("expect_eps_Ln.txt")

m = pa.mat()
m.load("expect_m_Ln.txt")


Cv = np.array(Cv)
Cv2 = np.array(Cv2)

plt.plot(T, Cv, 'rx', label='Numerical data L=20')
plt.plot(T, Cv2, 'bx', label='Numerical data L=40')
smooth = gaussian_filter1d(Cv, 6, axis=0)
smooth2 = gaussian_filter1d(Cv2, 6, axis=0)
plt.plot(T, smooth, label='Gaussian filter L=20')
plt.plot(T, smooth2, label='Gaussian filter L=40')

#plt.plot(eps)

#y = pa.mat()
#y.load("expect_eps_Ln.txt")

#plt.plot(y)


#plt.vlines(1, Cv.min(), Cv.max()*1.01, 'r', ls='--')
plt.legend()
plt.show()


chi20 = pa.mat()
chi20.load("chi_L20.txt")

chi40 = pa.mat()
chi40.load("chi_L40.txt")

chi20 = np.array(chi20); chi40 = np.array(chi40)

plt.plot(T, chi20, 'rx', label='Numerical data L=20')
plt.plot(T, chi40, 'bx', label='Numerical data L=40')

smooth = gaussian_filter1d(chi20, 6, axis=0)
smooth2 = gaussian_filter1d(chi40, 6, axis=0)
plt.plot(T, smooth, label='Gaussian filter L=20')
plt.plot(T, smooth2, label='Gaussian filter L=40')
plt.show()"""