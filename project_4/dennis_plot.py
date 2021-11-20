import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

"""
Cv = pa.mat()
Cv.load("Cv_L20.txt")

T = pa.mat()
T.load("Temp_L20.txt")

eps = pa.mat()
eps.load("expect_eps_Ln.txt")

m = pa.mat()
m.load("expect_m_Ln.txt")
plt.plot(T, Cv)
#plt.plot(eps)"""

x = pa.mat()
x.load("Temp.txt")

#y = pa.mat()
#y.load("expect_eps_Ln.txt")

#plt.plot(y)


#plt.vlines(1, Cv.min(), Cv.max()*1.01, 'r', ls='--')
plt.show()