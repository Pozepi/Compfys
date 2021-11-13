import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

eps = pa.mat()
eps.load("eps.txt")
eps = np.array(eps[:,0])
x = np.linspace(0,len(eps), len(eps))
print(np.std(eps))

fix, ax = plt.subplots(ncols=2)
ax[0].hist(eps, bins = 25)
ax[1].plot(x, eps)
plt.show()
