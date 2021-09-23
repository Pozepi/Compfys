import matplotlib.pyplot as plt
import numpy as np

N     = np.array([5,  10,  15,  20,  25])
count = np.array([30, 146, 349, 636, 1035])

def dense(N):
    return N**2 - 0.5*N

plt.plot(N, count, label='Result')
plt.xscale('log')
plt.yscale('log')
print(np.gradient(np.log(count), np.log(N)))
plt.plot(N, dense(N), label='Model')
plt.xlabel('Matrix column/row size N')
plt.ylabel('Number of iterations')
plt.grid()
plt.legend()
plt.savefig('iterations.png')