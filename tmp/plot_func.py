import matplotlib.pyplot as plt


x = []
u = []


with open('relerror.txt', 'r') as file:
    for line in file.readlines():
        l = line.split()
        x.append(float(l[0]))
        u.append(float(l[1]))

plt.plot(x,u, label='1')

"""
x2 = []
u2 = []
with open('values2.txt', 'r') as file:
    for line in file.readlines():
        l = line.split()
        x2.append(float(l[0]))
        u2.append(float(l[1]))

plt.plot(x,u, label='1')
plt.plot(x2,u2, label='2')
plt.legend()
plt.show()
"""
plt.yscale('log')
plt.legend()
plt.show()