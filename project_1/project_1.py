import matplotlib.pyplot as plt


x = []
u = []

with open('values2.txt', 'r') as file:
    for line in file.readlines():
        l = line.split()
        x.append(float(l[0]))
        u.append(float(l[1]))

plt.plot(x,u)
plt.show()
