import matplotlib.pyplot as plt
import numpy as np

filename = 'forward_euler'
t = []
r = [[], [], []]
v = [[], [], []]

Q = np.array([t,r,v])
i = 0 # quantity that Q will touch
with open(filename+'.txt', 'r') as file:
    for line in file.readlines():
        l = np.array(line.split())

        if len(l)==0:
            # then we have reacha new quantity
            i += 1

        elif len(l)==1:
            Q[i].append(float(l[0]))

        elif len(l)>1:
            for j in range(len(l)//3): 
                Q[i][0].append(l[j])
                Q[i][1].append(l[j+1])
                Q[i][2].append(l[j+2])

t,r,v=Q
""" Plot the position """
fig, [ax1, ax2] = plt.subplots(figsize=(15,10), ncols=2)
ax1.set_xlabel('X position'); ax1.set_ylabel('Y position')
ax2.set_xlabel('Time [micro seconds]'); ax2.set_ylabel('Z position')
[axi.grid() for axi in [ax1, ax2]]

ax1.plot(r[0], r[1])
ax2.plot(t, r[2])
plt.show()


""" Plot the velocity over time """
fig, [ax1, ax2, ax3] = plt.subplots(figsize=(20,10), ncols=3)

[axi.set_xlabel('Time [micro seconds]') for axi in [ax1, ax2, ax3]]
ax1.set_ylabel('X velocity')
ax2.set_ylabel('Y velocity')
ax3.set_ylabel('Z velocity')



"""
def magic(filename):
    lambda_ = [[]]
    v       = [[]]
    i = 0
    flip = 1
    with open(filename+'.txt', 'r') as file:
        for line in file.readlines():
            l = line.split()

            
            if len(l) == 0:
                flip = 1
                lambda_.append([])
                v.append([])
                i += 1

            elif l != '' and flip == 1:
                lambda_[i].append(float(l[0]))
                flip = 0

            elif l != '' and flip == 0:
                v[i].append(float(l[0]))
    
    v = v[:-1]
    lambda_ = lambda_[:-1]

    lambda_, v = zip(*sorted(zip(lambda_, v))) # sort
    #print(a)
    N = len(v)
    h = 1/(N+1)
    a_py = -1/(h**2)
    d_py =  2/(h**2)

    lambda_py = []; [lambda_py.append(d_py + 2*a_py*np.cos(i*np.pi/(N+1))) for i in range(1,N)]
    v_py = []; 
    for i in range(1-1,N):
        v_py.append([])
        for j in range(1,N+1):
            v_py[i].append(np.sin((i+1)*j*np.pi/(N + 1)))

    lambda_py, v_py = zip(*sorted(zip(lambda_py, v_py)))

    x = []
    with open('linspace_'+filename+'.txt', 'r') as file:
        for line in file.readlines():
            l = line.split()
            #print(l)
            try: 
                x.append(float(l[0]))
            except: 
                pass
    
    # Ploott
    for i in range(3):
        plt.plot(x, np.subtract(v[i], v[i][0]), label='v['+str(i)+']')
        plt.plot(x, np.subtract(v_py[i], v_py[i][0]), label='Analytical v['+str(i)+']')
    
    plt.legend()
    plt.xlabel('Length x')
    plt.ylabel('Displacement')
    plt.grid()
    plt.savefig(filename+'.png')
    plt.clf()
    
    return 0
"""