import matplotlib.pyplot as plt
import numpy as np


""" Plot one particle"""

filename = 'RK4_one_particle'#'forward_euler'

def find_values(filename):
    t = []
    rx = []; ry = []; rz = []
    vx = []; vy = []; vz = []
    i = 0
    with open(filename+'.txt', 'r') as file:
        for line in file.readlines():
            l = np.array(line.split())

            if len(l)==0:
                # then we have reacha new quantity
                i += 1
            elif len(l)==1:
                t.append(float(l[0]))
            elif len(l)>1:
                if len(l) < 4:
                    if i%2 != 0:
                        vx.append(float(l[0]))
                        vy.append(float(l[1]))
                        vz.append(float(l[2]))
                    if i%2 == 0:
                        rx.append(float(l[0]))
                        ry.append(float(l[1]))
                        rz.append(float(l[2]))

    return t, rx, ry, rz, vx, vy, vz

    
#t, rx, ry, rz, vx, vy, vz = find_values(filename)

""" Plot the position """
def plot_xy(filename):
    fig, [ax1, ax2] = plt.subplots(figsize=(10,5), ncols=2)
    ax1.set_xlabel('X position'); ax1.set_ylabel('Y position')
    ax2.set_xlabel('Time [micro seconds]'); ax2.set_ylabel('Z position')
    [axi.grid() for axi in [ax1, ax2]]

    ax1.plot(rx[0], ry[0], 'rx', label='Start pos')
    ax1.plot(rx[-1], ry[-1], 'kx', label='End pos')
    ax1.plot(rx, ry)
    ax1.legend()
    ax2.plot(t, rz)
    plt.show()

def plot_compare(filename):
    fig, [ax1, ax2] = plt.subplots(figsize=(10,5), ncols=2)

    [[axi.set_xlabel('X position'), axi.set_ylabel('Y position'), axi.grid()] for axi in [ax1, ax2]]
    ax1.plot(rx[0], ry[0], 'rx', label='Start pos')
    ax1.plot(rx[-1], ry[-1], 'kx', label='End pos')
    ax1.plot(rx, ry)
    ax1.legend()
    ax1.set_title('Runge Kutta 4')

    t, rx, ry, rz, vx, vy, vz = find_values('forward_euler')
    ax2.plot(rx[0], ry[0], 'rx', label='Start pos')
    ax2.plot(rx[-1], ry[-1], 'kx', label='End pos')
    ax2.plot(rx, ry)
    ax2.legend()
    ax2.set_title('Forward Euler')
    plt.show()

def plot_vel():
    """ Plot the velocity over time """
    fig, [ax1, ax2, ax3] = plt.subplots(figsize=(20,10), ncols=3)

    [axi.set_xlabel('Time [micro seconds]') for axi in [ax1, ax2, ax3]]
    ax1.set_ylabel('X velocity')
    ax2.set_ylabel('Y velocity')
    ax3.set_ylabel('Z velocity')


def plot_multiple_particles(filename):
    t = []
    rx = []; ry = []; rz = []
    vx = []; vy = []; vz = []
    i = 0
    with open(filename+'.txt', 'r') as file:
        for line in file.readlines():
            # get data
            l = np.array(line.split())

            # Stage 1: The time array
            if len(l)==1:
                t.append(float(l[0]))

            # Stage 2: approaching the velocity
            elif len(l)==0:
                # then we have reacha new quantity
                i += 1

            # Stage 3: append velocity and position
            elif len(l)>1:
                particle_N = len(l)//3
                for k in range(0,len(l),3):
                    if i%2 != 0:
                        vx.append(float(l[k+0])) # vx1p1, vx1p2, vx1p3, ..., vx2p1, vx2p2, vx3p3, .... 
                        vy.append(float(l[k+1]))
                        vz.append(float(l[k+2]))
                    if i%2 == 0:
                        rx.append(float(l[k+0]))
                        ry.append(float(l[k+1]))
                        rz.append(float(l[k+2]))
    
    ax = plt.axes(projection='3d')
    ax.set_xlabel('X position')
    ax.set_ylabel('Y position')
    ax.set_zlabel('Z position')

    #ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens');

    for i in range(particle_N):
        rx_i = rx[i::particle_N]
        print(rx_i)
        ry_i = ry[i::particle_N]
        rz_i = rz[i::particle_N]
        #ax.plot(rx_i[0], ry_i[-1], 'x')
        ax.plot3D(rx_i, ry_i, rz_i, label='Particle %.i'%i)
    ax.legend()
    plt.show()

plot_multiple_particles('RK4_two_particles_interaction')