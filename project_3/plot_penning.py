import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation


""" Plot one particle"""

filename = 'RK4_one_particle'#'forward_euler'

def find_values(filename, particle_jump=0):
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
                if particle_jump*3 > len(l):
                    print('Particle jump index is too high for array with ', len(l)//3, ' particles')
                    raise ValueError
                if i%2 != 0:
                    vx.append(float(l[particle_jump*3 + 0]))
                    vy.append(float(l[particle_jump*3 + 1]))
                    vz.append(float(l[particle_jump*3 + 2]))
                if i%2 == 0:
                    rx.append(float(l[particle_jump*3 + 0]))
                    ry.append(float(l[particle_jump*3 + 1]))
                    rz.append(float(l[particle_jump*3 + 2]))

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

def plot_compare(filename1, filename2):
    t, rx, ry, rz, vx, vy, vz = find_values(filename1)
    fig, [ax1, ax2] = plt.subplots(figsize=(10,5), ncols=2)

    [[axi.set_xlabel('X position'), axi.set_ylabel('Y position'), axi.grid()] for axi in [ax1, ax2]]
    ax1.plot(rx, ry)
    ax1.plot(rx[0], ry[0], 'rx', label='Start pos')
    ax1.plot(rx[-1], ry[-1], 'kx', label='End pos')
    ax1.legend()
    ax1.set_title(filename1)

    t, rx, ry, rz, vx, vy, vz = find_values(filename2)
    ax2.plot(rx, ry)
    ax2.plot(rx[0], ry[0], 'rx', label='Start pos')
    ax2.plot(rx[-1], ry[-1], 'kx', label='End pos')
    ax2.legend()
    ax2.set_title(filename2)
    plt.show()

def plot_vel():
    """ Plot the velocity over time """
    fig, [ax1, ax2, ax3] = plt.subplots(figsize=(20,10), ncols=3)

    [axi.set_xlabel('Time [micro seconds]') for axi in [ax1, ax2, ax3]]
    ax1.set_ylabel('X velocity')
    ax2.set_ylabel('Y velocity')
    ax3.set_ylabel('Z velocity')


def plot_multiple_particles(filename, animate=False):
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
    
    rx_i = []
    ry_i = []
    rz_i = []
    for i in range(particle_N):
        #rx_i = rx[i::particle_N]
        #ry_i = ry[i::particle_N]
        #rz_i = rz[i::particle_N]
        #plot both particles
        #ax.plot(rx_i[0], ry_i[-1], 'x')
        #ax.plot3D(rx_i, ry_i, rz_i, label='Particle %.i'%i)
        rx_i.append(rx[i::particle_N])
        ry_i.append(ry[i::particle_N])
        rz_i.append(rz[i::particle_N])

    #line = []
    # set up first frame
    if not animate:
        fig, ax = plt.subplots(figsize=(10,5))
        [ax.plot(rx_ii, ry_ii) for rx_ii, ry_ii in zip(rx_i, ry_i)]
        ax.grid()
        ax.set_xlabel('X position [μm]')
        ax.set_ylabel('Y position [μm]')
        d = 1e4
        ax.set_xlim(-d, d)
        ax.set_ylim(-d, d)
        plt.show()

    elif animate:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X position [μm]')
        ax.set_ylabel('Y position [μm]')
        ax.set_zlabel('Z position [μm]')

        lim = 1e4
        ax.set_xlim3d([-lim, lim])
        ax.set_ylim3d([-lim, lim]) 
        ax.set_zlim3d([-lim, lim])
        line = [ax.plot(rx_ij[0], ry_ij[0], rz_ij[0], 'ko')[0] for rx_ij, ry_ij, rz_ij in zip(rx_i, ry_i, rz_i)]
        """
        lines = []; [lines.append([]) for i in range(particle_N)] # particle_N lines
        for j in range(particle_N):
            lines[j].plot()
        N = len(rx_i[0])

        def update(i, line):
            for j in range(particle_N):
                for i in range(len(rx_i)):
                    line[j].set_data(rx_i[:i], ry_i[:i], rz_i[:i])
        """

        def update(i):
            #[[line[i].set_data([rx_ij[:i], ry_ij[:i], rz_ij[:i]]) for i in range(particle_N)] for rx_ij, ry_ij, rz_ij in zip(rx_i, ry_i, rz_i)]
            [ax.plot(rx_ij[:i], ry_ij[:i], rz_ij[:i], 'r') for rx_ij, ry_ij, rz_ij in zip(rx_i, ry_i, rz_i)]
            #[[linei.set_data(rx_ij[i], ry_ij[i]) for rx_ij, ry_ij in zip(rx_i, ry_i)] for linei in line]
            [[linei.set_data_3d(rx_ij[i], ry_ij[i], rz_ij[i]) for rx_ij, ry_ij, rz_ij in zip(rx_i, ry_i, rz_i)] for linei in line]
            #[[linei.set_3d_properties(rz_ij[i]) for rz_ij in rz_i] for linei in line]

        ani = animation.FuncAnimation(fig, update, len(rx_i[0]), blit=False, interval=1)


        #FFwriter = animation.FFMpegWriter()
        #ani.save('matplot003.mp4', writer=FFwriter)
        plt.show()
    

q = 1
B0 = 9.65e1
V0 = 9.65e8
d  = 1e4
m  = 38.97

omega0 = q*B0/m
omegaz = np.sqrt(2*q*V0/(m*d**2))

omega_pluss = 0.5*(omega0 + np.sqrt(omega0**2 - 2*omegaz**2))
omega_minus = 0.5*(omega0 - np.sqrt(omega0**2 - 2*omegaz**2))

def x(t, v0, x0):
    # v0 initial velocity in y
    # x0 initial position in x
    A_pluss = (v0+x0*omega_minus)/(omega_minus - omega_pluss)
    A_minus = -(v0+x0*omega_pluss)/(omega_minus - omega_pluss)
    return A_pluss*np.cos(omega_pluss*t) + A_minus*np.cos(omega_minus*t)

def y(t, v0, x0):
    # v0 initial velocity in y
    # x0 initial position in x
    A_pluss = (v0+x0*omega_minus)/(omega_minus - omega_pluss)
    A_minus = -(v0+x0*omega_pluss)/(omega_minus - omega_pluss)
    return -A_minus*np.sin(omega_minus*t) - A_pluss*np.sin(omega_pluss*t)

def z(t, z0):
    return z0*np.cos(omegaz*t)

def plot_compare_analytical(filename1):
    t, rx, ry, rz, vx, vy, vz = find_values(filename1)
    t = np.array(t)
    fig, [ax1, ax2] = plt.subplots(figsize=(10,5), ncols=2)

    ax1.set_xlabel('X position [μm]'); ax1.set_ylabel('Y position [μm]'); ax1.grid()
    ax1.plot(rx, ry, label='Numerical')
    ax1.plot(x(t, vy[0], rx[0]), y(t, vy[0], rx[0]), label='Analytical', ls='--')
    ax1.legend()
    #ax1.set_title()

    ax2.set_xlabel('Time [s]'); ax2.set_ylabel('Z position [μm]'); ax2.grid()
    ax2.plot(t, rz, label='Numerical')
    ax2.plot(t, z(t, rz[0]), label='Analytical', ls='--')
    ax2.legend()
    plt.show()


#plot_compare('RK4_one_particle', 'euler_one_particle')
#plot_compare('RK4_two_particles_interaction', 'euler_two_particles_interaction')

#plot_compare_analytical('RK4_against_analytical')

plot_multiple_particles('RK4_random_test_seed_1')