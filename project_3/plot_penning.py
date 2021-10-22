import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import pyarma as pa
import os
import glob


figure_path = 'figures/'

if not os.path.isdir(figure_path):
    os.mkdir(figure_path)

# Collect correct path strings
path = 'output_files/'
single_particle_paths = []; [single_particle_paths.append(mypath+'/'+mypath.split('/')[-1]) for mypath in glob.glob(path+'1_*')] #os.listdir(path)
two_particle_paths = []; [two_particle_paths.append(mypath+'/'+mypath.split('/')[-1]) for mypath in glob.glob(path+'2_*')] #os.listdir(path)
n_particle_paths = []; [n_particle_paths.append(mypath+'/'+mypath.split('/')[-1]) for mypath in glob.glob(path+'n_*')] #os.listdir(path)




""" Analytical """
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
    



""" Extract one particle coords """
def find_values_single(filename):
    t = pa.mat()
    r = pa.mat()
    v = pa.mat()

    t.load(filename+'_t.txt')
    r.load(filename+'_pos.txt')
    v.load(filename+'_v.txt')
    return np.array(t)[:,0], np.array(r), np.array(v)



""" Compare numerical to analytical for a single particle"""
def plot_compare_analytical(filename1, plot=False):
    t, r, v = find_values_single(filename1)
    fig, [ax1, ax2] = plt.subplots(figsize=(10,5), ncols=2)

    ax1.set_xlabel('X position [μm]'); ax1.set_ylabel('Y position [μm]'); ax1.grid()
    ax1.plot(r[:,0], r[:,1], label='Numerical')
    ax1.plot(x(t, v[0,1], r[0,0]), y(t, v[0,1], r[0,0]), label='Analytical', ls='--')
    ax1.legend()
    #ax1.set_title()

    ax2.set_xlabel('Time [s]'); ax2.set_ylabel('Z position [μm]'); ax2.grid()
    ax2.plot(t, r[:,2], label='Numerical')
    ax2.plot(t, z(t, r[0,2]), label='Analytical', ls='--')
    ax2.legend()

    titlename = r''
    if filename1[-6:-4] == 'EU':
        titlename += r'Euler with '
    elif filename1[-6:-4] == 'RK':
        titlename += r'Runge Kutta 4 with '

    titlename += r'dt = $10^{-%.i}$'%(int(filename1[-1:]))

    ax1.set_title(titlename + ' in x-y plane')
    ax2.set_title(titlename + ' in z-t plane')
    plt.savefig(figure_path+filename1.split('/')[-1]+'_compare_analytic.pdf')
    if plot:
        plt.show()





""" Find the relative error """
def relerror(plot=False):
    fig, ax = plt.subplots(figsize=(11, 6), ncols=2)

    for pathi in single_particle_paths:
        # numerical values        
        myname = pathi.split('/')[-1]
        t_n, r_n, v_n = find_values_single(pathi)
        jumpo = 100
        t_n = t_n[::jumpo]
        r_n = r_n[::jumpo]
        v_n = v_n[::jumpo]
        # analytical values
        x_a = x(t_n, v_n[0,1], r_n[0,0])
        y_a = y(t_n, v_n[0,1], r_n[0,0])
        z_a = z(t_n, r_n[0,2])

        R_n = np.sqrt(r_n[:,0]**2 + r_n[:,1]**2 + r_n[:,2]**2)
        R_a = np.sqrt(x_a**2 + y_a**2 + z_a**2)
        relerror = abs(R_a - R_n)

        if pathi[-6:-4] == 'EU':
            labelname = 'dt = $10^{-%.i}$'%int(pathi[-1:])
            ax[0].plot(t_n, relerror, label=labelname)
        elif pathi[-7:-5] == 'RK':
            labelname = 'dt = $10^{-%.i}$'%int(pathi[-1:])
            ax[1].plot(t_n, relerror, label=labelname)

    ax[0].set_title('Euler')
    ax[1].set_title('Runge Kutta 4')
    [[axi.set_yscale('log'), axi.legend(), axi.grid()] for axi in ax]
    plt.savefig(figure_path+'relerror.pdf')
    if plot:
        plt.show()

""" Plot the position """
"""
def plot_xy_zt(filename1,plot=False):
    t, r, v = find_values_single(filename1)
    fig, [ax1, ax2] = plt.subplots(figsize=(10,5), ncols=2)
    ax1.set_xlabel('X position'); ax1.set_ylabel('Y position')
    ax2.set_xlabel('Time [micro seconds]'); ax2.set_ylabel('Z position')
    [axi.grid() for axi in [ax1, ax2]]

    ax1.plot(r[0,0], r[0,1], 'rx', label='Start pos')
    ax1.plot(r[-1,0], r[-1,1], 'kx', label='End pos')
    ax1.plot(r[:,0], r[:,1])
    ax1.legend()
    ax2.plot(t, r[:,2])
    plt.savefig(figure_path)
    if plot:
        plt.show()
"""



""" Input filenames here not path """
def plot_compare(filename1, filename2, plot=False):

    t, r, v = find_values_single(path+filename1+'/'+filename1)
    fig, [ax1, ax2] = plt.subplots(figsize=(10,5), ncols=2)

    [[axi.set_xlabel('X position'), axi.set_ylabel('Y position'), axi.grid()] for axi in [ax1, ax2]]
    ax1.plot(r[:,0], r[:,1])
    ax1.plot(r[0,0], r[0,1], 'rx', label='Start pos')
    ax1.plot(r[-1,0], r[-1,1], 'kx', label='End pos')
    ax1.legend()
    if filename1[5:7] == 'EU':
        ax1.set_title('Euler with dt = $10^{-%.i}$'%int(filename1[-1:]))
    elif filename1[5:7] == 'RK':
        ax1.set_title('Runge Kutta with dt = $10^{-%.i}$'%int(filename1[-1:]))

    t, r, v = find_values_single(path+filename2+'/'+filename2)
    ax2.plot(r[:,0], r[:,0])
    ax2.plot(r[0,0], r[0,1], 'rx', label='Start pos')
    ax2.plot(r[-1,0], r[-1,1], 'kx', label='End pos')
    ax2.legend()
    if filename2[5:7] == 'EU':
        ax2.set_title('Euler with dt = $10^{-%.i}$'%int(filename2[-1:]))
    elif filename2[5:7] == 'RK':
        ax2.set_title('Runge Kutta with dt = $10^{-%.i}$'%int(filename2[-1:]))
    name = figure_path+'_'+filename1[5:7]+filename1[-1:]+'_vs_'+filename2[5:7]+filename2[-1:]+'.pdf'
    plt.savefig(name)
    if plot:
        plt.show()



def plot_vel(filename,plot=False,save=False):
    """ Plot the velocity over time """
    t, r, v = find_values_single(filename)
    fig, [ax1, ax2, ax3] = plt.subplots(figsize=(20,10), ncols=3)
    axes = [ax1, ax2, ax3]

    ax1.plot(r[:,0], v[:,0])
    ax2.plot(r[:,1], v[:,1])
    ax3.plot(r[:,2], v[:,2])

    ax1.set_xlabel('X position [μm]'); ax1.set_ylabel('X velocity [m/s]')
    ax2.set_xlabel('Y position [μm]'); ax2.set_ylabel('Y velocity [m/s]')
    ax3.set_xlabel('Z position [μm]'); ax3.set_ylabel('Z velocity [m/s]')

    name = filename.split('/')[-1]
    [axi.grid() for axi in axes]
    if name[5:7] == 'EU':
        ax2.set_title('Euler with dt = $10^{-%.i}$'%int(name[-1:]))
    elif name[5:7] == 'RK':
        ax2.set_title('Runge Kutta with dt = $10^{-%.i}$'%int(name[-1:]))
    if plot:
        plt.show()
    if save:
        
        savename = figure_path+name[5:7]+filename[-1:]+'_phasespace.pdf'
        plt.savefig(savename)



def plot_multiple_particles(filename, animate=False, anisave=False, plot3d=False, saveplot3d=False,
    plot=False, saveplot=False):
    t, r, v = find_values_single(filename)
    name = filename.split('/')[-1]
    method = name[5:7]
    power = name[-1:]
    interaction = name[2:4]

    if name[-7:-4] == 'noz':
        z_name = 'no z'
    else:
        z_name = 'z start'
    if interaction == 'no':
        interaction = 'no interaction'
    elif interaction == 'in':
        interaction = 'with interaction'
    if method == 'EU':
        method = 'Euler'
    elif method == 'RK':
        method = 'Runge Kutta'

    name_title = method + ' ' + interaction + ' ' + z_name
    
    particle_N = len(r[0])//3
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
        rx_i.append(r[:,i*3])
        ry_i.append(r[:,i*3+1])
        rz_i.append(r[:,i*3+2])
    rx_i = np.array(rx_i)
    ry_i = np.array(ry_i)
    rz_i = np.array(rz_i)


    #line = []
    # set up first frame
    """
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
    """
    if plot or saveplot:
        fig,ax = plt.subplots()
        ax.set_xlabel('X position [μm]')
        ax.set_ylabel('Y position [μm]')
        ax.set_title(name_title)
        ax.grid()
        for i in range(particle_N):
            ax.plot(rx_i[i], ry_i[i]) #, label='dt = $10^{-%.i}$'%int(power))
        if plot:
            plt.show()
        elif saveplot:
            plt.savefig(figure_path+name+'.pdf')
    
    if plot3d or saveplot3d:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X position [μm]')
        ax.set_ylabel('Y position [μm]')
        ax.set_zlabel('Z position [μm]')
        ax.set_title(name_title)
        for i in range(particle_N):
            ax.plot(rx_i[i], ry_i[i], rz_i[i]) #, label='dt = $10^{-%.i}$'%int(power))
        if plot3d:
            plt.show()
        elif saveplot3d:
            plt.savefig(figure_path+name+'_3d.pdf')

        

    elif animate:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X position [μm]')
        ax.set_ylabel('Y position [μm]')
        ax.set_zlabel('Z position [μm]')

        lim = 1e2
        ax.set_xlim3d([-lim, lim])
        ax.set_ylim3d([-lim, lim]) 
        ax.set_zlim3d([-lim, lim])

        line = [ax.plot(rx_ij[0], ry_ij[0], rz_ij[0], 'o')[0] for rx_ij, ry_ij, rz_ij in zip(rx_i, ry_i, rz_i)]
        lines = [ax.plot(rx_ij[0:1], ry_ij[0:1], rz_ij[0:1])[0] for rx_ij, ry_ij, rz_ij in zip(rx_i, ry_i, rz_i)]
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
            tail_length = 1000
            i *= 10 
            if i - tail_length < 0:
                j = 0
            elif i - tail_length >= 0:
                j = i - tail_length
            #print(j)
            #[[line[i].set_data([rx_ij[:i], ry_ij[:i], rz_ij[:i]]) for i in range(particle_N)] for rx_ij, ry_ij, rz_ij in zip(rx_i, ry_i, rz_i)]
            #[ax.plot(rx_ij[j:i], ry_ij[j:i], rz_ij[j:i], 'r') for rx_ij, ry_ij, rz_ij in zip(rx_i, ry_i, rz_i)]
            #[[linei.set_data(rx_ij[i], ry_ij[i]) for rx_ij, ry_ij in zip(rx_i, ry_i)] for linei in line]
            
            #lines[0].set_data(rx_i[0][i], ry_i[0][i])
            #lines[0].set_3d_properties(rz_i[0][i])

            #line[-1].set_data_3d(rx_i[-1][i], ry_i[-1][i], rz_i[-1][i])

            [linei.set_data_3d(rx_ij[i], ry_ij[i], rz_ij[i]) for linei, rx_ij, ry_ij, rz_ij in zip(line, rx_i, ry_i, rz_i)]
            [linesi.set_data_3d(rx_ij[j:i], ry_ij[j:i], rz_ij[j:i]) for linesi, rx_ij, ry_ij, rz_ij in zip(lines, rx_i, ry_i, rz_i)]
            #[[linei.set_3d_properties(rz_ij[i]) for rz_ij in rz_i] for linei in line]

        ani = animation.FuncAnimation(fig, update, len(rx_i[0])//10, blit=False, interval=1)

        if anisave:
            FFwriter = animation.FFMpegWriter(fps=1000)
            ani.save('penning.mp4', writer=FFwriter)
        plt.show()
    
#[plot_compare_analytical(pathi) for pathi in single_particle_paths]
#relerror()
#[plot_vel(pathi, save=True) for pathi in single_particle_paths]

#plot_xy('RK4_one_particle')

[plot_multiple_particles(pathi, saveplot3d=True) for pathi in two_particle_paths]



#plot_compare('RK4_one_particle', 'euler_one_particle')
#plot_compare('RK4_two_particles_interaction', 'euler_two_particles_interaction')
#plot_compare_analytical('RK4_against_analytical')
#plot_multiple_particles('RK4_random_test' ,animate=True)
#plot_multiple_particles('RK4_two_particles_interaction_z_start' ,animate=True)