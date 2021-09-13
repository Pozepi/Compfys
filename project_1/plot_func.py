import matplotlib.pyplot as plt
import numpy as np

print('Please input your number and hit enter \n \n')

print('[1] Problem 2\n')
print('[2] Problem 7\n')
print('[3] Problem 8\n')
print('[4] Problem 8.c\n')
print('[5] Problem 10\n')
problem_number = input()

if problem_number == '1':
    """Problem 2 plotting"""
    filename = 'values1'

    x = []
    u = []
    with open(filename+'.txt', 'r') as file:
        for line in file.readlines():
            l = line.split()
            x.append(float(l[0]))
            u.append(float(l[1]))

    plt.xlabel('x values')
    plt.ylabel('Poisson u(x)')
    plt.plot(x, u)
    plt.title('Analytical Poisson solution')
    plt.show()

elif problem_number == '2':
    """Problem 7 plotting goes here"""

    filenames = 'u_values_N_'
    
    x = [[], [], [], []]
    u = [[], [], [], []]
    
    # Load numerical values
    for i in range(4):
        with open(filenames+str(i)+'.txt', 'r') as file:
            for line in file.readlines():
                l = line.split()
                x[i].append(float(l[0]))
                u[i].append(float(l[1]))
    
    # Load analytical values
    x_a = []
    u_a = []
    with open('values1.txt', 'r') as file:
        for line in file.readlines():
            l = line.split()
            x_a.append(float(l[0]))
            u_a.append(float(l[1]))

    for i in range(4):
        plt.plot(x[i], u[i], label='N = 10^'+str(i+1))


    plt.plot(x_a, u_a, label='Analytical')
    plt.title('Numerical Poisson solution')

    plt.xlabel('x values')
    plt.ylabel('Poisson u(x)')
    plt.legend()
    plt.show()

elif problem_number == '3':
    """Problem 8 goes here: error and relative error"""
    x_e = [[], [], [], []]
    u_e = [[], [], [], []]

    x_r = [[], [], [], []]
    u_r = [[], [], [], []]

    filenames_e = 'error_values_N_'
    filenames_r = 'rel_error_values_N_'

    # Load numerical values
    for i in range(4):
        # open the error values
        with open(filenames_e+str(i)+'.txt', 'r') as file:
            for line in file.readlines():
                l = line.split()
                x_e[i].append(float(l[0]))
                u_e[i].append(float(l[1]))
                #x_e_val = l[0]
                #u_e_val = l[1]
                #x_e[i].append(x_e_val)
                #u_e[i].append(u_e_val)
        # open the rel error values
        with open(filenames_r+str(i)+'.txt','r') as file:
            for line in file.readlines():
                l = line.split()
                x_r[i].append(float(l[0]))
                u_r[i].append(float(l[1]))

    fig, ax = plt.subplots(ncols=2)
    

    for i in range(4):
        ax[0].plot(x_e[i], u_e[i], label='Error N = 10^'+str(i+1),)
        ax[1].plot(x_r[i], u_r[i], label='Rel Error N = 10^'+str(i+1))
        #plt.plot(x_e[i], u_e[i], label='Error N = 10^'+str(i+1),)
        #plt.plot(x_r[i], u_r[i], label='Rel Error N = 10^'+str(i+1))

    ax[0].set_ylabel('Error $\log_{10}(\Delta_i)$')
    ax[0].set_xlabel('x values')
    ax[0].set_yscale('log')
    ax[0].set_xticks(np.linspace(0,1,10))
    ax[0].set_title('Absolute Error')
    ax[0].legend()

    ax[1].set_ylabel('Relative Error $\log_{10}(\epsilon_i)$')
    ax[1].set_xlabel('x values')
    ax[1].set_yscale('log')
    ax[1].set_xticks(np.linspace(0,1,10))
    ax[1].set_title('Relative Error')
    ax[1].legend()

    plt.show()
    
elif problem_number == '4':
    """Problem 8.c goes here: error and relative error"""

    filename = 'max_eps'

    x = []
    u = []
    with open(filename+'.txt', 'r') as file:
        for line in file.readlines():
            l = line.split()
            x.append(float(l[0]))
            u.append(float(l[1]))

    plt.plot(x, u)
    plt.title('Max relative error')
    plt.xlabel('Array size N')
    plt.ylabel(r'log$_{10}($max($\epsilon_i$))')
    #plt.yscale('log')
    plt.xscale('log')
    plt.show()

elif problem_number == '5':
    """Time comparison"""
    filename0 = 'time_general'

    x0 = []
    y0 = []
    with open(filename0+'.txt', 'r') as file:
        for line in file.readlines():
            l = line.split()
            x0.append(float(l[0]))
            y0.append(float(l[1]))

    filename0 = 'time_special'

    x1 = []
    y1 = []
    with open(filename0+'.txt', 'r') as file:
        for line in file.readlines():
            l = line.split()
            x1.append(float(l[0]))
            y1.append(float(l[1]))

    plt.plot(x0, y0, 'r.', label='General algorithm time')
    plt.plot(x1, y1, 'b.', label='Special algorithm time')
    plt.xscale('log')
    plt.legend()
    plt.xlabel('Array size $N$')
    plt.ylabel('Runtime (s)')
    plt.title('Runtime for a given array size $N$')
    plt.show()
