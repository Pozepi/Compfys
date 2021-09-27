import matplotlib.pyplot as plt
import numpy as np

"""
38.4166
0.23053
0.387868
0.422061
0.322253
0.120131
-0.120131
-0.322253
-0.422061
-0.387868
-0.23053

445.583
-0.23053
0.387868
-0.422061
0.322253
-0.120131
-0.120131
0.322253
-0.422061
0.387868
-0.23053
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

    a, b = zip(*sorted(zip(lambda_, v))) # sort
    #print(a)

    x = []
    with open('linspace_'+filename+'.txt', 'r') as file:
        for line in file.readlines():
            l = line.split()
            print(l)
            try: 
                x.append(float(l[0]))
            except: 
                pass

    # Ploott
    for i in range(3):
        plt.plot(x, v[i], label='v['+str(i)+']')
    
    plt.legend()
    plt.xlabel('Length x')
    plt.ylabel('Displacement')
    plt.grid()
    plt.savefig(filename+'.png')
    plt.clf()
    
    return 0

magic('N_10')
magic('N_100')
# a, b = zip(*sorted(zip(a,b)))
