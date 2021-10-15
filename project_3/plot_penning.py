import matplotlib.pyplot as plt
import numpy as np

filename = 'RK4'#'forward_euler'

"""
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
                if i%2 != 0:
                    vx.append(float(l[0]))
                    vy.append(float(l[1]))
                    vz.append(float(l[2]))
                if i%2 == 0:
                    rx.append(float(l[0]))
                    ry.append(float(l[1]))
                    rz.append(float(l[2]))

    return t, rx, ry, rz, vx, vy, vz
"""


def find_values(filename):
    t = []
    r = [[], [], []]
    v = [[], [], []]
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
                if i%2 != 0:
                    v[0].append(float(l[0]))
                    v[1].append(float(l[1]))
                    v[2].append(float(l[2]))
                if i%2 == 0:
                    r[0].append(float(l[0]))
                    r[1].append(float(l[1]))
                    r[2].append(float(l[2]))

    return t, rx, ry, rz, vx, vy, vz

t, rx, ry, rz, vx, vy, vz = find_values('RK4')
