#!/usr/bin/env python3

# This script compares the curvature of the velocity in the two layers
# to the curvature that is expected based upon the force and viscosities

import sys
import math
import struct
import numpy as np

force = 1e-6
tau1 = 0.55
tau2 = 1

error_criteria = 1e-1


def read_file(filename, data):
    with open(filename, 'rb') as file:
        for i in range(data.size):
            data.ravel()[i] = struct.unpack('=d', file.read(8))[0]


def load_vel(direc='data'):
    with open(direc+"/Header.mat", 'rb') as headerFile:
        lx = struct.unpack('=i', headerFile.read(4))[0]
        ly = struct.unpack('=i', headerFile.read(4))[0]
        lz = struct.unpack('=i', headerFile.read(4))[0]
        ndim = struct.unpack('=i', headerFile.read(4))[0]
        tend = struct.unpack('=i', headerFile.read(4))[0]
        tinc = struct.unpack('=i', headerFile.read(4))[0]

    times = np.arange(0, tend+1, tinc)
    vel = np.zeros((len(times), lx, ly, lz, ndim))
    for it, t in enumerate(times):
        read_file(direc+"/Velocity_t%li.mat"%t, vel[it])
    return vel[:,0,2:-2,0,0]


def compare_curvature(vel):
    vel1 = vel[:len(vel)//2]
    vel2 = vel[len(vel)//2:]
    curv1 = - np.diff(np.diff(vel1))
    curv2 = - np.diff(np.diff(vel2))
    visc1 = (tau1 - 0.5) / 3
    visc2 = (tau2 - 0.5) / 3
    curv1_theory = force / visc1
    curv2_theory = force / visc2
    relative_error1 = np.abs(curv1 - curv1_theory) / curv1_theory
    relative_error2 = np.abs(curv2 - curv2_theory) / curv2_theory
    return np.max([relative_error1, relative_error2])

def plot(vel):
    import matplotlib.pyplot as plt
    plt.plot(vel[-1])
    plt.show()
    

vel = load_vel('data')
max_error = compare_curvature(vel[-1])

if (max_error > error_criteria):
    print(f'Error: Max relative error in velocity profile too large ({max_error:g} > {error_criteria:g})')
    sys.exit(1)
