#!/usr/bin/env python3

import sys
import math
import struct
import numpy as np
import matplotlib.pyplot as plt

force = 1e-6
viscosity = 1/6
ly = 50
width = ly-4

error_criteria = 1e-2


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


def poiseuille_profile(y):
    return force / (2 * viscosity) * y * (width - y)


def compare_velocity(vel):
    y = np.arange(0.5, len(vel))
    vel_poiseuille = poiseuille_profile(y)
    relative_error = np.abs(vel - vel_poiseuille) / vel_poiseuille
    return np.max(relative_error)


vel = load_vel('data')
max_error = compare_velocity(vel[-1])

if (max_error > error_criteria):
    print(f'Error: Max relative error in velocity profile too large ({max_error:g} > {error_criteria:g})')
    sys.exit(1)
