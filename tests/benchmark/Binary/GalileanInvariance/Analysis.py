#!/usr/bin/env python3

import sys
import math
import struct
import numpy as np

radius = 20
vel0 = [1e-3, 1e-3]
tmax = 10000

error_criteria = 1e-2


def read_header(header):
    with open(header, 'rb') as headerFile:
        lx = struct.unpack('=i', headerFile.read(4))[0]
        ly = struct.unpack('=i', headerFile.read(4))[0]
        lz = struct.unpack('=i', headerFile.read(4))[0]
        ndim = struct.unpack('=i', headerFile.read(4))[0]
        tend = struct.unpack('=i', headerFile.read(4))[0]
        tinc = struct.unpack('=i', headerFile.read(4))[0]
    return {'lx':lx, 'ly':ly, 'lz':lz, 'ndim':ndim, 'tend':tend, 'tinc':tinc}


def read_file(filename, data):
    with open(filename, 'rb') as file:
        for i in range(data.size):
            data.ravel()[i] = struct.unpack('=d', file.read(8))[0]


def load_phi(direc='data'):
    info = read_header(direc+"/Header.mat")
    times = np.arange(0, info['tend']+1, info['tinc'])
    phi = np.zeros((len(times), info['lx'], info['ly'], info['lz']))
    for it, t in enumerate(times):
        read_file(direc+"/OrderParameter_t%li.mat"%t, phi[it])
    return phi[:,:,:,0]


def load_vel(direc='data'):
    info = read_header(direc+"/Header.mat")
    times = np.arange(0, info['tend']+1, info['tinc'])
    vel = np.zeros((len(times), info['lx'], info['ly'], info['lz'], info['ndim']))
    for it, t in enumerate(times):
        read_file(direc+"/Velocity_t%li.mat"%t, vel[it])
    return vel[:,:,:,0]


def plot(phi, vel=None, t=-1):
    import matplotlib.pyplot as plt
    plt.contourf(phi[t].T, cmap='Blues')
    plt.colorbar()
    plt.contour(phi[t].T, levels=0, colors='k')
    if (vel is not None): plt.quiver(vel[t,:,:,0], vel[t,:,:,1])
    plt.gca().set_aspect('equal')
    plt.show()


def compare_centre(phi):
    c = 0.5 * (phi[-1] + 1)
    lx = c.shape[0]
    ly = c.shape[1]
    x, y = np.mgrid[:lx, :ly]
    dx = np.sum(x*c)/np.sum(c) - lx//2
    dy = np.sum(y*c)/np.sum(c) - ly//2
    dx_theory = vel0[0]*tmax
    dy_theory = vel0[1]*tmax
    relative_errorx = np.abs(dx - dx_theory) / dx_theory
    relative_errory = np.abs(dy - dy_theory) / dy_theory
    return max(relative_errorx, relative_errory)


def compare_vel(vel):
    relative_errorx = np.abs(vel[-1,:,:,0] - vel0[0]) / vel0[0]
    relative_errory = np.abs(vel[-1,:,:,1] - vel0[1]) / vel0[1]
    return np.max([relative_errorx, relative_errory])


phi = load_phi('data')
vel = load_vel('data')
# plot(phi, vel)

pos_error = compare_centre(phi)
if (pos_error > error_criteria):
    print(f'Error: Incorrect change in position (relative error, {pos_error:g} > {error_criteria:g})')
    sys.exit(1)

vel_error = compare_vel(vel)
if (vel_error > error_criteria):
    print(f'Error: Incorrect velocities (relative error, {vel_error:g} > {error_criteria:g})')
    sys.exit(1)
