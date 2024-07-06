#!/usr/bin/env python3

# This script is used to visualise the output of the simulation

import math
import struct
import numpy as np
import matplotlib.pyplot as plt


def coord_k(k, ly, lz):
    """From a k value, determines its xk, yk, and zk."""    
    xk = math.floor(k/(ly*lz))
    yk = math.floor((k - xk*lz*ly)/lz)
    zk = k - xk*lz*ly - yk*lz
    return xk, yk, zk


def read_data(direc='data'):
    print('Reading data...')
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
        try:
            vel_file = open(direc+"/Velocity_t%li.mat"%t, 'rb')

            for k in range(lx*ly*lz):
                (xk,yk,zk) = coord_k(k,ly,lz)
                for i in range(ndim):
                    vel[it,xk,yk,zk,i] = struct.unpack('=d', vel_file.read(8))[0]
            
            vel_file.close()
        except FileNotFoundError:
            vel = vel[:it]
            break
    return vel


def plot(vel):
    print('Plotting...')
    fig, ax = plt.subplots(2, 1)
    ax[0].set_aspect('equal')
    ax[1].set_ylabel('Velocity profile')
    for i in range(len(vel)):
        vel2d = vel[i,:,:,0,0]
        ax[0].contourf(vel2d.T, cmap='Reds')
        points, = ax[1].plot(vel2d[0], 'k.')
        plt.pause(0.5)
        points.remove()
    plt.show()
    #plt.savefig("test.png", dpi=500, format='png')


vel = read_data('data')
plot(vel)
