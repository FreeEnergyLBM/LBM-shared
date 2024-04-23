#!/usr/bin/env python3

# This script is used to visualise the output of the simulation and measures the contact angle

import math
import struct
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


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
    phi = np.zeros((len(times), lx, ly, lz))
    vel = np.zeros((len(times), lx, ly, lz, ndim))

    for it, t in enumerate(times):
        try:
            phi_file = open(direc+"/OrderParameter_t%li.mat"%t, 'rb')
            vel_file = open(direc+"/Velocity_t%li.mat"%t, 'rb')

            for k in range(lx*ly*lz):
                (xk,yk,zk) = coord_k(k,ly,lz)
                phi[it,xk,yk,zk] = struct.unpack('=d', phi_file.read(8))[0]
                for i in range(ndim):
                    vel[it,xk,yk,zk,i] = struct.unpack('=d', vel_file.read(8))[0]
            
            phi_file.close()
            vel_file.close()
        except FileNotFoundError:
            phi = phi[:it]
            vel = vel[:it]
            break

    return phi, vel


def plot(phi):
    print('Plotting...')
    for i in range(len(phi)):
        phi2d = phi[i,:,:,0]
        plt.contourf(phi2d.T, cmap='Blues', zorder=i,levels=5)
        plt.contour(phi2d.T, levels=[0], colors='k', zorder=i)
        #plt.pause(0.5)
    plt.gca().set_aspect('equal')
    plt.show()


def measure_angle(phi, it=-1):
    hsolid = 2
    phi2d = phi[it,:,hsolid:-hsolid,0]
    lx, ly = phi2d.shape

    def phi_circ(xy, x0, y0, r0, w):
        x = xy // ly
        y = xy % ly
        r = np.sqrt((x-x0)**2 + (y-y0)**2)
        return np.tanh((r0-r)/w)

    xy = np.arange(lx*ly)
    popt, _ = curve_fit(phi_circ, xy, phi2d.ravel(), [lx/2,0,20,1])
    x0, y0, r0, w0 = popt
    theta = 180/np.pi * np.arccos((-0.5-y0)/r0)
    print(f'Measured contact angle: {theta:g} degrees')


phi, vel = read_data('data')
measure_angle(phi)
plot(phi)
