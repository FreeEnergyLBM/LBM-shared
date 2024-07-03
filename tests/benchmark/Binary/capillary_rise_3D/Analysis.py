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
        plt.contourf(phi2d.T, cmap='Blues', zorder=i, levels=5)
        plt.contour(phi2d.T, levels=[0], colors='k', zorder=i)
        #plt.pause(0.5)
    plt.gca().set_aspect('equal')
    plt.show()
    
def compute_numerical_capillary_rise(phi, it  = -1):
    print('Computing capillary rise...')
    # Finding y_1 at the latest time step
    tol = 0.5
    y_1_index = np.where(np.abs(phi[it, 0, 2:-2, 0]) < tol)[0][0] + 2
    y_1 = y_1_index
    # Finding y_2 at the latest time step
    y_2_index = np.where(np.abs(phi[it, int(phi.shape[1] / 2), 2:-2, int(phi.shape[1] / 2)]) < tol)[0][0] + 2
    y_2 = y_2_index
    h_num = y_2 - y_1
    print(f'y_1: {y_1}, y_2: {y_2}')
    print(f'Numerical capillary rise: {h_num} grid points')
    return h_num


def compute_analytical_capillary_rise(lx, A, theta, rhog):   
    kappa = 2 * A
    r = lx / 8
    IFT = np.sqrt(8/9*kappa*A)
    h_ana = 2 * IFT * np.cos(np.deg2rad(theta)) / (r * rhog)
    print(f'Analytical capillary rise: {h_ana:.2f} grid points')
    return h_ana


phi, vel = read_data('data')
plot(phi)
h_num = compute_numerical_capillary_rise(phi)
h_ana = compute_analytical_capillary_rise(lx = 40, A = 1.00E-03, theta = 30, rhog = 1.0E-05)
print(f'The relative difference is: {abs(h_ana - h_num) / h_ana * 100:.2f}%')

