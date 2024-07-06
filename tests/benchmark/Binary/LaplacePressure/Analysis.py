#!/usr/bin/env python3

import sys
import math
import struct
import numpy as np

a = 0.00015
kappa = 0.0003
surface_tension = np.sqrt(8/9*a*kappa)
radius = 20

error_criteria = 0.05



def read_file(filename, data):
    with open(filename, 'rb') as file:
        for i in range(data.size):
            data.ravel()[i] = struct.unpack('=d', file.read(8))[0]


def load_parameter(parameter, direc='data', time=None):
    with open(direc+"/Header.mat", 'rb') as headerFile:
        lx = struct.unpack('=i', headerFile.read(4))[0]
        ly = struct.unpack('=i', headerFile.read(4))[0]
        lz = struct.unpack('=i', headerFile.read(4))[0]
        ndim = struct.unpack('=i', headerFile.read(4))[0]
        tend = struct.unpack('=i', headerFile.read(4))[0]
        tinc = struct.unpack('=i', headerFile.read(4))[0]

    times = np.arange(0, tend+1, tinc)
    data = np.zeros((len(times), lx, ly, lz))
    for it, t in enumerate(times):
        try:
            filename = f"{direc}/{parameter}_t{t}.mat"
            read_file(filename, data[it])
        except FileNotFoundError:
            data = data[:it]
            times = times[:it]
            break
    if (time is not None): time.extend(times)
    return data


def compare_pressure(density):
    centre = np.array(density.shape) // 2
    outside = np.zeros_like(centre)
    pdiff = (density[tuple(centre)] - density[tuple(outside)]) / 3
    pdiff_theory = surface_tension / radius
    return abs(pdiff - pdiff_theory) / pdiff_theory

def plot_errors():
    import matplotlib.pyplot as plt
    time = []
    density = load_parameter('Density', time=time)
    errors = np.zeros(len(density))
    for i in range(len(density)):
        errors[i] = compare_pressure(density[i])
    plt.figure()
    plt.plot(time, errors)
    plt.semilogy()
    plt.xlabel("Timesteps")
    plt.ylabel("Pressure relative error")
    plt.show()


density = load_parameter('Density')
error = compare_pressure(density[-1])

# plot_errors()

if (error > error_criteria):
    print(f'Error: Relative error in pressure difference too large ({error:g} > {error_criteria:g})')
    sys.exit(1)
