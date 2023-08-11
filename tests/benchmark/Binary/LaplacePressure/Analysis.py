#!/usr/bin/env python3

import sys
import math
import struct
import numpy as np

surface_tension = np.sqrt(8*0.0003*0.00015/9)
radius = 20

error_criteria = 1e-2



def read_file(filename, data):
    with open(filename, 'rb') as file:
        for i in range(data.size):
            data.ravel()[i] = struct.unpack('=d', file.read(8))[0]


def load_density(direc='data'):
    with open(direc+"/Header.mat", 'rb') as headerFile:
        lx = struct.unpack('=i', headerFile.read(4))[0]
        ly = struct.unpack('=i', headerFile.read(4))[0]
        lz = struct.unpack('=i', headerFile.read(4))[0]
        ndim = struct.unpack('=i', headerFile.read(4))[0]
        tend = struct.unpack('=i', headerFile.read(4))[0]
        tinc = struct.unpack('=i', headerFile.read(4))[0]

    times = np.arange(0, tend+1, tinc)
    density = np.zeros((len(times), lx, ly, lz))
    for it, t in enumerate(times):
        read_file(direc+"/Density_t%li.mat"%t, density[it])
    return density[:,:,:,0]


def compare_pressure(density):
    centre = np.array(density.shape) // 2
    pdiff = (density[tuple(centre)] - density[0,0]) / 3
    pdiff_theory = surface_tension / radius
    return abs(pdiff - pdiff_theory) / pdiff_theory


density = load_density('data')
error = compare_pressure(density[-1])

if (error > error_criteria):
    print(f'Error: Relative error in pressure difference too large ({error:g} > {error_criteria:g})')
    sys.exit(1)
