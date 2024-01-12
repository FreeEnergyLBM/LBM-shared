#!/usr/bin/env python3

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

force = 1e-4 / 49 # pdiff / length
viscosity = 1/6

error_criteria = 1e-2


def read_file(filename):
    def read_field(f, lx, ly, lz, data):
        next(f)
        structure, name, data_type = f.readline().lower().split()
        if (structure=='scalars'):
            data[name] = np.zeros(lx*ly*lz, dtype=data_type)
            next(f)
            for i in range(lx*ly*lz):
                data[name][i] = f.readline()
            data[name] = data[name].reshape((lx,ly,lz))
        elif (structure=='vectors'):
            data[name] = np.zeros((lx*ly*lz, 3), dtype=data_type)
            for i in range(lx*ly*lz):
                data[name][i] = f.readline().split()
            data[name] = data[name].reshape((lx,ly,lz,3))

    data = {}
    with open(filename, 'r') as f:
        # Header
        for _ in range(4): next(f)
        lz, ly, lx = map(int, f.readline().split()[1:])
        for _ in range(3): next(f)
        # Fields
        while (True):
            try:
                read_field(f, lx, ly, lz, data)
            except StopIteration:
                break
        return data


def load_data(direc='data'):
    extract_number = lambda filename: int(re.search(r'\d+', filename).group())
    files = sorted(glob.glob(direc+'/data_*.vtk'), key=extract_number)
    data = read_file(files[-1])
    return data


def poiseuille_profile(width):
    y = np.arange(0.5, width)
    return force / (2 * viscosity) * y * (width - y)


def compare_velocity(vel):
    vel_poiseuille = poiseuille_profile(len(vel))
    relative_error = np.abs(vel - vel_poiseuille) / vel_poiseuille
    return np.max(relative_error)


data = load_data()
vel_slice = data['velocity'][0,2:-2,0,0]
max_error = compare_velocity(vel_slice)

if (max_error > error_criteria):
    print(f'Error: Max relative error in velocity profile too large ({max_error:g} > {error_criteria:g})')
    sys.exit(1)
