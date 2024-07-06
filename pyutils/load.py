import numpy as np
import struct


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
