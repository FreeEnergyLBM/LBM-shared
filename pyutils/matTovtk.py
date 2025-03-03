import os
import struct
import numpy as np

def read_header(header_path):
    with open(header_path, 'rb') as headerFile:
        lx = struct.unpack('=i', headerFile.read(4))[0]
        ly = struct.unpack('=i', headerFile.read(4))[0]
        lz = struct.unpack('=i', headerFile.read(4))[0]
        ndim = struct.unpack('=i', headerFile.read(4))[0]
    return lx, ly, lz, ndim

def read_mat_file(file_path, lx, ly, lz, spatial_dim):
    if 'Velocity' in file_path:
        n_components_file = spatial_dim
        n_components_vtk = 3
    else:
        n_components_file = 1
        n_components_vtk = 1

    data = np.zeros((lx, ly, lz, n_components_vtk))
    expected_size = lx * ly * lz * n_components_file * 8
    actual_size = os.path.getsize(file_path)
    
    if actual_size != expected_size:
        raise ValueError(f"Size mismatch in {os.path.basename(file_path)}")

    with open(file_path, 'rb') as f:
        for k in range(lx * ly * lz):
            x = k // (ly * lz)
            y = (k % (ly * lz)) // lz
            z = k % lz
            if 'Velocity' in file_path:
                for i in range(n_components_file):
                    data[x, y, z, i] = struct.unpack('=d', f.read(8))[0]
            else:
                data[x, y, z, 0] = struct.unpack('=d', f.read(8))[0]
    return data

def write_vtk(data, lx, ly, lz, property_name, time, output_dir):
    vtk_filename = f"{property_name}_t{time}.vtk"
    output_path = os.path.join(output_dir, vtk_filename)
    
    with open(output_path, 'w') as vtk_file:
        vtk_file.write("# vtk DataFile Version 3.0\n")
        vtk_file.write(f"{property_name} at timestep {time}\n")
        vtk_file.write("ASCII\n")
        vtk_file.write("DATASET STRUCTURED_POINTS\n")
        vtk_file.write(f"DIMENSIONS {lx} {ly} {lz}\n")
        vtk_file.write("ORIGIN 0 0 0\n")
        vtk_file.write("SPACING 1 1 1\n")
        vtk_file.write(f"POINT_DATA {lx * ly * lz}\n\n")

        if 'Velocity' in property_name:
            vtk_file.write("VECTORS Velocity float\n")
            for z_coord in range(lz):
                for y_coord in range(ly):
                    for x_coord in range(lx):
                        vx = data[x_coord, y_coord, z_coord, 0]
                        vy = data[x_coord, y_coord, z_coord, 1]
                        vz = data[x_coord, y_coord, z_coord, 2]
                        vtk_file.write(f"{vx:.8f} {vy:.8f} {vz:.8f}\n")
        else:
            vtk_file.write(f"SCALARS {property_name} float\n")
            vtk_file.write("LOOKUP_TABLE default\n")
            for z_coord in range(lz):
                for y_coord in range(ly):
                    for x_coord in range(lx):
                        vtk_file.write(f"{data[x_coord, y_coord, z_coord, 0]:.8f}\n")

def process_directory(directory, current_dir):
    header_path = os.path.join(directory, "Header.mat")
    lx, ly, lz, spatial_dim = read_header(header_path)
    
    # Create output directory
    vtk_dir = os.path.join(current_dir, "data_vtk")
    os.makedirs(vtk_dir, exist_ok=True)  # Create if doesn't exist
    
    for filename in os.listdir(directory):
        if filename.endswith(".mat") and filename != "Header.mat":
            if '_t' not in filename:
                continue

            parts = filename.rsplit('_t', 1)
            property_name = parts[0]
            time = parts[1].split('.mat')[0]
            file_path = os.path.join(directory, filename)
            
            data = read_mat_file(file_path, lx, ly, lz, spatial_dim)
            write_vtk(data, lx, ly, lz, property_name, time, vtk_dir)
            print(f"Generated: {os.path.join(vtk_dir, f'{property_name}_t{time}.vtk')}")

# Example usage
directory = "data"
current_dir = os.path.dirname(os.path.abspath(__file__))
process_directory(directory, current_dir)