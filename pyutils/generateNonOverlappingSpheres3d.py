# Description: This script generates a 3D image with non-overlapping spheres.
# The spheres are placed randomly in the domain with a minimum distance between them.
# The radii of the spheres are sampled from a truncated log-normal distribution.
# The image slices are saved as TIFF files in the specified output folder.
# The frequency distribution of the placed radii is also plotted.

import os
import numpy as np
import matplotlib.pyplot as plt
from imageio import imwrite

# Enable LaTeX-style fonts
plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family": "cmr10",
    "axes.linewidth": 1.8,  # Thicker axes lines
    "axes.formatter.use_mathtext": True,  # Use mathtext for scientific notation
    
})

def generate_spheres_3d(Nx=200, Ny=200, Nz=200, m=200, min_distance=5, min_radius=10, max_radius=40, x_offset=2, y_offset=2, z_offset=2, output_folder='geometry3d', max_trials=1000):
    # Calculate the new dimensions with offsets
    Nx_new = Nx + 2 * x_offset
    Ny_new = Ny + 2 * y_offset
    Nz_new = Nz + 2 * z_offset

    # Create an empty binary 3D image with the new dimensions
    img = 255 * np.ones((Nz_new, Ny_new, Nx_new), dtype=np.uint8)

    # Set the offset regions to black, leaving a 2-pixel margin
    img[:, :, :x_offset] = 0  # Left offset
    img[:, :, -x_offset:] = 0  # Right offset
    img[:, :y_offset, :] = 0  # Front offset
    img[:, -y_offset:, :] = 0  # Back offset
    img[:z_offset, :, :] = 0  # Bottom offset
    img[-z_offset:, :, :] = 0  # Top offset

    # Generate m random radii from a truncated log-normal distribution
    mu = 0  # Mean of the log-normal distribution
    sigma = 1  # Standard deviation of the log-normal distribution
    radii = []  # Store all radii
    radii_placed = []  # Store radii of spheres that were successfully placed

    for _ in range(m):
        while True:
            r = np.random.lognormal(mean=mu, sigma=sigma)
            if min_radius <= r <= max_radius:
                radii.append(r)
                break

    # Sort radii in descending order to place larger spheres first
    radii = sorted(radii, reverse=True)

    # Initialize an array to store sphere positions
    sphere_positions = []

    # Generate sphere positions with trial mechanism
    for r in radii:
        placed = False
        trials = 0
        while not placed and trials < max_trials:
            x = x_offset + r + 2 + (Nx - 2 * r - 4) * np.random.rand()
            y = y_offset + r + 2 + (Ny - 2 * r - 4) * np.random.rand()
            z = z_offset + r + 2 + (Nz - 2 * r - 4) * np.random.rand()

            if not sphere_positions:
                # If no spheres yet, place the first sphere
                sphere_positions.append((x, y, z, r))
                placed = True
            else:
                # Check if the sphere overlaps with existing spheres
                if all(np.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2) >= (r + cr + min_distance) for cx, cy, cz, cr in sphere_positions):
                    # If no overlap, add the sphere
                    sphere_positions.append((x, y, z, r))
                    radii_placed.append(r)
                    placed = True
            trials += 1

        if not placed:
            print(f"Skipping sphere with radius {r:.2f} after {max_trials} trials.")

    # Draw the spheres on the 3D image
    X, Y, Z = np.meshgrid(np.arange(Nx_new), np.arange(Ny_new), np.arange(Nz_new), indexing='ij')
    for cx, cy, cz, r in sphere_positions:
        mask = (X - cx)**2 + (Y - cy)**2 + (Z - cz)**2 <= r**2
        img[mask] = 0

    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    for file in os.listdir(output_folder):
        os.remove(os.path.join(output_folder, file))

    # Save each slice of the 3D image as a TIFF file
    for z in range(Nz_new):
        filepath = os.path.join(output_folder, f'slice_{z:03d}.tif')
        imwrite(filepath, img[z, :, :])

    # Plot the frequency distribution of the radii
    plt.figure()
    plt.hist(radii_placed, bins=np.arange(min_radius, max_radius + 1), edgecolor='black')
    plt.title('Frequency Distribution of Sphere Radii')
    plt.xlabel('Radius (px)')
    plt.ylabel('Frequency (-)')
    plt.show()

# Main script
if __name__ == "__main__":
    Nx = 200  # Original width of the domain
    Ny = 200  # Original height of the domain
    Nz = 200  # Original depth of the domain
    m = 100  # Number of spheres
    min_distance = 5  # Minimum distance between spheres
    min_radius = 10  # Minimum radius of spheres
    max_radius = 30  # Maximum radius of spheres
    x_offset = 2  # Offset in x-direction
    y_offset = 2  # Offset in y-direction
    z_offset = 2  # Offset in z-direction
    output_folder = "geometry3d"  # Output folder for image slices

    generate_spheres_3d(Nx, Ny, Nz, m, min_distance, min_radius, max_radius, x_offset, y_offset, z_offset, output_folder)
