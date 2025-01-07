# Description: This script generates a 2D binary image with non-overlapping circles.
# The circles are placed randomly in the domain with a minimum distance between them.
# The radii of the circles are sampled from a truncated log-normal distribution.
# The image is saved as a TIFF file in the specified output folder.
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
})

def generate_circles_2d(Nx = 200, Ny = 200, m = 200, min_distance = 5, min_radius = 10, max_radius = 40, x_offset = 3, y_offset = 3, output_folder = 'geometry2d', max_trials=1000):
    # Calculate the new dimensions with offsets
    Nx_new = Nx + 2 * x_offset
    Ny_new = Ny + 2 * y_offset

    # Create an empty binary image with the new dimensions
    img = 255 * np.ones((Ny_new, Nx_new), dtype=np.uint8)

    # Set the offset regions to black, leaving a 2-pixel margin
    img[:, :x_offset] = 0  # Left offset
    img[:, -x_offset:] = 0  # Right offset
    img[:y_offset, :] = 0  # Top offset
    img[-y_offset:, :] = 0  # Bottom offset

    # Generate m random radii from a truncated log-normal distribution
    mu = 0  # Mean of the log-normal distribution
    sigma = 1  # Standard deviation of the log-normal distribution
    radii = [] # Store all radii
    radii_placed = [] # Store radii of circles that were successfully placed

    for _ in range(m):
        while True:
            r = np.random.lognormal(mean=mu, sigma=sigma)
            if min_radius <= r <= max_radius:
                radii.append(r)
                break

    # Sort radii in descending order to place larger circles first
    radii = sorted(radii, reverse=True)

    # Initialize an array to store circle positions
    circle_positions = []

    # Generate circle positions with trial mechanism
    for r in radii:
        placed = False
        trials = 0
        while not placed and trials < max_trials:
            x = x_offset + r + 2 + (Nx - 2 * r - 4) * np.random.rand()
            y = y_offset + r + 2 + (Ny - 2 * r - 4) * np.random.rand()

            if not circle_positions:
                # If no circles yet, place the first circle
                circle_positions.append((x, y, r))
                placed = True
            else:
                # Check if the circle overlaps with existing circles
                if all(np.sqrt((x - cx)**2 + (y - cy)**2) >= (r + cr + min_distance) for cx, cy, cr in circle_positions):
                    # If no overlap, add the circle
                    circle_positions.append((x, y, r))
                    radii_placed.append(r)
                    placed = True
            trials += 1

        if not placed:
            print(f"Skipping circle with radius {r:.2f} after {max_trials} trials.")

    # Draw the circles on the image
    X, Y = np.meshgrid(np.arange(Nx_new), np.arange(Ny_new))
    for cx, cy, r in circle_positions:
        mask = (X - cx)**2 + (Y - cy)**2 <= r**2
        img[mask] = 0

    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    for file in os.listdir(output_folder):
        os.remove(os.path.join(output_folder, file))

    # Save the image as a TIFF file
    filepath = os.path.join(output_folder, 'circles_image_with_inlet_outlet.tif')
    imwrite(filepath, img)

    # Plot the frequency distribution of the radii
    plt.figure()
    plt.hist(radii_placed, bins=np.arange(min_radius, max_radius + 1), edgecolor='black')
    plt.title('Frequency Distribution of Circle Radii')
    plt.xlabel('Radius (px)')
    plt.ylabel('Frequency (-)')
    plt.show()
    
    

# Main script
if __name__ == "__main__":
    Nx = 200  # Original width of the domain
    Ny = 200  # Original height of the domain
    m = 100  # Number of spheres
    min_distance = 5  # Minimum distance between spheres
    min_radius = 10  # Minimum radius of spheres
    max_radius = 30  # Maximum radius of spheres
    x_offset = 2  # Offset in x-direction
    y_offset = 2  # Offset in y-direction
    output_folder = "geometry2d"  # Output folder for image slices

    generate_circles_2d(Nx, Ny, m, min_distance, min_radius, max_radius, x_offset, y_offset, output_folder)
