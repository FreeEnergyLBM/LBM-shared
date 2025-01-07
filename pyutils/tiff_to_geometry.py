import os
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function to modify the geometry data
def modify_geometry(geometry_data):
    modified_geometry = geometry_data.copy()
    ny, nx, nz = geometry_data.shape
    
    for k in range(nz):
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                if geometry_data[j, i, k] == 1:  # Boundary node
                    if nz == 1:  # 2D case
                        if np.sum(geometry_data[max(j-1, 0):min(j+2, ny), max(i-1, 0):min(i+2, nx), k] == 1) == 9:
                            modified_geometry[j, i, k] = -1
                    else:  # 3D case
                        if np.sum(geometry_data[max(j-1, 0):min(j+2, ny), max(i-1, 0):min(i+2, nx), max(k-1, 0):min(k+2, nz)] == 1) == 27:
                            modified_geometry[j, i, k] = -1
    return modified_geometry


# Main script to process the image sequence
def process_image_sequence(image_dir='ImageSeq', output_filename='geometry.dat'):
    # Get a list of all TIFF files in the directory
    image_files = [f for f in os.listdir(image_dir) if f.endswith('.tif')]
    num_images = len(image_files)

    if num_images == 0:
        raise ValueError("No TIFF images found in the directory.")

    # Read the first image to get its size
    first_image = Image.open(os.path.join(image_dir, image_files[0]))
    first_image = np.array(first_image)
    ny, nx = first_image.shape

    # Initialize the geometry data matrix
    geometry_data = np.zeros((ny, nx, num_images), dtype=int)

    # Read each image and store the data in the 3D matrix
    for i, filename in enumerate(image_files):
        img = Image.open(os.path.join(image_dir, filename))
        img = np.array(img)

        # Convert the image to binary (0 or 255)
        img = (img == 255).astype(int)

        # Invert the image so that solid nodes are 1 and fluids are 0
        img = np.logical_not(img).astype(int)

        # Store the image data in the 3D geometry matrix
        geometry_data[:, :, i] = img

    # Modify the geometry data
    modified_geometry = modify_geometry(geometry_data)

    # Write the modified data to a new file
    with open(output_filename, 'w') as f:
        for k in range(modified_geometry.shape[2]):
            for j in range(modified_geometry.shape[0]):
                for i in range(modified_geometry.shape[1]):
                    f.write(f"{modified_geometry[j, i, k]}\n")

    print(f"geometry.dat file generated successfully.")

    # Visualize the modified geometry matrix in 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y, z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(num_images))

    # Create slices of the 3D data
    ax.scatter(x[modified_geometry == -1], y[modified_geometry == -1], z[modified_geometry == -1], c='b', marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Grains')
    plt.show()

# Run the script
if __name__ == "__main__":
    process_image_sequence(image_dir='geometry2d', output_filename='geometry.dat')
