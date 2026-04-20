import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import re

# Enable LaTeX-style fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]
})

# Path to the folder containing .npy files
folder_path = '.'

# Function to extract number from filename for sorting
def extract_number(filename):
    match = re.search(r'(\d+)', filename)
    return int(match.group(1)) if match else float('inf')

# Get all .npy files in the folder and sort by number in the filename
file_list = sorted([file for file in os.listdir(folder_path) if file.endswith('.npy')], key=extract_number)

# Define figure size and other constants
num_files = len(file_list)
fig_height = 6
fig_width = 3 * num_files  # Adjust width based on number of files
xStrech = 6

# Desired number of intervals for the colorbar and image
num_intervals = 11
vmax=55
boundariesC = np.linspace(0, vmax, num_intervals + 1)  # Define boundaries for intervals
boundariesI = np.linspace(0, 50, num_intervals + 1)  # Define boundaries for intervals


# Create a figure with subplots
fig, axes = plt.subplots(1, num_files, figsize=(fig_width, fig_height), constrained_layout=True)

for i, file_name in enumerate(file_list):
    # Load each .npy file
    img_array = np.load(os.path.join(folder_path, file_name))
    
    # Filter out values above 30% of the maximum
    vmax_ = np.max(img_array)
    filtered_img_array = np.where(img_array > vmax_ * 0.3, vmax_ * 0.3, img_array)
    
    # Rotate the image 180 degrees clockwise
    rotated_img_array = np.rot90(filtered_img_array, 2)
    
    # Calculate dimensions
    max_y = rotated_img_array.shape[0]
    min_y = 0
    max_x = rotated_img_array.shape[1] * xStrech
    min_x = 0
    
    # Define tick positions
    num_ticks = 5
    tick_yPositions = np.linspace(max_y - 1, min_y, num_ticks)
    yTicks = np.linspace(vmax, 0, num_ticks)
    xTicks = [0, 1.8, 1.8 + 0.65, 1.8 * 2 + 0.65]
    tick_xPositions = [(max_x / (1.8 * 2 + 0.65)) * n for n in xTicks]
    
    # Create a BoundaryNorm instance for both image and colorbar
    normC = colors.BoundaryNorm(boundariesC, plt.cm.viridis.N)
    normI = colors.BoundaryNorm(boundariesI, plt.cm.viridis.N)
    
    # Plotting on each subplot
    ax = axes[i]
    im = ax.imshow(rotated_img_array, cmap='plasma', extent=[min_x, max_x, min_y, max_y], norm=normI, interpolation='none')

    # Custom axis labels and ticks
    xlabel = r'Width [mm]'
    xlabel = re.sub(r'\[(.*)\]', r'[$\\mathrm{\1}$]', xlabel)  # Convert units to non-italic
    ax.set_xlabel(xlabel, fontsize=12)
    if i == 0:
        ylabel = r'Length [mm]'
        ylabel = re.sub(r'\[(.*)\]', r'[$\\mathrm{\1}$]', ylabel)  # Convert units to non-italic
        ax.set_ylabel(ylabel, fontsize=12)
    else:
        ax.set_yticklabels([])
    ax.set_yticks(tick_yPositions)
    ax.set_yticklabels([f'{tick:.2f}' for tick in yTicks], fontsize=10)
    ax.set_xticks(tick_xPositions)
    ax.set_xticklabels([f'{tick:.2f}' for tick in xTicks], fontsize=10)
    
    # Add vertical lines for membrane interfaces
    ax.axvline(x=1.8 * (max_x / (1.8 * 2 + 0.65)), color='b', linestyle='--', linewidth=1, label='membrane interface')
    ax.axvline(x=(1.8 + 0.65) * (max_x / (1.8 * 2 + 0.65)), color='b', linestyle='--', linewidth=1)
    
    # Title for each subplot (filename without .npy extension)
    title = file_name[:-4]  # Remove .npy extension
    title = re.sub(r'\[(.*)\]', r'[\\mathrm{\1}]', title)  # Convert units to non-italic
    ax.set_title(r'{}'.format(title), fontsize=14)

# Add colorbar to the right of all subplots
cbar = fig.colorbar(im, ax=axes, boundaries=boundariesC, ticks=boundariesC, location='right')
barTitle = r'Volumefraction [\%]'
barTitle = re.sub(r'\[(.*)\]', r'[$\\mathrm{\1}$]', barTitle) 
cbar.set_label(barTitle, fontsize=12)  # Set label for colorbar

# Add overall title
fig.suptitle(r'Volumefraction Distribution', fontsize=16)

# Adjust overall layout and display the figure
plt.subplots_adjust(wspace=0.4)  # Adjust space between subplots
plt.show()

