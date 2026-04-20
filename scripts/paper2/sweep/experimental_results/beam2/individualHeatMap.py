import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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
file_list = [file for file in os.listdir(folder_path) if file.endswith('.npy')]
file_list.sort(key=extract_number)

# Desired number of intervals for the colorbar and image
num_intervals = 10
vmax = 60
cmax = 100
boundariesC = np.linspace(0, cmax, num_intervals + 1)  # for colorbar ticks
boundariesI = np.linspace(0, cmax, num_intervals + 1)    # for normI

xStrech = 6  # (from your original code, mm per pixel in x)

for file_name in file_list:
    # ---------------------------
    # 1) Load and preprocess
    # ---------------------------
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

    # Define tick positions (same logic as your original)
    num_ticks = 5
    tick_yPositions = np.linspace(max_y - 1, min_y, num_ticks)
    yTicks = np.linspace(vmax, 0, num_ticks)
    #xTicks = [0, 1.8, 1.8 + 0.5, 1.8 * 2 + 0.5]
    xTicks = [0, 1.8 * 2 + 0.5]
    tick_xPositions = [(max_x / (1.8 * 2 + 0.5)) * n for n in xTicks]

    # Create colormap normalization
    normC = colors.BoundaryNorm(boundariesC, plt.cm.viridis.N)
    normI = colors.BoundaryNorm(boundariesI, plt.cm.viridis.N)

    # ---------------------------
    # 2) Create a figure (one per file)
    # ---------------------------
    fig, ax = plt.subplots(figsize=(3, 6), constrained_layout=True)

    # Show image
    im = ax.imshow(rotated_img_array, cmap='plasma',
                   extent=[min_x, max_x, min_y, max_y],
                   norm=normI, interpolation='none')

    # Custom axis labels and ticks
    xlabel = r'Width [mm]'
    xlabel = re.sub(r'\[(.*)\]', r'[$\\mathrm{\1}$]', xlabel)
    ax.set_xlabel(xlabel, fontsize=12)

    ylabel = r'Length [mm]'
    ylabel = re.sub(r'\[(.*)\]', r'[$\\mathrm{\1}$]', ylabel)
    ax.set_ylabel(ylabel, fontsize=12)

    ax.set_yticks(tick_yPositions)
    ax.set_yticklabels([f'{tick:.2f}' for tick in yTicks], fontsize=10)
    ax.set_xticks(tick_xPositions)
    ax.set_xticklabels([f'{tick:.2f}' for tick in xTicks], fontsize=10)

    # Add vertical lines for membrane interfaces
    ax.axvline(x=1.8 * (max_x / (1.8 * 2 + 0.65)), color='b', linestyle='--', linewidth=1, label='membrane interface')
    ax.axvline(x=(1.8 + 0.65) * (max_x / (1.8 * 2 + 0.65)), color='b', linestyle='--', linewidth=1)

    # Title for each figure (filename without .npy extension)
    #plt.show()
    # plt.close(fig)  # If you prefer not to show interactively, uncomment this
    base_title = os.path.splitext(file_name)[0]  # Remove the '.npy' extension

    # Extract jXYZ from the filename (e.g., j800)
    match_j = re.search(r'j(\d+)', base_title)
    j_val = match_j.group(1) if match_j else '???'

    # Extract QXYZ from the filename (e.g., Q110)
    match_Q = re.search(r'Q(\d+)', base_title)
    Q_val = match_Q.group(1) if match_Q else '???'

    # Extract the configuration (flat, channel, etc.)
    # Adjust this if you have different possible configurations
    match_cfg = re.search(r'(flat|channel|somethingElse)', base_title)
    cfg_val = match_cfg.group(1) if match_cfg else '???'

    # Build a "nice" title from these pieces
    title_str = (
        fr"i={j_val}, "#+r"[$\frac{mA}{cm^2}$],"   # i=800 [mA/cm^2]
        fr"Q={Q_val}, "#+r"[$\frac{mL}{min}$],"   # Q=110 [mL/min]
        fr"type={cfg_val}"
    )
    print(title_str)
    ax.set_title(title_str, fontsize=14)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, boundaries=boundariesC, ticks=boundariesC, location='right',shrink=0.6)
    barTitle = 'Volumefraction [%]'  # no \mathrm in brackets
    cbar.set_label(barTitle, fontsize=12)

    # Optionally add a main (super) title
    fig.suptitle(r'Volumefraction Distribution', fontsize=16)

    # Save the figure as .png, .svg, .pdf (using the file base name)
    out_name = os.path.splitext(file_name)[0]  # remove .npy extension
    fig.savefig(out_name + ".png", dpi=300)
    fig.savefig(out_name + ".svg", format='svg')
    fig.savefig(out_name + ".pdf", format='pdf')


