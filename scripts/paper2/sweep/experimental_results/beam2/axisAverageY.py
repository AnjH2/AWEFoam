import os
import re
import numpy as np
import matplotlib.pyplot as plt

# Enable LaTeX fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]
})

# Function to extract number from filename for sorting
def extract_number(filename):
    match = re.search(r'(\d+)', filename)
    return int(match.group(1)) if match else float('inf')

folder_path = '.'  # Folder with your .npy files

# Sort by numeric content in the file name
file_list = [f for f in os.listdir(folder_path) if f.endswith('.npy')]
file_list = sorted(file_list, key=extract_number)

x_max_mm = 4.10
y_max_mm = 60
i_total = 30

def slice_and_average_by_mm(array_2d, x_min_mm, x_max_mm, y_min_mm, y_max_mm, pixel_size):
    ny, nx = array_2d.shape

    # Convert mm -> pixel in X
    x_min_px = int(x_min_mm * pixel_size[0])
    x_max_px = int(x_max_mm * pixel_size[0])
    # Convert mm -> pixel in Y
    y_min_px = int(y_min_mm * pixel_size[1])
    y_max_px = int(y_max_mm * pixel_size[1])

    # Clamp to array bounds
    x_min_px = max(0, x_min_px)
    x_max_px = min(nx, x_max_px)
    y_min_px = max(0, y_min_px)
    y_max_px = min(ny, y_max_px)

    if x_min_px >= x_max_px or y_min_px >= y_max_px:
        print("Warning: Slice is empty.")
        return np.nan

    sub_array = array_2d[y_min_px:y_max_px, x_min_px:x_max_px]
    
    return np.nanmean(sub_array)

# Define two regions: (region_name, x_lo_mm, x_hi_mm)
regions = [
    ("cathode",   0.0,  1.8),
    ("anode", 2.3, 4.10)
]

region_colors = {
    "cathode":   "blue",
    "anode": "red"
}

# Ensure an output directory for text files
save_txt_dir = os.path.join(folder_path, "extractedData")
os.makedirs(save_txt_dir, exist_ok=True)

for file_name in file_list:
    img_array = np.load(os.path.join(folder_path, file_name))

    # Filter out values above 30% of the maximum
    vmax_ = np.max(img_array)
    filtered_img_array = np.where(img_array > vmax_ * 0.3, vmax_ * 0.3, img_array)

    # If needed, rotate:
    rotated_img_array = filtered_img_array

    max_y = rotated_img_array.shape[0]
    max_x = rotated_img_array.shape[1]
    pixel_size = [max_x / x_max_mm, max_y / y_max_mm]

    # region_data = { region_name: {"y": [...], "avg": [...]} }
    region_data = {}
    for (region_name, xlo_mm, xhi_mm) in regions:
        region_data[region_name] = {"y": [], "avg": []}

    # Loop over i_total y-slices
    for i in range(i_total):
        L1 = (i / i_total) * y_max_mm
        L2 = ((i + 1) / i_total) * y_max_mm

        # For each region, compute average
        for (region_name, xlo_mm, xhi_mm) in regions:
            avg_val = slice_and_average_by_mm(
                rotated_img_array,
                x_min_mm = xlo_mm,
                x_max_mm = xhi_mm,
                y_min_mm = L1,
                y_max_mm = L2,
                pixel_size = pixel_size
            )
            y_mid = 0.5 * (L1 + L2)

            region_data[region_name]["y"].append(y_mid)
            region_data[region_name]["avg"].append(avg_val)

    # PLOT
    fig, ax = plt.subplots(figsize=(4, 3))

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
    match_cfg = re.search(r'(flat|channel|conical)', base_title)
    cfg_val = match_cfg.group(1) if match_cfg else '???'

    # Build a "nice" title from these pieces
    title_str = (
        fr"i={j_val}, "#+r"[$\frac{mA}{cm^2}$],"   # i=800 [mA/cm^2]
        fr"Q={Q_val}, "#+r"[$\frac{mL}{min}$],"   # Q=110 [mL/min]
        fr"type={cfg_val}"
    )
    print(title_str)
    ax.set_title(title_str, fontsize=14)

    # Scatter plot each region in distinct color
    for (region_name, xlo_mm, xhi_mm) in regions:
        y_vals   = region_data[region_name]["y"]
        avg_vals = region_data[region_name]["avg"]
        color    = region_colors.get(region_name, "black")
        ax.scatter(y_vals, avg_vals, color=color, marker='o',
                   label=region_name, alpha=0.7)
        ax.set_ylim([0,70])
        # --------------------------------------------------------
        # Write out y_vals vs avg_vals to a text file
        # --------------------------------------------------------
        # Example: "myfile_anode.txt" or "myfile_cathode.txt"
        base_no_ext = file_name[:-4]  # remove ".npy"
        out_txt_name = f"{base_no_ext}_{region_name}.txt"
        out_txt_path = os.path.join(save_txt_dir, out_txt_name)

        with open(out_txt_path, 'w') as f:
            f.write("# y_mid[mm], AverageValue\n")
            for (yv, av) in zip(y_vals, avg_vals):
                f.write(f"{yv:.6g}\t{av:.6g}\n")  # e.g. tab-separated

    ax.set_xlabel("Y [mm]")
    ax.set_ylabel("Average Value")
    ax.legend()

    plt.tight_layout()
    plt.show()

