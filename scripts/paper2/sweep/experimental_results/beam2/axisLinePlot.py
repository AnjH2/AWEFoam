#!/usr/bin/env python3
"""
Pixel-by-pixel line extractor: x-axis = length (mm), y-axis = pixel value.

Author: chatGBT
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt

# ─── USER SETTINGS ─────────────────────────────────────────────────────────────
folder_path     = "."                 # directory with .npy files
save_dir        = "extractedData_lineX_b10"  # sub-folder for txt output
x_max_mm        = 4.45                # physical width of the frame
y_max_mm        = 60.0                # physical height of the frame

y_values_mm     = [10, 30, 50]   # heights you want (mm)
band_half_mm    = 10                  # 0 → single row; >0 → ±band_half_mm average
# ------------------------------------------------------------------------------
# leave the next three lines if you *ever* want an x-window instead of full width
# x_lo_mm, x_hi_mm = 0.0, x_max_mm      # physical window in X
# ------------------------------------------------------------------------------
# optional: LaTeX fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

def mm_to_px(mm, px_per_mm):
    return int(round(mm * px_per_mm))

def row_or_band(img, y_c_mm, half_mm, pix_sz):
    """Return a 2-D slice: shape (band_height_px, nx)."""
    ny, _ = img.shape
    if half_mm == 0.0:                       # single row
        y = mm_to_px(y_c_mm, pix_sz[1])
        y = np.clip(y, 0, ny - 1)
        return img[y:y+1, :]                 # keep 2-D dims
    # thin vertical band
    y1 = max(0,  mm_to_px(y_c_mm - half_mm, pix_sz[1]))
    y2 = min(ny, mm_to_px(y_c_mm + half_mm, pix_sz[1]))
    if y1 >= y2:
        return None
    return img[y1:y2, :]

# ensure output dir exists
os.makedirs(os.path.join(folder_path, save_dir), exist_ok=True)

# sort numerically if filenames contain integers
file_list = sorted([f for f in os.listdir(folder_path) if f.endswith(".npy")],
                   key=lambda s: int(re.search(r"(\d+)", s).group(1)) if re.search(r"(\d+)", s) else float("inf"))

for fname in file_list:
    img = np.load(os.path.join(folder_path, fname))

    # clip everything > 0.3 · max
    vmax = img.max()
    img  = np.where(img > 0.3 * vmax, 0.3 * vmax, img)

    ny, nx  = img.shape
    px_size = (nx / x_max_mm, ny / y_max_mm)   # px / mm in (x, y)

    # vector of x-centres in mm (one per column)
    x_vec_mm = (np.arange(nx) + 0.5) / px_size[0]

    # figure ---------------------------------------
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    base = fname[:-4]
    base1 = fname[:-17]
    title = re.sub(r"\[(.*?)\]", r"[\\mathrm{\1}]", base)
    ax.set_title(rf"{title}", fontsize=14)

    # one coloured line per requested height -------
    for y_mm in y_values_mm:
        band = row_or_band(img, y_mm, band_half_mm, px_size)
        if band is None:
            print(f"⚠  Empty band at y={y_mm} mm in {fname}")
            continue

        # average across band thickness; result = shape (nx,)
        vals = band.mean(axis=0)

        # save txt -----------------------------------------------------------
        out_name = f"{base1}_{y_mm}.txt"
        out_fp   = os.path.join(folder_path, save_dir, out_name)
        with open(out_fp, "w") as f:
            f.write("# x_mm\tvalue\n")
            for x, v in zip(x_vec_mm, vals):
                f.write(f"{x:.6g}\t{v:.6g}\n")

        # plot ---------------------------------------------------------------
        ax.plot(x_vec_mm, vals, "-", marker="", label=f"y = {y_mm} mm")

    ax.set_xlabel(r"$x$ [mm]")
    ax.set_ylabel("Gas volume fraction")
    ax.legend(title="Horizontal slices", loc="best", fontsize=8)
    ax.grid(True, ls="--", lw=0.4, alpha=0.4)
    plt.tight_layout()
    plt.show()

