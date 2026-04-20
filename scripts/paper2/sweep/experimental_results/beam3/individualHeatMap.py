
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import AutoMinorLocator

# ---------- CONFIG ----------
folder_path = "."
out_dir = "figures"
os.makedirs(out_dir, exist_ok=True)

# Pixel size (meters per pixel)  ⟶ 22.188 µm
PIXEL_SIZE_M = 22.188e-6
MM_PER_M = 1000.0

# Color scaling
clip_max = 0.5
norm = Normalize(vmin=0.0, vmax=clip_max)

# Figure sizing
BASE_WIDTH_IN = 4.0

# Colormap and colorbar
CMAP = "plasma"
COLORBAR_LABEL = "Volume fraction [%]"
NUM_CB_TICKS = 6      # colorbar ticks
NUM_Y_AX_TICKS = 6    # target number of y-axis ticks

# ---- Center-relative X ticks ----
# Provide your desired tick labels (in mm) RELATIVE TO THE CENTER (e.g., symmetric around 0).
X_TICKS_MM_REL = [-5,-2.5, 0.0, 2.5, 5]

# Manual center offset (in mm). Positive shifts the center to the RIGHT.
# Use this if the visual center is slightly off the geometric center.
CENTER_OFFSET_MM = 0.0

# Axis & grid style
SHOW_GRID = True
ORIGIN = "upper"  # change to "lower" if you prefer y to increase upward

# ---------- HELPERS ----------
def extract_number(filename: str):
    m = re.search(r"(\d+)", filename)
    return int(m.group(1)) if m else float("inf")

def nice_ticks(max_val, nticks=6):
    """Return 'nice' ticks from 0..max_val (inclusive) with ~nticks positions."""
    if max_val <= 0 or nticks < 2:
        return [0, max_val] if max_val > 0 else [0]
    raw_step = max_val / (nticks - 1)
    multiples = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50])
    exp = np.floor(np.log10(raw_step))
    candidates = multiples * (10 ** exp)
    step = candidates[np.argmin(np.abs(candidates - raw_step))]
    ticks = np.arange(0, max_val + 0.5 * step, step)
    if abs(ticks[-1] - max_val) > 0.25 * step:
        ticks = np.append(ticks, max_val)
    ticks = np.unique(np.round(ticks, 6))
    return ticks.tolist()

# ---------- FILES ----------
file_list = [f for f in os.listdir(folder_path) if f.endswith(".npy")]
file_list.sort(key=extract_number)
if not file_list:
    raise SystemExit("No .npy files found in the specified folder.")

# ---------- RENDER ----------
for fname in file_list:
    path = os.path.join(folder_path, fname)
    img = np.load(path)

    # Ensure 2D
    if img.ndim > 2:
        img = np.squeeze(img)
        if img.ndim > 2:
            img = img[..., 0]
    if img.ndim != 2:
        print(f"Skipping {fname}: not a 2D array after squeezing.")
        continue

    # Physical extent in mm from pixel size
    nrows, ncols = img.shape
    width_mm  = ncols * PIXEL_SIZE_M * MM_PER_M
    height_mm = nrows * PIXEL_SIZE_M * MM_PER_M
    extent = [0.0, width_mm, 0.0, height_mm]

    # Figure size preserves aspect
    aspect = height_mm / width_mm if width_mm > 0 else 1.0
    fig_h_in = BASE_WIDTH_IN * aspect
    fig, ax = plt.subplots(figsize=(BASE_WIDTH_IN, fig_h_in), constrained_layout=True)

    im = ax.imshow(
        img,
        cmap=CMAP,
        extent=extent,
        norm=norm,
        interpolation="none",
        origin=ORIGIN
    )

    # ----- AXES -----
    ax.set_xlabel("Width [mm]")
    ax.set_ylabel("Length [mm]")

    # Y ticks: auto "nice" ticks over full physical height
    ax.set_yticks(nice_ticks(height_mm, NUM_Y_AX_TICKS))

    # X ticks: user-provided RELATIVE-to-center list, placed at absolute positions
    center_mm = (width_mm / 2.0) + CENTER_OFFSET_MM
    # Convert relative list -> absolute positions, keep only those in bounds
    x_ticks_abs = [center_mm + rel for rel in X_TICKS_MM_REL
                   if 0.0 <= center_mm + rel <= width_mm]
    # Labels show the relative mm values (the list you provided)
    x_tick_labels = [f"{rel:.2f}" for rel in X_TICKS_MM_REL
                     if 0.0 <= center_mm + rel <= width_mm]
    ax.set_xticks(x_ticks_abs)
    ax.set_xticklabels(x_tick_labels)

    # optional: draw a vertical line at the (offset) center
    ax.axvline(center_mm, color="k", linestyle="--", linewidth=0.8, alpha=0.5)


# optional: minor ticks and grid
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if SHOW_GRID:
        ax.grid(True, which="major", linewidth=0.5, alpha=0.3)

    # --- Colorbar ---
    ticks = np.linspace(0, clip_max, NUM_CB_TICKS)
    cbar = fig.colorbar(im, ax=ax, ticks=ticks, fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel(COLORBAR_LABEL)

    # Save
    base = os.path.splitext(os.path.basename(fname))[0]
    fig.savefig(os.path.join(out_dir, base + ".png"), dpi=300)
    fig.savefig(os.path.join(out_dir, base + ".svg"), format="svg")
    fig.savefig(os.path.join(out_dir, base + ".pdf"), format="pdf")
    plt.close(fig)
    print(f"Saved: {os.path.join(out_dir, base + '.png')}")
