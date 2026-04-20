import os
import re
import numpy as np
import matplotlib.pyplot as plt

# ───────────────────────────────── CONFIG ────────────────────────────────
folder_path = "."                 # Folder with your .npy files
save_txt_dirname = "extractedData"

# Pixel size (defaults to Script B's dataset)
PIXEL_SIZE_M = 22.188e-6          # meters per pixel (22.188 µm/px)
MM_PER_M = 1000.0
MM_PER_PIXEL = PIXEL_SIZE_M * MM_PER_M   # e.g. 0.022188 mm/pixel
PX_PER_MM = 1.0 / MM_PER_PIXEL

# Vertical analysis range
y_max_mm_requested = 41.0         # clamp to the image height per file
i_total = 3                       # number of horizontal slices along y
y_min_mm_cutoff = 3.0             # >>> NEW <<< ignore the first N mm from the bottom (set 0.0 to disable)

# Regions in x (mm)
# When CENTER_BY_ZERO_RUN = False, these are ABSOLUTE x-positions in mm from the image left edge.
# When CENTER_BY_ZERO_RUN = True, these are RELATIVE to the detected center per slice.
regions = [
    ("cathode", -2.225, -0.425),
    ("anode",   0.425, 2.225),
]

region_colors = {
    "cathode": "blue",
    "anode":   "red",
}

# Intensity cutoff control
#   Set one of these to a number and the other to None.
#   - If CLIP_AT_FRAC is not None (e.g. 0.30), values above (frac * per-image max) are clipped.
#   - If CLIP_AT_ABS is not None (e.g. 0.12), values above that absolute level are clipped.
#   - If both are None, no clipping.
CLIP_AT_FRAC = 0.0
CLIP_AT_ABS  = None

# Auto-centering (optional): detect a central ~zero run, align its midpoint to x=0 per y-slice.
CENTER_BY_ZERO_RUN = True
ZERO_EPS_ABS = None            # absolute threshold; if None, use relative below
ZERO_EPS_REL = 0.01            # as fraction of line max (e.g. 1%)
MIN_RUN_PX = 5                 # min consecutive "zero" pixels
PREFER_NEAR_CENTER = True
CENTER_WINDOW_MM = None        # e.g., 2.0 to require run center within ±2 mm of geometric center

# Plot styling
USE_LATEX = True
if USE_LATEX:
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"]
    })

# ─────────────────────────────── HELPERS ─────────────────────────────────
def extract_number(filename):
    m = re.search(r'(\d+)', filename)
    return int(m.group(1)) if m else float('inf')

def mm_to_px(mm):
    return int(round(mm * PX_PER_MM))

def slice_and_average_by_mm(array_2d, x_min_mm, x_max_mm, y_min_mm, y_max_mm,
                            px_per_mm_x, px_per_mm_y):
    """
    Convert the given mm-window to pixel indices using the TRUE pixel sizes,
    clamp to array bounds, and return the NaN-safe mean of the sub-array.
    """
    ny, nx = array_2d.shape
    x_min_px = int(x_min_mm * px_per_mm_x)
    x_max_px = int(x_max_mm * px_per_mm_x)
    y_min_px = int(y_min_mm * px_per_mm_y)
    y_max_px = int(y_max_mm * px_per_mm_y)

    x_min_px = max(0, x_min_px); x_max_px = min(nx, x_max_px)
    y_min_px = max(0, y_min_px); y_max_px = min(ny, y_max_px)

    if x_min_px >= x_max_px or y_min_px >= y_max_px:
        return np.nan

    sub = array_2d[y_min_px:y_max_px, x_min_px:x_max_px]
    return np.nanmean(sub)

def find_zero_center(vals_1d, x_abs_mm):
    """
    Detect runs of ~zero values and pick the one near the center (or longest).
    Return (center_mm, (i0, i1, im)) or (None, None) if not found.
    """
    v = np.asarray(vals_1d)
    v = np.nan_to_num(v, nan=0.0)

    if ZERO_EPS_ABS is not None:
        thr = float(ZERO_EPS_ABS)
    else:
        vmax = np.max(v) if np.any(np.isfinite(v)) else 0.0
        thr = ZERO_EPS_REL * vmax

    zero_mask = v <= thr

    runs = []
    i = 0
    n = len(zero_mask)
    while i < n:
        if zero_mask[i]:
            j = i + 1
            while j < n and zero_mask[j]:
                j += 1
            if j - i >= MIN_RUN_PX:
                mid = (i + j - 1) // 2
                runs.append((i, j, mid))
            i = j
        else:
            i += 1

    if not runs:
        return None, None

    if CENTER_WINDOW_MM is not None:
        geom_center_mm = 0.5 * (x_abs_mm[0] + x_abs_mm[-1])
        runs = [r for r in runs if abs(x_abs_mm[r[2]] - geom_center_mm) <= CENTER_WINDOW_MM]
        if not runs:
            return None, None

    if PREFER_NEAR_CENTER:
        mid_idx = n // 2
        runs.sort(key=lambda r: abs(r[2] - mid_idx))
    else:
        runs.sort(key=lambda r: r[1] - r[0], reverse=True)

    i0, i1, im = runs[0]
    return float(x_abs_mm[im]), (i0, i1, im)

# ────────────────────────────── MAIN LOGIC ───────────────────────────────
file_list = [f for f in os.listdir(folder_path) if f.endswith('.npy')]
file_list = sorted(file_list, key=extract_number)
if not file_list:
    raise SystemExit("No .npy files found in the specified folder.")

save_txt_dir = os.path.join(folder_path, save_txt_dirname)
os.makedirs(save_txt_dir, exist_ok=True)

for file_name in file_list:
    arr = np.load(os.path.join(folder_path, file_name))

    # Squeeze to 2D if needed (like Script B)
    if arr.ndim > 2:
        arr = np.squeeze(arr)
        if arr.ndim > 2:
            arr = arr[..., 0]
    if arr.ndim != 2:
        print(f"Skipping {file_name}: not a 2D array after squeezing.")
        continue

    # Intensity clipping
    if CLIP_AT_ABS is not None and np.isfinite(CLIP_AT_ABS):
        cap = float(CLIP_AT_ABS)
        filtered = np.where(arr > cap, cap, arr)
    else:
        if CLIP_AT_FRAC is not None and CLIP_AT_FRAC > 0:
            vmax_ = np.nanmax(arr)
            if not np.isfinite(vmax_) or vmax_ <= 0:
                print(f"Skipping {file_name}: invalid max value.")
                continue
            filtered = np.where(arr > vmax_ * CLIP_AT_FRAC, vmax_ * CLIP_AT_FRAC, arr)
        else:
            filtered = arr

    rotated = filtered  # keep hook if rotation is ever needed

    ny, nx = rotated.shape
    width_mm  = nx * MM_PER_PIXEL
    height_mm = ny * MM_PER_PIXEL

    # >>> NEW <<< Compute effective vertical window [y_start_mm, y_stop_mm]
    y_start_mm = max(0.0, float(y_min_mm_cutoff))
    y_stop_mm  = min(float(y_max_mm_requested), height_mm)

    if y_start_mm >= y_stop_mm:
        print(f"Warning [{file_name}]: vertical window empty "
              f"(start={y_start_mm:.3f} mm, stop={y_stop_mm:.3f} mm). Skipping file.")
        continue

    y_span_eff_mm = y_stop_mm - y_start_mm
    print(f"[{file_name}] analyzing vertical window: y = {y_start_mm:.2f} .. {y_stop_mm:.2f} mm "
          f"(span {y_span_eff_mm:.2f} mm) split into {i_total} slices")

    # Precompute absolute x (mm) at pixel centers
    x_abs_mm = (np.arange(nx) + 0.5) * MM_PER_PIXEL

    # Warn if any absolute region (if not centering) exceeds width
    if not CENTER_BY_ZERO_RUN:
        for rn, xlo_mm, xhi_mm in regions:
            if xhi_mm > width_mm or xlo_mm < 0:
                print(f"Warning [{file_name}]: region '{rn}' [{xlo_mm}, {xhi_mm}] mm exceeds "
                      f"image width 0..{width_mm:.3g} mm. It will be clamped.")

    # Collectors
    region_data = {rn: {"y": [], "avg": []} for (rn, _, _) in regions}

    # Build i_total slices along y within the effective window
    for i in range(i_total):
        # >>> CHANGED <<< use y_start_mm/y_span_eff_mm instead of 0..y_span_mm
        L1 = y_start_mm + (i / i_total) * y_span_eff_mm
        L2 = y_start_mm + ((i + 1) / i_total) * y_span_eff_mm
        y_mid = 0.5 * (L1 + L2)

        # Per-slice center from band mean profile
        y1_px = max(0, mm_to_px(L1))
        y2_px = min(ny, mm_to_px(L2))
        if y1_px >= y2_px:
            center_mm = None
        else:
            band = rotated[y1_px:y2_px, :]
            line_vals = np.nanmean(band, axis=0)

            center_mm, _run = (None, None)
            if CENTER_BY_ZERO_RUN:
                c_mm, _run = find_zero_center(line_vals, x_abs_mm)
                if c_mm is None:
                    c_mm = 0.5 * (x_abs_mm[0] + x_abs_mm[-1])
                center_mm = c_mm

        for (rn, xlo_mm, xhi_mm) in regions:
            if CENTER_BY_ZERO_RUN and center_mm is not None:
                # relative → absolute
                x_min_mm_abs = center_mm + xlo_mm
                x_max_mm_abs = center_mm + xhi_mm
                x_min_mm_rel = xlo_mm
                x_max_mm_rel = xhi_mm
            else:
                # absolute coordinates
                x_min_mm_abs = xlo_mm
                x_max_mm_abs = xhi_mm
                x_min_mm_rel = None
                x_max_mm_rel = None

            # Print both coordinate systems
            if CENTER_BY_ZERO_RUN and center_mm is not None:
                print(f"[{file_name}] slice {i+1}/{i_total}, region '{rn}': "
                      f"rel = ({x_min_mm_rel:.3f} .. {x_max_mm_rel:.3f}) mm, "
                      f"abs = ({x_min_mm_abs:.3f} .. {x_max_mm_abs:.3f}) mm, "
                      f"center = {center_mm:.3f} mm, "
                      f"y = {L1:.2f}..{L2:.2f} mm")
            else:
                print(f"[{file_name}] slice {i+1}/{i_total}, region '{rn}': "
                      f"abs = ({x_min_mm_abs:.3f} .. {x_max_mm_abs:.3f}) mm, "
                      f"y = {L1:.2f}..{L2:.2f} mm")

            avg_val = slice_and_average_by_mm(
                rotated,
                x_min_mm=x_min_mm_abs, x_max_mm=x_max_mm_abs,
                y_min_mm=L1,           y_max_mm=L2,
                px_per_mm_x=PX_PER_MM, px_per_mm_y=PX_PER_MM
            )
            region_data[rn]["y"].append(y_mid)
            region_data[rn]["avg"].append(avg_val)

    # ── Plot
    fig, ax = plt.subplots(figsize=(4.8, 3.4))

    # Title with bracket->\mathrm substitution
    title_str = file_name[:-4]
    title_str = re.sub(r'\[(.*)\]', r'[\\mathrm{\1}]', title_str)
    ax.set_title(r'{}'.format(title_str), fontsize=14)

    # Scatter each region
    for (rn, _, _) in regions:
        y_vals   = region_data[rn]["y"]
        avg_vals = region_data[rn]["avg"]
        color    = region_colors.get(rn, "black")
        ax.scatter(y_vals, avg_vals, color=color, marker='o', label=rn, alpha=0.8)

        # Write y_mid vs avg to txt
        base_no_ext = file_name[:-4]
        out_txt_name = f"{base_no_ext}_{rn}.txt"
        out_txt_path = os.path.join(save_txt_dir, out_txt_name)
        with open(out_txt_path, 'w') as f:
            f.write("# y_mid_mm\taverage_value\n")
            for (yv, av) in zip(y_vals, avg_vals):
                f.write(f"{yv:.6g}\t{(np.nan if not np.isfinite(av) else av):.6g}\n")

    ax.set_xlabel("Y [mm]")
    ax.set_ylabel("Average Value")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend()
    plt.tight_layout()
    plt.show()

