"""
Line extractor with auto-centering:
- Pixel size: 22.188 µm
- For each y-slice, find the central zero-run and align its midpoint to x = 0
- X ticks are user-specified in mm, relative to this per-line center
- Exports both absolute and relative x for downstream use
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ─── USER SETTINGS ─────────────────────────────────────────────────────────────
folder_path      = "."
save_dir         = "extractedData_lineX_centered"

# Horizontal slices (mm) and band thickness (±)
y_values_mm      = [9.5, 22, 35]
band_half_mm     = 12.8/2

# Pixel size (true): 22.188 µm
PIXEL_SIZE_M     = 22.188e-6
MM_PER_M         = 1000.0
MM_PER_PIXEL     = PIXEL_SIZE_M * MM_PER_M      # 0.022188 mm/px
PX_PER_MM        = 1.0 / MM_PER_PIXEL           # ~45.065 px/mm

# Intensity handling (optional)
CLIP_AT_FRAC     = 0.30     # clip at 30% of per-image max; set 0 to disable

# ---- Auto-centering by zero-run detection ------------------------------------
ZERO_EPS_ABS     = None     # absolute threshold for "zero"; if None, use relative
ZERO_EPS_REL     = 0.01     # relative to the line max (e.g., 1% of max)
MIN_RUN_PX       = 5        # require at least this many consecutive "zeros"
PREFER_NEAR_CENTER = True   # among runs, choose the one whose midpoint is closest to line center
# Optional: restrict the chosen zero-run to a window around the geometric center (mm)
CENTER_WINDOW_MM = None     # e.g., 2.0 to require the run center within ±2 mm of the geometric center

# X-axis ticks you want to see (in mm, relative to the detected center)
X_TICKS_MM_REL   = [-4.3,-2.225,-0.425, 0.0,0.425, 2.225, 4.3]

# Plot / style
USE_LATEX        = True
FIG_W_IN, FIG_H_IN = 4.8, 3.2
SHOW_GRID        = True

# ─── FONT SETUP ────────────────────────────────────────────────────────────────
if USE_LATEX:
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })

# ─── HELPERS ───────────────────────────────────────────────────────────────────
def extract_int_for_sort(s: str) -> float:
    m = re.search(r"(\d+)", s)
    return int(m.group(1)) if m else float("inf")

def mm_to_px(mm: float) -> int:
    return int(round(mm * PX_PER_MM))

def row_or_band(img: np.ndarray, y_c_mm: float, half_mm: float):
    """Return a 2-D slice (band_height_px, nx). Clamps to image bounds."""
    ny, _ = img.shape
    if half_mm == 0.0:
        y = np.clip(mm_to_px(y_c_mm), 0, ny - 1)
        return img[y:y+1, :]
    y1 = max(0, mm_to_px(y_c_mm - half_mm))
    y2 = min(ny, mm_to_px(y_c_mm + half_mm))
    if y1 >= y2:
        return None
    return img[y1:y2, :]

def find_zero_center(vals: np.ndarray, x_abs_mm: np.ndarray):
    """
    Detect runs of ~zero values and pick the one near the center.
    Return center_mm (absolute), (i_start, i_end, i_mid) or (None, None).
    """
    v = np.asarray(vals)
    v = np.nan_to_num(v, nan=0.0)

    # Threshold for "zero"
    if ZERO_EPS_ABS is not None:
        thr = float(ZERO_EPS_ABS)
    else:
        vmax = np.max(v) if np.any(np.isfinite(v)) else 0.0
        thr = ZERO_EPS_REL * vmax

    zero_mask = v <= thr

    # Find consecutive True runs
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

    # Optional window around geometric center
    if CENTER_WINDOW_MM is not None:
        geom_center_mm = 0.5 * (x_abs_mm[0] + x_abs_mm[-1])
        runs = [r for r in runs if abs(x_abs_mm[r[2]] - geom_center_mm) <= CENTER_WINDOW_MM]
        if not runs:
            return None, None

# Choose by proximity to center or by longest length
    if PREFER_NEAR_CENTER:
        mid_idx = n // 2
        runs.sort(key=lambda r: abs(r[2] - mid_idx))
    else:
        runs.sort(key=lambda r: r[1] - r[0], reverse=True)

    i0, i1, im = runs[0]
    return float(x_abs_mm[im]), (i0, i1, im)

# ─── MAIN ──────────────────────────────────────────────────────────────────────
os.makedirs(os.path.join(folder_path, save_dir), exist_ok=True)

file_list = sorted(
    [f for f in os.listdir(folder_path) if f.endswith(".npy")],
    key=extract_int_for_sort
)

for fname in file_list:
    img = np.load(os.path.join(folder_path, fname))
    if img.ndim > 2:
        img = np.squeeze(img)
        if img.ndim > 2:
            img = img[..., 0]
    if img.ndim != 2:
        print(f"⚠️ Skipping {fname}: not a 2D array after squeezing.")
        continue

    # Optional clip at fraction of per-image max
    if CLIP_AT_FRAC and CLIP_AT_FRAC > 0:
        vmax = np.nanmax(img)
        if not np.isfinite(vmax) or vmax <= 0:
            print(f"⚠️ Skipping {fname}: invalid max.")
            continue
        img = np.where(img > CLIP_AT_FRAC * vmax, CLIP_AT_FRAC * vmax, img)

    ny, nx   = img.shape
    width_mm = nx * MM_PER_PIXEL
    height_mm= ny * MM_PER_PIXEL

    # Column-center positions (absolute mm)
    x_abs_mm = (np.arange(nx) + 0.5) * MM_PER_PIXEL

    # Plot
    fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
    base_no_ext = os.path.splitext(fname)[0]
    title = re.sub(r"\[(.*?)\]", r"[\\mathrm{\1}]", base_no_ext)
    ax.set_title(rf"{title}", fontsize=14)

    x_min_rel, x_max_rel = np.inf, -np.inf

    for y_mm in y_values_mm:
        band = row_or_band(img, y_mm, band_half_mm)
        if band is None:
            print(f"⚠️ Empty band at y={y_mm} mm in {fname}")
            continue

        vals = band.mean(axis=0)

        # --- find center from zero-run ---
        center_mm, run = find_zero_center(vals, x_abs_mm)
        if center_mm is None:
            # Fallback to geometric center if no valid run detected
            center_mm = 0.5 * (x_abs_mm[0] + x_abs_mm[-1])

        # Relative x for plotting (center at 0)
        x_rel_mm = x_abs_mm - center_mm

        # Track overall rel-x limits for axis
        x_min_rel = min(x_min_rel, float(x_rel_mm[0]))
        x_max_rel = max(x_max_rel, float(x_rel_mm[-1]))

        # Save TXT: x_rel, value, x_abs, center_mm
        out_name = f"{base_no_ext}_y{y_mm:g}.txt"
        out_fp   = os.path.join(folder_path, save_dir, out_name)
        with open(out_fp, "w") as f:
            f.write("# x_rel_mm\tvalue\tx_abs_mm\tcenter_mm\n")
            for xr, xa, v in zip(x_rel_mm, x_abs_mm, vals):
                f.write(f"{xr:.6g}\t{v:.6g}\t{xa:.6g}\t{center_mm:.6g}\n")

        # Plot line
        lbl = f"y = {y_mm} mm"
        if run is not None:
            i0, i1, im = run
            ax.axvspan(x_abs_mm[i0]-center_mm, x_abs_mm[i1-1]-center_mm,
                       color="k", alpha=0.06, lw=0)  # visualize detected zero-run
        ax.plot(x_rel_mm, vals, "-", label=lbl)

    # Axes, ticks, center line at x=0
    ax.set_xlabel(r"$x$ (centered) [mm]")
    ax.set_ylabel("Gas volume fraction")
    ax.axvline(0.0, color="k", linestyle="--", linewidth=0.8, alpha=0.6)

    # Use your relative tick list directly
    if X_TICKS_MM_REL:
        ax.set_xticks(X_TICKS_MM_REL)
        ax.set_xticklabels([f"{t:.2f}" for t in X_TICKS_MM_REL])

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if SHOW_GRID:
        ax.grid(True, ls="--", lw=0.4, alpha=0.4)

    ax.legend(title="Horizontal slices", fontsize=8, loc="best")
    plt.tight_layout()
    plt.show()
