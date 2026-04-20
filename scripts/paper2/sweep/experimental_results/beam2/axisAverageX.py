#!/usr/bin/env python3
"""
Scan intensity profiles **along the X‑direction** while averaging over a user‑selectable Y band
and export figures and text files **relative to the diaphragm interface**.

Revision 12 May 2025 – *distance‑from‑diaphragm edition (bug‑fix 2)*
------------------------------------------------------------------
This hot‑patch removes a duplicate code block that still referenced
``cfg["intervals_mm"]`` when that key was absent, causing a ``KeyError``.
No functional changes beyond the fix.
"""

import os
import re
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Any

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

# --------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------

def extract_number(filename: str) -> int:
    m = re.search(r"(\d+)", filename)
    return int(m.group(1)) if m else float("inf")


def slice_and_average_by_mm(
    array_2d: np.ndarray,
    *,
    x_min_mm: float,
    x_max_mm: float,
    y_min_mm: float,
    y_max_mm: float,
    pixel_size: Tuple[float, float],
) -> float:
    ny, nx = array_2d.shape

    x_min_px = int(np.floor(x_min_mm * pixel_size[0]))
    x_max_px = int(np.ceil( x_max_mm * pixel_size[0]))
    y_min_px = int(np.floor(y_min_mm * pixel_size[1]))
    y_max_px = int(np.ceil( y_max_mm * pixel_size[1]))

    x_min_px = max(0, x_min_px)
    x_max_px = min(nx, x_max_px)
    y_min_px = max(0, y_min_px)
    y_max_px = min(ny, y_max_px)

    if x_max_px <= x_min_px:
        x_max_px = min(nx, x_min_px + 1)
    if y_max_px <= y_min_px:
        y_max_px = min(ny, y_min_px + 1)

    return np.nanmean(array_2d[y_min_px:y_max_px, x_min_px:x_max_px])

# --------------------------------------------------------------
# USER CONFIG – tweak as needed
# --------------------------------------------------------------
folder_path = Path(".")

x_max_mm_full = 4.10
y_max_mm_full = 60.0

y_band_min_mm = 10.0
y_band_max_mm = 60.0

regions: Dict[str, Dict[str, Any]] = {
    "cathode": {
        "x_start_mm": 0.0,
        "x_end_mm": 1.8,
        "n_intervals": 30,
        # "align_edge": "end",
    },
    "anode": {
        "x_start_mm": 2.3,
        "x_end_mm": 4.10,
        "n_intervals": 30,
        # "align_edge": "start",
    },
}

region_colors = {"cathode": "blue", "anode": "red"}

# --------------------------------------------------------------
# CLI
# --------------------------------------------------------------
parser = argparse.ArgumentParser(description="X‑scan extractor – distance from diaphragm")
parser.add_argument("--verbose", action="store_true")
parser.add_argument("--save", action="store_true")
parser.add_argument("--outdir", default="figures")
args = parser.parse_args()

outdir_path = Path(args.outdir)
if args.save:
    outdir_path.mkdir(exist_ok=True)

# --------------------------------------------------------------
# File lists and output dirs
# --------------------------------------------------------------
file_list = sorted([f for f in os.listdir(folder_path) if f.endswith(".npy")], key=extract_number)
save_txt_dir = folder_path / "extractedData"
save_txt_dir.mkdir(exist_ok=True)

# --------------------------------------------------------------
# Determine diaphragm reference edges once, from config
# --------------------------------------------------------------
edge_dict: Dict[str, Dict[str, float]] = {}
for name, cfg in regions.items():
    if "x_start_mm" in cfg:
        start = cfg["x_start_mm"]
    elif "intervals_mm" in cfg:
        start = cfg["intervals_mm"][0][0]
    else:
        raise KeyError(f"Region '{name}' needs 'x_start_mm' or 'intervals_mm'.")

    if "x_end_mm" in cfg:
        end = cfg["x_end_mm"]
    elif "intervals_mm" in cfg:
        end = cfg["intervals_mm"][-1][1]
    else:
        raise KeyError(f"Region '{name}' needs 'x_end_mm' or 'intervals_mm'.")

    edge_dict[name] = {"start": start, "end": end}

align_edges: Dict[str, str] = {n: cfg["align_edge"].lower() for n, cfg in regions.items() if "align_edge" in cfg}

if len(align_edges) < len(regions):
    edges = []
    names = list(regions.keys())
    for i, r1 in enumerate(names):
        for r2 in names[i + 1:]:
            for k1, x1 in edge_dict[r1].items():
                for k2, x2 in edge_dict[r2].items():
                    edges.append((abs(x1 - x2), r1, k1, r2, k2))
    _, r1, k1, r2, k2 = min(edges, key=lambda t: t[0])
    align_edges.setdefault(r1, k1)
    align_edges.setdefault(r2, k2)

if args.verbose:
    print("Chosen diaphragm reference edges:")
    for r, k in align_edges.items():
        print(f"  {r}: {k} (x = {edge_dict[r][k]:.3f} mm)")

x0_dict = {r: edge_dict[r][k] for r, k in align_edges.items()}

# --------------------------------------------------------------
# Main processing loop
# --------------------------------------------------------------
for file_name in file_list:
    img = np.load(folder_path / file_name)
    vmax = np.nanmax(img)
    #img = np.where(img > vmax * 0.3, vmax * 0.3, img)

    ny, nx = img.shape
    pixel_size = (nx / x_max_mm_full, ny / y_max_mm_full)
    if args.verbose:
        print(f"{file_name}: {pixel_size[0]:.2f} px/mm (x), {pixel_size[1]:.2f} px/mm (y)")

    region_data: Dict[str, Dict[str, List[float]]] = {}
    for name, cfg in regions.items():
        region_data[name] = {"x": [], "dist": [], "avg": []}

        if "intervals_mm" in cfg:
            intervals = cfg["intervals_mm"]
        else:
            n = cfg.get("n_intervals", 30)
            edges = np.linspace(cfg["x_start_mm"], cfg["x_end_mm"], n + 1)
            intervals = list(zip(edges[:-1], edges[1:]))

        x0 = x0_dict[name]
        for x_lo_mm, x_hi_mm in intervals:
            x_mid = 0.5 * (x_lo_mm + x_hi_mm)
            dist = abs(x_mid - x0)
            avg_val = slice_and_average_by_mm(img, x_min_mm=x_lo_mm, x_max_mm=x_hi_mm, y_min_mm=y_band_min_mm, y_max_mm=y_band_max_mm, pixel_size=pixel_size)
            region_data[name]["x"].append(x_mid)
            region_data[name]["dist"].append(dist)
            region_data[name]["avg"].append(avg_val)

    # plotting
    fig, ax = plt.subplots(figsize=(4, 3))
    stem = Path(file_name).stem
    j = re.search(r"j(\d+)", stem)
    Q = re.search(r"Q(\d+)", stem)
    typ = re.search(r"(flat|channel|conical)", stem)
    ax.set_title(rf"i={j.group(1) if j else '???'}, Q={Q.group(1) if Q else '???'}, type={typ.group(1) if typ else '???'}")

    for name, data in region_data.items():
        order = np.argsort(data["dist"])
        ax.plot(np.array(data["dist"])[order], np.array(data["avg"])[order], "o-", label=name, color=region_colors.get(name, "black"))

    ax.set_xlabel("Distance from diaphragm [mm]")
    ax.set_ylabel(f"Average ({y_band_min_mm}–{y_band_max_mm} mm)")
    ax.legend()
    ax.set_ybound(0,100)
    plt.tight_layout()

    if args.save:
        out_png = outdir_path / f"{stem}_dist.png"
        fig.savefig(out_png, dpi=300)
        if args.verbose:
            print(f"Saved {out_png}")
    else:
        plt.show()

    base = file_name[:-4]
    for name, data in region_data.items():
        with open(save_txt_dir / f"{base}_{name}.txt", "w") as f:
            f.write("# x_mid_mm\tdist_mm\taverage_value\n")
            for x_mid, dist, avg in zip(data["x"], data["dist"], data["avg"]):
                f.write(f"{x_mid:.6g}\t{dist:.6g}\t{avg:.6g}\n")

