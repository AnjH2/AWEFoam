#!/usr/bin/env python3
"""
Overlay previously exported profile text files (**extractedData/*.txt**) with
**zero command-line flags** – everything is configured in the *USER CONFIG*
block below.

 * New in this revision (smart labelling):
   Each trace label is parsed from the **base name** so you get friendly
   strings like ``J=800 mA cm⁻², Q=110 mL min⁻¹`` instead of the full file
   stem.

Simply run:
    python plot_extracted_profiles.py
and the script will pop up two windows (anode + cathode) or save PNGs,
according to the options you set at the top of the file.
"""

from pathlib import Path
from typing import Dict, List, Tuple
import re

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    }
)

# -----------------------------------------------------------------------------
# USER CONFIG
# -----------------------------------------------------------------------------
EXTRACTED_DIR = Path("extractedData")  # where *_anode.txt / *_cathode.txt live

# Leave BASES as None to auto-discover everything in EXTRACTED_DIR.
# Or give a list like ["j800_Q110_flat", "j900_Q110_flat"]
BASES = None  # type: List[str] | None

# Which regions to plot – choose any subset of {"anode", "cathode"}
PLOT_REGIONS = ["anode", "cathode"]

SAVE_PNG = True           # True → write PNGs instead of showing interactively
PNG_OUTDIR = Path("figures")

VERBOSE = True             # print what the script is doing
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def extract_number(s: str) -> int:
    """Return first integer in *s* for human-natural sorting."""
    m = re.search(r"(\d+)", s)
    return int(m.group(1)) if m else float("inf")


def discover_bases(dir_: Path) -> List[str]:
    """Find every unique base stem (without _anode/_cathode) in *dir_*."""
    bases = set()
    for txt in dir_.glob("*.txt"):
        stem = txt.stem
        if stem.endswith("_anode"):
            bases.add(stem[:-6])
        elif stem.endswith("_cathode"):
            bases.add(stem[:-8])
    return sorted(bases, key=extract_number)


def pretty_label(base: str) -> str:
    """Return user-friendly label like 'J=800 mA cm⁻², Q=110 mL min⁻¹'."""
    m_j = re.search(r"j(\d+)", base, re.I)
    m_Q = re.search(r"Q(\d+)", base, re.I)
    parts = []
    if m_j:
        parts.append(fr"J={int(m_j.group(1))} mA cm$^{{-2}}$")
    if m_Q:
        parts.append(fr"Q={int(m_Q.group(1))} mL min$^{{-1}}$")
    return ", ".join(parts) if parts else base


def load_region_data(dir_: Path, bases: List[str]) -> Dict[str, List[Tuple[np.ndarray, np.ndarray, str]]]:
    """Read txt files and return dist/avg arrays plus pretty label."""
    data: Dict[str, List[Tuple[np.ndarray, np.ndarray, str]]] = {"anode": [], "cathode": []}
    for base in bases:
        label = pretty_label(base)
        for region in ["anode", "cathode"]:
            path = dir_ / f"{base}_{region}.txt"
            if not path.exists():
                if VERBOSE:
                    print(f"[skip] {path} missing")
                continue
            arr = np.loadtxt(path, comments="#")
            if arr.size == 0:
                if VERBOSE:
                    print(f"[skip] {path} empty")
                continue
            dist = arr[:, 1]
            avg = arr[:, 2]
            data[region].append((dist, avg, label))
    return data


def plot_region(region: str, datasets: List[Tuple[np.ndarray, np.ndarray, str]]):
    if not datasets:
        if VERBOSE:
            print(f"No datasets for {region}")
        return
    fig, ax = plt.subplots(figsize=(4, 3))
    for dist, avg, label in datasets:
        order = np.argsort(dist)
        ax.plot(dist[order], avg[order], "o-", label=label)
    ax.set_xlabel("Distance from diaphragm [mm]")
    ax.set_ylabel("Average Saturation [\%]")
    ax.set_title(f"{region.capitalize()} Saturation In Electrode Thickness")
    ax.legend(fontsize=7)
    plt.tight_layout()
    if SAVE_PNG:
        PNG_OUTDIR.mkdir(exist_ok=True)
        png_path = PNG_OUTDIR / f"aggregate_{region}.png"
        svg_path = PNG_OUTDIR / f"aggregate_{region}.svg"
        fig.savefig(png_path, dpi=300)
        fig.savefig(svg_path, dpi=300)
        if VERBOSE:
            print("Saved", png_path)
    else:
        plt.show()

# -----------------------------------------------------------------------------
# Main driver
# -----------------------------------------------------------------------------
if not EXTRACTED_DIR.exists():
    raise FileNotFoundError(f"Directory '{EXTRACTED_DIR}' not found.")

if BASES is None:
    BASES = discover_bases(EXTRACTED_DIR)
else:
    # Allow full filenames; strip _region suffix and extension
    norm_bases: List[str] = []
    for b in BASES:
        stem = Path(b).stem
        if stem.endswith("_anode"):
            norm_bases.append(stem[:-6])
        elif stem.endswith("_cathode"):
            norm_bases.append(stem[:-8])
        else:
            norm_bases.append(stem)
    BASES = norm_bases

if VERBOSE:
    print("Bases to plot:", BASES)

data_by_region = load_region_data(EXTRACTED_DIR, BASES)

for region in PLOT_REGIONS:
    plot_region(region, data_by_region[region])

