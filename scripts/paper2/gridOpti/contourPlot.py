#!/usr/bin/env python3
"""
Parse an optimization stdout log, reconstruct a parameter grid, and generate a
contour plot of the aggregated objective value.

What this script does
---------------------
This script reads a `stdout.log` file produced by the grid/pairs optimization
workflow, extracts:

- electrode thickness `t_elec`
- channel thickness `t_chan`
- aggregated result `agg`

and converts those values into a 2D contour plot over the `(t_elec, t_chan)`
parameter space.

The script can:
- parse optimization progress lines from `stdout.log`
- combine duplicate parameter points using a selected policy
- snap values to a regular grid to avoid floating-point binning problems
- build a rectangular `Z(t_elec, t_chan)` grid
- generate filled contours and contour lines
- optionally overlay the original sample points
- optionally mark user-defined points
- compute and overlay an optimum curve for each `t_chan` row
- save the parsed points and optimum curve to CSV
- save the final figure as PDF, SVG, and PNG

Typical use case
----------------
Use this script after running `optimize_grid_simple.py`, when you want to inspect
the response surface of the aggregated objective and identify the best
electrode/channel-thickness combinations.

Expected input format
---------------------
The script expects log lines similar to:

    [23/81] te=1.4 tc=2.1 jobs=['12345', '12346'] agg=-2.3118 per_case=[...]

Only lines containing numeric values for:
- `te=...`
- `tc=...`
- `agg=...`

are used.

Important note:
---------------
In the current implementation, the parsed aggregate value is converted with:

    ag = abs(float(...))

so the contour plot uses the absolute value of `agg`, not the signed value.

Inputs
------
The script is controlled entirely through the `CONFIG` dictionary.

Important configuration keys:
- `LOG_PATH`
    Path to the optimization `stdout.log`
- `OUT_DIR`
    Output directory for plots and CSV files
- `OUT_STEM`
    Base name for all generated files
- `DUPLICATE_POLICY`
    How duplicate `(t_elec, t_chan)` points are combined:
    - `"last"`
    - `"mean"`
    - `"min"`
    - `"max"`
- `SNAP_STEP`, `SNAP_DECIMALS`
    Used to stabilize floating-point values such as `1.899999999` -> `1.9`
- `N_LEVELS`
    Number of filled contour levels
- `N_ISO_LEVELS`
    Number of contour-line levels
- `PLOT_OPTIMUM`
    If True, overlay one optimum point per `t_chan` row
- `FIT_OPTIMUM_WITH_POLY`
    If True, fit a polynomial in each row before finding the minimum
- `POLY_ORDER`
    Polynomial order used for row-wise fitting
- `MARKERS`
    Optional labeled points to overlay on the contour plot
- `USE_LATEX`
    Enable LaTeX text rendering in Matplotlib

Outputs
-------
The script writes files to:

    <OUT_DIR>/<OUT_STEM>_*

Main outputs:
- `<OUT_STEM>_points.csv`
    Parsed and de-duplicated `(t_elec, t_chan, agg)` points actually used
- `<OUT_STEM>_optimum.csv`
    Optional optimum curve as:
    - `t_chan`
    - `t_elec_star`
    - `z_min`
- `<OUT_STEM>.pdf`
- `<OUT_STEM>.svg`
- `<OUT_STEM>.png`

How to use
----------
1. Set `LOG_PATH` to the optimization log file.
2. Adjust contour, labeling, optimum, and marker settings in `CONFIG`.
3. Run:

       python contour_from_stdout.py

Basic example
-------------
This example reads an optimization log, combines duplicates by keeping the last
entry, fits a quadratic optimum per row, and writes the figure plus CSV outputs:

    CONFIG = {
        "LOG_PATH": "./opt_results/Grid_1769328570/stdout.log",
        "OUT_DIR": "./opt_results/contours",
        "OUT_STEM": "contour_from_stdout",
        "X_LABEL": r"$t_\\mathrm{el}\\;\\mathrm{[mm]}$",
        "Y_LABEL": r"$t_\\mathrm{ch}\\;\\mathrm{[mm]}$",
        "CBAR_LABEL": r"Cell potential $[\\mathrm{V}]$",
        "DUPLICATE_POLICY": "last",
        "USE_LATEX": False,
        "N_LEVELS": 10,
        "N_ISO_LEVELS": 10,
        "FIT_OPTIMUM_WITH_POLY": True,
        "POLY_ORDER": 2,
        "PLOT_OPTIMUM": True,
        "SAVE_OPTIMUM_CSV": True,
        "SNAP_STEP": 0.1,
        "SNAP_AT_PARSE": True,
    }

Then run:

    python contour_from_stdout.py

How the optimum curve is computed
---------------------------------
The script can determine one optimum `t_elec*` for each fixed `t_chan` row.

Two modes are available:
- discrete minimum:
    choose the minimum among the sampled grid points
- polynomial minimum:
    fit a polynomial `z(t_elec)` for each row, then minimize that fit

If polynomial fitting fails or there are too few valid points, the row is skipped.

Dependencies
------------
This script uses:
- numpy
- matplotlib

Important notes
---------------
- This script does not run simulations. It only parses an existing log file.
- Missing grid points remain as `NaN` in the reconstructed grid.
- Duplicate parameter points are merged before plotting.
- Marker positions can also be snapped to the same grid as the parsed data.
- The figure is saved in three formats: PDF, SVG, and PNG.
- LaTeX rendering requires a working TeX installation when enabled.

Assumptions
-----------
- The log file contains lines with numeric `te`, `tc`, and `agg` fields.
- The parsed points form a meaningful 2D parameter grid.
- Matplotlib contouring is appropriate for the amount and distribution of data.

Limitations
-----------
- The script assumes a rectangular grid when building `Z`.
- Sparse or irregular point sets may produce many `NaN` cells.
- In the current implementation, some contour-level settings are overridden later
  inside `make_contour`, so not every level-related config option behaves
  independently.
"""

import os
import re
import csv
from collections import defaultdict
from typing import List, Tuple, Optional

from matplotlib.ticker import FormatStrFormatter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe


# ==============================
# ===== OVERVIEW / CONFIG ======
# ==============================
CONFIG = {
    # Path to the stdout.log you want to parse
    "LOG_PATH": "./opt_results/Grid_1769328570/stdout.log",

    # Output directory (created if missing)
    "OUT_DIR": "./opt_results/contours",

    # Base filename for outputs (no extension)
    "OUT_STEM": "contour_from_stdout",

    # What to plot on Z: we already read "agg", so nothing to change here
    #"TITLE": r"Aggregate $\mathrm{areaAverage}(U_{\mathrm{Ne}})$ vs $t_\mathrm{elec} \times t_\mathrm{chan}$",
    "X_LABEL": r"$t_\mathrm{el}\;\mathrm{[mm]}$",
    "Y_LABEL": r"$t_\mathrm{ch}\;\mathrm{[mm]}$",
    "CBAR_LABEL": r"Cell potential $[{\mathrm{V}}]$",

    # Number of filled levels in contourf
    "N_LEVELS": 10,

    # Overlay the sample points as small dots?
    "SCATTER_POINTS": False,

    # If duplicate (te, tc) entries exist, choose how to combine: "last", "mean", "min", "max"
    "DUPLICATE_POLICY": "last",

    # Use LaTeX rendering for all text (requires a TeX installation)
    "USE_LATEX": False,

    # Optional manual tick spacing (None = auto)
    "XTICK_STEP": None,
    "YTICK_STEP": None,
    
        "cmap": "plasma",

    # --- contour level controls ---
    # one of: "auto", "linspace", "manual"
    "LEVELS_MODE": "linspace",
    # used when LEVELS_MODE == "linspace"
    "Z_MIN": None,        # None -> use finite min(Z)
    "Z_MAX": None,        # None -> use finite max(Z)

    # used when LEVELS_MODE == "manual"
    "LEVELS": None,       # e.g. [ -2.42, -2.40, -2.38, -2.36 ]
    # NEW: control iso-line count separately
    "N_ISO_LEVELS": 10,   # e.g., set to 8 to force 8 line levels

    # if True, show iso-lines only at these values (overrides line levels only)
    "ISO_LEVELS": None,   # e.g. [ -2.40, -2.35 ]

    # put markers above everything
    "MARKER_ZORDER": 10,
    
    # Overlay the per-row optimum (t_el that minimizes Z for each t_ch)
    "PLOT_OPTIMUM": True,
    # Save the optimal curve to CSV
    "SAVE_OPTIMUM_CSV": True,
        # Find optimum via polynomial fit (per t_ch row)
    "FIT_OPTIMUM_WITH_POLY": True,   # set False to fall back to discrete argmin
    "POLY_ORDER": 2,                 # 2 for quadratic, 3 for cubic, etc.
    "POLY_REQUIRE_POINTS": None,     # min finite points needed; None -> POLY_ORDER + 1
    "POLY_CLIP_TO_GRID": True,       # keep optimum inside [min(te), max(te)]
    "PLOT_OPTIMUM": True,            # overlay the optimal t_el(t_ch) curve
    "SAVE_OPTIMUM_CSV": True,        # save the curve to CSV
    # Optional: draw a smooth preview of each fitted row (usually off to avoid clutter)
    "PLOT_FIT_PREVIEW": False,       # if True, draws faint polylines for each row fit
    "FIT_PREVIEW_SAMPLES": 200,      # samples for preview in TE range
    # Snap te/tc to avoid 1.8999999999997 vs 1.9 bins
    "SNAP_STEP": 0.1,      # use this if your grid step is known (e.g. 0.1)
    "SNAP_DECIMALS": None, # else set e.g. 3 to round(x, 3). If SNAP_STEP is set, this is ignored.
    "SNAP_AT_PARSE": True, # snap immediately while parsing (recommended)
    
        "MARKERS": [
        {"te": 0.5, "tc": 0.5, "label": "ca","text_dx": 0.04, "text_dy": 0.04},
        {"te": 2, "tc": 0.5, "label": "cb","text_dx": 0.04, "text_dy": 0.04},
        {"te": 2.5, "tc": 0.5, "label": "cc","text_dx": -0.1, "text_dy": 0.04},
        {"te": 2, "tc": 1.5, "label": "bb","text_dx": 0.04, "text_dy": 0.04},
        {"te": 2, "tc": 2.5, "label": "ab","text_dx": 0.04, "text_dy": -0.04},
        
    ],
    # Default marker look (can be overridden per item)
    "MARKER_DEFAULTS": {
        "color": "deepskyblue",          # blue
        "marker": "*",         # cross
        "size": 100,            # scatter size
        "linewidth": 1.0,      # cross thickness
        "text_dx": 0.04,       # label offset in data units (x-direction)
        "text_dy": 0.02,       # label offset in data units (y-direction)
    },
    # Snap markers to your grid like data? (helps 1.899999999 -> 1.9)
    "SNAP_MARKERS": True,

"FIG_SIZE": (3.3*2, 2.645*2),
}
# ==============================


# ---------- Helpers ----------
def _snap_value(x: float, step=None, decimals=None) -> float:
    if step and step > 0:
        return round(x / step) * step
    if decimals is not None:
        return round(x, int(decimals))
    return float(f"{x:.12g}")
    
def argmin_per_row(Z: np.ndarray, te_vals: List[float], tc_vals: List[float]):
    """
    For each row j (fixed t_ch), find the column i (t_el) that minimizes Z[j, i].
    Returns a list of tuples (t_ch, t_el*, z_min) for rows that have any finite data.
    """
    out = []
    for j, tc in enumerate(tc_vals):
        row = Z[j, :]
        if not np.any(np.isfinite(row)):
            continue
        i = int(np.nanargmin(row))
        out.append((tc, te_vals[i], float(row[i])))
    return out

def _polyfit_argmin_1row(te_vals, z_row, order: int,
                         require_points: Optional[int],
                         clip_to_grid: bool):
    """
    Fit poly of given 'order' to (te, z) for one t_ch row and return (te*, z*),
    where te* minimizes the polynomial over the TE range (optionally clipped).

    Returns None if there aren't enough finite points or the fit fails.
    """
    import numpy as np

    te_vals = np.asarray(te_vals, float)
    z_row   = np.asarray(z_row, float)

    # finite mask
    m = np.isfinite(z_row)
    if not np.any(m):
        return None

    x = te_vals[m]
    y = z_row[m]

    # how many points required to fit this order?
    min_pts = (order + 1) if require_points is None else int(require_points)
    if x.size < min_pts:
        return None

    # To improve conditioning, scale x to zero mean and unit span
    xm = x.mean()
    xs = x.std() if x.std() > 0 else (x.max() - x.min())
    if not np.isfinite(xs) or xs == 0:
        # degenerate row
        i = int(np.nanargmin(z_row))
        return float(te_vals[i]), float(z_row[i])

    x_scaled = (x - xm) / xs

    try:
        # fit polynomial in scaled-x
        coeff_scaled = np.polyfit(x_scaled, y, order)  # high->low
    except Exception:
        return None

    # derivative roots in scaled-x
    p = np.poly1d(coeff_scaled)
    dp = np.polyder(p)
    roots_scaled = dp.r

    # consider only real roots
    roots_scaled = roots_scaled[np.isreal(roots_scaled)].real
    # de-scale to original x
    roots = roots_scaled * xs + xm

    # evaluate candidates: interior critical points + bounds
    if clip_to_grid:
        lo, hi = float(te_vals.min()), float(te_vals.max())
    else:
        # allow slightly outside the grid (but still evaluate sensibly)
        pad = 0.1 * (float(te_vals.max()) - float(te_vals.min()))
        lo, hi = float(te_vals.min()) - pad, float(te_vals.max()) + pad

    cand = [lo, hi]
    for r in roots:
        if clip_to_grid:
            if lo <= r <= hi:
                cand.append(r)
        else:
            cand.append(r)

    # best candidate by evaluating the polynomial (in scaled domain)
    def eval_poly(x0):
        xs0 = (x0 - xm) / xs
        return float(p(xs0))

    z_vals = []
    for c in cand:
        try:
            zc = eval_poly(c)
        except Exception:
            zc = np.nan
        z_vals.append(zc)

    if not np.any(np.isfinite(z_vals)):
        # fallback to discrete
        i = int(np.nanargmin(z_row))
        return float(te_vals[i]), float(z_row[i])

    i_best = int(np.nanargmin(z_vals))
    te_star = float(cand[i_best])
    z_star  = float(z_vals[i_best])

    # hard-clip if requested (e.g., a root on the boundary might nudge outside)
    if clip_to_grid:
        te_star = min(max(te_star, lo), hi)

    return te_star, z_star


def argmin_poly_per_row(Z: np.ndarray, te_vals: list, tc_vals: list,
                        order: int, require_points: Optional[int],
                        clip_to_grid: bool):
    """
    For each row j (fixed t_ch), fit a polynomial z(te) and return the argmin:
    list of (t_ch, t_el*, z_min). Skips rows that can't be fitted.
    """
    out = []
    for j, tc in enumerate(tc_vals):
        res = _polyfit_argmin_1row(
            te_vals, Z[j, :], order=order,
            require_points=require_points,
            clip_to_grid=clip_to_grid
        )
        if res is None:
            continue
        te_star, z_star = res
        out.append((float(tc), float(te_star), float(z_star)))
    return out
def save_optimum_csv(opt_points, out_csv: str):
    """
    Save (t_ch, t_el_star, z_min) as a CSV.
    """
    with open(out_csv, "w", newline="") as cf:
        w = csv.writer(cf)
        w.writerow(["t_chan", "t_elec_star", "z_min"])
        for tc, te_star, zmin in opt_points:
            w.writerow([tc, te_star, zmin])

def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path

def setup_matplotlib_latex(use_latex: bool):
    if use_latex:
        mpl.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            # Add any macros/packages you like:
            "text.latex.preamble": r"\usepackage{siunitx}\usepackage{amsmath}",
        })
    else:
        mpl.rcParams.update({
            "text.usetex": False,
            "font.family": "serif",
        })

def parse_stdout_log(path: str) -> List[Tuple[float, float, float]]:
    """
    Returns a list of (te, tc, agg).
    Accepts lines containing 'te=<float>' 'tc=<float>' and 'agg=<float>'.
    Ignores lines missing any of those.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Log not found: {path}")

    # Regexes tolerate scientific notation and signs
    rex_te = re.compile(r"\bte\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")
    rex_tc = re.compile(r"\btc\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")
    rex_ag = re.compile(r"\bagg\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")

    rows: List[Tuple[float, float, float]] = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            m_te = rex_te.search(s)
            m_tc = rex_tc.search(s)
            m_ag = rex_ag.search(s)
            if not (m_te and m_tc and m_ag):
                continue
            try:
                te = float(m_te.group(1))
                tc = float(m_tc.group(1))
                ag = abs(float(m_ag.group(1)))
                if CONFIG.get("SNAP_AT_PARSE", True):
                    te = _snap_value(te, CONFIG.get("SNAP_STEP"), CONFIG.get("SNAP_DECIMALS"))
                    tc = _snap_value(tc, CONFIG.get("SNAP_STEP"), CONFIG.get("SNAP_DECIMALS"))
                if np.isfinite(te) and np.isfinite(tc) and np.isfinite(ag):
                    rows.append((te, tc, ag))
            except Exception:
                continue
    return rows

def combine_duplicates(rows, policy: str):
    if not rows:
        return rows

    step  = CONFIG.get("SNAP_STEP")
    decs  = CONFIG.get("SNAP_DECIMALS")
    buckets = defaultdict(list)

    for te, tc, ag in rows:
        te_s = _snap_value(te, step, decs)
        tc_s = _snap_value(tc, step, decs)
        buckets[(te_s, tc_s)].append(ag)

    out = []
    for (te_s, tc_s), vals in buckets.items():
        if policy == "last":
            z = vals[-1]
        elif policy == "mean":
            z = float(np.mean(vals))
        elif policy == "min":
            z = float(np.min(vals))
        elif policy == "max":
            z = float(np.max(vals))
        else:
            z = vals[-1]
        out.append((te_s, tc_s, z))

    out.sort(key=lambda t: (t[1], t[0]))  # (tc, te)
    return out

def build_grid(rows: List[Tuple[float, float, float]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[float], List[float]]:
    """
    Build a rectangular grid Z over unique TE (x) and TC (y).
    Returns T_e, T_c, Z, te_vals(sorted), tc_vals(sorted)
    """
    if not rows:
        raise ValueError("No (te, tc, agg) rows to plot.")

    te_vals = sorted({te for te, _, _ in rows})
    tc_vals = sorted({tc for _, tc, _ in rows})

    te_index = {v: i for i, v in enumerate(te_vals)}
    tc_index = {v: j for j, v in enumerate(tc_vals)}

    Z = np.full((len(tc_vals), len(te_vals)), np.nan, dtype=float)
    for te, tc, z in rows:
        i = te_index[te]
        j = tc_index[tc]
        Z[j, i] = z  # rows = tc, cols = te

    T_e, T_c = np.meshgrid(te_vals, tc_vals, indexing="xy")
    return T_e, T_c, Z, te_vals, tc_vals

def save_points_csv(rows: List[Tuple[float, float, float]], out_csv: str):
    with open(out_csv, "w", newline="") as cf:
        w = csv.writer(cf)
        w.writerow(["t_elec", "t_chan", "agg"])
        w.writerows(rows)

def make_contour(T_e, T_c, Z, te_vals, tc_vals, cfg, out_stem: str,
                 opt_points=None, fit_preview_lines=None):
    fig, ax = plt.subplots(figsize=CONFIG.get("FIG_SIZE", (7.2, 5.6)))

    # --- choose contour levels ---
    Zfin = Z[np.isfinite(Z)]
    if Zfin.size == 0:
        raise ValueError("All Z are NaN; cannot contour.")

    levels_mode = cfg.get("LEVELS_MODE", "linspace")
    n_levels    = int(cfg.get("N_LEVELS", 12))

    if levels_mode == "manual" and cfg.get("LEVELS"):
        levels_filled = list(cfg["LEVELS"])
    else:
        zmin = cfg.get("Z_MIN", None)
        zmax = cfg.get("Z_MAX", None)
        if zmin is None: zmin = float(np.nanmin(Zfin))
        if zmax is None: zmax = float(np.nanmax(Zfin))
        # guard if min==max
        if not np.isfinite(zmin) or not np.isfinite(zmax) or zmin == zmax:
            zmin, zmax = float(np.nanmin(Zfin)), float(np.nanmax(Zfin)) + 1e-6
        levels_filled = np.linspace(zmin, zmax, n_levels)

    # iso-line levels (optional override)
    iso_levels = cfg.get("ISO_LEVELS", None)
    if iso_levels is None:
        levels_lines = levels_filled
    else:
        levels_lines = iso_levels

    # ----- choose filled levels (you already have N_LEVELS) -----
    Zfin = Z[np.isfinite(Z)]
    zmin = float(np.nanmin(Zfin))
    zmax = float(np.nanmax(Zfin))
    levels_filled = np.round(np.linspace(zmin, zmax, int(cfg.get("N_LEVELS", 12))),2)
    levels_filled[0]=round(zmin,3)
    levels_filled[-1]=round(zmax,3)
    print(levels_filled)
    # ----- choose line (iso) levels -----
    iso_levels_cfg = cfg.get("ISO_LEVELS", None)
    n_iso = cfg.get("N_ISO_LEVELS", None)

    if iso_levels_cfg is not None:
        levels_lines = iso_levels_cfg
    elif n_iso is not None:
        levels_lines = np.round(np.linspace(zmin, zmax, int(n_iso)),2)
    else:
        # default: same as filled levels
        levels_lines = levels_filled

    # Filled contours (lower zorder so lines sit on top)
    cs = ax.contourf(T_e, T_c, Z, levels=levels_filled, cmap=cfg["cmap"], zorder=1)
    cbar = fig.colorbar(cs, ax=ax)
    cbar.set_label(cfg["CBAR_LABEL"])
    cbar.update_ticks()
    # Line contours above fills
    lines = ax.contour(T_e, T_c, Z, levels=levels_lines, colors="k",
                       linewidths=0.7, zorder=5)
    ax.clabel(lines, inline=True, fontsize=9, fmt="%.3f", zorder=6)

    # Optional sample points
    if cfg["SCATTER_POINTS"]:
        ys, xs = np.where(np.isfinite(Z))
        ax.scatter(np.array(te_vals)[xs], np.array(tc_vals)[ys], s=12, alpha=0.7)

    # >>>>> NEW: draw the optimum curve (t_el* vs t_ch) <<<<<
    if cfg.get("PLOT_OPTIMUM", True) and opt_points:
        # Sort by channel thickness to get a clean polyline
        opt_points_sorted = sorted(opt_points, key=lambda t: t[0])  # (t_ch, t_el*, zmin)
        tc_opt = [p[0] for p in opt_points_sorted]
        te_opt = [p[1] for p in opt_points_sorted]

        # White line on top, small markers in black
        ax.plot(te_opt, tc_opt, linewidth=1, color="w", alpha=1,
                zorder=10, label="poly argmin")
        ax.plot(te_opt, tc_opt, "o", ms=3.5, color="k", alpha=1, zorder=11)

    # Optional: faint preview lines for each row fit
    if cfg.get("PLOT_FIT_PREVIEW", False) and fit_preview_lines:
        for tc, xs, ys in fit_preview_lines:
            ax.plot(xs, [tc]*len(xs), linewidth=0.8, alpha=0.4)
        ax.legend(loc="best", frameon=False)

    # Labels and ticks
    ax.set_xlabel(cfg["X_LABEL"])
    ax.set_ylabel(cfg["Y_LABEL"])

    if cfg["XTICK_STEP"]:
        tmin, tmax = min(te_vals), max(te_vals)
        ax.set_xticks(np.arange(tmin, tmax + 1e-9, cfg["XTICK_STEP"]))
    if cfg["YTICK_STEP"]:
        tmin, tmax = min(tc_vals), max(tc_vals)
        ax.set_yticks(np.arange(tmin, tmax + 1e-9, cfg["YTICK_STEP"]))

    fig.tight_layout()
    markers = cfg.get("MARKERS") or []
    if markers:
        dflt = cfg.get("MARKER_DEFAULTS", {})
        snap_markers = cfg.get("SNAP_MARKERS", False)

        # optional: reuse same snap step/decimals as data
        step = cfg.get("SNAP_STEP")
        decs = cfg.get("SNAP_DECIMALS")

        for m in markers:
            te = float(m["te"])
            tc = float(m["tc"])
            if snap_markers:
                te = _snap_value(te, step, decs)
                tc = _snap_value(tc, step, decs)

            label      = str(m.get("label", ""))
            color      = m.get("color", dflt.get("color", "C0"))
            marker     = m.get("marker", dflt.get("marker", "x"))
            size       = m.get("size", dflt.get("size", 70))
            lw         = m.get("linewidth", dflt.get("linewidth", 1.0))
            text_dx    = m.get("text_dx", dflt.get("text_dx", 0.04))
            text_dy    = m.get("text_dy", dflt.get("text_dy", 0.02))

            zmark = cfg.get("MARKER_ZORDER", 50)  # you can add this key in CONFIG

            # marker on top
            ax.scatter([te], [tc],
                marker=marker,
                s=size,
                c=[color],
                linewidths=lw,
                edgecolors="w",
                zorder=zmark)

            # label on top with white halo
            txt = ax.text(te + text_dx, tc + text_dy, label,
              color=color, fontsize=10, weight="bold",
              va="center", ha="left", zorder=zmark+1)
            txt.set_path_effects([pe.withStroke(linewidth=1, foreground="white")])

    for ext in ("pdf", "svg", "png"):
        of = f"{out_stem}.{ext}"
        fig.savefig(of, bbox_inches="tight", dpi=300 if ext == "png" else None)
        print(f"[PLOT] saved {of}")
    plt.close(fig)

# ---------- Main ----------
def main():
    setup_matplotlib_latex(CONFIG["USE_LATEX"])
    out_dir = ensure_dir(CONFIG["OUT_DIR"])

    rows_raw = parse_stdout_log(CONFIG["LOG_PATH"])
    if not rows_raw:
        raise SystemExit(f"No usable lines found in {CONFIG['LOG_PATH']}")

    rows = combine_duplicates(rows_raw, CONFIG["DUPLICATE_POLICY"])

    # Save the points we actually used (after duplicate policy)
    points_csv = os.path.join(out_dir, f"{CONFIG['OUT_STEM']}_points.csv")
    save_points_csv(rows, points_csv)
    print(f"[CSV] saved {points_csv}")

    # Build grid and plot
    # Build grid
    T_e, T_c, Z, te_vals, tc_vals = build_grid(rows)

    opt_points = None
    fit_preview_lines = None

    if CONFIG.get("FIT_OPTIMUM_WITH_POLY", True):
        # Per-row polynomial argmin
        opt_points = argmin_poly_per_row(
            Z, te_vals, tc_vals,
            order=int(CONFIG["POLY_ORDER"]),
            require_points=CONFIG["POLY_REQUIRE_POINTS"],
            clip_to_grid=bool(CONFIG["POLY_CLIP_TO_GRID"]),
        )

        # (Optional) build faint preview curves for each row
        if CONFIG.get("PLOT_FIT_PREVIEW", False):
            fit_preview_lines = []
            order = int(CONFIG["POLY_ORDER"])
            require_points = CONFIG["POLY_REQUIRE_POINTS"]
            clip_to_grid = bool(CONFIG["POLY_CLIP_TO_GRID"])

            # preview domain in TE range
            xs_plot = np.linspace(min(te_vals), max(te_vals), CONFIG["FIT_PREVIEW_SAMPLES"])
            for j, tc in enumerate(tc_vals):
                # reuse the same scaling/fit as in _polyfit_argmin_1row
                import numpy as np
                row = Z[j, :]
                m = np.isfinite(row)
                if not np.any(m):
                    continue
                x = np.asarray(te_vals, float)[m]
                y = np.asarray(row, float)[m]
                min_pts = (order + 1) if require_points is None else int(require_points)
                if x.size < min_pts:
                    continue
                xm = x.mean()
                xs = x.std() if x.std() > 0 else (x.max() - x.min())
                if not np.isfinite(xs) or xs == 0:
                    continue
                x_scaled = (x - xm) / xs
                try:
                    coeff_scaled = np.polyfit(x_scaled, y, order)
                except Exception:
                    continue
                p = np.poly1d(coeff_scaled)
                ys_plot = p((xs_plot - xm) / xs)
                fit_preview_lines.append((float(tc), xs_plot, ys_plot))
    else:
        # Fallback: discrete row-wise argmin
        opt_points = argmin_per_row(Z, te_vals, tc_vals)

    # Optionally save the optimum curve
    if CONFIG.get("SAVE_OPTIMUM_CSV", False) and opt_points:
        opt_csv = os.path.join(out_dir, f"{CONFIG['OUT_STEM']}_optimum.csv")
        save_optimum_csv(opt_points, opt_csv)
        print(f"[CSV] saved {opt_csv}")

    # Plot
    out_stem = os.path.join(out_dir, CONFIG["OUT_STEM"])
    make_contour(T_e, T_c, Z, te_vals, tc_vals, CONFIG, out_stem,
                 opt_points=(opt_points if CONFIG.get("PLOT_OPTIMUM", True) else None),
                 fit_preview_lines=fit_preview_lines)

    print("\nDone.")

if __name__ == "__main__":
    main()

