#!/usr/bin/env python3
"""
Plot saved homogenized fields from JSON output and combine multiple cases into one figure.

What this script does
---------------------
This script reads previously saved homogenized results and renders them as
publication-style 2D patch plots using Matplotlib. It is intended to work with
data produced by the homogenization script that writes files such as:

- `<stem>_meta.json`
- `<stem>_polygons.json`

For each configured scene, the script:
- loads one or more saved homogenized cases
- reconstructs the coarse-cell polygon geometry
- colors each polygon by its stored homogenized value
- optionally rescales values before plotting
- arranges multiple cases in one combined figure using user-defined translations
- adds case labels such as "a)", "b)", etc.
- optionally adds a colorbar and circular domain outline
- either shows the figure or saves it to an image file

Typical use case
----------------
Use this script when you already have homogenized field data saved to disk and
want to create comparison figures across multiple CFD/OpenFOAM cases in a single
layout for reports, papers, or presentations.

Inputs
------
The script expects saved homogenized output files generated earlier, using a
common file stem such as:

    stored/coarse/alpha.gas_Base

From this stem, the script reads:

- `stored/coarse/alpha.gas_Base_meta.json`
- `stored/coarse/alpha.gas_Base_polygons.json`

The main user input is the `CONFIGS` list, which defines one or more scenes to plot.

Main configuration objects
--------------------------
- `SavedCaseEntry`
    Defines one saved dataset and how it should appear in the final figure.

    Fields:
    - `stem`: file stem of the saved homogenized result
    - `label`: text label shown on the plot
    - `translate`: manual shift in plot coordinates [m]

- `SavedSceneConfig`
    Defines the full rendering settings for one output figure.

    Includes:
    - list of cases
    - colormap and color limits
    - number of color levels
    - edge visibility and line style
    - plot rotation and coordinate display
    - label placement
    - axis visibility
    - figure size and output path
    - colorbar formatting

Outputs
-------
For each scene in `CONFIGS`, the script either:

- saves a figure to the path given by `cfg.screenshot`, or
- opens an interactive Matplotlib window if `cfg.screenshot` is `None`

Typical output examples:
- `out/coarse/homogenized_alpha_gas_collection.png`
- `out/coarse/homogenized_Ie_collection.png`

How to use
----------
1. Make sure the homogenized result files already exist on disk.
2. Set the base folder and layout variables, for example:
       Folder = "coarse"
3. Define one or more `SavedSceneConfig` objects in `CONFIGS`.
4. For each scene, list the saved cases to include using `SavedCaseEntry`.
5. Run the script:
       python your_script_name.py

Basic example
-------------
The example below plots two previously saved `alpha.gas` cases into one figure
and writes the result to disk:

    CONFIGS = [
        SavedSceneConfig(
            title="",
            cases=[
                SavedCaseEntry("stored/coarse/alpha.gas_Base", "a)", translate=(-0.09, 0.0)),
                SavedCaseEntry("stored/coarse/alpha.gas_noButtons", "b)", translate=(0.02, 0.0)),
            ],
            cmap="coolwarm",
            clim=(0.0, 0.8),
            n_color_levels=8,
            screenshot="out/coarse/homogenized_alpha_gas_collection.png",
            colorbar=True,
            colorbar_label=r"$\\langle \\mathrm{s_G}\\rangle^{\\mathrm{s}}\\;[-]$",
        )
    ]

Then run:

    python your_script_name.py

Dependencies
------------
This script uses:
- numpy
- matplotlib

It also relies on JSON files written by the homogenization workflow.

Important notes
---------------
- This script does not compute homogenization itself. It only plots results that
  were saved earlier.
- Geometry is reconstructed from `<stem>_polygons.json`.
- Plot values are taken from the saved `"value"` field in the polygon data.
- `translate` is applied in plot coordinates and is useful for building multi-panel
  figure layouts manually.
- If `rotate_x_up=True`, the displayed axes are rotated so that:
      plot X = -y
      plot Y =  x
- Coordinates remain stored internally in meters, even when displayed as mm.
- `value_scale` and `value_offset` affect only plotting, not the saved source data.
- This script uses hard-coded scene definitions in `CONFIGS`; it does not provide
  a command-line interface.

Assumptions
-----------
- Each case stem has matching `_meta.json` and `_polygons.json` files.
- The saved polygon data are compatible with the plotting functions here.
- LaTeX rendering is enabled in Matplotlib (`text.usetex=True`), so a working
  LaTeX installation may be required on your system.
"""
from __future__ import annotations

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon, Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize, BoundaryNorm
from matplotlib.ticker import FormatStrFormatter


# ============================================================
# CONFIG
# ============================================================

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times",
    "text.latex.preamble": r"\usepackage{amsmath}",
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 10,
})

@dataclass
class SavedCaseEntry:
    stem: str                      # e.g. "stored/alpha_gas_Base"
    label: str
    translate: Tuple[float, float] = (0.0, 0.0)   # translation in plot coordinates [m]


@dataclass
class SavedSceneConfig:
    title: str
    cases: List[SavedCaseEntry]
    
    # value scaling before plotting
    value_scale: float = 1.0
    value_offset: float = 0.0

    # rendering
    cmap: str = "viridis"
    clim: Optional[Tuple[float, float]] = None
    n_color_levels: int = 10
    show_edges: bool = True
    edgecolor: str = "black"
    linewidth: float = 0.25

    # geometry / view
    coords_in_mm: bool = True
    rotate_x_up: bool = True      # plot coordinates become X=-y, Y=x
    show_case_circle: bool = False
    circle_linestyle: str = "--"
    circle_linewidth: float = 0.8
    circle_color: str = "black"

    # labels
    label_font_size: int = 11
    label_offset: Tuple[float, float] = (0.0, 0.0)   # in plot coordinates [m]
    label_color: str = "black"

    # axes
    axis_off: bool = True
    equal_aspect: bool = True
    xlim: Optional[Tuple[float, float]] = None
    ylim: Optional[Tuple[float, float]] = None

    # output
    window_size: Tuple[int, int] = (700, 700)   # pixels
    screenshot: Optional[str] = None
    dpi: int = 300

    # colorbar
    colorbar: bool = True
    colorbar_vertical: bool = True
    colorbar_label: str = ""
    colorbar_ticks: Optional[List[float]] = None
    colorbar_fmt: str = "%0.2f"
    colorbar_fraction: float = 0.046
    colorbar_pad: float = 0.04
    colorbar_shrink: float = 0.9


dyy = 20e-3
dyI = -110e-3
Folder="coarse" #"coarse","covered"

CONFIGS = [
    SavedSceneConfig(
        title="",
        cases=[
            SavedCaseEntry("stored/"+Folder+"/alpha.gas_Base",      "a)", translate=(dyI +dyy,       0)),
            #SavedCaseEntry("stored/"+Folder+"/alpha.gas_LRIM",      "b)", translate=( dyy,       0)),
            #SavedCaseEntry("stored/"+Folder+"/alpha.gas_AO",        "c)", translate=(dyI +dyy,-110e-3)),
            #SavedCaseEntry("stored/"+Folder+"/alpha.gas_DXDY",      "d)", translate=( dyy,-110e-3)),
            #SavedCaseEntry("stored/"+Folder+"/alpha.gas_R45",       "e)", translate=(dyI + dyy, -220e-3)),
            SavedCaseEntry("stored/"+Folder+"/alpha.gas_noButtons", "b)", translate=( dyy,0)),
        ],
        cmap="coolwarm",
        clim=(0.0, 0.8),
        n_color_levels=8,
        show_edges=True,
        edgecolor="black",
        linewidth=0.4,
        coords_in_mm=True,
        rotate_x_up=True,
        show_case_circle=False,
        label_font_size=11,
        label_offset=(-42e-3, 42e-3),
        label_color="black",
        axis_off=True,
        equal_aspect=True,
        window_size=(2079, 2953),
        screenshot="out/"+Folder+"/homogenized_alpha_gas_collection.png",
        dpi=300,
        colorbar=True,
        colorbar_vertical=True,
        colorbar_label=r"$\langle \mathrm{s_G}\rangle^{\mathrm{s}}\;[-]$",
        colorbar_ticks=None,   # e.g. [0.0, 0.2, 0.4, 0.6, 0.8]
        colorbar_fmt="%0.1f",
        colorbar_fraction=0.2,
        colorbar_pad=0.0,
        colorbar_shrink=0.2,
    ),
    SavedSceneConfig(
        title="",
        cases=[
            SavedCaseEntry("stored/"+Folder+"/IeMag_Base",      "a)", translate=(dyI +dyy,       0)),
            #SavedCaseEntry("stored/"+Folder+"/IeMag_LRIM",      "b)", translate=( dyy,       0)),
            #SavedCaseEntry("stored/"+Folder+"/IeMag_AO",        "c)", translate=(dyI +dyy,-110e-3)),
            #SavedCaseEntry("stored/"+Folder+"/IeMag_DXDY",      "d)", translate=( dyy,-110e-3)),
            #SavedCaseEntry("stored/"+Folder+"/IeMag_R45",       "e)", translate=(dyI + dyy, -220e-3)),
            SavedCaseEntry("stored/"+Folder+"/IeMag_noButtons", "b)", translate=( dyy,0)),
        ],
        cmap="coolwarm",
        clim=(220, 420),
        n_color_levels=10,
        show_edges=True,
        edgecolor="black",
        linewidth=0.4,
        coords_in_mm=True,
        rotate_x_up=True,
        show_case_circle=False,
        label_font_size=11,
        label_offset=(-42e-3, 42e-3),
        label_color="black",
        axis_off=True,
        equal_aspect=True,
        window_size=(2079, 2953),
        screenshot="out/"+Folder+"/homogenized_Ie_collection.png",
        dpi=300,
        colorbar=True,
        colorbar_vertical=True,
        colorbar_label=r"$|\langle\mathbf{I}_{\mathrm{e}}\rangle^{\mathrm{s}}|\;[\mathrm{mA/cm^2}]$",
        colorbar_ticks=None,   # e.g. [0.0, 0.2, 0.4, 0.6, 0.8]
        colorbar_fmt="%0.0f",
        colorbar_fraction=0.2,
        colorbar_pad=0.0,
        colorbar_shrink=0.2,
        
        value_scale=0.1,
    ),
        SavedSceneConfig(
        title="",
        cases=[
            SavedCaseEntry("stored/"+Folder+"/UmagNew_Base",      "a)", translate=(dyI +dyy,       0)),
            #SavedCaseEntry("stored/"+Folder+"/UmagNew_LRIM",      "b)", translate=( dyy,       0)),
            #SavedCaseEntry("stored/"+Folder+"/UmagNew_AO",        "c)", translate=(dyI +dyy,-110e-3)),
            #SavedCaseEntry("stored/"+Folder+"/UmagNew_DXDY",      "d)", translate=( dyy,-110e-3)),
            #SavedCaseEntry("stored/"+Folder+"/UmagNew_R45",       "e)", translate=(dyI + dyy, -220e-3)),
            SavedCaseEntry("stored/"+Folder+"/UmagNew_noButtons", "b)", translate=( dyy,0)),
        ],
        cmap="coolwarm",
        clim=(0.0, 0.20),
        n_color_levels=8,
        show_edges=True,
        edgecolor="black",
        linewidth=0.4,
        coords_in_mm=True,
        rotate_x_up=True,
        show_case_circle=False,
        label_font_size=11,
        label_offset=(-42e-3, 42e-3),
        label_color="black",
        axis_off=True,
        equal_aspect=True,
        window_size=(2079, 2953),
        screenshot="out/"+Folder+"/homogenized_UNew_collection.png",
        dpi=300,
        colorbar=True,
        colorbar_vertical=True,
        colorbar_label=r"$|\langle\mathbf{u}\rangle^{\mathrm{s}}|\;[\mathrm{m/s}]$",
        colorbar_ticks=None,   # e.g. [0.0, 0.2, 0.4, 0.6, 0.8]
        colorbar_fmt="%0.2f",
        colorbar_fraction=0.2,
        colorbar_pad=0.0,
        colorbar_shrink=0.2,
    ),
    SavedSceneConfig(
        title="",
        cases=[
            SavedCaseEntry("stored/"+Folder+"/Ux_Base",      "a)", translate=(dyI +dyy,       0)),
            #SavedCaseEntry("stored/"+Folder+"/Ux_LRIM",      "b)", translate=( dyy,       0)),
            #SavedCaseEntry("stored/"+Folder+"/Ux_AO",        "c)", translate=(dyI +dyy,-110e-3)),
            #SavedCaseEntry("stored/"+Folder+"/Ux_DXDY",      "d)", translate=( dyy,-110e-3)),
            #SavedCaseEntry("stored/"+Folder+"/Ux_R45",       "e)", translate=(dyI + dyy, -220e-3)),
            SavedCaseEntry("stored/"+Folder+"/Ux_noButtons", "b)", translate=( dyy,0)),
        ],
        cmap="coolwarm",
        #clim=(0.0, 0.2),
        n_color_levels=8,
        show_edges=True,
        edgecolor="black",
        linewidth=0.4,
        coords_in_mm=True,
        rotate_x_up=True,
        show_case_circle=False,
        label_font_size=11,
        label_offset=(-42e-3, 42e-3),
        label_color="black",
        axis_off=True,
        equal_aspect=True,
        window_size=(2079, 2953),
        screenshot="out/"+Folder+"/homogenized_Ux_collection.png",
        dpi=300,
        colorbar=True,
        colorbar_vertical=True,
        colorbar_label=r"$\langle\mathrm{u_x}\rangle^{\mathrm{s}}\;[\mathrm{m/s}]$",
        colorbar_ticks=None,   # e.g. [0.0, 0.2, 0.4, 0.6, 0.8]
        colorbar_fmt="%0.3f",
        colorbar_fraction=0.2,
        colorbar_pad=0.0,
        colorbar_shrink=0.2,
    ),
    
    
]


# ============================================================
# IMPLEMENTATION
# ============================================================

def load_saved_case(stem: str | Path):
    stem = Path(stem)

    meta_path = stem.with_name(stem.name + "_meta.json")
    poly_path = stem.with_name(stem.name + "_polygons.json")

    with open(meta_path, "r", encoding="utf-8") as f:
        meta = json.load(f)

    with open(poly_path, "r", encoding="utf-8") as f:
        polys = json.load(f)

    return meta, polys


def transform_xy(x: np.ndarray, y: np.ndarray, rotate_x_up: bool):
    """
    If rotate_x_up=True:
        plot X = -y
        plot Y =  x
    """
    if rotate_x_up:
        return -y, x
    return x, y


def make_case_patches(
    poly_items: list,
    rotate_x_up: bool,
    translate: Tuple[float, float],
    value_scale: float = 1.0,
    value_offset: float = 0.0,
):
    """
    Build matplotlib patches for one saved case.
    Translation is applied in plot coordinates [m].
    Values are scaled before plotting:
        value_plot = value_scale * value_saved + value_offset
    """
    dx, dy = translate

    patches = []
    colors = []

    xmin = +np.inf
    xmax = -np.inf
    ymin = +np.inf
    ymax = -np.inf

    for item in poly_items:
        value_saved = float(item["value"])
        value_plot = value_scale * value_saved + value_offset

        for part in item["parts"]:
            exterior = np.asarray(part["exterior"], dtype=float)
            x = exterior[:, 0]
            y = exterior[:, 1]

            xp, yp = transform_xy(x, y, rotate_x_up)
            xp = xp + dx
            yp = yp + dy

            xy_plot = np.column_stack([xp, yp])
            patches.append(MplPolygon(xy_plot, closed=True))
            colors.append(value_plot)

            xmin = min(xmin, float(np.min(xp)))
            xmax = max(xmax, float(np.max(xp)))
            ymin = min(ymin, float(np.min(yp)))
            ymax = max(ymax, float(np.max(yp)))

    bounds = (xmin, xmax, ymin, ymax)
    return patches, np.asarray(colors, dtype=float), bounds


def make_case_label_position(
    bounds: Tuple[float, float, float, float],
    offset: Tuple[float, float],
):
    xmin, xmax, ymin, ymax = bounds
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    ox, oy = offset
    return cx + ox, cy + oy


def make_case_circle(meta: dict, rotate_x_up: bool, translate: Tuple[float, float]):
    """
    Create a matplotlib Circle artist in plot coordinates [m].
    """
    x0, y0 = meta.get("origin", [0.0, 0.0])
    r = float(meta["radius"])

    x0 = float(x0)
    y0 = float(y0)

    xp, yp = transform_xy(np.asarray([x0]), np.asarray([y0]), rotate_x_up)
    xp = float(xp[0] + translate[0])
    yp = float(yp[0] + translate[1])

    return xp, yp, r


def show_saved_scene(cfg: SavedSceneConfig):
    if not cfg.cases:
        raise ValueError("No cases configured.")

    # ------------------------------------------------------------
    # load all cases first
    # ------------------------------------------------------------
    loaded = []
    all_vals = []

    for entry in cfg.cases:
        meta, poly_items = load_saved_case(entry.stem)
        case_patches, case_colors, case_bounds = make_case_patches(
            poly_items=poly_items,
            rotate_x_up=cfg.rotate_x_up,
            translate=entry.translate,
            value_scale=cfg.value_scale,
            value_offset=cfg.value_offset,
        )
        loaded.append((entry, meta, poly_items, case_patches, case_colors, case_bounds))
        all_vals.append(case_colors)

    all_vals = np.concatenate(all_vals)
    clim = cfg.clim if cfg.clim is not None else (float(np.nanmin(all_vals)), float(np.nanmax(all_vals)))

    # ------------------------------------------------------------
    # figure / axes
    # ------------------------------------------------------------
    figsize = (cfg.window_size[0] / cfg.dpi, cfg.window_size[1] / cfg.dpi)
    fig, ax = plt.subplots(figsize=figsize, dpi=cfg.dpi)

    # ------------------------------------------------------------
    # colormap / normalization
    # ------------------------------------------------------------
    cmap = plt.get_cmap(cfg.cmap, cfg.n_color_levels)

    if cfg.n_color_levels and cfg.n_color_levels > 1:
        boundaries = np.linspace(clim[0], clim[1], cfg.n_color_levels + 1)
        norm = BoundaryNorm(boundaries, cmap.N, clip=True)
    else:
        norm = Normalize(vmin=clim[0], vmax=clim[1])

    # ------------------------------------------------------------
    # draw cases
    # ------------------------------------------------------------
    global_xmin = +np.inf
    global_xmax = -np.inf
    global_ymin = +np.inf
    global_ymax = -np.inf

    mappable = None

    for entry, meta, poly_items, case_patches, case_colors, case_bounds in loaded:
        pc = PatchCollection(
            case_patches,
            cmap=cmap,
            norm=norm,
            edgecolor=cfg.edgecolor if cfg.show_edges else "none",
            linewidth=cfg.linewidth if cfg.show_edges else 0.0,
        )
        pc.set_array(case_colors)
        ax.add_collection(pc)
        mappable = pc

        xmin, xmax, ymin, ymax = case_bounds
        global_xmin = min(global_xmin, xmin)
        global_xmax = max(global_xmax, xmax)
        global_ymin = min(global_ymin, ymin)
        global_ymax = max(global_ymax, ymax)

        if cfg.show_case_circle:
            cx, cy, r = make_case_circle(meta, cfg.rotate_x_up, entry.translate)
            ax.add_patch(
                Circle(
                    (cx, cy),
                    r,
                    fill=False,
                    linestyle=cfg.circle_linestyle,
                    linewidth=cfg.circle_linewidth,
                    edgecolor=cfg.circle_color,
                )
            )

        lx, ly = make_case_label_position(case_bounds, cfg.label_offset)
        ax.text(
            lx,
            ly,
            entry.label,
            fontsize=cfg.label_font_size,
            color=cfg.label_color,
            ha="center",
            va="center",
        )

    # ------------------------------------------------------------
    # axes formatting
    # ------------------------------------------------------------
    if cfg.equal_aspect:
        ax.set_aspect("equal", adjustable="box")

    pad = 0.03 * max(global_xmax - global_xmin, global_ymax - global_ymin)
    if cfg.xlim is not None:
        ax.set_xlim(*cfg.xlim)
    else:
        ax.set_xlim(global_xmin - pad, global_xmax + pad)

    if cfg.ylim is not None:
        ax.set_ylim(*cfg.ylim)
    else:
        ax.set_ylim(global_ymin - pad, global_ymax + pad)

    unit = "mm" if cfg.coords_in_mm else "m"
    scale = 1e3 if cfg.coords_in_mm else 1.0

    def _rescale_axis_to_unit():
        # only relabel visually; data remain in meters
        from matplotlib.ticker import FuncFormatter
        fmt = FuncFormatter(lambda val, pos: f"{val * scale:g}")
        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)

    _rescale_axis_to_unit()

    if cfg.rotate_x_up:
        ax.set_xlabel(f"-y [{unit}]")
        ax.set_ylabel(f"x [{unit}]")
    else:
        ax.set_xlabel(f"x [{unit}]")
        ax.set_ylabel(f"y [{unit}]")

    if cfg.axis_off:
        ax.axis("off")

    ax.set_title(cfg.title)

    # ------------------------------------------------------------
    # colorbar
    # ------------------------------------------------------------
    if cfg.colorbar and mappable is not None:
        cbar = fig.colorbar(
            mappable,
            ax=ax,
            orientation="vertical" if cfg.colorbar_vertical else "horizontal",
            fraction=cfg.colorbar_fraction,
            pad=cfg.colorbar_pad,
            shrink=cfg.colorbar_shrink,
        )
        

        if cfg.colorbar_ticks is not None:
            cbar.set_ticks(cfg.colorbar_ticks)
            cbar.set_ticklabels([cfg.colorbar_fmt % t for t in cfg.colorbar_ticks])
        
        if cfg.colorbar_vertical:
            cbar.ax.yaxis.set_major_formatter(FormatStrFormatter(cfg.colorbar_fmt))
            cbar.ax.tick_params(axis="y", length=0)
        else:
            cbar.ax.xaxis.set_major_formatter(FormatStrFormatter(cfg.colorbar_fmt))
            cbar.ax.tick_params(axis="x", length=0)
        # Remove colorbar border / outline
        cbar.outline.set_visible(False)
        cbar.outline.set_linewidth(0.0)

        # Put title above the colorbar as horizontal text
        if cfg.colorbar_label:
            if cfg.colorbar_vertical:
                cbar.ax.set_title(cfg.colorbar_label, pad=8, fontsize=11)
            else:
                cbar.set_label(cfg.colorbar_label)
        cbar.ax.minorticks_off()
    fig.tight_layout()

    # ------------------------------------------------------------
    # output
    # ------------------------------------------------------------
    if cfg.screenshot:
        out = Path(cfg.screenshot)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=cfg.dpi, bbox_inches="tight")
        print(f"Wrote: {out}")
        plt.close(fig)
    else:
        plt.show()


def main() -> int:
    for cfg in CONFIGS:
        show_saved_scene(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
