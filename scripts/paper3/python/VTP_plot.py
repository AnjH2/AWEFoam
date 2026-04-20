#!/usr/bin/env python3
"""
Render scalar or vector-derived fields from VTP meshes and export comparison figures.

What this script does
---------------------
This script loads one or more VTP surface meshes, extracts a scalar field for
plotting, renders the meshes with PyVista, and then composes the final figure
with Matplotlib so the colorbar and overall layout can be controlled precisely.

It is designed for side-by-side comparison figures where multiple meshes are:
- loaded from `.vtp` files
- colored by a scalar field or a vector-derived quantity
- translated in space to build a custom layout
- optionally rotated before plotting
- labeled individually (for example "a)", "b)", "c)")
- exported as publication-style PNG figures

The script supports:
- plotting a scalar field directly, such as `alpha.gas`
- plotting a vector magnitude, such as `|U|`
- plotting a single vector component, such as `U_x` or `U_y`
- applying simple scaling and offset to the plotted values
- controlling view direction, zoom, roll, and camera pan
- adding a discrete colorbar with custom formatting

Typical use case
----------------
Use this script when you want to create comparison images directly from raw CFD
surface output files, without first homogenizing or post-saving the data in a
separate format.

Inputs
------
The script expects:
- one or more `.vtp` mesh files
- mesh fields already present in point data or cell data
- one or more `SceneConfig` entries in `CONFIGS`

Typical fields used in this script:
- `alpha.gas`
- `U`
- `Ie`

Main configuration objects
--------------------------
- `MeshEntry`
    Defines one mesh to include in the figure.

    Fields:
    - `file`: path to the `.vtp` file
    - `label`: text label shown on the image
    - `translate`: translation applied to the mesh in world coordinates [m]
    - `label_offset`: optional manual label offset in world coordinates [m]
    - `rotate_x_180`: if True, rotate the mesh by 180 degrees around x before plotting

- `SceneConfig`
    Defines one full figure to render.

    Main options include:
    - `field_name`: name of the field to plot
    - `use_magnitude`: use vector magnitude if the field is vector-valued
    - `component`: use one vector component instead of the full vector magnitude
    - `meshes`: list of `MeshEntry` objects
    - colormap and color limits
    - label styling
    - camera/view settings
    - window size and output path
    - colorbar placement and formatting

Field selection rules
---------------------
The plotting field is selected using:
- `field_name`
- `use_magnitude`
- `component`

Rules:
- If the field is scalar, it is plotted directly.
- If the field is vector-valued:
  - set `use_magnitude=True` to plot its magnitude
  - or set `component=0/1/2` to plot one component
- Do not leave both `use_magnitude=False` and `component=None` for a vector field,
  because the script will not know what scalar quantity to plot.

Outputs
-------
For each scene in `CONFIGS`, the script either:

- saves an image to `cfg.screenshot`, or
- opens an interactive figure window if `cfg.screenshot` is `None`

Typical output examples:
- `out/xy_sG_cathode_cellSizes.png`
- `out/xy_U_cathode_H.png`
- `out/xy_Ie_cathode_inOut.png`

How to use
----------
1. Add one or more meshes using `MeshEntry`.
2. Define the field to visualize in a `SceneConfig`.
3. Set layout, camera, and colorbar options.
4. Add the scene to `CONFIGS`.
5. Run the script:

       python your_script_name.py

Basic example
-------------
The example below renders two VTP meshes colored by gas volume fraction and saves
the result to disk:

    CONFIGS = [
        SceneConfig(
            title="",
            field_name="alpha.gas",
            use_magnitude=False,
            component=None,
            meshes=[
                MeshEntry(
                    file="/path/to/case1/plane03.vtp",
                    label="a)",
                    translate=(0.0, 0.16, 0.0),
                    label_offset=(0.035, 0.040, 0.0),
                    rotate_x_180=True,
                ),
                MeshEntry(
                    file="/path/to/case2/plane03.vtp",
                    label="b)",
                    translate=(0.0, 0.0, 0.0),
                    label_offset=(0.070, 0.080, 0.0),
                    rotate_x_180=True,
                ),
            ],
            cmap="coolwarm",
            clim=(0.0, 1.0),
            view="xy",
            camera_zoom=1.4,
            camera_roll_deg=90,
            screenshot="out/xy_sG_comparison.png",
            colorbar_title=r"$\\mathrm{s_G}\\ [\\mathrm{-}]$",
        )
    ]

Then run:

    python your_script_name.py

Dependencies
------------
This script uses:
- numpy
- pyvista
- matplotlib

Install them with pip if needed.

Important notes
---------------
- This script reads raw `.vtp` meshes directly.
- If the requested field exists in point data, point values are plotted.
- If the requested field exists in cell data, cell values are plotted.
- `field_scale` and `field_offset` modify only the plotted values.
- `translate` is applied to the actual mesh position before rendering.
- Labels are placed in screen space based on a world-space reference point.
- The figure is rendered in two steps:
  1. PyVista renders the mesh scene
  2. Matplotlib adds the colorbar and writes the final image
- The colorbar is discrete, based on `n_color_levels`.
- This script uses hard-coded configuration in `CONFIGS`; it does not provide a
  command-line interface.

Assumptions
-----------
- The `.vtp` files exist and can be read by PyVista.
- The requested field name is present in each mesh used in a given scene.
- If the field is vector-valued, the user provides either `use_magnitude=True`
  or a valid `component`.
- Off-screen rendering with PyVista is available on the system.

Typical use cases in this script
--------------------------------
Examples already configured here include:
- gas volume fraction: `alpha.gas`
- velocity magnitude: `|U|`
- current-density magnitude: `|Ie|`
- velocity components such as `U_x` and `U_y`
"""
from __future__ import annotations
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib import cm, colors, ticker


# ============================================================
# CONFIG
# ============================================================

@dataclass
class MeshEntry:
    file: str
    label: str
    translate: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    label_offset: Optional[Tuple[float, float, float]] = None
    rotate_x_180: bool = False


@dataclass
class SceneConfig:
    title: str
    field_name: str
    use_magnitude: bool = False
    component: Optional[int] = None

    meshes: List[MeshEntry] = field(default_factory=list)

    # rendering
    cmap: str = "viridis"
    show_edges: bool = False
    lighting: bool = False
    clim: Optional[Tuple[float, float]] = None

    # labels
    label_font_size: int = 12
    label_offset: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    label_font_family: str = "times"   # "arial", "courier", "times"
    label_color: str = "black"

    # camera
    view: str = "xy"   # "xy", "xz", "yz", or "iso"
    camera_zoom: float = 1.0
    camera_roll_deg: float = 0.0
    camera_pan: Tuple[float, float] = (0.0, 0.0)  # (right, up) in world units along camera plane

    # output
    window_size: Tuple[int, int] = (1400, 900)
    screenshot: Optional[str] = None

    # colorbar
    n_color_levels: int = 10
    colorbar_vertical: bool = True
    colorbar_position_x: float = 0.82
    colorbar_position_y: float = 0.10
    colorbar_width: float = 0.10
    colorbar_height: float = 0.80
    colorbar_n_labels: int = 5
    colorbar_fmt: str = "%0.1f"

    colorbar_label_font_size: int = 11
    colorbar_title_font_size: int = 11

    # legacy / backward-compatible title fields
    colorbar_Title_position_x: float = 510
    colorbar_Title_position_y: float = 700
    colorbar_Title_dx: int = 7
    colorbar_Title_dy: int = -8
    colorbar_Title_dx2: int = 7
    colorbar_Title_dy2: int = -8
    colorbar_Title: str = ""
    colorbar_Title_sub: str = ""
    colorbar_Title2: str = ""
    colorbar_fake_title_main_size: int = 11
    colorbar_fake_title_sub_size: int = 7

    colorbar_ticks: Optional[List[float]] = None
    label_rotation_deg: float = 0.0
    title_rotation_deg: float = 0.0
    colorbar_tick_pad_px: int = 8
    colorbar_text_side: str = "succeed"    # "precede" or "succeed"

    colorbar_title: Optional[str] = None

    field_scale: float = 1.0
    field_offset: float = 0.0

    axes_viewport: Tuple[float, float, float, float] = (0.78, 0.02, 0.98, 0.22)
    show_axes_widget: bool = True


dyy = 30e-3
dyI = -110e-3
ctx = 480
cty = 540

CONFIGS = [
    SceneConfig(
        title="",
        field_name="alpha.gas",
        use_magnitude=False,
        component=None,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_3_12/surfaces/15/plane03.vtp",
                label="a)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_3_12/surfaces/15/plane03.vtp",
                label="b)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(0, 1),
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_sG_cathode_cellSizes.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\mathrm{s_G}\ [\mathrm{-}]$",
        colorbar_fmt="%0.1f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),

    SceneConfig(
        title="",
        field_name="U",
        use_magnitude=True,
        component=None,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_3_12/surfaces/15/plane03.vtp",
                label="c)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_3_12/surfaces/15/plane03.vtp",
                label="d)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(0, 0.4),
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_U_cathode_cellSizes.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\vert\mathbf{u}\vert\ [\mathrm{m/s}]$",
        colorbar_fmt="%0.2f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),
SceneConfig(
        title="",
        field_name="Ie",
        use_magnitude=True,
        component=None,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_3_12/surfaces/15/plane0.vtp",
                label="a)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_3_12/surfaces/15/plane0.vtp",
                label="b)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(150, 450),
        field_scale=0.1,
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_Ie_cathode_cellSizes.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\left\vert \mathbf{I}_{\mathrm{e}} \right\vert\ [\mathrm{mA/cm^2}]$",
        colorbar_fmt="%0.0f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),

    SceneConfig(
        title="",
        field_name="alpha.gas",
        use_magnitude=False,
        component=None,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_3_6/surfaces/15/plane03.vtp",
                label="a)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_3_6/surfaces/11.1/plane03.vtp",
                label="b)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(0, 1),
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_sG_cathode_inOut.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\mathrm{s_G}\ [\mathrm{-}]$",
        colorbar_fmt="%0.1f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),

    SceneConfig(
        title="",
        field_name="U",
        use_magnitude=True,
        component=None,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_3_6/surfaces/15/plane03.vtp",
                label="a)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_3_6/surfaces/11.1/plane03.vtp",
                label="b)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(0, 0.8),
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_U_cathode_inOut.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\vert\mathbf{u}\vert\ [\mathrm{m/s}]$",
        colorbar_fmt="%0.2f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),

    SceneConfig(
        title="",
        field_name="alpha.gas",
        use_magnitude=False,
        component=None,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_1_12/surfaces/15/plane03.vtp",
                label="a)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_1_12/surfaces/15/plane03.vtp",
                label="b)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(0, 1),
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_sG_cathode_H.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\mathrm{s_G}\ [\mathrm{-}]$",
        colorbar_fmt="%0.1f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),
SceneConfig(
        title="",
        field_name="Ie",
        use_magnitude=True,
        component=None,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_1_12/surfaces/15/plane0.vtp",
                label="a)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_1_12/surfaces/15/plane0.vtp",
                label="b)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(150, 450),
        field_scale=0.1,
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_Ie_cathode_H.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\left\vert \mathbf{I}_{\mathrm{e}} \right\vert\ [\mathrm{mA/cm^2}]$",
        colorbar_fmt="%0.0f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),
    SceneConfig(
        title="",
        field_name="U",
        use_magnitude=True,
        component=None,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_1_12/surfaces/15/plane03.vtp",
                label="c)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_1_12/surfaces/15/plane03.vtp",
                label="d)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(0, 1),
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_U_cathode_H.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\vert\mathbf{u}\vert\ [\mathrm{m/s}]$",
        colorbar_fmt="%0.2f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),

    SceneConfig(
        title="",
        field_name="U",
        use_magnitude=False,
        component=0,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_1_12/surfaces/15/plane03.vtp",
                label="a)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_1_12/surfaces/15/plane03.vtp",
                label="b)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(-0.05, 1),
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_Ux_cathode_H.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=12,
        n_color_levels=10,
        colorbar_title=r"$\mathrm{u_x}\ [\mathrm{m/s}]$",
        colorbar_fmt="%0.2f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),

    SceneConfig(
        title="",
        field_name="U",
        use_magnitude=False,
        component=1,
        meshes=[
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/100_1_12/surfaces/15/plane03.vtp",
                label="a)",
                translate=(0, 160e-3, 0),
                label_offset=(35e-3, 40e-3, 0),
                rotate_x_180=True,
            ),
            MeshEntry(
                file="/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/200_1_12/surfaces/15/plane03.vtp",
                label="b)",
                translate=(0, 0, 0),
                label_offset=(70e-3, 80e-3, 0),
                rotate_x_180=True,
            ),
        ],
        cmap="coolwarm",
        show_edges=False,
        lighting=False,
        clim=(-0.5, 0.5),
        label_font_size=7,
        view="xy",
        camera_zoom=1.4,
        camera_roll_deg=90,
        window_size=(585, 400),
        label_font_family="times",
        label_color="black",
        screenshot="out/xy_Uy_cathode_H.png",
        colorbar_vertical=True,
        colorbar_position_x=0.85,
        colorbar_position_y=0.25,
        colorbar_width=0.03,
        colorbar_height=0.5,
        colorbar_n_labels=11,
        n_color_levels=10,
        colorbar_title=r"$\mathrm{u_y}\ [\mathrm{m/s}]$",
        colorbar_fmt="%0.2f",
        label_rotation_deg=0.0,
        title_rotation_deg=0.0,
        colorbar_tick_pad_px=5,
        colorbar_label_font_size=9,
        colorbar_title_font_size=9,
        axes_viewport=(0.85, 0.0, 1.05, 0.2),
        camera_pan=(30e-3, 0),
    ),
]


# ============================================================
# IMPLEMENTATION
# ============================================================

def mpl_font_family(name: str) -> str:
    name = (name or "").lower()
    if name == "times":
        return "serif"
    if name == "courier":
        return "monospace"
    return "sans-serif"


def get_colorbar_title(cfg: SceneConfig) -> Optional[str]:
    if cfg.colorbar_title:
        return cfg.colorbar_title
    if cfg.colorbar_Title:
        return cfg.colorbar_Title
    return None


def pan_camera_in_view_plane(plotter: pv.Plotter, dx: float, dy: float):
    cam = plotter.camera
    pos = np.array(cam.position, dtype=np.float64)
    foc = np.array(cam.focal_point, dtype=np.float64)
    up = np.array(cam.up, dtype=np.float64)

    view_dir = foc - pos
    view_dir /= np.linalg.norm(view_dir)

    right = np.cross(view_dir, up)
    right /= np.linalg.norm(right)

    up2 = np.cross(right, view_dir)
    up2 /= np.linalg.norm(up2)

    shift = dx * right + dy * up2
    cam.position = tuple(pos + shift)
    cam.focal_point = tuple(foc + shift)


def world_to_display(plotter: pv.Plotter, xyz: tuple[float, float, float]) -> tuple[float, float]:
    ren = plotter.renderer
    ren.SetWorldPoint(float(xyz[0]), float(xyz[1]), float(xyz[2]), 1.0)
    ren.WorldToDisplay()
    x, y, _ = ren.GetDisplayPoint()
    return float(x), float(y)


def add_centered_screen_label(
    plotter: pv.Plotter,
    world_point: tuple[float, float, float],
    text: str,
    font_size: int,
    font_family: str = "times",
    color: str = "black",
):
    x, y = world_to_display(plotter, world_point)

    w = 0.6 * font_size * max(1, len(text))
    h = 1.2 * font_size

    pos = (int(x - w / 2), int(y - h / 2))
    plotter.add_text(
        text,
        position=pos,
        font_size=font_size,
        font=font_family,
        color=color,
    )


def get_scalar_array(mesh: pv.DataSet, field_name: str, use_magnitude: bool, component: Optional[int]):
    if field_name in mesh.point_data:
        arr = np.asarray(mesh.point_data[field_name], dtype=np.float64)
        assoc = "point"
    elif field_name in mesh.cell_data:
        arr = np.asarray(mesh.cell_data[field_name], dtype=np.float64)
        assoc = "cell"
    else:
        raise KeyError(
            f"Field '{field_name}' not found.\n"
            f"Point arrays: {list(mesh.point_data.keys())}\n"
            f"Cell arrays : {list(mesh.cell_data.keys())}"
        )

    if arr.ndim == 1:
        vals = arr
    elif arr.ndim == 2 and arr.shape[1] == 3:
        if use_magnitude:
            vals = np.linalg.norm(arr, axis=1)
        elif component is not None:
            vals = arr[:, int(component)]
        else:
            raise ValueError(
                f"Field '{field_name}' is vector-valued. Use use_magnitude=True or component=0/1/2."
            )
    else:
        raise ValueError(f"Unsupported array shape for '{field_name}': {arr.shape}")

    return assoc, vals


def apply_view(plotter: pv.Plotter, view: str):
    if view == "xy":
        plotter.view_xy()
    elif view == "xz":
        plotter.view_xz()
    elif view == "yz":
        plotter.view_yz()
    elif view == "iso":
        plotter.view_isometric()
    else:
        raise ValueError(f"Unknown view '{view}'")


def make_label_position(mesh: pv.PolyData, offset: Tuple[float, float, float]) -> Tuple[float, float, float]:
    b = mesh.bounds
    center = np.array([
        0.5 * (b[0] + b[1]),
        0.5 * (b[2] + b[3]),
        0.5 * (b[4] + b[5]),
    ], dtype=np.float64)
    return tuple(center + np.asarray(offset, dtype=np.float64))


def make_discrete_cmap_and_norm(cfg: SceneConfig, clim: Tuple[float, float]):
    boundaries = np.linspace(clim[0], clim[1], cfg.n_color_levels + 1)
    cmap = plt.get_cmap(cfg.cmap, cfg.n_color_levels)
    norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    return cmap, norm, boundaries


def make_colorbar_ticks(cfg: SceneConfig, clim: Tuple[float, float]):
    if cfg.colorbar_ticks is not None:
        return list(cfg.colorbar_ticks)

    if cfg.colorbar_n_labels is None or cfg.colorbar_n_labels < 2:
        return None

    return np.linspace(clim[0], clim[1], cfg.colorbar_n_labels)


def render_pyvista_scene(cfg: SceneConfig) -> tuple[np.ndarray, Tuple[float, float]]:
    if not cfg.meshes:
        raise ValueError("No meshes configured.")

    loaded = []
    all_vals = []

    for entry in cfg.meshes:
        mesh = pv.read(entry.file)
        if not isinstance(mesh, pv.PolyData):
            mesh = mesh.extract_surface()

        assoc, vals = get_scalar_array(mesh, cfg.field_name, cfg.use_magnitude, cfg.component)

        vals = cfg.field_scale * vals + cfg.field_offset

        mesh = mesh.copy()
        if assoc == "point":
            mesh.point_data["__plot_scalar__"] = vals
        else:
            mesh.cell_data["__plot_scalar__"] = vals

        if entry.rotate_x_180:
            mesh.rotate_x(180.0, inplace=True)

        mesh.translate(entry.translate, inplace=True)

        loaded.append((entry, mesh))
        all_vals.append(vals)

    all_vals = np.concatenate(all_vals)
    clim = cfg.clim if cfg.clim is not None else (
        float(np.nanmin(all_vals)),
        float(np.nanmax(all_vals)),
    )

    p = pv.Plotter(
        window_size=cfg.window_size,
        off_screen=True,
    )

    for entry, mesh in loaded:
        p.add_mesh(
            mesh,
            scalars="__plot_scalar__",
            clim=clim,
            cmap=cfg.cmap,
            n_colors=cfg.n_color_levels,
            interpolate_before_map=True,
            show_edges=cfg.show_edges,
            lighting=cfg.lighting,
            smooth_shading=True,
            show_scalar_bar=False,
        )

    p.add_title(cfg.title, font_size=11)

    if cfg.show_axes_widget:
        p.add_axes(viewport=cfg.axes_viewport)

    apply_view(p, cfg.view)
    p.reset_camera()
    p.camera.Zoom(cfg.camera_zoom)
    p.camera.Roll(cfg.camera_roll_deg)
    pan_camera_in_view_plane(p, cfg.camera_pan[0], cfg.camera_pan[1])

    p.render()

    for entry, mesh in loaded:
        this_label_offset = entry.label_offset if entry.label_offset is not None else cfg.label_offset
        label_pos = make_label_position(mesh, this_label_offset)
        add_centered_screen_label(
            p,
            label_pos,
            entry.label,
            font_size=cfg.label_font_size,
            font_family=cfg.label_font_family,
            color=cfg.label_color,
        )

    p.render()
    img = p.screenshot(return_img=True)
    p.close()

    return img, clim


def compose_with_matplotlib(scene_img: np.ndarray, cfg: SceneConfig, clim: Tuple[float, float]):
    width_px, height_px = cfg.window_size
    dpi = 100.0

    fig = plt.figure(figsize=(width_px / dpi, height_px / dpi), dpi=dpi)
    ax_img = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax_img.imshow(scene_img)
    ax_img.axis("off")

    cax = fig.add_axes([
        cfg.colorbar_position_x,
        cfg.colorbar_position_y,
        cfg.colorbar_width,
        cfg.colorbar_height,
    ])

    cmap, norm, boundaries = make_discrete_cmap_and_norm(cfg, clim)
    ticks = make_colorbar_ticks(cfg, clim)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    orientation = "vertical" if cfg.colorbar_vertical else "horizontal"

    cbar = fig.colorbar(
        sm,
        cax=cax,
        orientation=orientation,
        boundaries=boundaries,
        ticks=ticks,
        spacing="uniform",
        drawedges=False,
    )

    # cleaner colorbar styling
    cbar.outline.set_visible(False)
    cbar.ax.minorticks_off()
    for spine in cbar.ax.spines.values():
        spine.set_visible(False)

    fontfamily = mpl_font_family(cfg.label_font_family)
    fmt = ticker.FormatStrFormatter(cfg.colorbar_fmt)

    if cfg.colorbar_vertical:
        cbar.ax.yaxis.set_major_formatter(fmt)

        if cfg.colorbar_text_side == "precede":
            cbar.ax.yaxis.set_ticks_position("left")
        else:
            cbar.ax.yaxis.set_ticks_position("right")

        cbar.ax.tick_params(
            axis="y",
            which="both",
            labelsize=cfg.colorbar_label_font_size,
            pad=cfg.colorbar_tick_pad_px,
            length=0,
        )

        for ticklab in cbar.ax.get_yticklabels():
            ticklab.set_rotation(cfg.label_rotation_deg)
            ticklab.set_fontfamily(fontfamily)

    else:
        cbar.ax.xaxis.set_major_formatter(fmt)

        if cfg.colorbar_text_side == "precede":
            cbar.ax.xaxis.set_ticks_position("bottom")
        else:
            cbar.ax.xaxis.set_ticks_position("top")

        cbar.ax.tick_params(
            axis="x",
            which="both",
            labelsize=cfg.colorbar_label_font_size,
            pad=cfg.colorbar_tick_pad_px,
            length=0,
        )

        for ticklab in cbar.ax.get_xticklabels():
            ticklab.set_rotation(cfg.label_rotation_deg)
            ticklab.set_fontfamily(fontfamily)

    cb_title = get_colorbar_title(cfg)
    if cb_title:
        title_obj = cbar.ax.set_title(
            cb_title,
            fontsize=cfg.colorbar_title_font_size,
            pad=6,
            fontfamily=fontfamily,
        )
        title_obj.set_rotation(cfg.title_rotation_deg)

    return fig


def show_scene(cfg: SceneConfig):
    scene_img, clim = render_pyvista_scene(cfg)
    fig = compose_with_matplotlib(scene_img, cfg, clim)

    if cfg.screenshot:
        out = Path(cfg.screenshot)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=fig.dpi, bbox_inches=None, pad_inches=0)
        plt.close(fig)
        print(f"Wrote: {out}")
    else:
        plt.show()


def main() -> int:
    for cfg in CONFIGS:
        show_scene(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
