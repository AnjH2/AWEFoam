#!/usr/bin/env python3
"""
Homogenize surface fields from VTP meshes onto a coarse diamond lattice.

What this script does
---------------------
This script reads one or more VTP surface meshes, extracts scalar or vector-based
fields, and homogenizes them onto a coarse diamond-shaped lattice inside a circular
domain. It supports both direct fields already present in the mesh and derived
mass-flux fields computed from velocity and gas fraction data.

Main features:
- Read VTP meshes with PyVista
- Use existing scalar/vector fields from the mesh
- Compute derived mass flux from `U` and `alpha.gas`
- Homogenize data onto a coarse diamond grid
- Support two averaging methods:
  - `covered`: normalize by valid covered fine-cell area
  - `coarse`: normalize by clipped coarse-cell area
- Save results to CSV/JSON for later plotting or post-processing
- Optionally plot the homogenized result as a heatmap

Inputs
------
The script expects:
- One or more `.vtp` files
- Mesh fields such as:
  - `U` for velocity
  - `alpha.gas` for gas volume fraction
  - other scalar/vector fields such as `Ie`
- User-defined case, field, and homogenization settings in the `__main__` block

Main configuration objects
--------------------------
- `CaseSpec`:
    Defines a case name and path to a mesh file.
- `FieldSpec`:
    Defines which field to homogenize and how to interpret it.
- `HomogenizeConfig`:
    Defines grid spacing, domain radius, origin, masking, and coverage thresholds.

Outputs
-------
For each case/field/method combination, the script writes:

- `<stem>_cells.csv`
    Table of homogenized coarse-cell values and coverage information
- `<stem>_meta.json`
    Metadata describing field settings, grid parameters, and normalization method
- `<stem>_polygons.json`
    Polygon geometry of the coarse cells for exact reconstruction of the heatmap
    (only if `save_polygons=True`)

How to use
----------
1. Define the input cases in `CASES` or another list of `CaseSpec`.
2. Define the fields to process in `FIELDS` as `FieldSpec` objects.
3. Set homogenization parameters in `HCFG`.
4. Call:
       run_for_methods(...)
   or:
       run_homogenization_batch(...)
5. Run the script with Python:
       python your_script_name.py

Basic example
-------------
The example below homogenizes the magnitude of the velocity field `U`
for one case using the `coarse` method:

    CASES = [
        CaseSpec("Base", "/path/to/plane3.vtp"),
    ]

    FIELDS = [
        FieldSpec(out_name="Umag", kind="raw", field_name="U", magnitude=True),
    ]

    HCFG = HomogenizeConfig(
        spacing_x=10e-3,
        spacing_y=10e-3,
        radius=50e-3,
        origin=(0.0, 0.0),
        min_domain_fraction=0.1,
        min_data_coverage=0.1,
        require_center_inside_circle=False,
    )

    run_for_methods(
        CASES,
        FIELDS,
        HCFG,
        methods=["coarse"],
        out_root="stored",
        save_polygons=True,
        verbose=True,
    )

Dependencies
------------
This script uses:
- numpy
- pyvista
- shapely
- matplotlib

Install them with pip if needed.

Important notes
---------------
- Coordinates are assumed to be in SI units [m].
- If a field exists only as point data, it is converted to cell data before use.
- For vector fields, you must choose either:
  - one component (`component=0/1/2`), or
  - full magnitude (`magnitude=True`)
- `post_component_magnitude=True` means:
  homogenize the vector first, then compute the magnitude of the averaged vector.
- This script currently uses hard-coded configuration in the `__main__` block.
  It does not yet provide a command-line interface with arguments.

Typical use case
----------------
Use this script when you want to compare multiple CFD/OpenFOAM cases on a common
coarse spatial grid, while preserving area-weighted averages of scalar fields,
vector components, or derived mass-flux quantities.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import numpy as np
import pyvista as pv

from dataclasses import dataclass, replace
from time import perf_counter

import matplotlib.pyplot as plt

from shapely.geometry import Polygon, Point
from shapely.strtree import STRtree
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

@dataclass(frozen=True)
class CaseSpec:
    name: str
    path: str | Path


@dataclass(frozen=True)
class FieldSpec:
    """
    kind = "raw"       -> use a field already present in the VTP
    kind = "mass_flux" -> compute derived mass flux from U and alpha.gas
    """
    out_name: str
    kind: str = "raw"

    # for raw fields
    field_name: str | None = None
    component: int | None = None
    magnitude: bool = False
    post_component_magnitude: bool = False

    # for derived mass flux
    velocity_name: str = "U"
    alpha_name: str = "alpha.gas"
    rhoL: float | None = None
    rhoG: float | None = None
    mode: str = "gas"   # "gas", "liquid", "mixture"

    def validate(self) -> None:
        if self.magnitude and self.component is not None:
            raise ValueError(f"{self.out_name}: choose either component or magnitude")

        if self.post_component_magnitude and self.component is not None:
            raise ValueError(
                f"{self.out_name}: post_component_magnitude cannot be combined with component"
            )

        if self.post_component_magnitude and self.magnitude:
            raise ValueError(
                f"{self.out_name}: use either magnitude=True OR post_component_magnitude=True, not both"
            )

        if self.kind == "raw":
            if self.field_name is None:
                raise ValueError(f"{self.out_name}: raw field requires field_name")

        elif self.kind == "mass_flux":
            if self.rhoL is None or self.rhoG is None:
                raise ValueError(f"{self.out_name}: mass_flux field requires rhoL and rhoG")

        else:
            raise ValueError(f"{self.out_name}: unknown kind '{self.kind}'")


@dataclass(frozen=True)
class HomogenizeConfig:
    method: str = "coarse"   # "coarse" or "covered"
    spacing_x: float = 10e-3
    spacing_y: float = 10e-3
    radius: float = 50e-3
    origin: tuple[float, float] = (0.0, 0.0)
    include_mask_name: str | None = None
    include_threshold: float = 0.5
    circle_resolution: int = 128
    min_domain_fraction: float = 0.1
    min_data_coverage: float = 0.1
    require_center_inside_circle: bool = False

def log(msg: str, verbose: bool = True) -> None:
    if verbose:
        print(msg, flush=True)


def load_mesh(path: str | Path) -> pv.PolyData:
    mesh = pv.read(path)
    if not isinstance(mesh, pv.PolyData):
        mesh = mesh.extract_surface()
    return mesh


def get_homogenizer(method: str):
    method = method.lower()
    if method == "coarse":
        return homogenize_on_cone_grid_coarse_area
    if method == "covered":
        return homogenize_on_cone_grid
    raise ValueError("method must be 'coarse' or 'covered'")
def prepare_field_on_mesh(
    mesh: pv.PolyData,
    field: FieldSpec,
    *,
    verbose: bool = True,
) -> tuple[pv.PolyData, str, int | None, bool]:
    """
    Return:
        mesh_for_homogenization,
        field_name_to_use,
        component,
        magnitude
    """
    field.validate()

    if field.kind == "raw":
        log(
            f"    field={field.field_name}, component={field.component}, magnitude={field.magnitude}",
            verbose,
        )
        return mesh, field.field_name, field.component, field.magnitude

    if field.kind == "mass_flux":
        log(
            f"    derived mass flux: mode={field.mode}, component={field.component}, magnitude={field.magnitude}",
            verbose,
        )
        mesh_cell, values = compute_mass_flux_scalar(
            mesh,
            velocity_name=field.velocity_name,
            alpha_name=field.alpha_name,
            rhoL=field.rhoL,
            rhoG=field.rhoG,
            mode=field.mode,
            component=field.component,
            magnitude=field.magnitude,
        )

        # attach derived scalar to a copy so the original mesh is untouched
        work = mesh_cell.copy()
        derived_name = field.out_name
        work.cell_data[derived_name] = values

        # after derivation it is now just a scalar field
        return work, derived_name, None, False

    raise ValueError(f"Unsupported field kind: {field.kind}")   
    
def homogenize_case_field(
    mesh: pv.PolyData,
    field: FieldSpec,
    cfg: HomogenizeConfig,
    *,
    verbose: bool = True,
) -> dict:
    t0 = perf_counter()
    field.validate()

    # ------------------------------------------------------------
    # 1) vector-post-magnitude path
    # ------------------------------------------------------------
    if field.post_component_magnitude:
        log(f"    field={field.out_name} using post-component magnitude", verbose)

        if field.kind == "raw":
            if field.field_name is None:
                raise ValueError(f"{field.out_name}: raw field requires field_name")

            result = homogenize_vector_on_cone_grid(
                mesh,
                field_name=field.field_name,
                normalization=cfg.method,
                spacing_x=cfg.spacing_x,
                spacing_y=cfg.spacing_y,
                radius=cfg.radius,
                origin=cfg.origin,
                include_mask_name=cfg.include_mask_name,
                include_threshold=cfg.include_threshold,
                circle_resolution=cfg.circle_resolution,
                min_domain_fraction=cfg.min_domain_fraction,
                min_data_coverage=cfg.min_data_coverage,
                require_center_inside_circle=cfg.require_center_inside_circle,
            )

        elif field.kind == "mass_flux":
            mesh_cell, jvec = compute_mass_flux_vector(
                mesh,
                velocity_name=field.velocity_name,
                alpha_name=field.alpha_name,
                rhoL=field.rhoL,
                rhoG=field.rhoG,
                mode=field.mode,
            )

            work = mesh_cell.copy()
            work.cell_data[field.out_name] = jvec

            result = homogenize_vector_on_cone_grid(
                work,
                field_name=field.out_name,
                normalization=cfg.method,
                spacing_x=cfg.spacing_x,
                spacing_y=cfg.spacing_y,
                radius=cfg.radius,
                origin=cfg.origin,
                include_mask_name=cfg.include_mask_name,
                include_threshold=cfg.include_threshold,
                circle_resolution=cfg.circle_resolution,
                min_domain_fraction=cfg.min_domain_fraction,
                min_data_coverage=cfg.min_data_coverage,
                require_center_inside_circle=cfg.require_center_inside_circle,
            )
        else:
            raise ValueError(f"Unsupported field kind: {field.kind}")

    # ------------------------------------------------------------
    # 2) ordinary scalar path
    # ------------------------------------------------------------
    else:
        work_mesh, field_name, component, magnitude = prepare_field_on_mesh(
            mesh,
            field,
            verbose=verbose,
        )

        homogenizer = get_homogenizer(cfg.method)

        result = homogenizer(
            work_mesh,
            field_name=field_name,
            component=component,
            magnitude=magnitude,
            spacing_x=cfg.spacing_x,
            spacing_y=cfg.spacing_y,
            radius=cfg.radius,
            origin=cfg.origin,
            include_mask_name=cfg.include_mask_name,
            include_threshold=cfg.include_threshold,
            circle_resolution=cfg.circle_resolution,
            min_domain_fraction=cfg.min_domain_fraction,
            min_data_coverage=cfg.min_data_coverage,
            require_center_inside_circle=cfg.require_center_inside_circle,
        )

    # ------------------------------------------------------------
    # metadata
    # ------------------------------------------------------------
    result["field_key"] = field.out_name
    result["field_kind"] = field.kind
    result["averaging_method"] = cfg.method
    result["post_component_magnitude"] = field.post_component_magnitude

    if field.kind == "mass_flux":
        result["mass_flux_mode"] = field.mode
        result["rhoL"] = field.rhoL
        result["rhoG"] = field.rhoG

    dt = perf_counter() - t0
    log(f"    done in {dt:.2f} s | cells kept: {len(result['records'])}", verbose)

    return result
 
 
def run_homogenization_batch(
    cases: list[CaseSpec],
    fields: list[FieldSpec],
    cfg: HomogenizeConfig,
    *,
    out_dir: str | Path = "stored",
    save_polygons: bool = True,
    verbose: bool = True,
) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    total_jobs = len(cases) * len(fields)
    job_no = 0
    batch_t0 = perf_counter()

    log(
        f"Starting batch: {len(cases)} cases × {len(fields)} fields "
        f"using method='{cfg.method}'",
        verbose,
    )

    for case in cases:
        case_t0 = perf_counter()
        log(f"\nLoading case '{case.name}'", verbose)
        log(f"  path: {case.path}", verbose)

        mesh = load_mesh(case.path)

        dt_case_load = perf_counter() - case_t0
        log(
            f"  mesh loaded: {mesh.n_cells:,d} cells, {mesh.n_points:,d} points "
            f"({dt_case_load:.2f} s)",
            verbose,
        )

        for field in fields:
            job_no += 1
            log(f"  [{job_no}/{total_jobs}] case='{case.name}' field='{field.out_name}'", verbose)

            result = homogenize_case_field(mesh, field, cfg, verbose=verbose)

            stem = out_dir / f"{field.out_name}_{case.name}"
            save_homogenized_result(result, stem, save_polygons=save_polygons)

            log(f"    saved -> {stem}_*", verbose)

    batch_dt = perf_counter() - batch_t0
    log(f"\nFinished batch in {batch_dt:.2f} s", verbose)
 
def run_for_methods(
    cases: list[CaseSpec],
    fields: list[FieldSpec],
    base_cfg: HomogenizeConfig,
    methods: list[str],
    *,
    out_root: str | Path = "stored",
    save_polygons: bool = True,
    verbose: bool = True,
) -> None:
    out_root = Path(out_root)

    for method in methods:
        cfg = replace(base_cfg, method=method)
        log(f"\n==============================", verbose)
        log(f"Running method: {method}", verbose)
        log(f"==============================", verbose)

        run_homogenization_batch(
            cases,
            fields,
            cfg,
            out_dir=out_root / method,
            save_polygons=save_polygons,
            verbose=verbose,
        )
     
def _iter_face_point_ids(polydata: pv.PolyData):
    """Yield point-id arrays for each polygon cell in a PolyData."""
    faces = np.asarray(polydata.faces, dtype=np.int64)
    i = 0
    while i < len(faces):
        n = int(faces[i])
        ids = faces[i + 1 : i + 1 + n]
        yield ids
        i += n + 1


def _extract_scalar_cell_field(
    mesh: pv.DataSet,
    field_name: str,
    component: int | None = None,
    magnitude: bool = False,
) -> tuple[pv.PolyData, np.ndarray]:
    """
    Return mesh as PolyData and a scalar cell array for the requested field.

    Rules:
    - scalar field: component must be None
    - vector field: choose either component=0/1/2 or magnitude=True
    """
    work = mesh
    if not isinstance(work, pv.PolyData):
        work = work.extract_surface()

    # Ensure field is available on cells
    if field_name in work.cell_data:
        pass
    elif field_name in work.point_data:
        work = work.point_data_to_cell_data()
    else:
        raise KeyError(
            f"Field '{field_name}' not found.\n"
            f"Point arrays: {list(work.point_data.keys())}\n"
            f"Cell arrays : {list(work.cell_data.keys())}"
        )

    arr = np.asarray(work.cell_data[field_name])

    if arr.ndim == 1:
        if component is not None or magnitude:
            raise ValueError(
                f"Field '{field_name}' is scalar, so component/magnitude is not applicable."
            )
        values = arr.astype(np.float64)

    elif arr.ndim == 2:
        if magnitude and component is not None:
            raise ValueError("Choose either component or magnitude, not both.")

        if magnitude:
            values = np.linalg.norm(arr.astype(np.float64), axis=1)
        else:
            if component is None:
                raise ValueError(
                    f"Field '{field_name}' is vector-valued. "
                    f"Set component=0/1/2 or magnitude=True."
                )
            values = np.asarray(arr[:, component], dtype=np.float64)

    else:
        raise ValueError(f"Unsupported array shape for '{field_name}': {arr.shape}")

    return work, values


def _build_fine_cell_polygons(
    mesh: pv.PolyData,
    values: np.ndarray,
    include_mask_name: str | None = None,
    include_threshold: float = 0.5,
):
    """
    Convert mesh surface cells to shapely polygons in the XY plane.
    Optionally keep only cells where include_mask_name > include_threshold.
    """
    pts_xy = np.asarray(mesh.points[:, :2], dtype=np.float64)

    use_mask = include_mask_name is not None
    if use_mask:
        if include_mask_name in mesh.cell_data:
            mask_vals = np.asarray(mesh.cell_data[include_mask_name], dtype=np.float64)
        elif include_mask_name in mesh.point_data:
            tmp = mesh.point_data_to_cell_data()
            mask_vals = np.asarray(tmp.cell_data[include_mask_name], dtype=np.float64)
        else:
            raise KeyError(
                f"Mask '{include_mask_name}' not found.\n"
                f"Point arrays: {list(mesh.point_data.keys())}\n"
                f"Cell arrays : {list(mesh.cell_data.keys())}"
            )

    fine_polys = []
    fine_vals = []

    for icell, ids in enumerate(_iter_face_point_ids(mesh)):
        if use_mask and not (mask_vals[icell] > include_threshold):
            continue

        poly_xy = pts_xy[ids]
        poly = Polygon(poly_xy)

        if poly.is_empty or (not poly.is_valid) or poly.area <= 0.0:
            continue

        val = values[icell]
        if not np.isfinite(val):
            continue

        fine_polys.append(poly)
        fine_vals.append(float(val))

    if not fine_polys:
        raise ValueError("No valid fine-cell polygons were built.")

    return fine_polys, np.asarray(fine_vals, dtype=np.float64)


def _query_tree_indices(tree: STRtree, geom, geom_id_to_index: dict[int, int]):
    """
    Support both common STRtree return styles:
    - array of indices
    - array/list of geometries
    """
    hits = tree.query(geom)
    if len(hits) == 0:
        return []

    first = hits[0]
    if isinstance(first, (int, np.integer)):
        return list(map(int, hits))

    return [geom_id_to_index[id(g)] for g in hits]


def _diamond_cell(center_x: float, center_y: float, spacing_x: float, spacing_y: float):
    """
    Diamond cell bounded by the four neighboring concave centers.
    """
    return Polygon(
        [
            (center_x - spacing_x, center_y),
            (center_x, center_y + spacing_y),
            (center_x + spacing_x, center_y),
            (center_x, center_y - spacing_y),
        ]
    )

def homogenize_on_cone_grid_coarse_area(
    mesh: pv.DataSet,
    field_name: str,
    *,
    component: int | None = None,
    magnitude: bool = False,
    spacing_x: float = 10e-3,
    spacing_y: float = 10e-3,
    radius: float = 50e-3,
    origin: tuple[float, float] = (0.0, 0.0),
    include_mask_name: str | None = None,
    include_threshold: float = 0.5,
    circle_resolution: int = 128,
    min_domain_fraction: float = 0.0,
    min_data_coverage: float = 0.0,
    require_center_inside_circle: bool = False,
):
    """
    Homogenize onto the coarse diamond lattice, but normalize with the
    clipped coarse-cell area instead of the covered fine-cell area.

    Returned value:
        value = (sum over fine-cell overlaps of field * overlap_area)
                / (area of clipped coarse diamond)

    So uncovered parts of the coarse cell count as zero contribution.

    Extra controls:
    ----------------
    min_domain_fraction :
        Minimum fraction of the FULL diamond area that must lie inside the circle.

    min_data_coverage :
        Minimum fraction of the clipped coarse cell that must be covered by valid
        fine-mesh polygons in order to keep the coarse cell.

    require_center_inside_circle :
        If True, only keep coarse cells whose center lies inside the circle.
    """
    mesh_poly, fine_values = _extract_scalar_cell_field(
        mesh, field_name, component=component, magnitude=magnitude
    )

    fine_polys, fine_vals = _build_fine_cell_polygons(
        mesh_poly,
        fine_values,
        include_mask_name=include_mask_name,
        include_threshold=include_threshold,
    )

    tree = STRtree(fine_polys)
    geom_id_to_index = {id(g): i for i, g in enumerate(fine_polys)}

    x0, y0 = origin
    circle = Point(x0, y0).buffer(radius, quad_segs=circle_resolution)

    nx_max = int(np.ceil(radius / spacing_x)) + 2
    ny_max = int(np.ceil(radius / spacing_y)) + 2

    records = []

    for nx in range(-nx_max, nx_max + 1):
        for ny in range(-ny_max, ny_max + 1):
            # Only convex cone centers
            if (nx + ny) % 2 != 0:
                continue

            cx = x0 + nx * spacing_x
            cy = y0 + ny * spacing_y

            if require_center_inside_circle:
                if (cx - x0) ** 2 + (cy - y0) ** 2 > radius**2 + 1e-15:
                    continue

            diamond = _diamond_cell(cx, cy, spacing_x, spacing_y)
            full_area = diamond.area

            clipped_cell = diamond.intersection(circle)
            clipped_area = clipped_cell.area

            if clipped_cell.is_empty or clipped_area <= 0.0:
                continue

            domain_fraction = clipped_area / full_area
            if domain_fraction < min_domain_fraction:
                continue

            cand_ids = _query_tree_indices(tree, clipped_cell, geom_id_to_index)

            covered_area = 0.0
            value_area_sum = 0.0

            for idx in cand_ids:
                inter = fine_polys[idx].intersection(clipped_cell)
                a = inter.area
                if a <= 0.0:
                    continue

                covered_area += a
                value_area_sum += fine_vals[idx] * a

            data_coverage = covered_area / clipped_area
            if data_coverage < min_data_coverage:
                continue

            # New normalization:
            # divide by clipped coarse-cell area, not covered fine-cell area
            value = value_area_sum / clipped_area

            records.append(
                {
                    "nx": nx,
                    "ny": ny,
                    "x": cx,
                    "y": cy,
                    "area": clipped_area,              # denominator area
                    "covered_area": covered_area,      # actual fine-cell coverage
                    "value": value,
                    "value_area_sum": value_area_sum,  # numerator
                    "polygon": clipped_cell,
                    "domain_fraction": domain_fraction,
                    "data_coverage": data_coverage,
                }
            )

    return {
        "records": records,
        "field_name": field_name,
        "component": component,
        "magnitude": magnitude,
        "spacing_x": spacing_x,
        "spacing_y": spacing_y,
        "radius": radius,
        "origin": origin,
        "include_mask_name": include_mask_name,
        "min_domain_fraction": min_domain_fraction,
        "min_data_coverage": min_data_coverage,
        "require_center_inside_circle": require_center_inside_circle,
        "normalization": "coarse_clipped_area",
    }
    
def homogenize_on_cone_grid(
    mesh: pv.DataSet,
    field_name: str,
    *,
    component: int | None = None,
    magnitude: bool = False,
    spacing_x: float = 10e-3,
    spacing_y: float = 10e-3,
    radius: float = 50e-3,
    origin: tuple[float, float] = (0.0, 0.0),
    include_mask_name: str | None = None,
    include_threshold: float = 0.5,
    circle_resolution: int = 128,
    min_domain_fraction: float = 0.0,
    min_data_coverage: float = 0.0,
    require_center_inside_circle: bool = False,
):
    """
    Homogenize onto the coarse diamond lattice.

    Extra controls:
    ----------------
    min_domain_fraction :
        Minimum fraction of the FULL diamond area that must lie inside the circle.
        Example: 0.25 means keep a cell only if at least 25% of the coarse cell is
        inside the circular domain.

    min_data_coverage :
        Minimum fraction of the CLIPPED coarse cell that must be covered by valid
        fine-mesh polygons.

    require_center_inside_circle :
        If True, only keep coarse cells whose center lies inside the circle.
        If False, keep any coarse cell with enough clipped area.
    """
    mesh_poly, fine_values = _extract_scalar_cell_field(
        mesh, field_name, component=component, magnitude=magnitude
    )

    fine_polys, fine_vals = _build_fine_cell_polygons(
        mesh_poly,
        fine_values,
        include_mask_name=include_mask_name,
        include_threshold=include_threshold,
    )

    tree = STRtree(fine_polys)
    geom_id_to_index = {id(g): i for i, g in enumerate(fine_polys)}

    x0, y0 = origin
    circle = Point(x0, y0).buffer(radius, quad_segs=circle_resolution)

    nx_max = int(np.ceil(radius / spacing_x)) + 2
    ny_max = int(np.ceil(radius / spacing_y)) + 2

    records = []

    for nx in range(-nx_max, nx_max + 1):
        for ny in range(-ny_max, ny_max + 1):
            if (nx + ny) % 2 != 0:
                continue

            cx = x0 + nx * spacing_x
            cy = y0 + ny * spacing_y

            if require_center_inside_circle:
                if (cx - x0) ** 2 + (cy - y0) ** 2 > radius**2 + 1e-15:
                    continue

            diamond = _diamond_cell(cx, cy, spacing_x, spacing_y)
            full_area = diamond.area

            clipped_cell = diamond.intersection(circle)
            clipped_area = clipped_cell.area

            if clipped_cell.is_empty or clipped_area <= 0.0:
                continue

            domain_fraction = clipped_area / full_area
            if domain_fraction < min_domain_fraction:
                continue

            cand_ids = _query_tree_indices(tree, clipped_cell, geom_id_to_index)

            area_sum = 0.0
            value_area_sum = 0.0

            for idx in cand_ids:
                inter = fine_polys[idx].intersection(clipped_cell)
                a = inter.area
                if a <= 0.0:
                    continue

                area_sum += a
                value_area_sum += fine_vals[idx] * a

            if area_sum <= 0.0:
                continue

            data_coverage = area_sum / clipped_area
            if data_coverage < min_data_coverage:
                continue

            records.append(
                {
                    "nx": nx,
                    "ny": ny,
                    "x": cx,
                    "y": cy,
                    "area": area_sum,
                    "value": value_area_sum / area_sum,
                    "polygon": clipped_cell,
                    "domain_fraction": domain_fraction,
                    "data_coverage": data_coverage,
                }
            )

    return {
        "records": records,
        "field_name": field_name,
        "component": component,
        "magnitude": magnitude,
        "spacing_x": spacing_x,
        "spacing_y": spacing_y,
        "radius": radius,
        "origin": origin,
        "include_mask_name": include_mask_name,
        "min_domain_fraction": min_domain_fraction,
        "min_data_coverage": min_data_coverage,
        "require_center_inside_circle": require_center_inside_circle,
        "normalization": "covered_area",
    }


def plot_homogenized_heatmap(
    result: dict,
    *,
    cmap: str = "viridis",
    show_centers: bool = True,
    annotate_indices: bool = False,
    coords_in_mm: bool = True,
    title: str | None = None,
    rotate_x_up: bool = True,
    savepath: str | None = None,
    dpi: int = 300,
    show: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
):
    """
    Plot the homogenized field on the coarse diamond lattice.

    rotate_x_up=True:
        Rotate the plot by +90 degrees so the original +x direction points up.
        Plot coordinates become:
            X_plot = -y
            Y_plot =  x

    vmin, vmax:
        Optional colorbar/data range limits.
    """
    records = result["records"]
    if not records:
        raise ValueError("No homogenized cells to plot.")

    scale = 1e3 if coords_in_mm else 1.0
    unit = "mm" if coords_in_mm else "m"

    def transform_xy(x, y):
        if rotate_x_up:
            return -y * scale, x * scale
        return x * scale, y * scale

    patches = []
    colors = []
    centers_x = []
    centers_y = []

    fig, ax = plt.subplots(figsize=(7, 6))

    for rec in records:
        geom = rec["polygon"]
        geoms = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)

        for g in geoms:
            xy = np.asarray(g.exterior.coords, dtype=np.float64)
            x = xy[:, 0]
            y = xy[:, 1]

            if rotate_x_up:
                xy_plot = np.column_stack([-y * scale, x * scale])
            else:
                xy_plot = np.column_stack([x * scale, y * scale])

            patches.append(MplPolygon(xy_plot, closed=True))
            colors.append(rec["value"])

        cx, cy = transform_xy(rec["x"], rec["y"])
        centers_x.append(cx)
        centers_y.append(cy)

        if annotate_indices:
            ax.text(
                cx,
                cy,
                f"({rec['nx']},{rec['ny']})",
                ha="center",
                va="center",
                fontsize=7,
                color="white",
            )

    pc = PatchCollection(
        patches,
        cmap=cmap,
        edgecolor="black",
        linewidth=0.6,
    )

    colors = np.asarray(colors, dtype=np.float64)
    pc.set_array(colors)

    if vmin is not None or vmax is not None:
        pc.set_clim(vmin=vmin, vmax=vmax)

    ax.add_collection(pc)

    if show_centers:
        ax.scatter(centers_x, centers_y, s=12, c="black", zorder=3)

    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    x0, y0 = result["origin"]
    r = result["radius"]

    xc = x0 + r * np.cos(theta)
    yc = y0 + r * np.sin(theta)

    if rotate_x_up:
        ax.plot(-yc * scale, xc * scale, "k--", lw=1.0)
    else:
        ax.plot(xc * scale, yc * scale, "k--", lw=1.0)

    margin = 1.05 * r * scale
    x0p, y0p = transform_xy(x0, y0)
    ax.set_xlim(x0p - margin, x0p + margin)
    ax.set_ylim(y0p - margin, y0p + margin)
    ax.set_aspect("equal", adjustable="box")

    if rotate_x_up:
        ax.set_xlabel(f"-y [{unit}]")
        ax.set_ylabel(f"x [{unit}]")
    else:
        ax.set_xlabel(f"x [{unit}]")
        ax.set_ylabel(f"y [{unit}]")

    field_label = result["field_name"]
    if result["component"] is not None:
        field_label += f"[{result['component']}]"
    if result["magnitude"]:
        field_label += " magnitude"

    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label(f"Area-averaged {field_label}")

    if title is None:
        title = f"Homogenized field: {field_label}"
        if result["include_mask_name"] is not None:
            title += f" | mask={result['include_mask_name']}"
        if rotate_x_up:
            title += " | rotated 90°"

    ax.set_title(title)
    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)
def save_homogenized_result(
    result: dict,
    stem: str | Path,
    save_polygons: bool = True,
) -> None:
    """
    Save homogenized result to files for later plotting in another script.

    Outputs
    -------
    <stem>_cells.csv
        One row per coarse cell with scalar values and bookkeeping.

    <stem>_meta.json
        Metadata about field name, spacing, radius, normalization, etc.

    <stem>_polygons.json   (optional)
        Polygon coordinates for each coarse cell, so the exact heatmap geometry
        can be reconstructed later.

    Notes
    -----
    - Coordinates are saved in SI units [m].
    - Polygon coordinates are saved as nested lists.
    - Works for both covered-area and coarse-area normalization results.
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)

    records = result.get("records", [])
    if not records:
        raise ValueError("No records found in result.")

    # ------------------------------------------------------------
    # 1) metadata
    # ------------------------------------------------------------
    meta = {}
    for key, val in result.items():
        if key == "records":
            continue
        if isinstance(val, (str, int, float, bool)) or val is None:
            meta[key] = val
        elif isinstance(val, tuple):
            meta[key] = list(val)
        else:
            # keep metadata json-safe without overcomplicating it
            try:
                json.dumps(val)
                meta[key] = val
            except TypeError:
                meta[key] = str(val)

    meta_path = stem.with_name(stem.name + "_meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    # ------------------------------------------------------------
    # 2) cell table
    # ------------------------------------------------------------
    # Collect all scalar keys that may exist across records
    preferred_cols = [
        "nx",
        "ny",
        "x",
        "y",
        "value",
        "area",
        "covered_area",
        "value_area_sum",
        "domain_fraction",
        "data_coverage",
    ]

    extra_cols = []
    for rec in records:
        for k in rec.keys():
            if k not in preferred_cols and k != "polygon" and k not in extra_cols:
                extra_cols.append(k)

    fieldnames = preferred_cols + extra_cols

    cells_path = stem.with_name(stem.name + "_cells.csv")
    with open(cells_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for rec in records:
            row = {}
            for key in fieldnames:
                val = rec.get(key, "")
                if isinstance(val, (np.floating, np.integer)):
                    val = val.item()
                row[key] = val
            writer.writerow(row)

    # ------------------------------------------------------------
    # 3) polygons
    # ------------------------------------------------------------
    if save_polygons:
        poly_items = []

        for rec in records:
            geom = rec["polygon"]

            if geom.geom_type == "Polygon":
                geoms = [geom]
            else:
                geoms = list(geom.geoms)

            parts = []
            for g in geoms:
                exterior = np.asarray(g.exterior.coords, dtype=float).tolist()
                holes = [np.asarray(r.coords, dtype=float).tolist() for r in g.interiors]
                parts.append(
                    {
                        "exterior": exterior,
                        "holes": holes,
                    }
                )

            poly_items.append(
                {
                    "nx": rec["nx"],
                    "ny": rec["ny"],
                    "x": float(rec["x"]),
                    "y": float(rec["y"]),
                    "value": float(rec["value"]),
                    "parts": parts,
                }
            )

        poly_path = stem.with_name(stem.name + "_polygons.json")
        with open(poly_path, "w", encoding="utf-8") as f:
            json.dump(poly_items, f, indent=2)
            
def compute_mass_flux_scalar(
    mesh: pv.DataSet,
    *,
    velocity_name: str = "U",
    alpha_name: str = "alpha.gas",
    rhoL: float,
    rhoG: float,
    mode: str = "gas",
    component: int | None = None,
    magnitude: bool = False,
) -> tuple[pv.PolyData, np.ndarray]:
    """
    Compute a scalar mass-flux-like field from U and alpha.gas.

    Parameters
    ----------
    mesh : pv.DataSet
    velocity_name : str
        Vector velocity field, e.g. "U".
    alpha_name : str
        Scalar gas volume fraction field, e.g. "alpha.gas".
    rhoL, rhoG : float
        Constant liquid and gas densities [kg/m^3].
    mode : str
        "gas"     -> rhoG * alpha * U
        "liquid"  -> rhoL * (1-alpha) * U
        "mixture" -> (rhoG*alpha + rhoL*(1-alpha)) * U
    component : int | None
        For directional flux: 0, 1, or 2.
    magnitude : bool
        If True, use vector magnitude instead of one component.

    Returns
    -------
    mesh_poly : pv.PolyData
        Mesh with cell-compatible data.
    values : np.ndarray
        Scalar cell values [kg/(m^2 s)] if using directional component or magnitude of flux vector.
    """
    if magnitude and component is not None:
        raise ValueError("Choose either component or magnitude, not both.")

    work = mesh
    if not isinstance(work, pv.PolyData):
        work = work.extract_surface()

    # Ensure both arrays are cell data
    needed = [velocity_name, alpha_name]
    need_point_to_cell = False
    for name in needed:
        if name in work.cell_data:
            continue
        elif name in work.point_data:
            need_point_to_cell = True
        else:
            raise KeyError(
                f"Field '{name}' not found.\n"
                f"Point arrays: {list(work.point_data.keys())}\n"
                f"Cell arrays : {list(work.cell_data.keys())}"
            )

    if need_point_to_cell:
        work = work.point_data_to_cell_data()

    U = np.asarray(work.cell_data[velocity_name], dtype=np.float64)
    alpha = np.asarray(work.cell_data[alpha_name], dtype=np.float64).squeeze()

    if U.ndim != 2 or U.shape[1] != 3:
        raise ValueError(f"'{velocity_name}' must be a 3-component vector. Got {U.shape}")
    if alpha.ndim != 1 or alpha.shape[0] != U.shape[0]:
        raise ValueError(
            f"'{alpha_name}' must be scalar per cell and match U cell count. "
            f"Got {alpha.shape} vs {U.shape[0]}"
        )

    if mode == "gas":
        rho_eff = rhoG * alpha
    elif mode == "liquid":
        rho_eff = rhoL * (1.0 - alpha)
    elif mode == "mixture":
        rho_eff = rhoG * alpha + rhoL * (1.0 - alpha)
    else:
        raise ValueError("mode must be 'gas', 'liquid', or 'mixture'")

    j = rho_eff[:, None] * U   # vector mass flux [kg/(m^2 s)]

    if magnitude:
        values = np.linalg.norm(j, axis=1)
    else:
        if component is None:
            raise ValueError("Set component=0/1/2 or magnitude=True")
        values = j[:, component]

    return work, values
def _extract_vector_cell_field(
    mesh: pv.DataSet,
    field_name: str,
) -> tuple[pv.PolyData, np.ndarray]:
    """
    Return mesh as PolyData and a vector cell array of shape (n_cells, 3).
    """
    work = mesh
    if not isinstance(work, pv.PolyData):
        work = work.extract_surface()

    if field_name in work.cell_data:
        pass
    elif field_name in work.point_data:
        work = work.point_data_to_cell_data()
    else:
        raise KeyError(
            f"Field '{field_name}' not found.\n"
            f"Point arrays: {list(work.point_data.keys())}\n"
            f"Cell arrays : {list(work.cell_data.keys())}"
        )

    arr = np.asarray(work.cell_data[field_name], dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(
            f"Field '{field_name}' must be a 3-component vector. Got shape {arr.shape}."
        )

    return work, arr    
def homogenize_vector_on_cone_grid(
    mesh: pv.DataSet,
    field_name: str,
    *,
    normalization: str = "covered",   # "covered" or "coarse"
    spacing_x: float = 10e-3,
    spacing_y: float = 10e-3,
    radius: float = 50e-3,
    origin: tuple[float, float] = (0.0, 0.0),
    include_mask_name: str | None = None,
    include_threshold: float = 0.5,
    circle_resolution: int = 128,
    min_domain_fraction: float = 0.0,
    min_data_coverage: float = 0.0,
    require_center_inside_circle: bool = False,
):
    """
    Homogenize a 3-component vector field component-wise, then compute the
    magnitude of the homogenized vector.

    Returns records containing:
        value_x, value_y, value_z, value
    where:
        value = sqrt(value_x^2 + value_y^2 + value_z^2)
    """
    if normalization not in {"covered", "coarse"}:
        raise ValueError("normalization must be 'covered' or 'coarse'")

    mesh_poly, vec = _extract_vector_cell_field(mesh, field_name)

    fine_polys, _ = _build_fine_cell_polygons(
        mesh_poly,
        np.ones(mesh_poly.n_cells, dtype=np.float64),  # dummy values, geometry only
        include_mask_name=include_mask_name,
        include_threshold=include_threshold,
    )

    # Rebuild vector values only for the kept cells in the same order
    pts_xy = np.asarray(mesh_poly.points[:, :2], dtype=np.float64)

    use_mask = include_mask_name is not None
    if use_mask:
        if include_mask_name in mesh_poly.cell_data:
            mask_vals = np.asarray(mesh_poly.cell_data[include_mask_name], dtype=np.float64)
        elif include_mask_name in mesh_poly.point_data:
            tmp = mesh_poly.point_data_to_cell_data()
            mask_vals = np.asarray(tmp.cell_data[include_mask_name], dtype=np.float64)
        else:
            raise KeyError(f"Mask '{include_mask_name}' not found.")

    fine_vecs = []
    for icell, ids in enumerate(_iter_face_point_ids(mesh_poly)):
        if use_mask and not (mask_vals[icell] > include_threshold):
            continue

        poly_xy = pts_xy[ids]
        poly = Polygon(poly_xy)
        if poly.is_empty or (not poly.is_valid) or poly.area <= 0.0:
            continue

        v = vec[icell]
        if not np.all(np.isfinite(v)):
            continue

        fine_vecs.append(v.astype(np.float64))

    fine_vecs = np.asarray(fine_vecs, dtype=np.float64)

    if len(fine_polys) != len(fine_vecs):
        raise RuntimeError("Geometry/value mismatch while building vector cells.")

    tree = STRtree(fine_polys)
    geom_id_to_index = {id(g): i for i, g in enumerate(fine_polys)}

    x0, y0 = origin
    circle = Point(x0, y0).buffer(radius, quad_segs=circle_resolution)

    nx_max = int(np.ceil(radius / spacing_x)) + 2
    ny_max = int(np.ceil(radius / spacing_y)) + 2

    records = []

    for nx in range(-nx_max, nx_max + 1):
        for ny in range(-ny_max, ny_max + 1):
            if (nx + ny) % 2 != 0:
                continue

            cx = x0 + nx * spacing_x
            cy = y0 + ny * spacing_y

            if require_center_inside_circle:
                if (cx - x0) ** 2 + (cy - y0) ** 2 > radius**2 + 1e-15:
                    continue

            diamond = _diamond_cell(cx, cy, spacing_x, spacing_y)
            full_area = diamond.area

            clipped_cell = diamond.intersection(circle)
            clipped_area = clipped_cell.area

            if clipped_cell.is_empty or clipped_area <= 0.0:
                continue

            domain_fraction = clipped_area / full_area
            if domain_fraction < min_domain_fraction:
                continue

            cand_ids = _query_tree_indices(tree, clipped_cell, geom_id_to_index)

            covered_area = 0.0
            vec_area_sum = np.zeros(3, dtype=np.float64)

            for idx in cand_ids:
                inter = fine_polys[idx].intersection(clipped_cell)
                a = inter.area
                if a <= 0.0:
                    continue

                covered_area += a
                vec_area_sum += fine_vecs[idx] * a

            if covered_area <= 0.0:
                continue

            data_coverage = covered_area / clipped_area
            if data_coverage < min_data_coverage:
                continue

            denom = covered_area if normalization == "covered" else clipped_area
            vbar = vec_area_sum / denom
            vmag = np.linalg.norm(vbar)

            rec = {
                "nx": nx,
                "ny": ny,
                "x": cx,
                "y": cy,
                "area": denom,
                "covered_area": covered_area,
                "value_x": float(vbar[0]),
                "value_y": float(vbar[1]),
                "value_z": float(vbar[2]),
                "value": float(vmag),
                "polygon": clipped_cell,
                "domain_fraction": domain_fraction,
                "data_coverage": data_coverage,
            }

            if normalization == "coarse":
                rec["value_area_sum_x"] = float(vec_area_sum[0])
                rec["value_area_sum_y"] = float(vec_area_sum[1])
                rec["value_area_sum_z"] = float(vec_area_sum[2])

            records.append(rec)

    return {
        "records": records,
        "field_name": field_name,
        "spacing_x": spacing_x,
        "spacing_y": spacing_y,
        "radius": radius,
        "origin": origin,
        "include_mask_name": include_mask_name,
        "min_domain_fraction": min_domain_fraction,
        "min_data_coverage": min_data_coverage,
        "require_center_inside_circle": require_center_inside_circle,
        "normalization": "covered_area" if normalization == "covered" else "coarse_clipped_area",
        "vector_post_magnitude": True,
    }
    
def compute_mass_flux_vector(
    mesh: pv.DataSet,
    *,
    velocity_name: str = "U",
    alpha_name: str = "alpha.gas",
    rhoL: float,
    rhoG: float,
    mode: str = "gas",
) -> tuple[pv.PolyData, np.ndarray]:
    """
    Compute vector mass flux:
        gas     -> rhoG * alpha * U
        liquid  -> rhoL * (1-alpha) * U
        mixture -> (rhoG*alpha + rhoL*(1-alpha)) * U
    """
    work = mesh
    if not isinstance(work, pv.PolyData):
        work = work.extract_surface()

    needed = [velocity_name, alpha_name]
    need_point_to_cell = False
    for name in needed:
        if name in work.cell_data:
            continue
        elif name in work.point_data:
            need_point_to_cell = True
        else:
            raise KeyError(
                f"Field '{name}' not found.\n"
                f"Point arrays: {list(work.point_data.keys())}\n"
                f"Cell arrays : {list(work.cell_data.keys())}"
            )

    if need_point_to_cell:
        work = work.point_data_to_cell_data()

    U = np.asarray(work.cell_data[velocity_name], dtype=np.float64)
    alpha = np.asarray(work.cell_data[alpha_name], dtype=np.float64).squeeze()

    if U.ndim != 2 or U.shape[1] != 3:
        raise ValueError(f"'{velocity_name}' must be a 3-component vector. Got {U.shape}")
    if alpha.ndim != 1 or alpha.shape[0] != U.shape[0]:
        raise ValueError(
            f"'{alpha_name}' must be scalar per cell and match U cell count. "
            f"Got {alpha.shape} vs {U.shape[0]}"
        )

    if mode == "gas":
        rho_eff = rhoG * alpha
    elif mode == "liquid":
        rho_eff = rhoL * (1.0 - alpha)
    elif mode == "mixture":
        rho_eff = rhoG * alpha + rhoL * (1.0 - alpha)
    else:
        raise ValueError("mode must be 'gas', 'liquid', or 'mixture'")

    j = rho_eff[:, None] * U
    return work, j    
# ------------------------------------------------------------------
# Example usage
# ------------------------------------------------------------------
if __name__ == "__main__":
    CASES = [
        CaseSpec("Base", "/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/base/surface/10.5/plane3.vtp"),
        CaseSpec("noButtons", "/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/noDimples/surfaces/11.16/plane3.vtp"),
    ]
    CASESIe = [
        CaseSpec("Base", "/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/base/surface/10.5/plane0.vtp"),
        CaseSpec("LRIM", "/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/LRim/surfaces/11.4/plane0.vtp"),
        CaseSpec("DXDY", "/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/dxdy/surfaces/11.6/plane0.vtp"),
        CaseSpec("AO", "/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/AO/surfaces/11.3/plane0.vtp"),
        CaseSpec("R45", "/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/RA45/surfaces/11.1/plane0.vtp"),
        CaseSpec("noButtons", "/home/andreas/OpenFOAM/andreas-v2206/run/hisDriftFluxFoam/paper3/postProcessRaw/noDimples/surfaces/11.16/plane0.vtp"),
    ]

    FIELDS = [
        #FieldSpec(out_name="alpha.gas", kind="raw", field_name="alpha.gas"),
        #FieldSpec(out_name="Umag", kind="raw", field_name="U", magnitude=True),
        FieldSpec(out_name="UmagNew", kind="raw", field_name="U", post_component_magnitude=True), 
        #FieldSpec(out_name="Ux", kind="raw", field_name="U", component=0), 
        #FieldSpec(out_name="IeMag", kind="raw", field_name="Ie", magnitude=True),
        #FieldSpec(            out_name="JMix",            kind="mass_flux",            velocity_name="U",            alpha_name="alpha.gas",            rhoL=1231.0,            rhoG=0.126,            mode="mixture",            magnitude=True,        ),
        #FieldSpec(
        #    out_name="JGas",
        #    kind="mass_flux",
        #    velocity_name="U",
        #    alpha_name="alpha.gas",
        #    rhoL=1231.0,
        #    rhoG=0.126,
        #    mode="gas",
        #    magnitude=True,
        #),
    ]
    FIELDSIe = [
     FieldSpec(out_name="IeMag", kind="raw", field_name="Ie", magnitude=True),
    ]

    HCFG = HomogenizeConfig(
        spacing_x=10e-3,
        spacing_y=10e-3,
        radius=50e-3,
        origin=(0.0, 0.0),
        min_domain_fraction=0.1,
        min_data_coverage=0.1,
        require_center_inside_circle=False,
    )

    run_for_methods(
        CASES,
        FIELDS,
        HCFG,
        methods=["coarse"],#"covered","coarse"
        out_root="stored",
        save_polygons=True,
        verbose=True,
    )
    #
   # run_for_methods(
   #     CASESIe,
   #     FIELDSIe,
   #     HCFG,
   #     methods=["coarse"],
   #     out_root="stored",
   #     save_polygons=True,
   #     verbose=True,
   # )
    
