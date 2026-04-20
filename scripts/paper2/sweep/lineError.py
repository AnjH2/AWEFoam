"""
Compare simulation line data against experimental reference curves and compute
per-height, per-case, and overall cost metrics.

What this script does
---------------------
This module implements the comparison logic used by the optimization workflow.

Its main job is to:
- read experimental and simulation line data files
- align both datasets on a common x-axis
- optionally shift simulation x-values relative to experiment
- optionally exclude specified x-ranges from the comparison
- optionally filter points using a minimum alpha threshold
- interpolate simulation values onto experimental x-points
- compute error metrics such as RMSE, MAE, and R2
- aggregate results across heights and cases
- write plots and structured logs for each run

Typical use case
----------------
`run_experiment.py` calls `curveCompare(...)` after a batch has finished. This
returns a scalar cost used for ranking the parameter set, plus optional
breakdown information and saved figures.

Core functions
--------------
- `compare_curves_from_files(...)`
    Compare one experimental file against one simulation file at one height.

- `curveCompare(...)`
    High-level aggregation across:
    - multiple cases
    - multiple heights per case
    - weighted averaging across cases

- `save_comparison_plot(...)`
    Save a scatter plot of experimental vs interpolated simulation values.

- `save_case_summary(...)`
    Save a multi-panel summary plot for one case.

- `save_run_logs(...)`
    Write structured per-run outputs such as:
    - `OVERALL.txt`
    - `*_breakdown.json`
    - `per_case.csv`
    - `per_height.csv`

Metrics
-------
The module computes:
- `MSE`
- `MAE`
- `RMSE`
- `R2`
- `RMSE%`
- `MAE%`

The chosen optimization cost can be one of:
- `RMSE`
- `MSE`
- `MAE`
- `RMSE%`
- `MAE%`
- `1-R2`

Comparison procedure
--------------------
For each height:
1. read experimental data
2. read simulation data
3. convert simulation x from meters to millimeters
4. optionally flip the simulation x-axis
5. optionally shift simulation x relative to experiment
6. restrict to the overlapping x-range
7. interpolate simulation y onto experimental x-points
8. exclude user-specified x-ranges if configured
9. filter out low-alpha points if configured
10. compute metrics on the remaining points

Outputs
-------
This module can write:
- per-height scatter plots
- per-case summary figures
- run-level text, JSON, and CSV logs

These outputs are stored under a run-specific folder such as:

    <figures_root>/<run_tag>/

Important notes
---------------
- Experimental values are divided by 100 inside `_read_experimental(...)`.
- Simulation x-values are converted from meters to millimeters.
- `_read_simulation(...)` flips the x-axis by default.
- Missing cases can be penalized through `missing_case_penalty`.
- Exclusion ranges are defined on the experimental x-axis in millimeters.
- When `return_breakdown=True`, `curveCompare(...)` returns both the scalar cost
  and a detailed nested breakdown dictionary.

Assumptions
-----------
- Experimental files contain at least two columns: x and value.
- Simulation files contain at least two columns: x and value.
- File paths provided by the path provider are valid.
- Cases and heights are correctly matched to the experiment definition.

Limitations
-----------
- Path conventions are tightly coupled to the current project folder structure.
- Some helper functions still reflect older hard-coded workflows.
- Plotting and logging are mixed into the same module as the metric logic.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Dict, Any, List, Optional, Callable
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import re, os
import json
import csv
from datetime import datetime

# -------------------------------
# Core comparator (Option 1)
# -------------------------------

@dataclass
class ComparisonResult:
    x_eval_mm: np.ndarray           # experimental x used after overlap/filter
    y_exp: np.ndarray               # experimental y used after overlap/filter
    y_sim_interp: np.ndarray        # interpolated sim y on x_eval_mm
    metrics: Dict[str, float]
    x_range_mm: Tuple[float, float]
    # NEW: raw arrays (full files, sorted). sim x has shift applied.
    raw_x_exp_mm: np.ndarray
    raw_y_exp: np.ndarray
    raw_x_sim_mm: np.ndarray
    raw_y_sim: np.ndarray

def _read_experimental(
    path: str | Path,
    replace_nan_with: float | None = None,   # NEW: if set, replace non-finite y with this
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Experimental: header line(s) starting with '#', columns: x_mm, value
    Applies: y /= 100.0
    If replace_nan_with is not None:
        - drop rows with non-finite x
        - replace non-finite y with replace_nan_with
      else:
        - drop rows where x or y are non-finite
    """
    path = Path(path)
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Experimental file must have at least two columns (x_mm, value): {path}")

    x_mm = data[:, 0].astype(float)
    y = data[:, 1].astype(float) / 100.0

    if replace_nan_with is None:
        # original behavior: drop rows with non-finite x OR y
        m = np.isfinite(x_mm) & np.isfinite(y)
        x_mm, y = x_mm[m], y[m]
    else:
        # agile behavior: keep rows with finite x; for y, replace non-finite with given value
        m_x = np.isfinite(x_mm)
        x_mm, y = x_mm[m_x], y[m_x]
        bad_y = ~np.isfinite(y)
        if np.any(bad_y):
            y = y.copy()
            y[bad_y] = float(replace_nan_with)

    if x_mm.size == 0:
        raise ValueError(f"Experimental file has no usable rows after cleaning: {path}")
    return x_mm, y


def _read_simulation(path: str | Path, flip: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    path = Path(path)
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Simulation file must have at least two columns (x_m, y): {path}")
    x_m = data[:, 0].astype(float)
    y = data[:, 1].astype(float)
    x_mm = x_m * 1000.0

    if flip:
        x_mm = (x_mm.max() + x_mm.min()) - x_mm
        # keep pairing correct: reverse the order so x is ascending
        order = np.argsort(x_mm)
        x_mm, y = x_mm[order], y[order]

    return x_mm, y

def _r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    
def _fmt_num(x: float) -> str:
    """Nice filename number: drop trailing .0; keep decimals if needed."""
    x = round(float(x), 6)     # avoid 299.999999 artifacts
    return f"{x:g}"            # 'g' trims trailing zeros and the decimal point
    
# remove _fmt_height, or keep but call _fmt_num inside it:
def _fmt_height(h: float) -> str:
    return _fmt_num(h)
    
def _fmt_tag_val(v):
    # compact, filesystem-safe
    s = f"{float(v):g}" if isinstance(v, (float, int)) else str(v)
    s = s.replace("-", "m").replace(".", "p")
    return s

def make_run_tag_from_x(x):
    labels = _labels_for_x()
    # if lengths mismatch, zip will stop at the shorter one (safe)
    parts = [f"{lab}{_fmt_tag_val(val)}" for lab, val in zip(labels, x)]
    return "_".join(parts)

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path
def _build_exclusion_keep_mask(x_mm: np.ndarray, ranges: Optional[List[Tuple[float, float]]]):
    """
    Returns a boolean mask of points to KEEP (i.e., not inside any excluded range).
    'ranges' is a list of (lo, hi) in mm on the experimental x-axis.
    """
    if not ranges:
        return np.ones_like(x_mm, dtype=bool)
    keep = np.ones_like(x_mm, dtype=bool)
    for lo, hi in ranges:
        if lo is None or hi is None:
            continue
        lof, hif = float(min(lo, hi)), float(max(lo, hi))
        keep &= ~((x_mm >= lof) & (x_mm <= hif))
    return keep    
def _labels_for_x():
    # order must match how you build x in wrapper(...) XXCHECKXX WHEN CHANGEING VARS
    return ["SG","SC","DO","C0Ne","C0Pe","CDNe","CDPe"]

def compare_curves_from_files(
    exp_path: str | Path,
    sim_path: str | Path,
    min_alpha: float | None = None,
    sim_to_exp_shift_mm: float = 0.0,
    exclude_ranges: Optional[List[Tuple[float, float]]] = None,  # NEW
) -> ComparisonResult:
    """
    Interpolate simulation onto experimental x-points (within overlap) and compute metrics.
    Applies:
      - experimental y /= 100
      - simulation x: meters -> mm
      - optional filter: drop points where y_sim_interp <= min_alpha
    Metrics: MSE, MAE, RMSE, R2, N
    """
    x_exp_mm, y_exp = _read_experimental(exp_path,replace_nan_with=min_alpha)
    x_sim_mm, y_sim = _read_simulation(sim_path)
    # Apply shift so sim x is expressed in the experimental-relative frame
    if sim_to_exp_shift_mm:
        x_sim_mm = x_sim_mm + float(sim_to_exp_shift_mm)
    # sort by x
    sE = np.argsort(x_exp_mm)
    x_exp_mm, y_exp = x_exp_mm[sE], y_exp[sE]
    sS = np.argsort(x_sim_mm)
    x_sim_mm, y_sim = x_sim_mm[sS], y_sim[sS]
    raw_x_exp_mm = x_exp_mm.copy()
    raw_y_exp    = y_exp.copy()
    raw_x_sim_mm = x_sim_mm.copy()      # AFTER shift
    raw_y_sim    = y_sim.copy()
    # overlap
    xmin = max(float(np.min(x_exp_mm)), float(np.min(x_sim_mm)))
    xmax = min(float(np.max(x_exp_mm)), float(np.max(x_sim_mm)))
    if xmin >= xmax:
        raise ValueError("No overlapping x-range between experimental and simulation data.")

    mask = (x_exp_mm >= xmin) & (x_exp_mm <= xmax)
    if not np.any(mask):
        raise ValueError("No experimental points fall within the overlapping x-range.")

    # 1) select exp in overlap and interpolate sim
    x_eval = x_exp_mm[mask]
    y_exp_eval = y_exp[mask]
    y_sim_interp = np.interp(x_eval, x_sim_mm, y_sim)
    
    n_exp_nan_before = int(np.sum(~np.isfinite(y_exp_eval)))
    print(n_exp_nan_before)
    if n_exp_nan_before > 0 and (min_alpha is not None):
        y_exp_eval = np.where(np.isfinite(y_exp_eval), y_exp_eval, float(min_alpha))
        
    # 2) exclude user ranges (on experimental x-axis)
    keep_excl = _build_exclusion_keep_mask(x_eval, exclude_ranges)
    x_eval       = x_eval[keep_excl]
    y_exp_eval   = y_exp_eval[keep_excl]
    y_sim_interp = y_sim_interp[keep_excl]
    
        
    # 3) min_alpha threshold on interpolated sim
    if min_alpha is not None:
        keep_alpha = y_sim_interp > float(min_alpha)
        x_eval       = x_eval[keep_alpha]
        y_exp_eval   = y_exp_eval[keep_alpha]
        y_sim_interp = y_sim_interp[keep_alpha]

    # 4) final non-finite cleanup AFTER all filters
    finite = np.isfinite(x_eval) & np.isfinite(y_exp_eval) & np.isfinite(y_sim_interp)
    n_before = int(finite.size)
    x_eval       = x_eval[finite]
    y_exp_eval   = y_exp_eval[finite]
    y_sim_interp = y_sim_interp[finite]
    n_after = int(x_eval.size)
    n_interp_nan_dropped = n_before - n_after

    if x_eval.size == 0:
        raise ValueError("No points left after filtering / overlap.")

    # 5) metrics
    diff = y_exp_eval - y_sim_interp
    mse = float(np.mean(diff**2))
    mae = float(np.mean(np.abs(diff)))
    rmse = float(np.sqrt(mse))
    r2 = float(_r2_score(y_exp_eval, y_sim_interp))

    mean_abs_yexp = float(np.mean(np.abs(y_exp_eval))) if y_exp_eval.size > 0 else float("nan")
    rmse_pct = _safe_pct(rmse, mean_abs_yexp)
    mae_pct  = _safe_pct(mae,  mean_abs_yexp)

    metrics = {
    "MSE": mse,
    "MAE": mae,
    "RMSE": rmse,
    "R2": r2,
    "N": int(x_eval.size),
    "N_overlap": int(np.count_nonzero(mask)),
    "N_excluded_ranges": int(np.count_nonzero(~keep_excl)) if exclude_ranges else 0,
    "N_interp_nan_dropped": int(n_interp_nan_dropped),
    "RMSE%": rmse_pct,
    "MAE%":  mae_pct,
    }

    return ComparisonResult(
        x_eval_mm=x_eval,
        y_exp=y_exp_eval,
        y_sim_interp=y_sim_interp,
        metrics=metrics,
        x_range_mm=(xmin, xmax),
        raw_x_exp_mm=raw_x_exp_mm,
        raw_y_exp=raw_y_exp,
        raw_x_sim_mm=raw_x_sim_mm,
        raw_y_sim=raw_y_sim,
    )

# -------------------------------
# Higher-level aggregation
# -------------------------------



def _find_existing_sim(job_ids: List[str], template: Callable[[str], str]) -> Optional[tuple[str, str]]:
    """
    Try all job_ids to construct a candidate path via template(job_id),
    return the newest existing (path, jid), or None if none exist.
    """
    candidates: List[tuple[str, str]] = []
    for jid in job_ids:
        p = template(jid)
        if os.path.exists(p):
            candidates.append((p, jid))
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    candidates.sort(key=lambda pj: os.path.getmtime(pj[0]), reverse=True)
    return candidates[0]

def _heights_and_paths(case_idx: int, case_entry: List[float], batch_index: int):
    """
    Returns (heights, exp_path_fn, sim_path_fn)
    sim_path_fn expects (jid, height)
    Adjust the file naming here if your paths are slightly different.
    """
    if case_idx <= 3:
        heights = [10, 30, 50]
        def exp_path(h):
            current = case_entry[1] * 0.1  # A/m² -> mA/cm²
            return ("/home/anj/experimental_results/beam1/extractedData_lineX_b10/"
                    f"{_fmt_num(current)}_{_fmt_height(h)}.txt")
        def sim_path(jid, h):
            return (f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/doubleStructuredSearch/"
                    f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/5/s{_fmt_height(h)}_alpha.gas.xy")
    else:
        heights = [9.5, 22, 35]
        def exp_path(h):
            return f"/home/anj/experimental_results/beam3/extractedData_lineX_centered/staggered_y{_fmt_height(h)}.txt"
        def sim_path(jid, h):
            # If you truly need (h-3) in filename, change _fmt_height(h) to _fmt_height(h-3).
            return (f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/doubleStructuredSearch/"
                    f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/10/s{_fmt_height(h-3)}_alpha.gas.xy")
    return heights, exp_path, sim_path

def _metric_to_cost_value(metrics: Dict[str, float], metric_name: str) -> Optional[float]:
    """
    Convert chosen metric into a scalar where smaller is better:
      - "RMSE", "MSE", "MAE": used as-is (lower=better)
      - "1-R2": transforms R2 into cost = 1 - R2 (so higher R2 => lower cost)
    Returns None if metric missing or non-finite.
    """
    val = None
    if metric_name in ("RMSE", "MSE", "MAE","RMSE%", "MAE%"):
        val = metrics.get(metric_name, None)
    elif metric_name == "1-R2":
        r2 = metrics.get("R2", None)
        if r2 is not None and np.isfinite(r2):
            val = 1.0 - r2
    else:
        raise ValueError(f"Unsupported metric_name: {metric_name}")

    if val is None or not np.isfinite(val):
        return None
    return float(val)
def _safe_pct(numer: float, denom: float) -> float:
    """Return numer/denom*100, guarding zero/NaN denom."""
    if denom is None:
        return float("nan")
    try:
        denom = float(denom)
    except Exception:
        return float("nan")
    return float(numer / denom * 100.0) if np.isfinite(denom) and abs(denom) > 0 else float("nan")
def curveCompare(
    CaseList_, job_ids, *,
    batch_index: int,
    min_alpha: float = 1e-4,
    case_weights=None,
    metric_name: str = "RMSE",
    per_case_sim_shift_mm=None,
    per_case_flip=None,                  # if you’re using flip globally
    per_case_exclude_ranges: Optional[List[Optional[List[Tuple[float, float]]]]] = None,  # NEW
    path_provider: Optional[
        Callable[[int, List[float], int],
                 Tuple[List[float], Callable[[float], str], Callable[[str, float], str]]]
    ] = None,  # NEW
    return_breakdown: bool = False,
    outdir_root: str = "./figures",
    run_tag: str | None = None,
    save_per_height_scatter: bool = False,
    write_logs: bool = True,                    # NEW
    x_params: Optional[List[float]] = None,     # NEW
    missing_case_penalty: float | None = None,
):
    """
    For each case:
      - decide heights & build file paths
      - pick a SIM file by trying each job_id and choosing an existing (newest)
      - compute metrics for each height (skipping missing/failed)
      - average chosen metric across heights -> case value (simple mean of available)
    Then:
      - weighted average across cases -> scalar cost

    Returns:
      cost  (and breakdown dict if return_breakdown=True)
    """
    K = len(CaseList_)
    run_tag = run_tag or f"OBatch{batch_index}"
    run_dir = ensure_dir(os.path.join(outdir_root, run_tag))
    if case_weights is None:
        weights = np.ones(K, dtype=float)
    else:
        weights = np.asarray(case_weights, dtype=float)
        if weights.shape != (K,):
            raise ValueError(f"case_weights must have length {K}")
    if per_case_sim_shift_mm is None:
        per_case_sim_shift_mm = [0.0] * K
    elif len(per_case_sim_shift_mm) != K:
        raise ValueError(f"per_case_sim_shift_mm must have length {K}")

    per_case_vals: List[Optional[float]] = []
    breakdown = {"per_case": [], "overall": {}}
    if per_case_exclude_ranges is None:
        per_case_exclude_ranges = [None] * K
    elif len(per_case_exclude_ranges) != K:
        raise ValueError(f"per_case_exclude_ranges must have length {K}")
    print("1***")
    for k, case_entry in enumerate(CaseList_):
        print("2***"+str(k))
        try:
            current_mAcm2 = _fmt_num(case_entry[1] * 0.1)
            case_dir_name = f"case_{k}_{current_mAcm2}"
        except Exception:
            case_dir_name = f"case_{k}"

        case_dir = ensure_dir(os.path.join(run_dir, case_dir_name))
        if path_provider is None:
            heights, exp_path_fn, sim_path_fn = _heights_and_paths(k, case_entry, batch_index)
        else:
            heights, exp_path_fn, sim_path_fn = path_provider(k, case_entry, batch_index)

        height_records = []
        per_height_cost_vals: List[float] = []
        fig_items = []
        for h in heights:
            shift_mm = float(per_case_sim_shift_mm[k])
            print(f"[shift] case {k}, h={h}: applying {shift_mm} mm")
            exp_path = exp_path_fn(h)
            print("3***"+exp_path)
            if not os.path.exists(exp_path):
                height_records.append({
                    "height": h, "status": "missing_exp", "exp_path": exp_path, "sim_path": None, "metrics": None
                })
                continue

            # find an existing SIM file among job_ids
            found = _find_existing_sim(job_ids, lambda jid: sim_path_fn(jid, h))
            
            
                        
            if found is None:
                height_records.append({
                    "height": h, "status": "missing_sim", "exp_path": exp_path, "sim_path": None, "metrics": None
                })
                continue
            sim_path, sim_jid = found  
            print("4***"+sim_path) 
            try:
                print("5***")
                shift_mm = float(per_case_sim_shift_mm[k])
                ranges = per_case_exclude_ranges[k]
                print(ranges)
                res = compare_curves_from_files(exp_path, sim_path, min_alpha=min_alpha,sim_to_exp_shift_mm=shift_mm, exclude_ranges=ranges)
                print(metric_name)   
                print(res.metrics)             
                cost_val = _metric_to_cost_value(res.metrics, metric_name)
                print(cost_val)
                height_records.append({
                    "height": h,
                    "status": "ok",
                    "exp_path": exp_path,
                    "sim_path": sim_path,
                    "jid": sim_jid,
                    "metrics": res.metrics,
                    "x_range_mm": res.x_range_mm,
                })
                
                if save_per_height_scatter:
                    try:
                        save_comparison_plot(res, batch_index, k, sim_jid, h, outdir=case_dir)
                    except Exception as pe:
                        print(f"[WARN] Plot failed for case {k}, height {h}: {pe}")
                # collect for summary
                fig_items.append({"height": h, "jid": sim_jid, "res": res, "exclude_ranges": ranges})
                if cost_val is not None:
                    per_height_cost_vals.append(cost_val)
            except Exception as e:
                height_records.append({
                    "height": h, "status": f"error: {e}", "exp_path": exp_path, "sim_path": sim_path, "metrics": None
                })

        # mean across heights that yielded a finite value
        case_val = float(np.mean(per_height_cost_vals)) if len(per_height_cost_vals) > 0 else None
        per_case_vals.append(case_val)

        breakdown["per_case"].append({
            "case_index": k,
            "case_inputs": case_entry,
            "heights": heights,
            "metric_name": metric_name,
            "case_value": case_val,
            "per_height": height_records,
            "case_weights": weights.tolist()
        })
        if fig_items:
            try:
                save_case_summary(fig_items, batch_index, k, outdir=case_dir)
            except Exception as fe:
                print(f"[WARN] Case summary plot failed for case {k}: {fe}")

    # --- weighted average across cases ---
    # Option A: Penalize missing cases so they still contribute to the overall cost
    num_penalized = 0
    if missing_case_penalty is not None:
        penalty = float(missing_case_penalty)
        per_case_vals = [
            (penalty if v is None else float(v))
            for v in per_case_vals
        ]
        num_penalized = sum(1 for v in per_case_vals if v == penalty)

        # now compute with all finite values (skip any accidental NaNs)
        vals = np.array(per_case_vals, dtype=float)
        mask_ok = np.isfinite(vals)
    else:
        # original behavior: only use cases that produced a numeric value
        vals = np.array([v for v in per_case_vals if v is not None], dtype=float)
        mask_ok = np.ones_like(vals, dtype=bool)  # all 'vals' are ok in this branch

    if vals.size == 0 or not np.any(mask_ok):
        raise RuntimeError("No usable metrics found for any case (missing files or failures).")

    # align weights with the vector we’re averaging
    if missing_case_penalty is not None:
        weights_vec = np.asarray(weights, dtype=float)
        weights_ok = weights_vec[mask_ok]
        vals_ok = vals[mask_ok]
    else:
        # we filtered out None cases to build 'vals'. Filter weights the same way:
        weights_vec = np.asarray([w for v, w in zip(per_case_vals, weights) if v is not None], dtype=float)
        weights_ok = weights_vec[mask_ok]
        vals_ok = vals[mask_ok]

    # fallback if weights are all zero/invalid
    if not np.any(weights_ok > 0):
        weights_ok = np.ones_like(vals_ok)

    cost = float(np.dot(vals_ok, weights_ok) / np.sum(weights_ok))

    breakdown["overall"] = {
        "metric_name": metric_name,
        "min_alpha": float(min_alpha),
        "weighted_cost": cost,
        "per_case_values": [None if v is None else float(v) for v in per_case_vals],
        "effective_weights_used": weights_ok.tolist(),
        "missing_case_penalty": None if missing_case_penalty is None else float(missing_case_penalty),
        "num_penalized_cases": int(num_penalized),
    }

    if write_logs:
        try:
            save_run_logs(run_dir, run_tag, batch_index, breakdown, CaseList_, x_params)
        except Exception as le:
            print(f"[WARN] Failed to write run logs: {le}")
    return (cost, breakdown) if return_breakdown else cost

def print_curve_info(info):
    print("\n=== Curve Comparison Breakdown ===")
    print(f"Overall cost ({info['overall']['metric_name']}): {info['overall']['weighted_cost']:.6g}")
    print(f"min_alpha: {info['overall']['min_alpha']}")
    if "effective_weights_used" in info["overall"]:
        print(f"Effective weights used: {info['overall']['effective_weights_used']}")
    print("")

    for case in info["per_case"]:
        idx = case["case_index"]
        inputs = case["case_inputs"]
        case_val = case["case_value"]
        print(f"Case {idx} inputs={inputs}")
        if case_val is None:
            print("  Case value: None (no valid heights)")
        else:
            print(f"  Case value ({case['metric_name']}): {case_val:.6g}")

        for hrec in case["per_height"]:
            h = hrec["height"]
            status = hrec["status"]
            if status == "ok":
                m = hrec["metrics"]
                extra = ""
                if "N_excluded_ranges" in m:
                    extra = f", N_excl={m['N_excluded_ranges']}"
                print(f"    Height {h:>5}: ok | RMSE={m['RMSE']:.4g}, MSE={m['MSE']:.4g}, "
                      f"MAE={m['MAE']:.4g}, R2={m['R2']:.4g}, N={m['N']}{extra}")
        print("")
def save_comparison_plot(res, batch_index, case_idx, jid, height, outdir="./figures"):
    ensure_dir(outdir)
    tag = f"OBatch{batch_index}_{case_idx}-{jid}_h{height}"
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(res.y_exp, res.y_sim_interp, s=15, alpha=0.7, label="points")
    lo = min(np.min(res.y_exp), np.min(res.y_sim_interp))
    hi = max(np.max(res.y_exp), np.max(res.y_sim_interp))
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = 0.0, 1.0
    ax.plot([lo, hi], [lo, hi], "--", lw=1, label="ideal y=x")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Experimental"); ax.set_ylabel("Simulation (interp)")
    ax.set_title(tag); ax.legend(loc="lower right", fontsize=8, frameon=False)
    for ext in ("svg","pdf"):
        fig.savefig(os.path.join(outdir, f"{tag}_scatter.{ext}"), bbox_inches="tight")
    plt.close(fig)

def save_case_summary(fig_items, batch_index, case_idx, outdir="./figures"):
    ensure_dir(outdir)
    heights_sorted = sorted(fig_items, key=lambda d: float(d["height"]))
    jids = list({d["jid"] for d in heights_sorted})
    jid_tag = jids[0] if len(jids) == 1 else "mixedJIDs"
    tag = f"OBatch{batch_index}_{case_idx}-{jid_tag}_summary"

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for r in range(2):
        for c in range(3):
            axes[r, c].axis("off")

    # --- Top row: scatter (Exp vs Sim-interp). No x-range shading here ---
    for col, item in enumerate(heights_sorted[:3]):
        ax = axes[0, col]; ax.axis("on")
        res = item["res"]; h = item["height"]

        ax.scatter(res.y_exp, res.y_sim_interp, s=10, alpha=0.7, label="points")
        lo = min(np.min(res.y_exp), np.min(res.y_sim_interp))
        hi = max(np.max(res.y_exp), np.max(res.y_sim_interp))
        if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
            lo, hi = 0.0, 1.0
        ax.plot([lo, hi], [lo, hi], "--", lw=1, label="ideal y=x")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)

        m = res.metrics
        ax.set_title(f"h={h} | RMSE={m['RMSE']:.3g}, MAE={m['MAE']:.3g}, R2={m['R2']:.3g}, N={m['N']}", fontsize=10)
        ax.set_xlabel("Exp")
        if col == 0: ax.set_ylabel("Sim (interp)")
        ax.legend(loc="lower right", fontsize=8, frameon=False)

    # --- Bottom row: raw curves (x vs y), shade excluded x-ranges on EXP axis ---
    for col, item in enumerate(heights_sorted[:3]):
        ax = axes[1, col]; ax.axis("on")
        res = item["res"]; h = item["height"]

        ax.plot(res.raw_x_exp_mm, res.raw_y_exp, lw=1.2, label="Exp (raw)")
        ax.plot(res.raw_x_sim_mm, res.raw_y_sim, lw=1.2, label="Sim (raw, shifted)")

        # Shade excluded ranges (if any) on the x-axis
        ranges = item.get("exclude_ranges", None)
        if ranges:
            for lo, hi in ranges:
                lof, hif = float(min(lo, hi)), float(max(lo, hi))
                ax.axvspan(lof, hif, color="0.85", alpha=0.5, lw=0)

        ax.set_xlabel("x (mm)")
        if col == 0: ax.set_ylabel("value")
        ax.set_title(f"h={h} | raw curves", fontsize=10)
        ax.legend(loc="best", fontsize=8, frameon=False)

    fig.suptitle(f"Case summary: OBatch{batch_index}_{case_idx}  (jid: {jid_tag})", fontsize=12, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    for ext in ("svg","pdf"):
        fig.savefig(os.path.join(outdir, f"{tag}.{ext}"), bbox_inches="tight")
    plt.close(fig)
def save_run_logs(
    run_dir: str,
    run_tag: str,
    batch_index: int,
    info: dict,
    CaseList: List[List[float]],
    x_params: Optional[List[float]] = None
):
    """
    Write per-batch logs into run_dir:
      - OVERALL.txt (human readable)
      - breakdown.json (full structure)
      - per_case.csv
      - per_height.csv
    """
    ensure_dir(run_dir)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    overall_txt = os.path.join(run_dir, "OVERALL.txt")
    breakdown_json = os.path.join(run_dir, f"{run_tag}_breakdown.json")
    per_case_csv = os.path.join(run_dir, "per_case.csv")
    per_height_csv = os.path.join(run_dir, "per_height.csv")

    metric = info["overall"]["metric_name"]
    cost = info["overall"]["weighted_cost"]
    min_alpha = info["overall"]["min_alpha"]
    eff_w = info["overall"].get("effective_weights_used", [])
    case_weights = info["overall"].get("case_weights", [])

    # ---- OVERALL.txt (human readable) ----
    with open(overall_txt, "w") as f:
        f.write(f"Timestamp     : {ts}\n")
        f.write(f"Run tag       : {run_tag}\n")
        f.write(f"Batch index   : {batch_index}\n")
        if x_params is not None:
            labels = _labels_for_x()
            pairs = ", ".join(f"{lab}={_fmt_num(val)}" for lab, val in zip(labels, x_params or []))
            f.write(f"Params (x)    : {pairs}\n")
        f.write(f"Metric        : {metric}\n")
        f.write(f"Min alpha     : {min_alpha}\n")
        f.write(f"Weighted cost : {cost:.6g}\n")
        f.write(f"Case weights  : {case_weights}\n")
        f.write(f"Eff. weights  : {eff_w}\n")
        f.write("\nPer-case values:\n")
        for pc in info["per_case"]:
            f.write(f"  case {pc['case_index']}: {pc['case_value']}\n")

    # ---- breakdown.json (full) ----
    with open(breakdown_json, "w") as jf:
        json.dump(info, jf, indent=2)

    # ---- per_case.csv ----
    with open(per_case_csv, "w", newline="") as cf:
        w = csv.writer(cf)
        w.writerow(["case_index","inputs","case_value","metric_name","weight"])
        for idx, pc in enumerate(info["per_case"]):
            weight = case_weights[idx] if idx < len(case_weights) else ""
            w.writerow([pc["case_index"], repr(pc["case_inputs"]), pc["case_value"], pc["metric_name"], weight])

    # ---- per_height.csv ----
    with open(per_height_csv, "w", newline="") as hf:
        w = csv.writer(hf)
        w.writerow([
            "case_index","height","status","RMSE","RMSE%","MSE","MAE","MAE%","R2","N",
            "N_excluded_ranges","x_min_mm","x_max_mm","exp_path","sim_path","jid"
        ])
        for pc in info["per_case"]:
            case_idx = pc["case_index"]
            for hrec in pc["per_height"]:
                status = hrec["status"]
                m = hrec.get("metrics", {}) or {}
                x_rng = hrec.get("x_range_mm", (None, None))
                w.writerow([
                    case_idx,
                    hrec["height"],
                    status,
                    m.get("RMSE",""),
                    m.get("RMSE%",""),
                    m.get("MSE",""),
                    m.get("MAE",""),
                    m.get("MAE%",""),
                    m.get("R2",""),
                    m.get("N",""),
                    m.get("N_excluded_ranges",""),
                    x_rng[0], x_rng[1],
                    hrec.get("exp_path",""),
                    hrec.get("sim_path",""),
                    hrec.get("jid",""),
                ])
