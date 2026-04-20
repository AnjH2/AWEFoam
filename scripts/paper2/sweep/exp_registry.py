"""
Experiment registry for simulation-to-experiment comparisons.

What this script does
---------------------
This module defines named experiment configurations used by `run_experiment.py`.

Each experiment configuration describes:
- how to locate the experimental data files
- how to locate the matching simulation result files
- which heights should be compared
- optional per-case x-shifts
- optional per-case exclusion ranges
- default comparison metric
- default case weights
- default minimum alpha threshold

The registry allows the main driver script to switch between different
experimental setups using a short key such as `beam1`, `beam2`, or `beam1_2`.

Main concepts
-------------
An experiment config is a dictionary with keys such as:
- `name`
- `path_provider`
- `weights`
- `shifts_mm`
- `exclude`
- `metric`
- `min_alpha`

The most important item is `path_provider`, a callable:

    (case_idx, case_entry, batch_index) -> (heights, exp_path_fn, sim_path_fn)

where:
- `heights` is a list of comparison heights
- `exp_path_fn(h)` returns the experimental file path for height `h`
- `sim_path_fn(jid, h)` returns the simulation file path for job `jid` and height `h`

Main entry point
----------------
    load_experiment(name, K_cases)

This:
- validates the experiment key
- constructs the experiment dictionary
- broadcasts scalar settings such as `weights` or `shifts_mm` to length `K_cases`
- returns a normalized config dictionary

Available experiments
---------------------
The `EXPERIMENTS` dictionary maps short names to experiment factories, for example:
- `beam1`
- `beam2`
- `beam3`
- `beam1_2`
- `beam1_2S`

Important notes
---------------
- Paths are highly project-specific and mostly hard-coded.
- Some experiment definitions assume different simulation line folders such as
  `line/10`, `line/15`, or `samples/...`.
- Scalar values for `weights`, `shifts_mm`, and `exclude` can be broadcast to all
  cases by `load_experiment(...)`.

Assumptions
-----------
- Experimental files exist at the configured locations.
- Simulation outputs follow the expected `OBatch...` directory structure.
- The number of cases in the chosen experiment matches the active workflow.
"""
from __future__ import annotations
from typing import List, Tuple, Callable, Optional, Dict, Any
import os

# ---- helpers -------------------------------------------------------

def repeat_or_list(val_or_list, K):
    """Return a list of length K. If val_or_list is a scalar, repeat it."""
    if isinstance(val_or_list, (list, tuple)):
        if len(val_or_list) != K:
            raise ValueError(f"Expected list of length {K}, got {len(val_or_list)}")
        return list(val_or_list)
    return [val_or_list] * K

# ---- default provider (your current behavior) ----------------------

def default_path_provider() -> Callable:
    """
    Returns a function (case_idx, case_entry, batch_index) -> (heights, exp_path_fn, sim_path_fn)
    Mirrors the hardcoded logic from lineError._heights_and_paths.
    """
    def provider(case_idx: int, case_entry: List[float], batch_index: int):
        if case_idx <= 3:
            heights = [10, 30, 50]
            def exp_path(h):
                current = case_entry[1] * 0.1  # A/m² -> mA/cm²
                return ("/home/anj/experimental_results/beam1/extractedData_lineX_b10/"
                        f"{current:g}_{int(h) if float(h).is_integer() else h}.txt")
            def sim_path(jid, h):
                return (f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/"
                        f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/5/s{int(h) if float(h).is_integer() else h}_alpha.gas.xy")
        else:
            heights = [9.5, 22, 35]
            def exp_path(h):
                return f"/home/anj/experimental_results/beam3/extractedData_lineX_centered/staggered_y{int(h) if float(h).is_integer() else h}.txt"
            def sim_path(jid, h):
                return (f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/"
                        f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/10/s{int(h-3) if float(h-3).is_integer() else (h-3)}_alpha.gas.xy")
        return heights, exp_path, sim_path
    return provider

# ---- example experiment: beam1 only --------------------------------

def exp_beam1_b10() -> Dict[str, Any]:
    """
    Compare cases 0..3 to beam1 b=10 data (heights 10,30,50),
    case 4 to beam3 centered (9.5,22,35). Adjust as needed.
    """
    return {
        "name": "beam1_b10",
        "path_provider": default_path_provider(),
        # per-case arrays (length = number of cases)
        "weights":       None,                         # -> equal
        "shifts_mm":     [0.0,0.0], # example
        "exclude": [
            [(1.8, 2.65)], [(1.8, 2.65)]
        ],
        "metric": "RMSE%",
        "min_alpha": 1e-4,
    }

# ---- example experiment: fully custom paths ------------------------

def exp_beam3_centered_mixed(exp_root="/home/anj/experimental_results") -> Dict[str, Any]:
    """
    Example where every case uses beam3 centered data at heights [9.5,22,35],
    and sim files live under line/10 with s{h-3}_alpha.gas.xy names.
    """
    def provider(case_idx: int, case_entry: List[float], batch_index: int):
        heights = [9.5, 22, 35]
        def exp_path(h):
            return f"{exp_root}/beam3/extractedData_lineX_centered/staggered_y{int(h) if float(h).is_integer() else h}.txt"
        def sim_path(jid, h):
            tag = int(h-3) if float(h-3).is_integer() else (h-3)
            return (f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/"
                    f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/10/s{tag}_alpha.gas.xy")
        return heights, exp_path, sim_path

    return {
        "name": "beam3_centered_mixed",
        "path_provider": provider,
        "weights": None,
        "shifts_mm": 0.0,     # one scalar -> broadcast to all cases
        "exclude": None,
        "metric": "RMSE%",
        "min_alpha": 1e-4,
    }
def _h_token(h: float) -> str:
    # "10" instead of "10.0"; "9.5" stays as "9.5"
    try:
        hf = float(h)
        return str(int(hf)) if hf.is_integer() else f"{hf:g}"
    except Exception:
        return str(h)

def exp_beam2() -> dict:
    def path_provider(case_idx: int, case_entry: list[float], batch_index: int):
        heights = [10, 30, 50]

        def exp_path(h):
            current = float(case_entry[1]) * 0.1  # A/m² -> mA/cm²
            return (
                f"/home/anj/experimental_results/beam2/extractedData_lineX_b10/"
                f"{current:g}_{_h_token(h)}.txt"
            )

        def sim_path(jid, h):
            return (
                f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/"
                f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/7/"
                f"s{_h_token(h)}_alpha.gas.xy"
            )

        return heights, exp_path, sim_path

    return {
        "name": "beam2",
        "path_provider": path_provider,
        "weights": 1,
        "shifts_mm": -1.8,
        "exclude": [[(1.8, 2.65)]],
        "metric": "RMSE%",
        "min_alpha": 1e-4,
    }
    
def exp_beam1_beam2() -> Callable:
    """
    Returns a function (case_idx, case_entry, batch_index) -> (heights, exp_path_fn, sim_path_fn)
    Mirrors the hardcoded logic from lineError._heights_and_paths.
    """
    
    def provider(case_idx: int, case_entry: List[float], batch_index: int):
        if case_idx <= 1:
            heights = [10, 30, 50]
            def exp_path(h):
                current = case_entry[1] * 0.1  # A/m² -> mA/cm²
                return ("/home/anj/experimental_results/beam1/extractedData_lineX_b10/"
                        f"{current:g}_{_h_token(h)}.txt")
            def sim_path(jid, h):
                return (f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/"
                        f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/10/"
                        f"s{_h_token(h)}_alpha.gas.xy"
                        )
        else:
            heights = [10, 30, 50]
            def exp_path(h):
                current = float(case_entry[1]) * 0.1  # A/m² -> mA/cm²
                return (
                    f"/home/anj/experimental_results/beam2/extractedData_lineX_b10/"
                    f"{current:g}_{_h_token(h)}.txt"
                )

            def sim_path(jid, h):
                return (
                    f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/"
                    f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/15/"
                    f"s{_h_token(h)}_alpha.gas.xy"
                )

        return heights, exp_path, sim_path
    return {
        "name": "beam1_beam2",
        "path_provider": provider,
        # per-case arrays (length = number of cases)
        "weights":       None,                         # -> equal
        "shifts_mm":     [ 0.0, 0.0, -1.8,-1.8], # example
        "exclude": [
              [(1.8, 2.65)], [(1.8, 2.65)],[(1.8, 2.65)],[(1.8, 2.65)]
        ],
        "metric": "RMSE%",
        "min_alpha": 1e-4,
    }   
def exp_beam1_beam2_samples() -> Callable:
    """
    Returns a function (case_idx, case_entry, batch_index) -> (heights, exp_path_fn, sim_path_fn)
    Mirrors the hardcoded logic from lineError._heights_and_paths.
    """
    
    def provider(case_idx: int, case_entry: List[float], batch_index: int):
        if case_idx <= 1:
            heights = [10, 30, 50]
            def exp_path(h):
                current = case_entry[1] * 0.1  # A/m² -> mA/cm²
                return ("/home/anj/experimental_results/beam1/extractedData_lineX_b10/"
                        f"{current:g}_{_h_token(h)}.txt")
            def sim_path(jid, h):
                return (f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/samples/"
                        f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/10/"
                        f"s{_h_token(h)}_alpha.gas.xy"
                        )
        else:
            heights = [10, 30, 50]
            def exp_path(h):
                current = float(case_entry[1]) * 0.1  # A/m² -> mA/cm²
                return (
                    f"/home/anj/experimental_results/beam2/extractedData_lineX_b10/"
                    f"{current:g}_{_h_token(h)}.txt"
                )

            def sim_path(jid, h):
                return (
                    f"/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/samples/"
                    f"OBatch{batch_index}_{case_idx}-{jid}/postProcessing/line/15/"
                    f"s{_h_token(h)}_alpha.gas.xy"
                )

        return heights, exp_path, sim_path
    return {
        "name": "beam1_beam2_samples",
        "path_provider": provider,
        # per-case arrays (length = number of cases)
        "weights":       None,                         # -> equal
        "shifts_mm":     [ 0.0, 0.0, -1.8,-1.8], # example
        "exclude": [
              [(1.8, 2.65)], [(1.8, 2.65)],[(1.8, 2.65)],[(1.8, 2.65)]
        ],
        "metric": "RMSE%",
        "min_alpha": 1e-4,
    }  
# ---- registry lookup -----------------------------------------------

EXPERIMENTS = {
    "beam1": exp_beam1_b10,
    "beam3": exp_beam3_centered_mixed,
    "beam2": exp_beam2,
    "beam1_2": exp_beam1_beam2,
    "beam1_2S": exp_beam1_beam2_samples
}

def load_experiment(name: str, K_cases: int) -> Dict[str, Any]:
    if name not in EXPERIMENTS:
        raise KeyError(f"Unknown experiment '{name}'. Available: {list(EXPERIMENTS)}")
    cfg = EXPERIMENTS[name]()
    # normalize per-case arrays to length K_cases
    for key in ("weights", "shifts_mm", "exclude"):
        if key in cfg and cfg[key] is not None:
            cfg[key] = repeat_or_list(cfg[key], K_cases)
    return cfg

