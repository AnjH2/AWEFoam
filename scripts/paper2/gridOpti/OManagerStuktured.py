#!/usr/bin/env python3
"""
Run a parameter sweep over electrode and channel thickness, submit simulation jobs,
wait for completion, extract current-density results, and write one summary CSV.

What this script does
---------------------
This script automates a batch study for two geometric parameters:

- electrode thickness (`t_elec`)
- channel thickness (`t_chan`)

For each requested parameter combination, it:
- prepares one or more simulation cases using `batch_PrepAndRun`
- submits the jobs to SLURM
- waits until the submitted jobs for that parameter point are finished
- reads the final area-averaged `UNe` value from each case
- aggregates the case results into one value for that parameter point
- writes all results to a single CSV file
- writes a run log with job IDs and per-case values

The script supports two task-generation modes:
- `grid`: build a full 2D grid from `TE_LIST` × `TC_LIST`
- `pairs`: read explicit `(t_elec, t_chan)` pairs from a CSV file

All user input is configured inside the `CONFIG` dictionary. This script does not
use command-line arguments.

Typical use case
----------------
Use this script when you want to explore how simulation output changes with
electrode thickness and channel thickness, and you want the full run process
(preparation, submission, waiting, extraction, and CSV export) handled automatically.

Inputs
------
The script expects:

- a valid `CONFIG` dictionary
- working imports for:
  - `batch_PrepAndRun`
  - `lineError.make_run_tag_from_x`
- accessible SLURM commands such as `sbatch` and `squeue`
- simulation templates referenced by `TEMPLATES`
- a simulation root directory where `OBatch...` folders appear
- simulation output files at:

      OBatch{batch_index}_{case_idx}-{jobid}/postProcessing/planeNeAvg/0/surfaceFieldValue.dat

Key configuration entries
-------------------------
- `TASK_MODE`
    `"grid"` or `"pairs"`

- `TE_LIST`, `TC_LIST`
    Parameter grids used when `TASK_MODE == "grid"`

- `PAIRS_CSV`
    CSV file with columns `t_elec,t_chan` used when `TASK_MODE == "pairs"`

- `CASE_LIST`
    List of case definitions passed to `batch_PrepAndRun`

- `TEMPLATES`
    Template paths passed to `batch_PrepAndRun`
    Must match `CASE_LIST` length in the current setup

- `SIM_ROOT`
    Root directory where output case folders are created

- `SLURM_USERNAME`
    Username used when polling `squeue`

- `AGGREGATE`
    How per-case values are combined into one value per parameter point:
    - `"mean"`
    - `"max"`
    - `"min"`
    - integer case index, for example `0`

- `MAX_CONCURRENCY`
    Maximum number of parameter points processed at the same time

- `OUT_ROOT`
    Root directory for result folders

Outputs
-------
For each run, the script creates a timestamped or user-defined output folder:

    <OUT_ROOT>/<RUN_TAG>/

Inside that folder it writes:

- `results.csv`
    One row per parameter point with:
    - `t_elec`
    - `t_chan`
    - aggregate mode
    - aggregate value
    - one column per case result

- `stdout.log`
    Run metadata, progress messages, job IDs, and per-case values

How to use
----------
1. Edit the `CONFIG` dictionary.
2. Choose whether to use:
   - `TASK_MODE = "grid"`
   - or `TASK_MODE = "pairs"`
3. Set the template paths, case list, SLURM username, and simulation root.
4. Run the script:

       python optimize_grid_simple.py

Basic example
-------------
This example runs explicit parameter pairs from a CSV file and averages the case
results for each pair:

    CONFIG = {
        "TASK_MODE": "pairs",
        "PAIRS_CSV": "./opt_results/contours/contour_from_stdout_failed_points.csv",
        "SNAP_STEP": 0,
        "SNAP_DECIMALS": None,
        "CASE_LIST": [
            [91, 3000],
        ],
        "TEMPLATES": [
            "sBatchFolder/batchRunSemi2.orig",
        ],
        "SIM_ROOT": "/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/gridOpti",
        "SLURM_USERNAME": "anj",
        "AGGREGATE": "mean",
        "OUT_ROOT": "./opt_results",
        "POLL_INTERVAL": 60,
        "MAX_CONCURRENCY": 48,
        "RUN_TAG": None,
    }

Then run:

    python optimize_grid_simple.py

How results are read
--------------------
For each submitted case, the script reads:

    postProcessing/planeNeAvg/0/surfaceFieldValue.dat

and extracts the last numeric value from the second column, interpreted here as
the final area-averaged `UNe` value.

Important notes
---------------
- This script submits real jobs to SLURM.
- It waits for submitted jobs by repeatedly checking `squeue`.
- Parallelism is applied at the parameter-point level, not inside one parameter point.
- Submission is protected by `SUBMIT_LOCK` because `batch_PrepAndRun` edits shared
  files in `sBatchFolder`.
- In `pairs` mode, duplicate `(t_elec, t_chan)` pairs are removed.
- In `grid` mode, `batch_index` is assigned in row-major order over `(TE_LIST, TC_LIST)`.
- If a result file is missing or unreadable, that case value becomes `None`.
- Aggregation ignores missing and non-finite values for `"mean"`, `"max"`, and `"min"`.
- If an integer case index is used for `AGGREGATE`, that exact per-case value is returned.

Assumptions
-----------
- `batch_PrepAndRun.batch_PrepAndRun(...)` returns job IDs in the same order as `CASE_LIST`.
- The expected `OBatch...` folder naming scheme is used by the submitted jobs.
- `surfaceFieldValue.dat` contains numeric data in the expected format.
- The user has permission to submit and query SLURM jobs.

Limitations
-----------
- No command-line interface is provided.
- No retry logic is implemented for failed submissions or missing outputs.
- No plotting is performed in this script.
- Error handling is minimal; failed futures are logged, but not automatically rerun.
"""

import os
import csv
import time
import subprocess
from typing import List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

import numpy as np

# ====== YOUR MODULES ======
import batch_PrepAndRun
from lineError import make_run_tag_from_x

# ==============================
# ===== OVERVIEW / CONFIG ======
# ==============================
CONFIG = {
    # --- Grid (used when TASK_MODE="grid") ---
    "TE_LIST": list(np.arange(0.5, 2.5001, 0.1)),   # electrode thickness grid [mm]
    "TC_LIST": list(np.arange(0.5, 2.5001, 0.1)),   # channel thickness grid [mm]

    # --- Alternate mode: run explicit pairs from a CSV (t_elec,t_chan) ---
    # Set TASK_MODE to "pairs" and point PAIRS_CSV to your file.
    "TASK_MODE": "pairs",                 # "grid" or "pairs"
    "PAIRS_CSV": "./opt_results/contours/contour_from_stdout_failed_points.csv",          # only used if TASK_MODE == "pairs"

    # Optional snapping (applies to PAIRS_CSV rows)
    "SNAP_STEP": 0,                    # if >0, snap x -> round(x/step)*step
    "SNAP_DECIMALS": None,               # else round(x, decimals). Ignored if SNAP_STEP is set.

    # --- Cases (each entry whatever your templates expect) ---
    "CASE_LIST": [
        [91, 3000]
    ],

    # --- Template for each case (same length as CASE_LIST) ---
    "TEMPLATES": [
        "sBatchFolder/batchRunSemi2.orig"
    ],

    # --- Where OBatch* folders appear (so we can read results) ---
    "SIM_ROOT": "/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/gridOpti",

    # --- SLURM user (for squeue polling) ---
    "SLURM_USERNAME": "anj",

    # --- How to aggregate per-case current into a single value per (te, tc):
    #       "mean", "max", "min", or a specific case index (int, 0-based)
    "AGGREGATE": "mean",

    # --- Output root; a run folder will be created inside this with a timestamp tag ---
    "OUT_ROOT": "./opt_results",

    # --- Polling interval (seconds) while waiting for jobs to finish (per grid point) ---
    "POLL_INTERVAL": 60,

    # --- Parallelism: how many grid points to run at the same time ---
    "MAX_CONCURRENCY": 48,

    # --- Optional static run tag (None -> auto timestamp) ---
    "RUN_TAG": None,
}
# ==============================


# ------------- Helpers -------------
def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path

def run_cmd(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

def check_any_jobs_running(username: str, job_ids: List[str]) -> bool:
    res = run_cmd(["squeue", "-u", username, "-o", "%A"])
    running = set(res.stdout.splitlines())
    return any(j in running for j in job_ids)

def wait_for_jobs(username: str, job_ids: List[str], interval: int):
    while True:
        if not check_any_jobs_running(username, job_ids):
            break
        time.sleep(interval)

def read_planeNeAvg_UNe(file_path: str) -> Optional[float]:
    """
    Parse planeNeAvg/0/surfaceFieldValue.dat and return the last numeric
    areaAverage(UNe) value (2nd column) from the last data row (ignores comments).
    """
    if not os.path.exists(file_path):
        return None
    last_val = None
    try:
        with open(file_path, "r") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#") or s.startswith("//"):
                    continue
                toks = s.split()
                if len(toks) < 2:
                    continue
                try:
                    val = float(toks[1])  # column 1 = areaAverage(UNe)
                except Exception:
                    continue
                if np.isfinite(val):
                    last_val = val
        return last_val
    except Exception:
        return None

def aggregate(values: List[Optional[float]], mode) -> Optional[float]:
    vals = [v for v in values if v is not None and np.isfinite(v)]
    if isinstance(mode, str):
        if mode == "mean":
            return float(np.mean(vals)) if vals else None
        if mode == "max":
            return float(np.max(vals)) if vals else None
        if mode == "min":
            return float(np.min(vals)) if vals else None
        try:
            mode = int(mode)
        except Exception:
            return None
    if isinstance(mode, int) and 0 <= mode < len(values):
        return values[mode]
    return None

def _snap_value(x: float, step=None, decimals=None) -> float:
    if step and step > 0:
        return round(x / step) * step
    if decimals is not None:
        return round(x, int(decimals))
    return float(f"{x:.12g}")

# -------- Task builders --------
SUBMIT_LOCK = threading.Lock()   # serialize sbatch file edits inside batch_PrepAndRun

def build_tasks_grid(te_list: List[float], tc_list: List[float]) -> List[Tuple[int, int, int, float, float, str]]:
    """
    Returns tasks:
      (batch_index, ii, jj, te, tc, run_tag)
    batch_index = 0..N-1, row-major over (te, tc).
    """
    tasks = []
    batch_index = 0
    for ii, te in enumerate(te_list):
        for jj, tc in enumerate(tc_list):
            x = [te, tc]
            run_tag = make_run_tag_from_x(x)
            tasks.append((batch_index, ii, jj, te, tc, run_tag))
            batch_index += 1
    return tasks

def load_pairs_csv(path: str, snap_step=None, snap_decimals=None) -> List[Tuple[float, float]]:
    """
    Read a CSV with columns: t_elec,t_chan -> list of (te, tc) unique pairs.
    Header required; extra columns ignored.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"PAIRS_CSV not found: {path}")

    pairs = []
    with open(path, "r", newline="") as cf:
        rdr = csv.DictReader(cf)
        if "t_elec" not in rdr.fieldnames or "t_chan" not in rdr.fieldnames:
            raise ValueError(f"CSV must contain columns 't_elec' and 't_chan' (got {rdr.fieldnames})")
        for row in rdr:
            try:
                te = float(row["t_elec"])
                tc = float(row["t_chan"])
            except Exception:
                continue
            te = _snap_value(te, snap_step, snap_decimals)
            tc = _snap_value(tc, snap_step, snap_decimals)
            pairs.append((te, tc))

    # de-duplicate and stable-sort by (tc, te)
    pairs = sorted(set(pairs), key=lambda t: (t[1], t[0]))
    return pairs

def build_tasks_pairs(pairs: List[Tuple[float, float]]) -> List[Tuple[int, int, int, float, float, str]]:
    """
    Build tasks from explicit (te, tc) pairs.
    ii/jj are just ordinal indices (0..), not tied to a full grid.
    """
    tasks = []
    for bidx, (te, tc) in enumerate(pairs):
        run_tag = make_run_tag_from_x([te, tc])
        # use ii=bidx, jj=0 (not used downstream except logging)
        tasks.append((bidx, bidx, 0, te, tc, run_tag))
    return tasks

def run_one_point(
    batch_index: int,
    ii: int,
    jj: int,
    te: float,
    tc: float,
    run_tag: str,
    *,
    case_list,
    templates,
    sim_root: str,
    user: str,
    poll_interval: int,
    aggregate_mode,
) -> dict:
    """
    Submit, wait, collect per-case currents, aggregate, return a result dict.
    """
    # --- Submit (serialize to avoid sbatch file conflicts) ---
    with SUBMIT_LOCK:
        job_ids = batch_PrepAndRun.batch_PrepAndRun(
            case_list, batch_index, [te, tc], templates=templates
        )

    # --- Wait only for these jobs ---
    wait_for_jobs(user, job_ids, interval=poll_interval)

    # --- Collect per-case values ---
    per_case_vals: List[Optional[float]] = []
    for case_idx in range(len(case_list)):
        jid = job_ids[case_idx]
        fpath = os.path.join(
            sim_root,
            f"OBatch{batch_index}_{case_idx}-{jid}",
            "postProcessing", "planeNeAvg", "0", "surfaceFieldValue.dat",
        )
        val = read_planeNeAvg_UNe(fpath)
        per_case_vals.append(val)

    agg = aggregate(per_case_vals, aggregate_mode)
    return {
        "batch_index": batch_index,
        "ii": ii,
        "jj": jj,
        "te": te,
        "tc": tc,
        "job_ids": job_ids,
        "per_case_vals": per_case_vals,
        "agg": agg,
        "run_tag": run_tag,
    }


# ------------- Main -------------
def main():
    TE = CONFIG["TE_LIST"]
    TC = CONFIG["TC_LIST"]
    CASE_LIST = CONFIG["CASE_LIST"]
    TEMPLATES = CONFIG["TEMPLATES"]
    SIM_ROOT = CONFIG["SIM_ROOT"]
    USER = CONFIG["SLURM_USERNAME"]
    MODE = CONFIG["AGGREGATE"]
    OUT_ROOT = CONFIG["OUT_ROOT"]
    POLL = CONFIG["POLL_INTERVAL"]
    MAX_CONCURRENCY = int(CONFIG.get("MAX_CONCURRENCY", 4))
    TASK_MODE = CONFIG.get("TASK_MODE", "grid")
    SNAP_STEP = CONFIG.get("SNAP_STEP")
    SNAP_DECIMALS = CONFIG.get("SNAP_DECIMALS")

    if len(TEMPLATES) != len(CASE_LIST):
        raise SystemExit(f"TEMPLATES length ({len(TEMPLATES)}) must match CASE_LIST length ({len(CASE_LIST)}).")

    auto_tag = f"Grid_{int(time.time())}"
    run_tag = CONFIG["RUN_TAG"] or auto_tag
    outdir = ensure_dir(os.path.join(OUT_ROOT, run_tag))
    log_path = os.path.join(outdir, "stdout.log")
    csv_path = os.path.join(outdir, "results.csv")

    # Build task list
    if TASK_MODE == "grid":
        tasks = build_tasks_grid(TE, TC)
        print(f"[PLAN] TASK_MODE=grid  |TE|={len(TE)}, |TC|={len(TC)}  -> planned jobs = {len(tasks)}")
    elif TASK_MODE == "pairs":
        pairs_csv = CONFIG.get("PAIRS_CSV")
        pairs = load_pairs_csv(pairs_csv, SNAP_STEP, SNAP_DECIMALS)
        tasks = build_tasks_pairs(pairs)
        print(f"[PLAN] TASK_MODE=pairs |pairs|={len(pairs)}          -> planned jobs = {len(tasks)}")
    else:
        raise SystemExit(f"Unknown TASK_MODE: {TASK_MODE} (use 'grid' or 'pairs').")

    total = len(tasks)
    if total == 0:
        print("[PLAN] No tasks produced. Check CONFIG.")
        # still create empty outputs for traceability
        ensure_dir(outdir)
        with open(csv_path, "w", newline="") as cf:
            w = csv.writer(cf)
            w.writerow(["t_elec", "t_chan", "aggregate_mode", "aggregate_value"] +
                       [f"case{i}_UNe" for i in range(len(CASE_LIST))])
        with open(log_path, "w") as lf:
            lf.write(f"RUN TAG: {run_tag}\n")
            lf.write(f"MODE: {MODE}\n")
            lf.write(f"TASK_MODE: {TASK_MODE}\n")
            lf.write("No tasks.\n")
        print("\n=== DONE ===")
        print(f"CSV   : {csv_path}")
        print(f"Log   : {log_path}")
        return

    # CSV header
    with open(csv_path, "w", newline="") as cf:
        w = csv.writer(cf)
        w.writerow(["t_elec", "t_chan", "aggregate_mode", "aggregate_value"]
                   + [f"case{i}_UNe" for i in range(len(CASE_LIST))])

    with open(log_path, "w") as lf:
        lf.write(f"RUN TAG: {run_tag}\n")
        lf.write(f"MODE: {MODE}\n")
        lf.write(f"TASK_MODE: {TASK_MODE}\n")
        lf.write(f"MAX_CONCURRENCY: {MAX_CONCURRENCY}\n\n")

    # Run in parallel
    results: List[dict] = []
    done = 0
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENCY) as ex:
        futures = []
        for (batch_index, ii, jj, te, tc, rt) in tasks:
            fut = ex.submit(
                run_one_point,
                batch_index, ii, jj, te, tc, rt,
                case_list=CONFIG["CASE_LIST"],
                templates=CONFIG["TEMPLATES"],
                sim_root=SIM_ROOT,
                user=USER,
                poll_interval=POLL,
                aggregate_mode=MODE,
            )
            futures.append(fut)

        for fut in as_completed(futures):
            try:
                res = fut.result()
                results.append(res)
                done += 1
                print(f"[{done}/{total}] te={res['te']} tc={res['tc']} -> agg={res['agg']}")
                with open(log_path, "a") as lf:
                    lf.write(f"[{done}/{total}] te={res['te']} tc={res['tc']} "
                             f"jobs={res['job_ids']} agg={res['agg']} per_case={res['per_case_vals']}\n")
            except Exception as e:
                done += 1
                print(f"[{done}/{total}] ERROR: {e}")
                with open(log_path, "a") as lf:
                    lf.write(f"[{done}/{total}] ERROR: {e}\n")

    # Write CSV in a stable order by batch_index
    results.sort(key=lambda d: d["batch_index"])
    with open(csv_path, "a", newline="") as cf:
        w = csv.writer(cf)
        for r in results:
            w.writerow([r["te"], r["tc"], MODE, "" if r["agg"] is None else r["agg"]] +
                       r["per_case_vals"])

    print("\n=== DONE ===")
    print(f"CSV   : {csv_path}")
    print(f"Log   : {log_path}")

if __name__ == "__main__":
    main()

