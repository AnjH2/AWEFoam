#!/usr/bin/env python3
"""
Run a parameter sweep, submit simulation batches, compare each run against a
selected experiment, and record a scalar cost for every parameter combination.

What this script does
---------------------
This script is the main driver for an experiment-based optimization workflow.

For each parameter combination in the configured grid, it:
- prepares and submits one batch of simulation jobs
- waits until the submitted SLURM jobs finish
- compares simulation line data against experimental reference data
- computes a scalar cost using a chosen metric
- writes figures and per-run logs through `lineError.curveCompare`
- appends a compact summary line to a legacy optimization log

The experiment definition is selected at runtime with `--exp` and is loaded from
`exp_registry.py`.

Typical use case
----------------
Use this script when you want to scan a multi-parameter design space and rank
simulation runs by how well they match a chosen experimental dataset.

Inputs
------
This script uses two sources of input:

1. Command-line arguments
   - `--exp`
       Experiment key from `exp_registry.EXPERIMENTS`
   - `--metric`
       Optional metric override, for example:
       `RMSE%`, `MAE%`, `RMSE`, `MAE`, `1-R2`
   - `--weights`
       Optional comma-separated per-case weights
   - `--min-alpha`
       Optional minimum alpha threshold override
   - `--figures-root`
       Output root for figures and logs

2. Hard-coded parameter lists in this file
   - `SG_List`
   - `SC_List`
   - `DO_List`
   - `C0Ne_List`
   - `C0Pe_List`
   - `CDNe_List`
   - `CDPe_List`

These are combined into the full parameter grid.

Main workflow
-------------
For each parameter vector `x = [SG, SC, DO, C0Ne, C0Pe, CDNe, CDPe]`, the script:
- creates a filesystem-safe run tag
- submits one batch using `batch_PrepAndRun.batch_PrepAndRun(...)`
- waits for the returned job IDs to complete
- calls `curveCompare(...)` from `lineError.py`
- receives a scalar cost and detailed breakdown
- writes a one-line entry to `optimization_log.txt`

Outputs
-------
Main outputs include:
- per-run figures and logs under `--figures-root/<run_tag>/`
- appended summary lines in `optimization_log.txt`
- printed terminal summary of completed batches

Per-run outputs are written by `lineError.py` and typically include:
- scatter plots
- per-case summary plots
- `OVERALL.txt`
- `*_breakdown.json`
- `per_case.csv`
- `per_height.csv`

How to use
----------
Example:

    python run_experiment.py --exp beam1 --metric RMSE% --figures-root ./figures

Optional example with custom weights:

    python run_experiment.py --exp beam2 --weights 1,2 --min-alpha 1e-4

Important notes
---------------
- Job submission is serialized with `SUBMIT_LOCK` because batch preparation edits
  shared files.
- Multiple batches can still be processed concurrently at the higher level.
- The parameter grid is defined in code, not through CLI arguments.
- Missing simulation cases are penalized inside `curveCompare(...)` using
  `missing_case_penalty=200.0`.
- The legacy text log format is preserved for compatibility with downstream scripts.

Assumptions
-----------
- SLURM commands such as `sbatch` and `squeue` are available.
- Template files and simulation folder structure match the expectations of
  `batch_PrepAndRun.py`.
- The selected experiment key exists in `exp_registry.py`.
- Simulation output paths produced by the experiment path provider are valid.

Limitations
-----------
- No resume or retry mechanism is implemented for failed batches.
- Most configuration is hard-coded in the script body.
- Cleanup of generated batch directories is currently disabled by default.
"""

import os, glob, shutil
import time
import argparse
import subprocess
from typing import List, Dict, Tuple, Optional

import matplotlib.pyplot as plt  # if you need it somewhere else
import numpy as np

import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import batch_PrepAndRun
from lineError import curveCompare, print_curve_info, make_run_tag_from_x
from exp_registry import load_experiment

# ---------------- CLI ----------------

def parse_args():
    ap = argparse.ArgumentParser(description="Run batch + compare against a chosen experiment.")
    ap.add_argument("--exp", required=True, help="Experiment key (see exp_registry.EXPERIMENTS)")
    ap.add_argument("--metric", default=None, help="Override metric (e.g., RMSE%%, MAE%%, RMSE, 1-R2)")
    ap.add_argument("--weights", default=None, help="Comma-separated per-case weights to override experiment default")
    ap.add_argument("--min-alpha", type=float, default=None, help="Override min_alpha threshold")
    ap.add_argument("--figures-root", default="./figures", help="Output root for figures/logs")
    return ap.parse_args()
    
def delete_batch_dirs(batch_index: int, base_dir: str):
    """
    Remove all directories named like:
      {base_dir}/OBatch{batch_index}_*
    """
    pattern = os.path.join(base_dir, f"OBatch{batch_index}_*")
    for d in glob.glob(pattern):
        if os.path.isdir(d):
            try:
                shutil.rmtree(d)
                print(f"[CLEAN] Removed {d}")
            except Exception as e:
                print(f"[CLEAN][ERROR] {d}: {e}")
    
def parse_weights_csv(s: str):
    return [float(tok.strip()) for tok in s.split(",") if tok.strip() != ""]

# ------------- (your existing helpers snipped for brevity) ----------

def check_jobs_running(job_ids):
    """Check if the jobs are running using squeue."""
    # Run squeue to get the list of running jobs
    result = subprocess.run(['squeue', '-u', username, '-o', '%A'], capture_output=True, text=True)
    
    # Get a list of all currently running job IDs from squeue output
    running_jobs = result.stdout.splitlines()

    # Check if each job ID in job_ids is in the list of running jobs
    running_status = {job_id: (job_id in running_jobs) for job_id in job_ids}
    
    return running_status

def check_jobs_until_complete(job_ids, check_interval=60, max_runtime=None):
    """
    Check if the jobs are running and wait until they complete.
    
    Args:
    - job_ids: List of job IDs to monitor.
    - check_interval: Time (in seconds) between checks. Default is 60 seconds.
    - max_runtime: Maximum runtime (in seconds) before timeout. Default is None (no limit).
    
    Returns:
    - True if all jobs completed within the max_runtime, False if timeout occurred.
    """
    start_time = time.time()  # Record the starting time
    
    while True:
        # Check the status of the jobs
        status = check_jobs_running(job_ids)
        
        # Print the current status of each job
        #for job_id, is_running in status.items():
       #     if is_running:
       #         print(f"Job {job_id} is still running.")
       #     else:
       #         print(f"Job {job_id} has completed or does not exist.")
        
        # If all jobs have completed, return True and break the loop
        if all(not is_running for is_running in status.values()):
            print("All jobs have completed.")
            return True
        
        # If max_runtime is specified, check if the time limit is reached
        elapsed_time = time.time() - start_time
        if max_runtime is not None and elapsed_time >= max_runtime:
            print(f"Max runtime of {max_runtime / 60:.2f} minutes reached. Exiting check.")
            return False
        
        # Wait for the specified interval before checking again
        time.sleep(check_interval)

username = "anj"
CaseList = [
    [25, 2000],
    [25, 3000],
]

log_file   = "optimization_log.txt"
LOG_LOCK   = threading.Lock()     # serialize writes to the legacy log
SUBMIT_LOCK = threading.Lock()    # serialize sbatch file edits inside batch_PrepAndRun
#SG=2.74864, SC=0.559858, DO=0, C0Ne=5.5, C0Pe=3.5, CDNe=40, CDPe=32
SG_List    = [5]
SC_List    = [0.5,1]#[0,1]
DO_List    = [0]#[1,2]#[1,2,4,8] 
C0Ne_List    = [5,5.5,6]#[3,6]#[3,6,9]
C0Pe_List    = [2,2.5,3]#[4,5]#[4,5,6]
CDNe_List    = [4,8]#[16,32,48]#[64]
CDPe_List  = [4,8]#[0.06,0.07]#[0.01,0.02,0.04]

templates_list = [
    "sBatchFolder/batchRun.orig",      # case 1
    "sBatchFolder/batchRun.orig",      # case 2
]

# ---------------- PARALLEL ORCHESTRATION ---------------------------

def build_param_grid():
    """Yield (batch_index, x_params, run_tag) for the full nested grid."""
    i = 0
    for v0 in SG_List:
        for v1 in SC_List:
            for v2 in DO_List:
                for v3 in C0Ne_List:
                    for v4 in C0Pe_List:
                        for v5 in CDNe_List:
                            for v6 in CDPe_List:
                                x = [v0, v1, v2, v3, v4, v5, v6]
                                run_tag = make_run_tag_from_x(x)
                                yield (i, x, run_tag)
                                i += 1

def run_one_batch(
    batch_index: int,
    x: List[float],
    run_tag: str,
    *,
    path_provider,
    per_case_shifts,
    per_case_exclude,
    metric_name: str,
    min_alpha: float,
    weights,
    figures_root: str,
):
    """Submit, wait, compare, log; returns (batch_index, run_tag, cost)."""

    # --- Submit (serialize this step to avoid sbatch file conflicts) ---
    with SUBMIT_LOCK:
        job_ids = batch_PrepAndRun.batch_PrepAndRun(
            CaseList, batch_index, x, templates=templates_list
        )

    # --- Wait for JUST these jobs ---
    check_jobs_until_complete(job_ids, check_interval=60)

    # --- Compare ---
    cost, info = curveCompare(
        CaseList, job_ids,
        batch_index=batch_index,
        min_alpha=min_alpha,
        case_weights=weights,
        metric_name=metric_name,
        per_case_sim_shift_mm=per_case_shifts,
        per_case_exclude_ranges=per_case_exclude,
        return_breakdown=True,
        outdir_root=figures_root,
        run_tag=run_tag,
        save_per_height_scatter=True,
        write_logs=True,
        x_params=x,
        path_provider=path_provider,
        missing_case_penalty=200.0,
    )

    # --- Print a short summary to stdout ---
    print(f"[DONE] RUN={run_tag} | COST ({metric_name}) = {cost:.6g}")
    # If you want the long, per-height breakdown:
    # print_curve_info(info)

    # --- Append to legacy optimization log (serialize) ---
    with LOG_LOCK:
        with open(log_file, "a") as f:
            f.write(f"{batch_index}; {x[0]}; {x[1]}; {x[2]}; {x[3]}; "
                    f"{x[4]}; {x[5]}; {x[6]}; {cost}\n")

    return batch_index, run_tag, cost

def main():
    args = parse_args()
    K = len(CaseList)

    # Load experiment config (paths, per-case defaults)
    exp_cfg         = load_experiment(args.exp, K_cases=K)
    path_provider   = exp_cfg["path_provider"]
    per_case_shifts = exp_cfg.get("shifts_mm", [0.0] * K)
    per_case_exclude= exp_cfg.get("exclude",   [None] * K)
    metric_name     = args.metric or exp_cfg.get("metric", "RMSE%")
    min_alpha       = exp_cfg.get("min_alpha", 1e-4) if args.min_alpha is None else args.min_alpha
    weights         = exp_cfg.get("weights", None)

    if args.weights:
        user_w = parse_weights_csv(args.weights)
        if len(user_w) == K:
            weights = user_w
        else:
            print(f"[WARN] Provided {len(user_w)} weights but there are {K} cases; using experiment defaults.")

    # Build the full list of batches we plan to run
    tasks = list(build_param_grid())
    if not tasks:
        print("No parameter combinations to run.")
        return

    # Concurrency control
    MAX_CONCURRENCY = 24  # submit/wait/evaluate up to 3 batches at a time

    # Optional: disable clean(i) during parallel runs (dangerous!)
    # You can add a final cleanup pass afterward if you need it.

    results = []
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENCY) as ex:
        futures = []
        for (batch_index, x, run_tag) in tasks:
            fut = ex.submit(
                run_one_batch,
                batch_index, x, run_tag,
                path_provider=path_provider,
                per_case_shifts=per_case_shifts,
                per_case_exclude=per_case_exclude,
                metric_name=metric_name,
                min_alpha=min_alpha,
                weights=weights,
                figures_root=args.figures_root,
            )
            futures.append(fut)

        for fut in as_completed(futures):
            try:
                results.append(fut.result())
            except Exception as e:
                print(f"[ERROR] A batch failed: {e}")

    # Summary of all completed runs
    if results:
        results.sort(key=lambda t: t[0])  # sort by batch_index
        print("\n=== Completed Batches ===")
        for idx, tag, cost in results:
            print(f"#{idx:03d}  {tag}  ->  cost={cost:.6g}")

    #delete_batch_dirs(
    #    batch_index=batch_index,
    #    base_dir="/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/channelStructuredSearch"
    #)

if __name__ == "__main__":
    main()
