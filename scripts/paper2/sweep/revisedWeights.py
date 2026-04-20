#!/usr/bin/env python3
"""
Recompute weighted overall RMSE% for previously completed runs using new case
weights, write a consolidated ranking log, and generate a timeline plot.

What this script does
---------------------
This script is a post-processing utility for runs that have already been
evaluated and logged by the main optimization workflow.

It scans run folders under a figures directory, reads each run's
`*_breakdown.json`, recomputes the overall weighted RMSE% using a new set of
case weights, and writes the results to a consolidated text log.

It can also generate a timeline plot showing:
- per-case RMSE% over run index
- recomputed weighted overall RMSE% over run index
- number of missing cases per run

Typical use case
----------------
Use this script after an optimization campaign when you want to change the case
weighting without rerunning simulations.

Inputs
------
Command-line arguments:
- `--locLog`
    Root folder containing run subfolders, default `./figures`
- `--weights`
    Required comma-separated case weights
- `--out`
    Consolidated output log path
- `--overwrite`
    Overwrite the output log and rewrite its header
- `--plot-out`
    Output path stem for timeline plots
- `--show`
    Show the plot interactively after saving

How run folders are detected
----------------------------
A folder is treated as a valid run directory if it contains at least one file
matching:

    *_breakdown.json

For each such folder, the script:
- loads the breakdown JSON
- computes mean RMSE% per case across valid heights
- identifies missing cases
- recomputes the weighted overall RMSE%
- optionally writes one consolidated output line if the run is complete enough

How parameters are recovered
----------------------------
The script tries to recover the parameter vector from:
1. the run folder name
2. `OVERALL.txt`
3. the run tag

The expected parameter-label order is:

    ["SG", "SC", "DO", "C0Ne", "C0Pe", "CDNe", "CDPe"]

Outputs
-------
Main outputs:
- revised consolidated log, default:
      `./revised_optiLog.txt`
- timeline plots, default:
      `./revised_optiPlot.svg`
      `./revised_optiPlot.pdf`

The consolidated log contains:
- one header with the chosen weights
- one line per accepted run:
      `SG; SC; DO; C0Ne; C0Pe; CDNe; CDPe; cost`

Only runs with complete parameters, no missing cases, and finite recomputed cost
are appended to the consolidated log.

Important notes
---------------
- This script does not rerun simulations.
- It uses only existing per-run JSON outputs.
- Per-case values are recomputed as the mean RMSE% across valid heights.
- Runs with missing cases are still included in the plot, but may be excluded
  from the consolidated log.
- Legacy suffixes such as `; inc_before=...` are stripped if found in the output log.

How to use
----------
Example:

    python revisedWeights.py --weights 1,2 --locLog ./figures --out ./revised_optiLog.txt

Example with plot display:

    python revisedWeights.py --weights 1,1 --plot-out ./revised_optiPlot --show

Assumptions
-----------
- Run folders contain valid `*_breakdown.json` files written by `lineError.py`.
- The case order in the JSON matches the intended weight order.
- Parameter labels in folder names follow the expected naming convention.

Limitations
-----------
- The script currently recomputes only weighted RMSE%, not arbitrary metrics.
- Parameter recovery depends on naming consistency in the run folder structure.
- Some CLI argument names still reflect older naming conventions.
"""

import os
import json
import argparse
import math
from statistics import mean
import matplotlib.pyplot as plt
import re

FIGURES_ROOT_DEFAULT = "./figures"
OUTPUT_LOG_DEFAULT   = "./revised_optiLog.txt"
PLOT_OUT_DEFAULT     = "./revised_optiPlot"

# The label order must match how you build x in your runs
#X_LABELS = ["aL","bP","cL","dP","CC","CD","R_Pe","R_Ne","nSat","ds","KR","KEL"]
X_LABELS = ["SG","SC","DO","C0Ne","C0Pe","CDNe","CDPe"]

# ------------------------
# CLI
# ------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Recompute overall RMSE%% with new case weights for all runs and plot timeline.")
    ap.add_argument("--locLog", default=FIGURES_ROOT_DEFAULT,
                    help="Root folder containing run subfolders (default: ./figures)")
    ap.add_argument("--weights", required=True,
                    help="Comma-separated list of case weights, e.g. '1,1,2,1,1'")
    ap.add_argument("--out", default=OUTPUT_LOG_DEFAULT,
                    help="Path to consolidated output log (default: ./revised_optiLog.txt)")
    ap.add_argument("--overwrite", action="store_true",
                    help="Overwrite the output log (write header again).")
    ap.add_argument("--plot-out", default=PLOT_OUT_DEFAULT,
                    help="Path stem for the plot (without extension). Saves .svg and .pdf.")
    ap.add_argument("--show", action="store_true",
                    help="Show the plot window after saving.")
    return ap.parse_args()

# ------------------------
# Utilities
# ------------------------
def _strip_trailing_inc_before(s: str) -> str:
    """
    Remove any trailing '; inc_before=...' (or variations of spacing/case).
    We only strip that suffix; everything before stays intact.
    """
    # case-insensitive, strip from the *last* '; inc_before=...' to end
    import re
    return re.sub(r";\s*inc_before\s*=\s*[^;\n\r#]*\s*$", "", s, flags=re.IGNORECASE)
    
def parse_params_from_folder(run_dir, labels=X_LABELS):
    """
    Parse parameter values from the run folder name without splitting on '_'.
    Works with labels that contain underscores (e.g., 'R_Pe', 'R_Ne').
    Accepts both plain decimal values and compact encoding (m -> '-', p -> '.').
    Returns a list of values in the same order as `labels`, or None if any missing.
    """
    name = os.path.basename(os.path.normpath(run_dir))
    lab2val = {}

    # Allowed characters after label: digits, optional sign (via 'm' or '-'),
    # decimal (via 'p' or '.'), possibly more digits.
    # We anchor matches to underscores or string boundaries to avoid partial matches.
    for lab in labels:
        # Look for: (^|_)<lab><value>(_|$)
        # Value = one or more of [0-9mp.-]
        pattern = rf'(?:^|_){re.escape(lab)}([0-9mp\.\-]+)(?:_|$)'
        m = re.search(pattern, name)
        if not m:
            lab2val[lab] = None
            continue
        raw = m.group(1)

        # Try plain float first
        try:
            lab2val[lab] = float(raw)
            continue
        except Exception:
            pass

        # Fall back to compact encoding: m -> '-', p -> '.'
        dec = raw.replace('m', '-').replace('p', '.')
        try:
            lab2val[lab] = float(dec)
        except Exception:
            lab2val[lab] = None

    vals = [lab2val.get(lab, None) for lab in labels]
    if any(v is None for v in vals):
        return None
    return vals
def parse_weights(s):
    try:
        w = [float(tok.strip()) for tok in s.split(",") if tok.strip() != ""]
        if not w:
            raise ValueError("no weights provided")
        return w
    except Exception as e:
        raise SystemExit(f"Failed to parse --weights: {e}")

def is_finite(x):
    try:
        return x is not None and math.isfinite(float(x))
    except:
        return False

def find_run_dirs(root):
    if not os.path.isdir(root):
        return []
    run_dirs = []
    for name in sorted(os.listdir(root)):
        full = os.path.join(root, name)
        if not os.path.isdir(full):
            continue
        has_json = any(fn.endswith("_breakdown.json") for fn in os.listdir(full))
        if has_json:
            run_dirs.append(full)
    return run_dirs

def load_json_from_run_dir(run_dir):
    json_files = sorted([fn for fn in os.listdir(run_dir) if fn.endswith("_breakdown.json")])
    if not json_files:
        raise FileNotFoundError(f"No *_breakdown.json in {run_dir}")
    path = os.path.join(run_dir, json_files[0])
    with open(path, "r") as f:
        return json.load(f), os.path.basename(path).replace("_breakdown.json","")

def parse_params_from_overall(run_dir):
    """Load param vector x from OVERALL.txt if present."""
    overall_path = os.path.join(run_dir, "OVERALL.txt")
    if not os.path.exists(overall_path):
        return None
    try:
        with open(overall_path, "r") as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("Params (x)"):
                after = line.split(":",1)[1]
                pairs = [p.strip() for p in after.split(",")]
                lab2val = {}
                for p in pairs:
                    if "=" in p:
                        lab, val = p.split("=",1)
                        lab2val[lab.strip()] = float(val.strip())
                xs = [lab2val.get(lab, None) for lab in X_LABELS]
                return xs if not any(v is None for v in xs) else None
    except Exception:
        return None
    return None

def parse_params_from_run_tag(run_tag):
    segs = run_tag.split("_")
    lab2str = {}
    for seg in segs:
        for lab in sorted(X_LABELS, key=len, reverse=True):
            if seg.startswith(lab):
                lab2str[lab] = seg[len(lab):]
                break
    vals = []
    for lab in X_LABELS:
        sval = lab2str.get(lab, None)
        if not sval:
            vals.append(None)
            continue
        sval = sval.replace("m","-").replace("p",".")
        try:
            vals.append(float(sval))
        except:
            vals.append(None)
    return vals

def per_case_rmsep_list(data):
    """Return list per case: mean RMSE% across that case's valid heights (or None)."""
    per_case = data.get("per_case", [])
    out = []
    for case in per_case:
        vals = []
        for hrec in case.get("per_height", []):
            if hrec.get("status") == "ok":
                m = (hrec.get("metrics") or {})
                v = m.get("RMSE%")
                if is_finite(v):
                    vals.append(float(v))
        out.append(mean(vals) if vals else None)
    return out

def missing_case_indices(case_vals):
    """Indices of cases with no usable value."""
    return [i for i, v in enumerate(case_vals) if not is_finite(v)]

def recompute_rmsep_cost_from_case_vals(case_vals, weights):
    """Given per-case values and weights, compute weighted overall RMSE%."""
    mask = [v is not None and is_finite(v) for v in case_vals]
    vals_ok = [case_vals[i] for i, ok in enumerate(mask) if ok]
    if not vals_ok:
        return float("nan")
    if len(weights) != len(case_vals):
        w = [1.0]*len(case_vals)
    else:
        w = list(weights)
    w_ok = [w[i] for i, ok in enumerate(mask) if ok]
    if not any(wi > 0 for wi in w_ok):
        w_ok = [1.0]*len(vals_ok)
    denom = sum(w_ok)
    return sum(v*wi for v,wi in zip(vals_ok, w_ok)) / denom

def write_header_once(out_path, weights, overwrite=False):
    if overwrite or not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
        with open(out_path, "w") as f:
            f.write("# weights; " + "; ".join(str(w) for w in weights) + "\n")
            f.write("# columns: " + "; ".join(X_LABELS) + "; cost(RMSE%)\n")
            f.write("# a comment line will follow each entry: '# missing_cases: [...]'\n")

def append_result_line(out_path, x_vals, cost):
    xs = [(f"{v:g}" if v is not None else "") for v in (x_vals or [])]
    xs = (xs + [""]*len(X_LABELS))[:len(X_LABELS)]
    line = "; ".join(xs) + f"; {cost:.12g}\n"
    # extra safety: ensure no legacy 'inc_before' sneaks in
    line = _strip_trailing_inc_before(line)
    with open(out_path, "a") as f:
        f.write(line)
    return line.strip()
def clean_existing_outfile(out_path):
    if not os.path.exists(out_path):
        return
    with open(out_path, "r") as f:
        lines = f.readlines()
    cleaned = [_strip_trailing_inc_before(L) for L in lines]
    # only rewrite if something changed
    if cleaned != lines:
        with open(out_path, "w") as f:
            f.writelines(cleaned)
        print(f"[CLEAN] Stripped 'inc_before' tails in {out_path}")    
def params_complete(x_vals):
    """True if we have all labels and none are None."""
    return (
        x_vals is not None
        and len(x_vals) >= len(X_LABELS)
        and all(v is not None for v in x_vals[:len(X_LABELS)])
    )

def has_no_missing_cases(miss_idxs):
    """True if no cases were missing for this run."""
    return not miss_idxs  # empty list
# ------------------------
# Plotting
# ------------------------

def plot_timeline(per_run_case_vals, per_run_weighted, per_run_missing_counts, weights, out_stem, show=False):
    """
    per_run_case_vals: list of lists (runs x cases) values or None
    per_run_weighted:  list of floats (runs)
    per_run_missing_counts: list of ints (runs)
    """
    import numpy as np
    R = len(per_run_case_vals)
    if R == 0:
        print("[WARN] No runs to plot.")
        return

    # Determine number of cases (max length seen)
    K = max((len(v) for v in per_run_case_vals), default=0)
    xs = np.arange(R)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, gridspec_kw={"height_ratios":[3,1]})

    # --- Top: per-case + weighted overall
    for k in range(K):
        ys = []
        for r in range(R):
            vals = per_run_case_vals[r]
            y = vals[k] if k < len(vals) else None
            ys.append(np.nan if (y is None or not is_finite(y)) else float(y))
        ax1.plot(xs, ys, linewidth=1.5, marker="", label=f"Case {k}")

    ysW = [np.nan if not is_finite(v) else float(v) for v in per_run_weighted]
    ax1.plot(xs, ysW, linewidth=2.5, linestyle="--", marker="o", label="Weighted overall")

    ax1.set_ylabel("RMSE (%)")
    ax1.grid(True, linestyle=":", linewidth=0.8)

    wtxt = ", ".join(f"w{k}={w:g}" for k, w in enumerate(weights))
    ax1.set_title(f"RMSE% per case and weighted overall\nweights = [{wtxt}]")
    ax1.legend(loc="best", ncol=2, fontsize=9, frameon=False)

    # --- Bottom: missing count per run
    ax2.bar(xs, per_run_missing_counts, width=0.8, alpha=0.6)
    ax2.set_xlabel("Run index")
    ax2.set_ylabel("# missing cases")
    ax2.set_ylim(bottom=0)
    ax2.grid(True, axis="y", linestyle=":", linewidth=0.8)

    fig.tight_layout()
    for ext in ("svg", "pdf"):
        fname = f"{out_stem}.{ext}"
        fig.savefig(fname, bbox_inches="tight")
        print(f"[PLOT] saved {fname}")
    if show:
        plt.show()
    plt.close(fig)

# ------------------------
# Main
# ------------------------

def main():
    args = parse_args()
    weights = parse_weights(args.weights)

    run_dirs = find_run_dirs(args.locLog)
    if not run_dirs:
        print(f"No runs found under {args.figures_root}")
        return

    write_header_once(args.out, weights, overwrite=args.overwrite)
    clean_existing_outfile(args.out)  # optional but handy

    results = []                 # (cost, run_tag, x_vals, missing_idxs)
    per_run_case_vals = []       # per run: [case0, case1, ...]
    per_run_weighted = []        # weighted overall per run
    per_run_missing_counts = []  # number of missing cases per run

    for idx, run_dir in enumerate(run_dirs):
        try:
            data, run_tag = load_json_from_run_dir(run_dir)

            # parameters (for log row)
            x_vals = parse_params_from_folder(run_dir)
            #print(x_vals)
            # per-case means (RMSE% across heights)
            case_vals = per_case_rmsep_list(data)
            miss_idxs = missing_case_indices(case_vals)
            cost = recompute_rmsep_cost_from_case_vals(case_vals, weights)

            # append consolidated log line only if perfect:
            strict_ok = params_complete(x_vals) and has_no_missing_cases(miss_idxs)
            if strict_ok and is_finite(cost):
                line = append_result_line(args.out, x_vals, cost)
                results.append((cost, run_tag, x_vals, miss_idxs))
                print(f"[OK]  {os.path.basename(run_dir)} -> {line}")
            else:
                reason_bits = []
                if not params_complete(x_vals):
                    reason_bits.append("incomplete params")
                if not has_no_missing_cases(miss_idxs):
                    reason_bits.append(f"missing cases {miss_idxs}")
                if not is_finite(cost):
                    reason_bits.append("non-finite cost")
                reason = "; ".join(reason_bits) if reason_bits else "not strict-ok"
                print(f"[SKIP] {os.path.basename(run_dir)} ({reason})")

            # still collect for plotting (optional)
            per_run_case_vals.append(case_vals)
            per_run_weighted.append(cost if is_finite(cost) else float("nan"))
            per_run_missing_counts.append(len(miss_idxs))

            if not is_finite(cost):
                print(f"[WARN] {os.path.basename(run_dir)}: no usable RMSE%%; skipped.")
                continue
            

        except Exception as e:
            print(f"[WARN] Failed to process {run_dir}: {e}")

    # Top 3 best fits
    if results:
        results.sort(key=lambda r: r[0])  # smaller is better

        def _fmt(v):
            try:
                return f"{float(v):g}"
            except Exception:
                return "NA"

        print("\n=== Top 30 Fits (Lowest RMSE%) ===")
        for rank, (cost, run_tag, x, miss_idxs) in enumerate(results[:30], 1):
            x = x or [None] * len(X_LABELS)
            x = (list(x) + [None]*len(X_LABELS))[:len(X_LABELS)]
            param_str = ", ".join(f"{lab}={_fmt(v)}" for lab, v in zip(X_LABELS, x))
            miss_txt = f" | MISS cases: {miss_idxs}" if miss_idxs else ""
            print(f"{rank}. RMSE%={cost:.4f} | {param_str} | {run_tag}{miss_txt}")
    else:
        print("\nNo valid results to summarize.")

    # Plot
    if per_run_case_vals:
        plot_timeline(per_run_case_vals, per_run_weighted, per_run_missing_counts, weights, args.plot_out, show=args.show)

    print(f"\nDone. Updated {args.out} with {len(per_run_weighted)} run(s).")

if __name__ == "__main__":
    main()

