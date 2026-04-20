#!/usr/bin/env python3
"""
Fix previously exported line-scan .txt files by reversing the value column
for specified datasets, keeping x in increasing order.
"""

import os
import glob
import numpy as np

# ─── USER SETTINGS ─────────────────────────────────────────────────────────────
folder_path       = "."                          # root dir containing save_dir
save_dir          = "extractedData_lineX_b10"    # where the .txt files live
reversed_list_path= "reversed_files.txt"         # same list used above
# ------------------------------------------------------------------------------

def load_reversed_set(list_path):
    if list_path is None or not os.path.isfile(list_path):
        return set()
    names = set()
    with open(list_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            names.add(os.path.splitext(os.path.basename(s))[0])
    return names

REVERSED = load_reversed_set(reversed_list_path)

txt_dir = os.path.join(folder_path, save_dir)
all_txt = glob.glob(os.path.join(txt_dir, "*.txt"))

# Build a map: dataset base -> list of its slice files
by_base = {}
for fp in all_txt:
    name = os.path.basename(fp)
    # Expect pattern like "{base1}_{y}.txt" — base may include underscores
    if "_" not in name:
        continue
    base_part = name.rsplit("_", 1)[0]  # everything before the last underscore
    by_base.setdefault(base_part, []).append(fp)

num_changed = 0
for base_part, files in by_base.items():
    # match against reversed basenames (allow prefix truncation: we match if any reversed name is a prefix of base_part)
    should_fix = any(base_part.startswith(rname) for rname in REVERSED)
    if not should_fix:
        continue

    for fp in files:
        with open(fp, "r") as f:
            lines = [ln.rstrip("\n") for ln in f]
        if not lines:
            continue

        header = lines[0]
        data = [ln for ln in lines[1:] if ln.strip()]

        # parse columns
        xs, vs = [], []
        for ln in data:
            parts = ln.split()
            if len(parts) < 2:
                continue
            xs.append(float(parts[0]))
            vs.append(float(parts[1]))

        if len(vs) == 0:
            continue

        vs = vs[::-1]  # reverse values only; xs remain increasing
        with open(fp, "w") as f:
            f.write(header + "\n")
            for x, v in zip(xs, vs):
                f.write(f"{x:.6g}\t{v:.6g}\n")

        num_changed += 1
        print(f"fixed: {fp}")

print(f"Done. Files rewritten: {num_changed}")

