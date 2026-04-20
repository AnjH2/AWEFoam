"""
Prepare batch shell scripts from templates, submit them to SLURM, and return job IDs.

What this script does
---------------------
This module generates runnable shell scripts for a batch of simulation cases,
replaces placeholder tokens inside template `.orig` files, registers those
scripts in a shared submission list, submits them with `sbatch`, and returns the
resulting SLURM job IDs.

It is used by `run_experiment.py` to launch one simulation batch for each
parameter vector.

Main responsibilities
---------------------
- resolve which template file to use for each case
- fill placeholders in the chosen template
- write generated `batchRun*.sh` files into `sBatchFolder/`
- append new entries to `sBatchFolder/sbatchToRun.txt`
- submit pending scripts from the shared queue file
- return submitted job IDs

Main function
-------------
    batch_PrepAndRun(CaseList, iteration, varList, templates=None)

Parameters:
- `CaseList`
    List of case definitions, for example temperature/current pairs
- `iteration`
    Batch index used in naming generated folders and scripts
- `varList`
    Parameter vector in the current workflow:
    `[SG, SC, DO, C0Ne, C0Pe, CDNe, CDPe]`
- `templates`
    Optional template selector. Supported forms:
    - one string
    - tuple `(default, special)`
    - list of strings of length `K`
    - callable `(i, case, iteration, x) -> template path`

Outputs
-------
Returns:
- list of submitted SLURM job IDs as strings

Side effects:
- writes new shell scripts into `sBatchFolder/`
- modifies the shared queue file `sBatchFolder/sbatchToRun.txt`
- submits jobs via `sbatch`

Placeholder replacement
-----------------------
The current script replaces placeholders such as:
- `TEMPERATURE`
- `_JGOAL_`
- `SGNE_`, `SGPE_`
- `SNE_`, `SPE_`
- `DO_`
- `CC_NE_`, `CC_PE_`
- `CD_NE_`, `CD_PE_`
- `HR`
- `CName`
- `NPRO`

Important notes
---------------
- This module edits shared files and is not safe for parallel calls unless the
  caller protects it with a lock.
- Template resolution is handled by `_resolve_template(...)`.
- Submission state is tracked through a hard-coded queue file path.
- The current implementation assumes a specific project directory layout.

Assumptions
-----------
- All required template placeholders exist in the template files.
- `sbatch` is available and returns job IDs in its standard output.
- The queue file exists and follows the expected `path:flag` format.
"""
import numpy as np
import sys
import os
import subprocess

def run_command(command, cwd=None):
    print(f"Running command: {command}")
    result = subprocess.run(command, shell=True, cwd=cwd,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running command: {command}\n{result.stderr}")
    return result.stdout + result.stderr

def sed(find, replaceWith, filePath, newFilePath=None):
    if newFilePath is not None:
        return f"sed -e 's|{find}|{replaceWith}|g' {filePath} > {newFilePath}"
    else:
        return f"sed -i 's|{find}|{replaceWith}|g' {filePath}"

# ---------------- NEW: flexible template resolution ---------------- #

def _resolve_template(orig_root: str, i: int, K: int, case, iteration, x, templates):
    """
    Decide which .orig template to use for case index i.
    templates can be:
      - str: same template for all cases
      - tuple (default, special): special used only for last case (i == K-1)
      - list[str] length K: per-case paths
      - callable fn(i, case, iteration, x) -> template path (relative or absolute)
    Returns absolute path.
    """
    # default (backward compatible) behavior if None:
    if templates is None:
        # old hardcoded rule: last case uses batchRunSemi.orig
        rel = "sBatchFolder/batchRunSemi.orig" if i == (K - 1) else "sBatchFolder/batchRun.orig"
        return os.path.join(orig_root, rel)

    # single string -> use for all
    if isinstance(templates, str):
        rel = templates
        return rel if os.path.isabs(rel) else os.path.join(orig_root, rel)

    # tuple -> keep previous "last-case special"
    if isinstance(templates, tuple) and len(templates) == 2:
        default_t, special_t = templates
        rel = special_t if i == (K - 1) else default_t
        return rel if os.path.isabs(rel) else os.path.join(orig_root, rel)

    # list per case
    if isinstance(templates, (list, tuple)) and len(templates) == K and all(isinstance(t, str) for t in templates):
        rel = templates[i]
        return rel if os.path.isabs(rel) else os.path.join(orig_root, rel)

    # callable
    if callable(templates):
        rel = templates(i, case, iteration, x)
        return rel if os.path.isabs(rel) else os.path.join(orig_root, rel)

    raise ValueError(
        "Unsupported 'templates' argument. "
        "Use: str | (str,str) | list[str] len=K | callable(i,case,iteration,x) -> str"
    )

# ------------------------------------------------------------------ #

def batch_PrepAndRun(CaseList, iteration, varList, templates=None):
    """
    Prepare and submit a batch of runs.
    - CaseList: list of cases
    - iteration: batch index (used in names)
    - varList: your x params [aL, bP, cL, dP, CC_Ne, CD_Ne, CC_Pe, CD_Pe]
    - templates: (optional) template selection (see _resolve_template)
    """
    orig = os.path.abspath(os.getcwd())
    NPro = 2
    TotalTime = 4
    K = len(CaseList)       # number of cases
    NPara = K #+ 1           # number of parallel runs started (kept as in your code) # it submits 1 to many jobs1 before

    FName = 'OBatch' + str(iteration) + '_'

    print('Number of cases to run ' + str(K) + '\n')

    with open(orig + '/sBatchFolder/sbatchToRun.txt', 'r') as file:
        data = file.readlines()
    j0 = len(data)

    for i in range(K):
        print('Preparing case ' + str(j0 + i) + ' with ' + str(CaseList[i]))

        # --- pick template for this case --- (NEW)
        BatchOrig = _resolve_template(orig, i, K, CaseList[i], iteration, varList, templates)

        BatchDict = orig + "/sBatchFolder/batchRun" + str(j0 + i) + ".sh"
        os.system('echo ' + BatchDict)
        os.system('echo ' + str(CaseList[i][0] + 273.15))

        run_command(sed('TEMPERATURE', str(CaseList[i][0] + 273.15), BatchOrig, BatchDict))
        run_command(sed('_JGOAL_', str(CaseList[i][1]), BatchDict))
        run_command(sed('SGNE_',    str(varList[0]), BatchDict))
        run_command(sed('SGPE_',    str(varList[0]), BatchDict))
        run_command(sed('SNE_',     str(varList[1]), BatchDict))
        run_command(sed('SPE_',     str(varList[1]), BatchDict))
        run_command(sed('DO_',      str(varList[2]), BatchDict))
        run_command(sed('CC_NE_',   str(varList[3]), BatchDict))
        run_command(sed('CC_PE_',   str(varList[4]), BatchDict))
        run_command(sed('CD_NE_',   str(varList[5]), BatchDict))
        run_command(sed('CD_PE_',   str(varList[6]), BatchDict))
        run_command(sed('HR', str(TotalTime), BatchDict))
        run_command(sed('CName', FName + str(i), BatchDict))
        run_command(sed('NPRO', str(NPro), BatchDict))

        with open(orig + '/sBatchFolder/sbatchToRun.txt', 'a') as file:
            file.write(BatchDict + ':' + str(0) + "\n")

    FOLDER = os.path.expanduser(
        '/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/lastStride/sBatchFolder/sbatchToRun.txt'
    )
    with open(FOLDER, 'r') as file:
        data = file.readlines()
    jobs = []
    for h in range(NPara):
        CASE = "NOCASE"
        j = 0
        for i in data:
            if (i.split(':')[1] == "0\n"):
                CASE = i.split(':')[0]
                break
            j = j + 1

        if CASE != "NOCASE":
            data[j] = CASE + ':1\n'
            with open(FOLDER, 'w') as file:
                file.writelines(data)
            jobInfo = run_command('sbatch ' + CASE)
            jobs.append(jobInfo)
            jobs = [line.split()[-1].strip() for line in jobs]
    os.system('echo jobs runing:')
    print(jobs)
    return jobs

