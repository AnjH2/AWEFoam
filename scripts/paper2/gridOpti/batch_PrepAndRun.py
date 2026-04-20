"""
Prepare and submit batches of simulation cases from template shell scripts.

What this script does
---------------------
This module builds runnable batch scripts from one or more template `.orig` files,
replaces placeholder values using `sed`, appends the generated scripts to a shared
submission list, and submits jobs with `sbatch`.

It is designed to support parameter studies where each batch corresponds to one
parameter point and each batch may contain multiple related simulation cases.

Main responsibilities
---------------------
- choose the correct template file for each case
- compute grading parameters using `CalGrading`
- replace placeholders in shell-script templates
- generate `batchRun*.sh` files
- append them to `sBatchFolder/sbatchToRun.txt`
- submit jobs to SLURM with `sbatch`
- return the submitted job IDs

Inputs
------
Main function:

    batch_PrepAndRun(CaseList, iteration, varList, templates=None)

Parameters:
- `CaseList`
    List of case definitions. Each entry is passed through the template replacement logic.
- `iteration`
    Batch index used in generated folder and case names.
- `varList`
    Parameter list. In the current use, this is:
        [t_elec, t_chan]
- `templates`
    Optional template selector. Supported forms:
    - one string for all cases
    - tuple `(default, special)`
    - list of strings, one per case
    - callable `fn(i, case, iteration, x) -> path`

Outputs
-------
Returns:
- `jobs`: list of submitted SLURM job IDs as strings

It also creates or updates:
- generated batch shell scripts in `sBatchFolder/`
- the shared file `sBatchFolder/sbatchToRun.txt`

How template selection works
----------------------------
Template selection is handled by `_resolve_template(...)`.

Supported forms:
- `None`
    Use built-in default behavior
- `"path/to/template.orig"`
    Same template for all cases
- `("default.orig", "special.orig")`
    Use special template only for the last case
- `["case0.orig", "case1.orig", ...]`
    Explicit per-case template list
- callable
    Return a template path dynamically per case

Important notes
---------------
- This module modifies shared files in `sBatchFolder`, so concurrent calls should
  be protected externally with a lock.
- The current implementation assumes a specific folder structure relative to the
  current working directory.
- Placeholder replacement is done with shell `sed` commands.
- Job submission is done with `sbatch`.
- Some path handling is hard-coded, especially the expanded path to `sbatchToRun.txt`.

Assumptions
-----------
- The required template placeholders exist in the `.orig` files.
- `CalGrading.AllFunctions` is available and returns valid grading parameters.
- `sbatch` is available on the system.
- The generated shell scripts are valid and runnable.
"""
import numpy as np
import sys
import os
import subprocess
import CalGrading as CG

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
    TotalTime = 24
    K = len(CaseList)       # number of cases
    NPara = K + 1           # number of parallel runs started (kept as in your code)
    
    tDi=0.85 ## have to be the same everywhere
    tE=varList[0];
    tCh=varList[1];
    widthA=0.018*2;
    cratio=1.2;
    all_funcs=CG.AllFunctions();
    nDi,rDi=all_funcs.calculate('calculate_nr_and_ratio_from_cratio_widthA_length', cratio, widthA/2,tDi/2)
    nE,rE=all_funcs.calculate('calculate_nr_and_ratio_from_cratio_widthA_length', cratio, widthA,tE/2)
    nCh,rCh=all_funcs.calculate('calculate_nr_and_ratio_from_cratio_widthA_length', cratio, widthA,tCh/2)
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
        run_command(sed('_TE_', str(varList[0]), BatchDict))
        run_command(sed('_TC_', str(varList[1]), BatchDict))
        
        run_command(sed('_nDi_', str(nDi*2), BatchDict))
        run_command(sed('_rDi_', str(rDi), BatchDict))
        run_command(sed('_nE_', str(nE*2), BatchDict))
        run_command(sed('_rE_', str(rE), BatchDict))
        run_command(sed('_nCh_', str(nCh*2), BatchDict))
        run_command(sed('_rCh_', str(rCh), BatchDict))
        
        run_command(sed('HR', str(TotalTime), BatchDict))
        run_command(sed('CName', FName + str(i), BatchDict))
        run_command(sed('NPRO', str(NPro), BatchDict))

        with open(orig + '/sBatchFolder/sbatchToRun.txt', 'a') as file:
            file.write(BatchDict + ':' + str(0) + "\n")

    FOLDER = os.path.expanduser(
        '/home/anj/OpenFOAM/anj-v2206/run/alkaWEFoam/gridOpti/sBatchFolder/sbatchToRun.txt'
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

