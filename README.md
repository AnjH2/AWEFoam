# AlkaWEFoam — Paper 1 Version

This branch contains the version of **AlkaWEFoam** associated with the first publication using this solver. The branch is maintained to improve usability, tutorials, documentation, and reproducibility for the original model version.

## Overview

AlkaWEFoam is an OpenFOAM-based solver developed for simulations of alkaline water electrolysis. The solver includes two-phase mixture momentum transport using a drift-flux closure model, potential distribution in porous electrodes and electrolyte, species transport, and interfacial mass transfer.

This branch represents the first published version of the solver. Later development versions may differ in model formulation, implementation details, available options, or tutorial structure.

The solver version in this branch is associated with the following publication:

**Multiphysics simulations of alkaline water electrolyzer cells — A sensitivity study on the effect of two-phase flow modeling**
2025
DOI: `10.1016/j.electacta.2025.147148`

## OpenFOAM version

This version is intended for:

```text
OpenFOAM v2206
```

The solver has also been compiled with OpenFOAM v2506. However, the publication results were obtained using OpenFOAM v2206.

## Solver location

The main solver is located in:

```text
applications/solvers/2206/alkaWEFoamP1
```

## Tutorials

Tutorial cases are located in:

```text
tutorials/A.Jacobsen_2025
```

Each tutorial contains the necessary OpenFOAM case files and is intended to demonstrate a specific solver capability or model setup.

## Compilation

Source the OpenFOAM environment first, then compile the solver and required libraries.

Example:

```bash
source /path/to/OpenFOAM-v2206/etc/bashrc
cd applications/solvers/2206/alkaWEFoam
./Allwmake
```

## Running a tutorial

Enter one of the tutorial case directories, then run the case using the provided script.

Example:

```bash
cd tutorials/A.Jacobsen_2025/<tutorial-case>
./allRun
```

To clean a tutorial case:

```bash
./allClean
```

Note that Linux file names are case-sensitive. If the scripts in a tutorial are named `Allrun` or `Allclean`, use those names instead.

## Reproducibility

This branch may receive documentation updates, tutorial improvements, and minor usability fixes.

The `paper1-version` branch is intended as a maintained and user-friendly version of the original solver line. For strict reproducibility, use the release or tag corresponding to the exact code version used for the publication, if available.

