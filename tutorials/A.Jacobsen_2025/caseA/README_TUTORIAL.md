# Tutorial: <case name>

## Purpose

Briefly describe what this tutorial demonstrates.


> this case is for BASELINE case A of the 2025 paper.
New relative velocity models and more, have been implemented thus the results might not be replicable. 
please note that not all initial fields are written, 
    they have to be copied from 0/ to the startTime

only a single potential is solved for, how to do a sweep is shown in 0.orig/UNe

---

## Case overview

| Item | Description |
|---|---|
| Geometry | `<2D without channel>` |
| Operating mode | `<potentiostatic >` |
| Temperature | `<364.5>` |
| Pressure | `<1atm>` |
| Electrolyte | `<30 wt% KOH>` |
| Mesh | `<blockmesh` |

---

## Models used

List only the important run-time selections for this case.

```text
Free Relative velocity:   Schillings
    Phase diffusion:      Schillings
Pore Relative velocity: stokes
    Phase diffusion:    constant
Coverage:            linearVogt
Viscosity:           Beckermann
```



## Running the tutorial

### Serial

```bash
./Allrun
```



