# Tutorial: <case name>

## Purpose

Briefly describe what this tutorial demonstrates.


> This case returns the the point BB in the 2026 paper by A.Jacobsen.

---

## Case overview

| Item | Description |
|---|---|
| Geometry | `<2D with channel>` |
| Operating mode | `<galvanostatic / PID controled potential >` |
| Temperature | `<364K>` |
| Pressure | `<1atm>` |
| Electrolyte | `<30wt% KOH>` |
| Mesh | `<blockmesh` |

---

## Models used

List only the important run-time selections for this case.

```text
Free Relative velocity:   Schillings
    Phase diffusion:     Schillings
Pore Relative velocity: stuckBubbles
    Base Relative velocity: stokes
    Base diffusion:     solidPressure
Coverage:            linearVerticalVogt
Viscosity:           Beckermann
```



## Running the tutorial

### Serial

```bash
./Allrun
```



