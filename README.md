# alkaWEFoam

`alkaWEFoam` is an OpenFOAM-based multiphysics solver developed for numerical modelling of alkaline water electrolysis.

The solver combines two-phase flow, gas transport, porous-media effects, species transport, electrochemical reactions, gas-liquid mass transfer, electrode bubble coverage, and phase-relative transport. Several physical closures are run-time selectable so that individual models can be exchanged without modifying the main solver.

> **Important:** the custom boundary-condition library used by the cases must be compiled in addition to the solver and its supporting libraries. The boundary-condition source is maintained separately from the solver source in this branch.

---

## Main capabilities

The solver currently includes:

- incompressible gas-liquid two-phase transport;
- porous electrodes, separator/membrane, and free-flow regions;
- drift-flux / relative-velocity modelling for the gas phase;
- separate models for:
  - buoyancy- or flow-induced relative drift;
  - phase diffusion / hydrodynamic dispersion;
- species transport in the electrolyte;
- electrochemical source terms based on Butler-Volmer kinetics;
- electrolyte and solid-phase potential equations;
- gas-liquid mass transfer and water-vapour treatment;
- electrode bubble-coverage models;
- concentration-loss models;
- gas-mixture and effective mixture-viscosity models;
- run-time selectable closures for several of the above model families.

The implementation is associated with the alkaline-electrolysis modelling work documented in the accompanying publications.

---

## Solver structure

The main executable is built around the OpenFOAM PIMPLE framework. A typical time step contains the following operations:

1. update water partial pressure and mixture properties;
2. update relative-velocity and phase-diffusion models;
3. update bubble coverage;
4. solve the gas-volume-fraction equation;
5. update density, viscosity, and porous drag;
6. solve the mixture momentum and pressure equations;
7. update effective species diffusivities and electrochemical properties;
8. solve electrode/electrolyte potentials;
9. solve species transport equations;
10. write current, phase, transport, and electrochemical fields.

The exact correction order is controlled by the solver implementation and the iteration settings in `fvSolution`.

---

## Model families

### Relative velocity

The relative-velocity framework describes non-diffusive motion of the dispersed gas phase relative to the mixture.

Available models in this branch include:

- `noDrift`
- `stokes`
- `DahlkildU`
- `SchillingsU`
- `stuckBubbles`

`stuckBubbles` combines a mobile-bubble base model with an attached/static gas fraction determined from a coverage model.

### Phase diffusion

Phase diffusion is separated from the non-diffusive relative velocity and represented through a diffusion/dispersion tensor.

Available models include:

- `noDiffusion`
- `constant`
- `solidPressure`
- `Schillings`
- `Dahlkild`
- `BearScheidegger`

This separation allows buoyancy/drift and saturation-gradient-driven transport to be selected independently.

### Bubble coverage

Coverage models determine the electrode surface fraction covered by gas bubbles and are used by the electrochemical and interfacial-transfer models.

Available models include:

- `zeroCoverage`
- `linearCoverage`
- `vogt`
- `linearVogt`
- `verticalVogt`
- `linearVerticalVogt`

### Other supporting models

The solver also contains dedicated classes for:

- porous properties;
- species properties;
- mixture viscosity;
- mass and species transfer;
- concentration losses;
- reaction properties;
- interfacial mass-transfer correlations.

---

## Requirements

A compatible OpenFOAM installation is required.

Before compiling, load the OpenFOAM environment in the usual way for the installed version, for example:

```bash
source <path-to-OpenFOAM>/etc/bashrc
```
The solver was written for v2206, but have been running without problems on 2506.


## Compilation

The supporting libraries should be compiled before the solver.

From each library source directory, use the appropriate OpenFOAM build command, typically:

```bash
wmake libso
```

Then compile the solver itself from its application directory:

```bash
./Allwmake
```

### Custom boundary condition

The tutorial cases also depend on a custom boundary condition that is **not compiled automatically with the solver**.

Compile the boundary-condition library separately before running the cases:

```bash
cd /src/BC/pidUNePlaneFvPatchScalarField
wmake
```

The resulting library must be available to OpenFOAM at run time, normally through `$FOAM_USER_LIBBIN`, and the corresponding case must load it through its OpenFOAM configuration where required.

The exact boundary-condition library name and source path are not documented here because they were not included with the files used to prepare this README.

---

## Case setup

A case generally follows the normal OpenFOAM directory structure:

```text
case/
├── 0/          # initial and boundary fields
├── constant/   # physical properties and model selection
└── system/     # numerical schemes, solvers, and run control
```

In addition to the standard flow fields, a case may contain fields for:

- gas volume fraction;
- electrolyte species concentrations;
- electrode and electrolyte potentials;
- porous-region indicators and properties;
- temperature;
- electrochemical source terms.

The exact required fields depend on the selected models.

Model coefficients and run-time selections are supplied through the corresponding dictionaries in `constant/` and `system/`.

---



## Numerical configuration

Important solver controls are located in `system/fvSolution`.

The solver-specific correction settings include:

- species corrections;
- potential corrections;
- outer chemical corrections.

Standard OpenFOAM PIMPLE, linear-solver, relaxation, and phase-fraction controls must also be defined consistently with the selected case.

Transport discretisation is configured in `system/fvSchemes`.

---

## Output

Depending on the selected models, the solver can write fields including:

- velocity and pressure;
- gas volume fraction;
- mixture and phase properties;
- species concentrations;
- electrode and electrolyte potentials;
- current distributions;
- bubble coverage;
- gas-liquid mass-transfer rates;
- relative velocity;
- phase-diffusion coefficients;
- porous-media properties.

Not every field is written by every model.

---

## Tutorials

Tutorial cases should contain their own README describing:

- the purpose of the case;
- geometry and mesh;
- selected physical models;
- important boundary and initial conditions;
- required custom libraries;
- execution commands;
- expected qualitative or quantitative behaviour.

A reusable tutorial README template is included as:

```text
README_TUTORIAL.md
```

---

## Publications

The solver is associated with the alkaline water electrolysis modelling work supplied with this repository/branch. Consult the accompanying publications for physical-model derivations, assumptions, and validation cases.

When publishing results produced with this solver, cite the relevant associated publication(s).

---

## Development status

This is research software. Model availability and dictionary entries may change between branches or solver iterations.

Before transferring an old case to a newer branch, check:

- model names;
- dictionary keywords;
- required fields;
- custom libraries;
- boundary conditions;
- numerical settings.

A case that runs with an earlier solver revision should not be assumed to be directly compatible with this branch.

---

## License

Individual source files retain their existing OpenFOAM/GPL license headers. Refer to the source files and repository license information for the applicable licensing terms.
