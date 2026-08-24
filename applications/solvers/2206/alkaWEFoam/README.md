# alkaWEFoam

`alkaWEFoam` is a research OpenFOAM solver for multiphysics simulation of alkaline water electrolysis.  
It extends the drift-flux mixture formulation with electrochemical reactions, species transport, gas evolution, electrode surface coverage, porous-electrode effects, and gas/liquid mass transfer.

The code was developed for research use and is intended primarily to support reproducibility and further model development.

## Model overview

The solver couples:

- gas/liquid two-phase flow using a drift-flux mixture formulation;
- transport of dissolved and gaseous species;
- electrolyte and electrode potential equations;
- Butler-Volmer reaction kinetics;
- porous-electrode properties;
- electrode surface coverage by gas;
- concentration-loss corrections;
- attached- and detached-bubble mass transfer;
- water-vapour effects.

The main solver uses a PIMPLE loop to couple the hydrodynamic and electrochemical parts of the model.

## Main solver structure

The central executable is `alkaWEFoam.C`. The main sequence is approximately:

```text
Create mesh, fields and physical models
        |
        v
Update water partial pressure
        |
        v
PIMPLE loop
    |
    +-- Update relative velocity
    +-- Update electrode gas coverage
    +-- Solve gas volume fraction
    +-- Update mixture properties
    +-- Solve momentum equation
    +-- Solve pressure equation
    +-- Update turbulence
    +-- Update species diffusivities and reaction properties
    |
    +-- Electrochemical correction loop
            |
            +-- Update concentration-loss model
            +-- Solve electrode/electrolyte potentials
            +-- Solve species transport equations
```

## Important equation files

| File | Purpose |
| --- | --- |
| `alkaWEFoam.C` | Main solver, time loop and coupling sequence |
| `createFields.H` | Creates fields and runtime-selectable physical models |
| `alphaEqnSubCycle.H` | Gas/liquid volume-fraction transport |
| `UEqn.H` | Mixture momentum equation |
| `pEqn.H` | Pressure correction and conservative flux update |
| `potentialEqn.H` | Electrode and electrolyte potential equations with Butler-Volmer kinetics |
| `CiEqn.H` | Species transport, electrochemical source terms and gas transfer |

## Mass and species transfer

Mass transfer is handled through the runtime-selectable base class:

```text
massAndSpeciesTransferModel
        |
        +-- mixedSaturation
```

`massAndSpeciesTransferModel` provides the common fields and interfaces for:

- electrochemical species production (`Psi_BV`);
- water and pure-water vapour pressure;
- species mass-transfer rates;
- wall/bubble mass-transfer source terms.

`mixedSaturation` combines two transfer mechanisms:

1. **Attached bubbles** - transfer at gas-covered electrode surfaces.
2. **Detached bubbles** - transfer between dissolved species and bubbles in the bulk electrolyte.

Saturation concentrations are evaluated using Henry-law-type relations, while transfer coefficients are supplied through the Sherwood-number models.

## Electrode surface coverage

Electrode gas coverage is represented by the runtime-selectable `coverageModel`.

Available models include:

| Model | Description |
| --- | --- |
| `zeroCoverage` | Sets gas coverage to zero |
| `linearCoverage` | Coverage based on local gas volume fraction |
| `vogt` | Current-density, temperature and pressure dependent coverage correlation |
| `linearVogt` | Uses the larger value from the gas-fraction and Vogt-type correlations |

The resulting coverage field `theta` reduces the electrochemically active electrode area through factors such as `(1 - theta)`.

## Electrochemistry

`potentialEqn.H` solves the potentials in:

- the negative electrode (`UNe`);
- the positive electrode (`UPe`);
- the electrolyte (`Ue`).

The Butler-Volmer expressions are linearized to allow the potential equations to be solved implicitly.  
The resulting transfer-current field `J` is subsequently used in the species source terms and coverage models.

## Species transport

`CiEqn.H` solves the transported species concentrations.

The equations include:

- transient accumulation;
- advection with the liquid-phase flux;
- diffusion;
- electrochemical production/consumption;
- gas-liquid mass transfer;
- optional water-vapour transfer.

The current implementation uses the species ordering defined by the corresponding species/reaction property classes.

## Runtime-selectable models

Several physical submodels use the standard OpenFOAM runtime-selection mechanism.  
The selected model and its coefficients are therefore supplied through the case dictionaries rather than being hard-coded into the solver.

Examples include:

- mass and species transfer models;
- electrode coverage models;
- Sherwood-number models;
- relative-velocity models;
- concentration-loss models.

## Building and running

The solver follows the standard OpenFOAM build and case structure.  
Compile the required libraries and solver using the repository `Make/files` and `Make/options`, then run `alkaWEFoam` from a prepared OpenFOAM case.

Exact dictionary entries and boundary conditions depend on the physical case being simulated.

## Repository notes

This repository contains research code. Some classes and solver components originate from or are derived from OpenFOAM and retain their original copyright and GPL license headers.

When modifying the model, the most important coupling points are:

- `potentialEqn.H` for electrochemical kinetics and current distribution;
- `CiEqn.H` for species source terms;
- `massAndSpeciesTransferModel` / `mixedSaturation` for gas-liquid transfer;
- `coverageModel` and its derived classes for electrode gas coverage;
- `alphaEqnSubCycle.H`, `UEqn.H` and `pEqn.H` for hydrodynamics.

## Reference

The model formulation and validation are described in the associated 2025 publication supplied with this code release.

If this repository is used in academic work, please cite the accompanying publication.

## License

The source files derived from OpenFOAM are distributed under the GNU General Public License as stated in their individual file headers. See the repository license and source-file headers for details.
