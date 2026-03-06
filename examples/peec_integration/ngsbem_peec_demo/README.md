# ngsbem_peec_demo — PEEC Impedance Extraction with ngsbem

PEEC impedance extraction using ngsbem (NGSolve BEM), with material coupling
via EFIE + SIBC for both magnetic cores and conducting shields.

## Notebooks

| Notebook | Description |
|----------|-------------|
| `peec_impedance.ipynb` | Full PEEC Loop-Star workflow: mesh, BEM assembly, 3-mode impedance sweep, material coupling |

## Scripts

| Script | Requires Radia? | Description |
|--------|-----------------|-------------|
| `compute_L_final.py` | No | Loop inductance of rectangular frame via ngbem HelmholtzSL |
| `benchmark_4case.py` | **Yes** | 4-case benchmark: Air / Ferrite / Shield / Both (PEEC filament + ngbem) |

## Modules

| Module | Description |
|--------|-------------|
| `ngbem_peec.py` | PEEC Loop-Star solver (`NGBEMPEECSolver`) with 3 modes: mqs, full, stabilized |
| `ngbem_eddy.py` | Eddy current solver (`ShieldBEMSIBC`): EFIE + SIBC for conducting/magnetic bodies |
| `ngbem_interface.py` | Bridge: ngbem matrices to Radia PEEC topology |
| `ngbem_coupled.py` | Coupled PEEC + Radia MMM solver (lazy Radia import) |

## Subdirectories

| Directory | Description |
|-----------|-------------|
| `fasthenry-like/` | Demos using Radia PEEC filament model (requires `peec_matrices.pyd`) |
| `ngbem/` | Pure ngbem demos (no Radia dependency) |

## Prerequisites

```bash
pip install ngsolve==6.2.2405
pip install ngsolve-ngsbem
pip install matplotlib numpy scipy
```

## Key concept: PEEC Loop-Star = Stabilized EFIE

The `HDivSurface` x `SurfaceL2` product space used for PEEC Loop-Star decomposition
is identical to the one in Weggler's stabilized EFIE formulation. This provides:

- O(1) condition number from DC to RF (no low-frequency breakdown)
- Natural separation of inductive (Loop) and capacitive (Star) DOFs
- Same matrices serve both circuit extraction (PEEC) and scattering (EFIE)

Reference: [Maxwell_DtN_Stabilized.ipynb](https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb)
