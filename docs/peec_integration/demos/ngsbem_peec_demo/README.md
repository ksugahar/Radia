# ngsbem_peec_demo — PEEC Impedance Extraction with ngsbem

PEEC impedance extraction using ngsbem (NGSolve BEM), with material coupling
via EFIE + SIBC for both magnetic cores and conducting shields.

## Result-Saved Notebooks

The old notebook draft was retired because it had no saved execution results
and depended on helper scripts that now live in `validation_test`. The public
result-bearing entry points are:

| Notebook | Description |
|----------|-------------|
| `../../public_demo.ipynb` | Human-facing PEEC demo selection and representative source excerpts |
| `../../post_examples_migration.ipynb` | Migration ledger for docs demos, validation scripts, and distilled deletions |

## Scripts

| Script | Requires Radia? | Description |
|--------|-----------------|-------------|
| `validation_test/peec_integration/ngsbem_peec_demo/compute_L_final.py` | No | Loop inductance of rectangular frame via ngbem HelmholtzSL |
| `validation_test/peec_integration/ngsbem_peec_demo/benchmark_4case.py` | **Yes** | 4-case benchmark: Air / Ferrite / Shield / Both (PEEC filament + ngbem) |

## Modules

| Module | Description |
|--------|-------------|
| `validation_test/peec_integration/ngsbem_peec_demo/ngbem_peec.py` | PEEC Loop-Star solver (`NGBEMPEECSolver`) with 3 modes: mqs, full, stabilized |
| `validation_test/peec_integration/ngsbem_peec_demo/ngbem_eddy.py` | Eddy current solver (`ShieldBEMSIBC`): EFIE + SIBC for conducting/magnetic bodies |
| `validation_test/peec_integration/ngsbem_peec_demo/ngbem_interface.py` | Bridge: ngbem matrices to Radia PEEC topology |
| `validation_test/peec_integration/ngsbem_peec_demo/ngbem_coupled.py` | Coupled PEEC + Radia MMM solver (lazy Radia import) |

## Subdirectories

| Directory | Description |
|-----------|-------------|
| `validation_test/peec_integration/ngsbem_peec_demo/fasthenry-like/` | Demos using Radia PEEC filament model (requires `peec_matrices.pyd`) |
| `validation_test/peec_integration/ngsbem_peec_demo/ngbem/` | Pure ngbem demos (no Radia dependency) |

## Prerequisites

```bash
pip install ngsolve>=6.2.2601
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
