# ngsbem_peec_demo — PEEC Impedance Extraction Notebook

Jupyter notebook demonstrating PEEC impedance extraction using ngsbem (NGSolve BEM),
with coupling to magnetic cores (Radia MMM) and conducting shields (BEM+SIBC).

## Notebook

**[peec_impedance.ipynb](peec_impedance.ipynb)** — Full PEEC Loop-Star workflow:

1. Surface mesh generation (Netgen OCC)
2. BEM matrix assembly (inductance L, potential coefficients P)
3. MQS impedance sweep (R + jwL, 10 Hz -- 1 MHz)
4. Full Loop-Star impedance with Schur complement
5. Vector potential from solved currents
6. High-order convergence study (p = 0, 1, 2)
7. Ferrite core coupling via Radia MMM (inductance increase)
8. Conducting shield coupling via BEM+SIBC (inductance decrease, Lenz's law)
9. Connection to Weggler's stabilized EFIE (condition number comparison,
   combined ferrite/shield coupling plot)

## Prerequisites

- NGSolve with ngsbem (Sections 1--6, 9)
- Radia (`pip install radia`) for Sections 7--8 (optional, gracefully skipped if missing)

## Files

| File | Description |
|------|-------------|
| `peec_impedance.ipynb` | Main notebook |
| `ngbem_peec.py` | PEEC Loop-Star solver (HDivSurface + SurfaceL2) |
| `ngbem_coupled.py` | Coupled PEEC + core solver (Radia MMM / FEM-BEM) |
| `ngbem_eddy.py` | Eddy current solvers (FEM-BEM, BEM+SIBC) |
| `ngbem_interface.py` | Bridge: ngbem matrices to Radia PEEC topology |

All Python modules are self-contained copies from `src/radia/` for portability.

## Running

```bash
jupyter notebook peec_impedance.ipynb
```

## Key concept: PEEC Loop-Star = Stabilized EFIE

The `HDivSurface` x `SurfaceL2` product space used for PEEC Loop-Star decomposition
is identical to the one in Weggler's stabilized EFIE formulation. This provides:

- O(1) condition number from DC to RF (no low-frequency breakdown)
- Natural separation of inductive (Loop) and capacitive (Star) DOFs
- Same matrices serve both circuit extraction (PEEC) and scattering (EFIE)

Reference: [Maxwell_DtN_Stabilized.ipynb](https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb)

## See also

- Radia source modules: `src/radia/ngbem_*.py`
- Demo scripts: [`../ngbem/`](../ngbem/)
- Validation: [`../verification/validate_stabilized_bem.py`](../verification/validate_stabilized_bem.py)
