# BEM Inductance Extraction (Hodge Decomposition)

Self-inductance extraction using Hodge decomposition + `ngsolve.bem.LaplaceSL`.

## Method (EFIE-based)

Solves for the current distribution via Hodge decomposition on a genus-1 surface:

1. **Divergence matrix** D: `HDivSurface` -> `SurfaceL2`
2. **Vertex-edge incidence** C (graph curl constraint)
3. **Harmonic subspace**: `ker([D; C^T M_J])` (2D for genus-1 torus)
4. **Generalized eigenvalue**: `SL_harm * c = lambda * M_harm * c`
5. **Inductance**: `L = mu_0 * lambda * R / a` (3D torus scaling)

Reference: Lucy Weggler's ngsbem framework (Section 2 of `low_freq_efie_ngbem_applications.ipynb`).

## Key Settings

- **`use_fmm=False`**: Reproducible results, faster dense extraction (Joachim Schoeberl, 2026-03-22)
- **`ToDense().NumPy()`**: Optimized dense extraction
- **`Glue(torus.faces)`**: Required for OCC surface-only mesh (correct Euler characteristic)
- **`mesh.Curve(p)`**: Capped at p=2 for Cubit meshes (Python callback bottleneck)

## Verification

Neumann formula: `L = mu_0 * R * (ln(8R/a) - 2)` = 149.67 nH (R=50mm, a=5mm)

| Mesh | p | n_J | L (BEM) | Error |
|------|---|-----|---------|-------|
| OCC cs=1.0 | 1 | 1,203 | 148.01 nH | -1.11% |
| OCC cs=1.5 | 1 | 2,859 | 149.40 nH | -0.18% |
| OCC cs=2.0 | 1 | 5,103 | 149.73 nH | -0.04% |
| Cubit | 2 | 4,611 | 149.21 nH | -0.31% |

## Files

| File | Description |
|------|-------------|
| `inductance_hodge.py` | OCC torus: Hodge decomposition (standalone, no Cubit) |
| `test_bem_inductance.py` | Cubit torus: Hodge decomposition + GMSH export |
| `inductance_torus.py` | Cubit model creation (torus with gap, source/sink blocks) |
| `taskmanager_bem_test.ipynb` | TaskManager reproducibility investigation |

## Usage

```bash
python inductance_hodge.py                     # OCC torus, convergence study
python test_bem_inductance.py                  # Cubit torus + GMSH export
python test_bem_inductance.py --order 1        # Cubit, flat mesh (p=1)
```
