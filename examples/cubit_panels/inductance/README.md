# BEM Inductance Extraction (Energy Method)

Self-inductance extraction using `ngsolve.bem.LaplaceSL` on conductor surface meshes.

## Method

Energy method (NOT the inverse formula `1/(e^T L^{-1} e)` which diverges with mesh refinement):

```
L = mu_0 * J^T * SL * J    (for I = 1)
```

where:
- `SL` is the LaplaceSL single-layer operator matrix on the conductor surface
- `J` is the toroidal surface current density for unit total current, projected onto `HDivSurface(order=0)`
- No matrix inversion required -- just a matrix-vector product

## Port Modes

| Cubit Blocks | Mode | Description |
|-------------|------|-------------|
| `source` + `sink` | **1-port** | Gap conductor (current in at source, out at sink) |
| No `sink` block | **loop** | Closed loop (e.g., complete torus) |

Block names are case-insensitive. The `source` and `sink` blocks contain surface elements (tri/quad) on the port faces.

## Key Settings (Joachim Schoeberl, 2026-03-22)

- **`use_fmm=False`**: Reproducible results, faster dense matrix extraction
- **`ToDense().NumPy()`**: Optimized dense extraction (replaces manual column-by-column)
- **`TaskManager`**: Must wrap BOTH setup and extraction, or neither

## Mesh Curving

Cubit provides 1st order mesh. Netgen fork (`export_NGSolveCurvedMesh`) curves to order `p` using ACIS geometry via `CallbackGeometry`. Higher `p` improves geometric accuracy (area, normals) on curved surfaces.

## Verification

Neumann formula for thin torus (R >> a): `L = mu_0 * R * (ln(8R/a) - 2)`

| Mesh | ndof | L (BEM) | L (Neumann) | Error |
|------|------|---------|-------------|-------|
| OCC cs=0.5 | 269 | 106.5 nH | 149.7 nH | -28.9% |
| OCC cs=1.0 | 1,790 | 142.1 nH | 149.7 nH | -5.0% |
| OCC cs=2.0 | 10,977 | 148.7 nH | 149.7 nH | -0.6% |
| Cubit (4611 DOF) | 4,611 | 148.7 nH | 149.7 nH | -0.7% |
| Cubit (8082 DOF) | 8,082 | 149.5 nH | 149.7 nH | -0.1% |

## Files

| File | Description |
|------|-------------|
| `test_bem_inductance.py` | Test script (creates torus, runs BEM, compares with Neumann) |
| `inductance_torus.py` | Cubit model creation (torus with gap, source/sink blocks) |
| `inductance_torus.cub5` | Pre-built Cubit model |
| `taskmanager_bem_test.ipynb` | TaskManager reproducibility investigation |

## Usage

```bash
python test_bem_inductance.py                    # tet mesh, p=1 + p=2 comparison
python test_bem_inductance.py --order 2          # tet mesh, p=2 only
python test_bem_inductance.py --mesh hex         # hex mesh
```
