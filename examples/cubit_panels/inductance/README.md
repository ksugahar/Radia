# BEM Inductance Extraction (Source/Sink Saddle Point EFIE)

Self-inductance extraction using `ngsolve.bem.LaplaceSL` with source/sink constrained EFIE.

## Method

Solves the constrained EFIE on the conductor surface:

```
[SL  D^T] [J] = [0]
[D   0  ] [p] = [g]
```

where:
- SL = LaplaceSL (single layer BEM operator on surface currents)
- D = divergence matrix (HDivSurface -> SurfaceL2)
- g = source/sink current injection (+1 at source, -1 at sink)
- J = surface current, p = Lagrange multiplier

Inductance: `L = mu_0 * J^T @ SL @ J`

B-distribution: direct Biot-Savart `B(x) = mu_0/(4pi) * sum J_e x r / |r|^3 * A_e`

## Key Settings

- **`use_fmm=False`**: Reproducible results, faster dense extraction
- **`ToDense().NumPy()`**: Optimized dense extraction (68x faster)
- **`SurfaceL2(order=0)`**: Constraint always order=0 (higher order causes rank deficiency)
- **`mesh.Curve(2)`**: Quadratic curving for geometry accuracy

## Verification

Neumann formula: `L = mu_0 * R * (ln(8R/a) - 2)` = 149.67 nH (R=50mm, a=5mm)

| Method | n_J | L (BEM) | Error |
|--------|-----|---------|-------|
| Source/sink (OCC, maxh=a) | 5,085 | 144.5 nH | -3.5% (gap) |
| Source/sink (Cubit hex) | 2,208 | 144.9 nH | -3.2% |
| Energy (OCC, closed) | 5,103 | 149.7 nH | -0.04% |

Note: -3.5% error is due to the 5-degree gap (not mesh error). Gap -> 0 converges to analytical.

## GMSH Visualization

Results output as combined `inductance.geo` merging:
- `inductance_B.msh` — volume |B| and B vector fields
- `inductance_J.msh` — surface |J| and J vector fields
- `inductance_coil.msh` — coil wireframe (1D line elements)

All views independently toggleable in GMSH Post-processing tree.
`.geo` sets `Mesh.NumSubEdges = 4` for curved Tri6 display.

## Files

| File | Description |
|------|-------------|
| `inductance_source_sink.py` | Source/sink EFIE (OCC gapped torus, standalone) |
| `inductance_hodge.py` | Hodge decomposition (OCC closed torus, legacy) |
| `test_bem_inductance.py` | Cubit torus: Hodge decomposition + GMSH export |
| `inductance_torus.py` | Cubit model creation (torus with gap) |
| `inductance_torus.cub5` | Pre-built Cubit model |

## Cubit Panel

The Cubit panel (`src/radia/panels/calc_inductance.py`) provides GUI access:
- Journal editor with default hex sweep torus
- Auto source/sink block detection
- Solve + Open GMSH button (combined J + B visualization)

## Usage

```bash
python inductance_source_sink.py         # OCC gapped torus (standalone)
python inductance_hodge.py               # OCC closed torus (legacy)
```
