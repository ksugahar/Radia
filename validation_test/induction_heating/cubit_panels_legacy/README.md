# BEM Inductance + ESIM Surface Impedance

This directory is the legacy validation corpus moved from the former IH
Cubit-panel examples location.  It is not a public examples tier; use it to
preserve old checks while reusable kernels migrate to `src/` and public
walkthroughs migrate to result-saved docs notebooks.

Self-inductance extraction using `ngsolve.bem.LaplaceSL` with source/sink constrained EFIE.
Optional workpiece surface impedance (ESIM/Dowell) for induction heating analysis.

## Method

### Inductance (BEM)

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

### Workpiece Surface Impedance (ESIM/Dowell)

When a `workpiece` block is defined, the pipeline extends:

1. **BEM** on coil surface -> surface current J
2. **Biot-Savart** from J -> H at workpiece surface panels
3. **ESIM** or **Dowell** -> surface impedance Z_s(H), power loss P per panel
4. Integration -> total R, P, Q for workpiece

**ESIM** (Effective Surface Impedance Method): solves 1D cell problem
`rho * d^2H/dz^2 + j*omega*mu(|H|)*H = 0` per panel. Handles nonlinear BH curves.

**Dowell**: analytical formula `Z_s = (rho/a) * gamma*a * tanh(gamma*a)` for linear materials.

## Key Settings

- **`use_fmm=False`**: Reproducible results, faster dense extraction
- **`ToDense().NumPy()`**: Optimized dense extraction (68x faster)
- **`SurfaceL2(order=0)`**: Constraint always order=0 (higher order causes rank deficiency)
- **`mesh.Curve(2)`**: Quadratic curving for geometry accuracy

## Verification

### Inductance

Neumann formula: `L = mu_0 * R * (ln(8R/a) - 2)` = 149.67 nH (R=50mm, a=5mm)

| Method | n_J | L (BEM) | Error |
|--------|-----|---------|-------|
| Source/sink (OCC, maxh=a) | 5,085 | 144.5 nH | -3.5% (gap) |
| Source/sink (Cubit hex) | 2,208 | 144.9 nH | -3.2% |
| Energy (OCC, closed) | 5,103 | 149.7 nH | -0.04% |

Note: -3.5% error is due to the 5-degree gap (not mesh error). Gap -> 0 converges to analytical.

### ESIM

Verified against NGSolve H1 FEM (p=4) as independent method (`verify_esim.py`):

| Test | ESIM vs | Max error | Result |
|------|---------|-----------|--------|
| Linear Z_s | Analytical rho*gamma*tanh(gamma*a) | 0.25% | PASS |
| Linear Z_s | NGSolve H1 FEM (p=4) | 0.04% | PASS |
| Nonlinear Z_s (steel BH) | NGSolve H1 FEM + Picard | 0.78% | PASS |
| H(z) profile | Analytical / FEM | < 0.2% | PASS |

### Coupled BEM coil terminal Delta_L (2026-04-12)

`src/radia/bem_coupled_solver.py::CoupledBEMSolver` computes the coil
terminal inductance change due to a workpiece, with **per-DOF
back-reaction RHS** (`f_back[i] = int v_i.Trace() . A_wp dS_coil`
assembled via NGSolve LinearForm). The previous v1 implementation
used a scalar rescale and produced wrong-signed Delta_L; the new
version is sign-correct in both the Lenz screening and flux
concentration regimes.

**Validation**: independent FEM-Kelvin SIBC solve via
`calc_fem_kelvin.py` on the same `radia_model.vol`.

| Material | mu_r | f | L_BEM | L_FEM | diff |
|---|---|---|---|---|---|
| copper, sigma=5.8e7 | 1   | 50 kHz | **84.31 nH** | **84.56 nH** | **+0.29%** |
| steel,  sigma=2e6   | 100 | 50 kHz | 87.92 nH | 89.43 nH | +1.72% |

L_air (coil only) = 87.81 nH. Both methods report dL < 0 for copper
(Lenz screening) and dL > 0 for steel (flux concentration in the
ferromagnetic skin layer). The 0.3% agreement on copper is the
strongest validation we have for the coupled BEM solver.

Run the cross-check:
```bash
python validation_test/induction_heating/cubit_panels_legacy/compare_bem_coupled_vs_fem_kelvin.py
```

Background notes: `memory/bem_coupled_solver_existing.md`

## GMSH Visualization

Results output as combined `inductance.geo` merging:
- `inductance_B.msh` -- volume |B| and B vector fields
- `inductance_J.msh` -- surface |J| and J vector fields
- `inductance_coil.msh` -- coil wireframe (1D line elements)

All views independently toggleable in GMSH Post-processing tree.
`.geo` sets `Mesh.NumSubEdges = 4` for curved Tri6 display.

## Files

| File | Description |
|------|-------------|
| `inductance_source_sink.py` | Source/sink EFIE (OCC gapped torus, standalone) |
| `inductance_hodge.py` | Hodge decomposition (OCC closed torus, legacy) |
| `inductance_torus.py` | Cubit model creation (torus with gap) |
| `inductance_torus.cub5` | Pre-built Cubit model |
| `impedance_esim.py` | BEM coil + ESIM workpiece coupled analysis |
| `verify_esim.py` | ESIM verification against analytical + NGSolve FEM |
| **`compare_bem_coupled_vs_fem_kelvin.py`** | **Coupled BEM vs FEM-Kelvin SIBC cross-check (canonical regression record, 2026-04-12)** |

## Cubit Panel

The Cubit panel (`src/radia/panels/calc_inductance.py`) provides GUI access:
- Journal editor with default hex sweep torus
- Auto source/sink block detection
- **Workpiece block detection** -> ESIM/Dowell settings appear automatically
- Model selection: ESIM (nonlinear) / Dowell (analytical)
- Material selection: Steel / Copper / Aluminum (sigma auto-set)
- Result table: L, R, P, Q, skin depth, total Z
- Solve + Open GMSH button (combined J + B visualization)

### Cubit blocks for workpiece

```
block N add volume <workpiece_vid>
block N name "workpiece"
```

When the `workpiece` block exists, the panel shows additional settings:
- **Model**: ESIM / Dowell
- **Material**: Steel / Copper / Aluminum
- **Frequency** [Hz]
- **Sigma** [S/m] (auto-updated by material selection)
- **Half-thickness** [m] (slab model parameter)

## Usage

```bash
# Inductance only (no workpiece)
python inductance_source_sink.py

# BEM + ESIM coupled analysis
python impedance_esim.py                        # Steel, 50 kHz (default)
python impedance_esim.py --material copper      # Copper workpiece
python impedance_esim.py --material steel --sweep  # Frequency sweep 1kHz-1MHz
python impedance_esim.py --freq 100000 --wp-radius 0.015

# ESIM verification
python verify_esim.py           # All tests (analytical + NGSolve FEM)
python verify_esim.py --test 3  # Nonlinear only
```

## Tests

```bash
# Layer 1: Computation (CI, no Cubit)
pytest tests/panels/test_calc_inductance.py -v    # 8 tests, ~60s

# Layer 2: UI logic (CI, no Cubit, no Qt)
pytest tests/panels/test_panel_ui_logic.py -v     # 17 tests, ~0.5s

# Layer 2b: GUI widget (CI, no Cubit, requires PySide6)
pytest tests/panels/test_panel_gui.py -v          # 13 tests, ~1.5s

# Layer 3: Integration (local only, requires Cubit)
python tests/panels/test_panel_integration.py
```
