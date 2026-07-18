# ESIM (Effective Surface Impedance Method): CLI Usage Guide

User-facing companion to [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md).
This document explains how to invoke ESIM through Radia's three Layer-4
calc scripts and what flags each accepts.  For physics / discretisation
details refer to the mathematical analysis.

**Method origin.**  The "Karl iteration" mentioned throughout this
document is **Karl Hollaus's** outer Picard fixed-point loop for the
nonlinear scalar-potential ESIM formulation:

> K. Hollaus, M. Kaltenbacher, J. Schöberl, *"A Nonlinear Effective
> Surface Impedance in a Magnetic Scalar Potential Formulation,"*
> **IEEE Trans. Magn.**, 2025.  DOI:
> [10.1109/TMAG.2025.3613932](https://doi.org/10.1109/TMAG.2025.3613932).

External readers should cite this paper at first use of "Karl
iteration"; the in-source code uses the name as lab shorthand for the
canonical Hollaus-type Picard relaxation.

---

## 1. When to choose ESIM over linear SIBC

| Scenario | Use |
|---|---|
| Cu / Al / brass workpiece (`μ_r ≈ 1`) | `--impedance-model sibc` (linear Dowell) |
| Steel / ferrite workpiece, **mid-frequency**, `|H_t|` stays on one side of the BH knee | `--impedance-model sibc` with reasonable `--mu-r` |
| Steel workpiece traversing the BH knee (saturation pattern matters) | **`--impedance-model esim`** + `--bh-file <table.txt>` |
| Strong spatial saturation contrast across the workpiece surface | `--impedance-model esim --esim-per-panel` (BEM), or the local HCurl-VIM API in Section 3.4 |

Linear SIBC is closed-form and fast.  The scalar/per-panel CLI paths add an
outer Karl iteration and runtime one-dimensional cell solves.  The local
HCurl-VIM API can instead use the persistent LUT in Section 3.4; its online
outer iterations perform no cell solves.  Per-panel runtime ESIM without that
LUT multiplies cell-solve work by the number of sampled panels/DOFs.

---

## 2. BH-curve file format

A BH curve is required for `--impedance-model esim`.  The format is
two-column whitespace-separated `H [A/m]   B [T]`, ascending in H,
including the origin:

```
0.0     0.0
10.0    0.00126
50.0    0.00628
100.0   0.0125
500.0   0.061
1000.0  0.12
5000.0  0.55
10000.0 1.20
50000.0 1.85
100000.0 1.95
```

Optional third / fourth columns are ignored.  Comment lines start with
`#`.  Steel sample BH curves ship under
[`src/radia/panels/samples/`](../../src/radia/panels/samples/) (e.g.
`em_sample_bh.txt`).

---

## 3. Three production CLIs

### 3.1 `calc_inductance.py` (BEM-SIBC workpiece, PEEC or BEM-A coil)

The PEEC+BEM weak-coupling path.  Workpiece BIE solved on a scalar
potential; the SIBC enters as a complex Robin coefficient on the
boundary integral equation.

```bash
python src/radia/panels/calc_inductance.py \
    --coil-solver peec --coil-step coil.step \
    --vol workpiece.vol --wp-label sibc \
    --sigma 2e6 --mu-r 100 --half-thickness 0.005 \
    --frequency 100e3 --current 1.0 \
    --impedance-model esim --bh-file em_sample_bh.txt \
    --esim-max-iter 15 --esim-tol 1e-3 \
    --esim-relax 0.5
```

ESIM-specific flags:

| Flag | Default | Meaning |
|---|---|---|
| `--impedance-model esim` | `sibc` | Switch from linear Dowell to ESIM cell-problem |
| `--bh-file <path>` | (none) | Required when `--impedance-model esim` |
| `--esim-max-iter` | 15 | Outer Karl iteration cap |
| `--esim-tol` | 1e-3 | Outer convergence on `max\|dZ_s\|/\|Z_s\|` |
| `--esim-relax` | 0.5 | Karl damping (under-relaxation); lower if oscillation observed |
| `--esim-per-panel` | False | Per-DOF Z_s mode (BEM-A path only; raises if combined with `--wp-bem-backend hacapk`) |

The workpiece geometry is the `cylinder` cell-problem mode in all cases
(see § 3.2 of [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md)).
`--half-thickness` is the cylinder radius for the cell-problem; for
solid bulk workpieces, pass `min(R_wp, H_wp/2)`.

### 3.2 `calc_fem_kelvin.py` (PEEC coil + FEM workpiece + Kelvin)

The FEM-SIBC + Kelvin path; coil is a PEEC filament bundle (line
integral), workpiece is HCurl A-formulation with a Robin BC.

```bash
python src/radia/panels/calc_fem_kelvin.py \
    --vol workpiece.vol \
    --fes-order 1 \
    --material custom --sigma 2e6 --mu-r 100 \
    --half-thickness 0.005 \
    --frequency 100e3 --current 1.0 \
    --formulation total \
    --impedance esim --bh-file em_sample_bh.txt \
    --max-iter 15 \
    --solver pardiso \
    --peec-step coil.step --peec-sigma 5.8e7 \
    --peec-n-peri 16 --peec-nwinc 3 --peec-nhinc 3 \
    --require-kelvin
```

ESIM-specific flags (note: the flag is `--impedance`, not
`--impedance-model`, and the iteration cap is `--max-iter`, not
`--esim-max-iter`):

| Flag | Default | Meaning |
|---|---|---|
| `--impedance esim` | `sibc` | Switch from linear Dowell to ESIM |
| `--bh-file <path>` | (none) | Required when `--impedance esim` |
| `--max-iter` | 15 | Outer Karl iteration cap |

There is no `--esim-tol` equivalent here — the tolerance is hard-coded
at `1e-3` in the script.  Per-panel ESIM is not yet wired into the FEM
Robin coefficient.

### 3.3 `calc_fem_coilmesh.py` (Full FEM A-V + workpiece SIBC + Kelvin)

The most physically complete (and most expensive) path: coil is a
volumetric mesh with H1 source potential and HCurl A, workpiece SIBC
appears as a Robin BC.

```bash
python src/radia/panels/calc_fem_coilmesh.py \
    --vol workpiece.vol \
    --frequency 100e3 --current 1.0 \
    --coil-sigma 5.8e7 --sigma 2e6 --mu-r 100 \
    --half-thickness 0.005 \
    --fes-order 1 \
    --solver pardiso \
    --sibc-bnd sibc --source-bnd source --sink-bnd sink \
    --coil-mat coil \
    --impedance-model esim --bh-file em_sample_bh.txt \
    --esim-max-iter 15 --esim-tol 1e-3 \
    --require-kelvin
```

The flag set matches `calc_inductance.py` (note: this one uses
`--impedance-model` and `--esim-max-iter`, unlike `calc_fem_kelvin.py`).
Per-panel ESIM is not wired here either — the Robin term uses a scalar
`s/Z_s`.

### 3.4 Local ESIM on the HCurl Eddy-Bubble VIM

The topology-aware HCurl-VIM path has a separate local surface API. It keeps
one ESIM value per air-facing surface quadrature sample, assembles the complex
surface Gram without spatial averaging, and repeats the reduced VIM solve in a
Hollaus-type outer Picard iteration:

```python
import numpy as np
import radia.vim as vim

model = vim.LocalESIMSurfaceModel(
    bh_curve=np.loadtxt("steel_bh.txt")[:, :2],
    sigma=2.0e6,
    bins=12,
    n_nodes=100,
)
result = hybrid.solve_local_esim(
    model,
    frequency_hz=100_000.0,
    outer_tolerance=1.0e-3,
    outer_max_iterations=20,
    outer_relaxation=0.5,
)
assert result.converged
```

`hybrid` is the `TopologyAwareHybridVIM` returned by
`NgsolveTopologyAwareHybridVIM` or `NgsolveEddyBubbleHybridVIM`. The lower-level
entry is `SolveLocalESIMSurfaceVIM(system, surface_basis, rhs, model, frequency)`.
Because ESIM is nonlinear, both entries accept exactly one physical excitation;
multi-column port superposition is rejected.

The one-dimensional cell problem deliberately uses SciPy's sparse/tridiagonal
solvers. HCurl basis construction, VIM assembly, HACApK action, mixed-Galerkin
Schur reduction, and the reduced linear solve remain on the NGSolve/Radia native
path. `SurfaceImpedanceGram` records the surface mass used at assembly and is
rejected if passed to a different hybrid system.

For repeated design sweeps, precompute the cell problems once as a persistent
two-dimensional `Z_s(f, |H_t|)` LUT:

```python
frequencies = np.geomspace(1.0e3, 1.0e6, 13)
fields = np.geomspace(model.h_floor, 3.0e4, 96)
lut = vim.BuildLocalESIMSurfaceLUT(
    model,
    frequencies,
    fields,
    output_path="steel_local_esim.npz",
)
quality = vim.ValidateLocalESIMSurfaceLUT(model, lut)
print(quality["max_relative_error"])

# A later process performs no SciPy cell solve in the Karl iteration.
lut = vim.LocalESIMSurfaceLUT.load("steel_local_esim.npz")
production_model = model.with_lut(lut)
result = hybrid.solve_local_esim(
    production_model,
    frequency_hz=100_000.0,
)
assert all(row["cell_solve_count"] == 0 for row in result.history)
```

The NPZ archive is pickle-free and carries a SHA-256 signature of the BH curve,
conductivity, and cell discretization/tolerances. A mismatched model is rejected.
Interpolation of `log(|Z_s|)` and continuous impedance phase is bilinear in
`log(f)` and `log(|H_t|)`; extrapolation is also rejected, so the offline grid
must cover the full expected operating envelope.
Generated LUT files are runtime/cache artifacts and are not committed as Radia
source binaries.

A LUT is sufficient when the local ESIM cell response is single-valued in the
chosen state variables and a node-refinement check bounds interpolation error.
`ValidateLocalESIMSurfaceLUT` evaluates every field-interval midpoint at all
frequency nodes and frequency midpoints with direct cell solves. Select the
field-node count from that reported error, rather than treating 96 as universal.
Use separate tables when BH data or conductivity changes, including a material
temperature update. A two-dimensional LUT is not a substitute for additional
state variables required by hysteresis, rotational magnetization, or a genuinely
two/three-dimensional corner cell problem.

The converged Gram can be passed to `CoupledHDivHybridVIMSystem.solve_frequency`
for a fixed HDiv-MMM magnetic operator. A simultaneous outer iteration that
also updates the ordinary nonlinear HDiv constitutive operator is not yet a
one-call API; it must not be described as completed production functionality.

---

## 4. Reading the JSON output

The Karl iteration history is exposed in the result JSON for
diagnostic plotting:

```json
{
  "impedance_model": "esim",
  "esim_iterations": 6,
  "esim_converged": true,
  "esim_history": [
    {"iteration": 0, "Z_s_abs": 3.58e-2, "H_t_rms": 247.3, "dZ": 1.0, "t_solve": 0.21},
    {"iteration": 1, "Z_s_abs": 3.52e-2, "H_t_rms": 261.0, "dZ": 0.017, "t_solve": 0.20},
    {"iteration": 2, "Z_s_abs": 3.49e-2, "H_t_rms": 268.4, "dZ": 0.008, "t_solve": 0.20},
    {"iteration": 3, "Z_s_abs": 3.48e-2, "H_t_rms": 270.8, "dZ": 0.003, "t_solve": 0.20},
    {"iteration": 4, "Z_s_abs": 3.48e-2, "H_t_rms": 271.4, "dZ": 0.001, "t_solve": 0.20},
    {"iteration": 5, "Z_s_abs": 3.48e-2, "H_t_rms": 271.5, "dZ": 0.0003, "t_solve": 0.20}
  ],
  ...
}
```

When `--esim-per-panel` is used (BEM path only) the schema gains:

- `esim_per_panel: true`
- `esim_per_panel_Z_s_real: [...]`, `esim_per_panel_Z_s_imag: [...]`:
  per-DOF Z_s ndarray, listed in BEM DOF order
- `esim_per_panel_H_t: [...]`: per-DOF |H_t| at convergence
  (radia ≥ 4.55.x; same DOF order as `esim_per_panel_Z_s_real`).
  Use this for spatial visualisation (see
  [`EXAMPLES.md`](EXAMPLES.md) → `plot_zs_per_dof_map.py`).
- `Z_s_wp_real` / `Z_s_wp_imag`: area-weighted mean (for back-compat)
- Each `esim_history` entry adds `Z_s_abs_min / Z_s_abs_max`,
  `H_t_per_dof_mean / H_t_per_dof_max`, `dZ_max`

---

## 5. Choosing `--esim-relax`

The default `0.5` works for all production samples we ship.  Lower
values trade convergence speed for stability.  Guidance:

| Symptom | Suggested `--esim-relax` |
|---|---|
| Convergence in <8 iter, dZ monotone decreasing | leave at 0.5 |
| Oscillation: dZ goes up after first few iter | drop to 0.3 |
| Iteration count > 20 with monotone dZ at default | raise to 0.7 |
| Deep saturation (B saturates at `|H| > 30 kA/m`) and oscillation | drop to 0.2 + raise `--esim-max-iter` to 30 |

Anderson acceleration (planned, roadmap § 7 of
[`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md)) will make the
relaxation parameter less relevant.

---

## 6. Troubleshooting

| Error message | Cause | Fix |
|---|---|---|
| `--impedance-model esim requires --bh-file` | ESIM requested without a BH table | Pass `--bh-file <path>`; see § 2 |
| `BH curve is empty` / `not monotone in H` | Malformed BH-file | Verify two-column ASCII, ascending H, includes (0, 0) |
| `ESIM:NOT-CONVERGED after N iter` | Karl loop hit `max_iter` before `dZ < tol` | First inspect `esim_history` with [`plot_karl_history.py`](../ih_esim_benchmark/plot_karl_history.py).  If `Z_s_abs` / `H_t_rms` are plateaued and only the per-DOF `dZ_max` failed to drop, the run is usable (see [`IMPLEMENTATION.md`](IMPLEMENTATION.md) § 3.4).  Otherwise: raise `--esim-max-iter` or lower `--esim-relax`; check BH curve monotonicity. |
| `--esim-per-panel ... wp-bem-backend hacapk` | per-panel ESIM not yet supported on HACApK | Use `--wp-bem-backend intree-dense`, or fall back to scalar Z_s |
| `cell solver SCIPY_AVAILABLE False` | scipy not installed in the calc-side Python | `pip install scipy` |

---

## 7. Cross-references

- Mathematical formulation:
  [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md)
- Cell-problem implementation:
  [`src/radia/esim_cell_problem.py`](../../src/radia/esim_cell_problem.py)
- Integration tests:
  [`tests/test_esim_integration.py`](../../tests/test_esim_integration.py)
- Benchmark (Bessel linear baseline):
  [`docs/ih_esim_benchmark/`](../ih_esim_benchmark/)
- Material loader (BH file parser):
  [`src/radia/em_material.py`](../../src/radia/em_material.py)
- Radia IH notebook workbench (drives the CLIs):
  [`src/radia/panels/notebooks/radia_ih.ipynb`](../../src/radia/panels/notebooks/radia_ih.ipynb)
  plus [`src/radia/ih_design.py`](../../src/radia/ih_design.py)

---

**Document version**: 2026-05-30, written against radia v4.67.0+.
