# (ACA+)+TSVD least-norm solver — entry point

> **For the full documentation, see [docs/stream_function/](stream_function/README.md)**.
> That folder has structured sub-docs for theory, single-stroke chain construction,
> regularisation choices, surface deformation, API reference, benchmarks,
> NGSolve 6.2.2604+ `ngsolve.bem` integration, and the paper outline.
>
> This file is the SHORT entry point; treat it as the README that hasn't been
> split out yet.  All numbers, all design choices, all dead-ends are documented
> in detail in `docs/stream_function/`.

Accelerated, **kernel-agnostic** least-norm solver for field-synthesis /
inverse-source problems -- the numerical core of the *stream function method*
of coil design, generalised to any Radia source family.

- **Full docs**: [`docs/stream_function/`](stream_function/README.md)
- Module: [`src/radia/stream_function.py`](../src/radia/stream_function.py)
- C++ core: `src/core/rad_stream_function.{h,cpp}`
- Executed public notebooks: [`theory.ipynb`](stream_function/theory.ipynb),
  [`regularization.ipynb`](stream_function/regularization.ipynb), and
  [`deformation.ipynb`](stream_function/deformation.ipynb)
- Tests: `tests/test_stream_function.py`
- MCP knowledge: `aca_tsvd(topic=session_2026_05_30)` for the full session-log
  narrative; `aca_tsvd(topic=single_stroke)` for the chain construction detail.

## 1. Problem

Given `M` field (observation) points and `N` basis sources, the field-coupling
matrix `A` and the design system are

```
A phi = B            A in R^{M x N},  usually M < N  (underdetermined)
A(i,j) = (a field component) at observation i produced by basis source j
B(i)   = desired field at observation i
phi(j) = unknown strength of basis source j
```

The least-norm solution uses the **truncated-SVD (TSVD) pseudo-inverse**

```
A ~= U diag(S) V^T            (SVD)
phi = V diag(1/S) U^T B       (truncated to the first k singular triplets)
```

Truncation at `k` modes regularises the otherwise ill-posed inverse: small
singular values (which amplify noise) are dropped.  Sweeping `k` traces the
familiar L-curve of residual `||A phi - B||` vs solution norm `||phi||`.

## 2. Why (ACA+)+TSVD

Forming the dense `A` and taking a full SVD costs `O(N M^2)`.  For smooth
kernels `A` is numerically **low rank**, so we instead:

1. **ACA+** (Adaptive Cross Approximation, plus pivoting) factors
   `A ~= C D^T` with `C in R^{M x k_aca}`, `D in R^{N x k_aca}`,
   `k_aca << min(M,N)`, evaluating only `O(k_aca (M+N))` entries of `A` --
   never the full `M N`.
2. **Recompression to a truncated SVD** of the small factors -- the standard
   "SVD of a low-rank product" (peer review JIAM-2026-36).  QR each tall-skinny
   factor `C = Qc Rc`, `D = Qd Rd`, then take ONE small SVD of the
   `k_aca x k_aca` core `Rc Rd^T = U_bar Sigma V_bar^T`.  This gives the exact
   SVD of the ACA approximation:

   ```
   A ~= C D^T = Qc (Rc Rd^T) Qd^T = (Qc U_bar) Sigma (Qd V_bar)^T
   U = Qc U_bar,   S = Sigma,   V = Qd V_bar
   ```

   QR (backward-stable, cheaper than an SVD of the same factor) orthogonalises
   the factors and the only SVD is the tiny `k_aca x k_aca` core, so `U` and `V`
   come out **orthonormal** -- exactly what the folded regularisation
   (`RegularizedTSVD`) requires.

   Exactly **two** methods are kept (the manuscript Method 2/3 were removed; see
   `memory/aca_tsvd_qr_recompression.md`): `method="aca_qr_tsvd"` (this fast
   ACA + QR + TSVD path, the default) and `method="dense"` -- the plain/direct
   TSVD (materialise `A` and take its dense SVD), the **exact reference** the
   review noted the paper lacked.  NO backward compatibility: legacy `method=2/3`
   and the terse `"qr"`/`"aca"` now raise.

Net cost is roughly `(M/k_aca)^2` lower than the dense route.

**Measured** (`validation_test/stream_function/bench_aca_vs_dense.py`, smooth
`1/(1+alpha r^2)` kernel, `M = N/4`, same per-call kernel for both methods, LAB
2026-05-29).  The kernel is numerically low rank so `k_aca` stays ~constant
(~30) while `N` grows; the eval-count reduction `M*N -> ~k_aca(M+N)` therefore
widens with `N`, and -- since both methods call the same kernel -- shows up
almost one-for-one in wall-clock time:

| N | M | k_aca | kernel evals (naive -> ACA) | eval reduction | time (naive -> ACA) | speedup | dS/S |
|---|---|-------|-----------------------------|----------------|---------------------|---------|------|
| 256  | 64  | 31 | 16 384 -> 9 892    | 1.7x  | 37 ms -> 23 ms      | 1.6x  | 4.9e-9 |
| 512  | 128 | 33 | 65 536 -> 21 781   | 3.0x  | 141 ms -> 48 ms     | 2.9x  | 3.6e-9 |
| 1024 | 256 | 34 | 262 144 -> 45 588  | 5.8x  | 554 ms -> 96 ms     | 5.8x  | 2.0e-9 |
| 2048 | 512 | 36 | 1 048 576 -> 98 469 | 10.6x | 2217 ms -> 204 ms   | 10.8x | 2.0e-9 |

Singular values match the dense SVD to ~2e-9.  Notes:
- The speedup **grows with N** (it is `~N/(5 k_aca)` for `M=N/4`); at fixed
  `k_aca` the dense route is `O(N M^2)` while ACA is `O(k_aca(M+N))` evals.
- The win is largest when the kernel is **expensive** (Biot-Savart / fixed-magnet
  field, `radia.Fld`) -- every avoided `A(i,j)` is an avoided field evaluation.
- The matrix entry is supplied through a Python callback, so there is a
  per-call overhead; because both methods pay it equally, the *ratio* is
  eval-count-driven and a C++ kernel lowers both absolute times by the same
  factor without changing the speedup.

## 3. Kernel-agnostic design

The solver **embeds no field kernel**.  The matrix entry `A(i,j)` is supplied
by the caller as a callback, so the same machinery serves every Radia source
family using Radia's *already-implemented* field computation:

| Source family | Radia kernel |
|---------------|--------------|
| coils (thin wires) | Biot-Savart (`ObjFlmCur`, `ObjArcCur`) |
| permanent magnets / soft iron | Radia field evaluation (`ObjRecMag`, HDiv-VIM solved materials, ...) |

ACA+ itself is **delegated to the in-repo HACApK C library**
(`src/ext/HACApK/cHACApK_acaplus`) -- the single source of truth for ACA+ in
Radia. There is no second ACA+ implementation. The Biot-Savart kernel feeds
HACApK through the `HACApK_set_entry_func` override (default behaviour, the
in-repo material interaction matrix, is unchanged when the override is null). This module's only
numerical algorithm is the TSVD recompression, which HACApK does not provide.

## 4. API

```python
from radia.stream_function import (
    aca_tsvd, pseudo_inverse_solve, solve, radia_field_kernel, StreamTSVD,
)
```

### `aca_tsvd(M, N, entry, modes=None, kmax=None, aca_eps=1e-4, method="aca_qr_tsvd") -> StreamTSVD`

(ACA+)+TSVD of the `M x N` matrix whose entries are returned by
`entry(i, j) -> float` (0-based `i in [0,M)`, `j in [0,N)`). `entry` is called
on demand by ACA+, not over the full grid.

- `modes`   -- singular triplets to return (clamped to `k_aca`); default `kmax`.
- `kmax`    -- maximum ACA+ rank; default `min(M, N)`.
- `aca_eps` -- ACA+ stopping tolerance (absolute pivot threshold).
- `method`  -- `"aca_qr_tsvd"` (default) = ACA + QR + TSVD (the fast path above);
  `"dense"` / `"tsvd"` = the direct dense TSVD (exact reference, small problems).
  NO backward compatibility: legacy `2`/`3` and the terse `"qr"`/`"aca"` raise
  (the manuscript Method 2/3 were removed).

Returns a `StreamTSVD` with `U (M,modes)`, `S (modes,)`, `V (N,modes)` (row-major
NumPy arrays), `k_aca`, and `method`.

### `pseudo_inverse_solve(result, B, k_mode=None) -> phi`

`phi = V diag(1/S) U^T B` using the first `k_mode` modes (default `result.modes`).
`B` has length `M`; returns `phi` of length `N`.

### `solve(M, N, entry, B, modes=None, k_mode=None, ...) -> (phi, result)`

Convenience: `aca_tsvd` then `pseudo_inverse_solve`.

### `radia_field_kernel(obs_points, sources, component=2, field="b") -> entry`

Builds the `entry(i, j)` callback from a list of Radia object handles via
`radia.Fld`: `A(i,j) = (component of field) at obs_points[i] from sources[j]`.
Works for any source family. `component` is `0/1/2` for `x/y/z`; `field` is the
Radia field id (`"b"`, `"h"`, `"a"`, ...).

### Example

```python
import numpy as np, radia as rad
from radia.stream_function import aca_tsvd, pseudo_inverse_solve, radia_field_kernel

loops = [rad.ObjFlmCur(square_loop(c), 1.0) for c in centers]   # N current loops
obs   = ...                                                     # (M,3) field points

entry = radia_field_kernel(obs, loops, component=2)             # A(i,j) = Bz_i(loop_j)
res   = aca_tsvd(len(obs), len(loops), entry, modes=20)
phi   = pseudo_inverse_solve(res, B_target, k_mode=10)          # loop currents
```

### Headless CLI (Stage 2): `calc_stream_coil.py`

The cylindrical Gz gradient-coil demo
([`docs/stream_function/demo_coil_design_gz.py`](stream_function/demo_coil_design_gz.py))
is also packaged as a headless Stage-2 application script,
[`src/radia/panels/calc_stream_coil.py`](../src/radia/panels/calc_stream_coil.py).
It accepts argparse input, writes JSON on stdout, has no GUI, and never imports
`cubit`. It exposes every design knob as a CLI flag
(`--radius --length --gradient --dsv --n-rings --n-wires --modes --aca-eps
--method`) and prints `{k_aca, modes, n_wires, fitted_dBdz,
gradient_nonlinearity, continuous_residual, wire_z, wire_I}`.

```bash
python src/radia/panels/calc_stream_coil.py --radius 0.15 --length 0.5 --gradient 1.0
```

It is locked by a golden test,
[`validation_test/panels/test_stream_coil_golden.py`](../validation_test/panels/test_stream_coil_golden.py),
which runs the CLI with defaults and asserts the discrete wire coil reproduces
the unit gradient (`fitted_dBdz` in `[0.9, 1.1]`) with on-axis
`gradient_nonlinearity < 0.05`. This remains a Stage-2-only utility (validated
CLI + golden band); promotion should add a DesignSpec-backed masked Simulink
block and register it in the application interface manifest.

## 5. Validation

- **vs the true dense matrix**: reconstruction `||A - U S V^T|| / ||A||` reaches
  machine precision when `k_aca = min(M,N)` and tracks `aca_eps` otherwise.
- **vs the Fortran reference** `coil_solver.f90` (`method_aca_tsvd_1/2`, a
  faithful port of the same HACApK ACA+): identical `k_aca` and
  `||S_f90 - S_radia|| / ||S_f90|| ~ 1e-15` for both methods
  (`tests/test_stream_function.py::test_matches_f90_reference`, LAB-only).
- **magnetic-material path**: `test_radia_field_kernel_magnets` factors the
  fixed-magnet field matrix of a permanent-magnet array to `< 1e-5`.

## 6. Optimisation layer: linear (TSVD) vs nonlinear (CMA-ES)

The (ACA+)+TSVD solve is the **linear** design layer: it finds the source
*amplitudes* `phi` (with fixed source directions/positions) that best match a
target field.  When the design variables enter the field **nonlinearly** --
magnetization *directions* (angles), magnet positions, coil-region geometry --
the problem is no longer a linear least-norm solve and a black-box optimiser is
the right tool.  This mirrors the "ACA stream function **+ CMA-ES**" workflow
(SA-25-020): the fast linear inner solve is (ACA+)+TSVD; the nonlinear outer
search is CMA-ES.

For CMA-ES, use **Optuna's `CmaEsSampler`** (do not re-implement it):

```python
import optuna
study = optuna.create_study(direction="minimize",
                            sampler=optuna.samplers.CmaEsSampler(seed=42))
study.optimize(objective, n_trials=200)     # objective builds Radia magnets,
                                            # evaluates the field, returns a scalar
```

The public gallery/catalog record for `demo_cmaes_magnet_design.py` captures the
16-magnetization-angle CMA-ES case for a uniform transverse field.
See the `radia-mcp` `optuna_*` tools for sampler choice, multi-objective
(NSGA-II), pruning, and lab BBO recipes.

## 7. References

- Sugahara Lab, "ACA-accelerated stream function method + CMA-ES", IEEJ Joint
  Technical Meeting on Static Apparatus / Rotating Machinery, **SA-25-020**.
- M. Bebendorf, "Approximation of boundary element matrices", *Numer. Math.* 86
  (2000) -- ACA.
- HACApK (ppOpen-HPC, MIT license): `src/ext/HACApK/`.
