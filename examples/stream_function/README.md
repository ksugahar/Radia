# (ACA+)+TSVD least-norm solver -- examples

Accelerated least-norm solver for the **stream function method** of coil
design, generalised to **any Radia source family**.

Given M field (observation) points and N basis sources, the field-coupling
matrix is

```
A phi = B            A in R^{M x N},  M < N  (underdetermined)
A(i,j) = (a field component) at observation i produced by basis source j
```

The design problem "find source strengths `phi` that produce a desired field
`B`" is solved by the **TSVD-regularised pseudo-inverse**
`phi = V diag(1/S) U^T B`, truncated to `k` modes.  Building the dense `A` and
doing a full SVD is `O(N M^2)`.  Instead we factor `A ~= C D^T` with **ACA+**
(rank `k_aca << min(M,N)`) and TSVD only the small factors -- about `(M/k)^2`
faster.

Production module: [`src/radia/stream_function.py`](../../src/radia/stream_function.py)
(`aca_tsvd`, `pseudo_inverse_solve`, `solve`, `radia_field_kernel`).
C++ core: `src/core/rad_stream_function.cpp`.  See
[`../../docs/stream_function.md`](../../docs/stream_function.md) for the method,
API, and design notes.

## Kernel-agnostic by design

The solver embeds **no field kernel**.  The matrix entry `A(i,j)` is supplied by
a callback, so the same machinery serves any Radia source using Radia's
*already-implemented* field computation:

| Source family | Radia kernel | Example |
|---------------|--------------|---------|
| coils (thin wires) | Biot-Savart (`ObjFlmCur`, `ObjArcCur`) | `demo_coil_field_synthesis.py` |
| permanent magnets / soft iron | MMM / MSC (`ObjRecMag`, `ObjHexahedron`) | `demo_magnet_array.py` |

`radia_field_kernel(obs, sources, component, field)` builds the callback from a
list of Radia object handles via `radia.Fld` -- no per-application Biot-Savart
code.

ACA+ itself is delegated to the in-repo **HACApK** C library
(`cHACApK_acaplus`), the single source of truth for ACA+ in Radia.

## Files

| File | What it shows |
|------|---------------|
| `demo_coil_field_synthesis.py` | Coil design: N filament loops, solve loop currents for a target axial-region field; ACA+ compression + TSVD L-curve. |
| `demo_magnet_array.py` | Same solver on a permanent-magnet array (MMM/MSC field) -- proves kernel-agnosticism. |
| `bench_aca_vs_dense.py` | `(ACA+)+TSVD` vs naive dense `numpy.linalg.svd`: time / memory / rank, written to `results_aca_vs_dense.json`. |
| `demo_cmaes_magnet_design.py` | The **nonlinear counterpart**: CMA-ES (Optuna `CmaEsSampler`) optimises the magnetization *directions* (angles) of a magnet array for a uniform transverse field. Linear amplitude design -> (ACA+)+TSVD; nonlinear direction design -> CMA-ES (the "+ CMA-ES" half of SA-25-020). Needs `optuna` (optional). |
| `demo_coil_design_gz.py` | **End-to-end coil design**: cylindrical z-gradient (Gz) coil via the stream function method. Target `Bz=Gz*z` -> azimuthal ring currents (ACA+TSVD) -> stream function `psi(z)` -> equal-current wire rings -> verified on-axis gradient linearity. The axisymmetric Gz problem reduces to a full-ring (1D `psi(z)`) basis. |
| `demo_sf_to_peec_gz.py` | **Full workflow, loop closed**: SF design -> **single-stroke** (one continuous wire) smooth helix with blended crossovers -> CAD STEP (build123d Spline + Frenet swept solid) -> PEEC (`L`, `R`) -> exact Biot-Savart field -> verify `Bz` vs the design `Gz*z`. `--with-peec` adds the STEP + PEEC stages (needs build123d, in `radia`). |
| `demo_coil_design_gx.py` | **Transverse gradient (Gx), the 2D case**: a non-axisymmetric target `Bz=Gx*x` gives a genuine 2D surface stream function `psi(phi,z)` (a "fingerprint" pattern) -> marching-squares contour -> wire loops; verified `Bz` matches `Gx*x` to ~0.8% over the DSV. Numpy Biot-Savart kernel (avoids the ObjFlmCur bug); single-stroke connection of the nested fingerprint loops is future work. |

## Running

```bash
python demo_coil_field_synthesis.py
python demo_magnet_array.py
python bench_aca_vs_dense.py
python demo_cmaes_magnet_design.py        # needs optuna (pip install optuna)
python demo_coil_design_gz.py             # end-to-end Gz gradient coil design
python demo_sf_to_peec_gz.py --with-peec  # full SF -> CAD(STEP) -> PEEC -> field
python demo_coil_design_gx.py             # transverse Gx gradient (2D surface SF)
```

Each script is standalone (no Cubit, no panel UI).  `matplotlib` is optional:
if installed, the demos save a PNG next to the script; otherwise they print an
ASCII summary only.  `demo_cmaes_magnet_design.py` additionally needs `optuna`
(it prints a friendly message and exits cleanly if optuna is missing).

## Expected results

- **`demo_coil_field_synthesis.py`** (M=25 field points, N=64 loops): smooth
  off-plane field is low rank, so ACA+ stops at `k_aca` well below `min(M,N)=25`.
  The TSVD L-curve shows the residual `||A phi - B|| / ||B||` dropping as modes
  are added, with the solution norm `||phi||` rising -- the usual
  regularisation trade-off.  A few modes already reproduce the target to ~1%.
- **`demo_magnet_array.py`** (N permanent magnets): the `(ACA+)+TSVD`
  factorization reconstructs the Radia MMM/MSC coupling matrix to `< 1e-5`
  relative, identical machinery, zero coil-specific code.
- **`bench_aca_vs_dense.py`**: for a smooth (low-rank) kernel, `(ACA+)+TSVD`
  matches the dense singular values to ~1e-12 while running markedly faster as
  `N` grows and `k_aca` stays small.
- **`demo_cmaes_magnet_design.py`**: a 16-dimensional continuous optimisation
  (one magnetization angle per pixel).  CMA-ES drives the relative field-match
  objective down by ~3x from its first trial, producing an approximately uniform
  transverse field with small cross-components.  The residual is set by the
  finite array (physical limit), not the optimiser.
- **`demo_coil_design_gz.py`**: ACA+ compresses the on-axis ring operator to
  `k_aca ~ 7`; the continuous ring-current solution reproduces `Bz=Gz*z` to
  `~4e-4`; the contoured ~32-wire coil achieves `dBz/dz` within ~0.5% of target
  with ~1.4% on-axis nonlinearity over the DSV -- a textbook generalised
  Maxwell-pair gradient coil.
- **`demo_sf_to_peec_gz.py`** (`--with-peec`, ~16 turns): the single-stroke
  conductor (~15 m, one continuous wire) reproduces the design gradient
  (`dBz/dz ~ 0.99`, ~2.6% nonlinearity); the helix sweeps to a clean STEP solid
  (Frenet frame + auto wire radius so turns don't self-intersect); PEEC returns
  `L ~ 38 uH`, `R ~ 16 mOhm` at 1 kHz.  Confirms the SF design survives the
  single-stroke manufacturing constraint.  (CoilBuilder is for planar
  racetrack/saddle coils; a solenoidal helix uses the smooth-helix + Spline path.)
- **`demo_coil_design_gx.py`** (transverse Gx, the 2D case): the
  non-axisymmetric target `Bz=Gx*x` produces a genuine 2D surface stream
  function `psi(phi,z)` (the classic "fingerprint" pattern).  Unlike the
  low-rank axisymmetric Gz problem, the transverse target fills the operator's
  rank -- ACA+ reaches `k_aca = min(M,N) = 123` (no compression here; the 2D
  problem is intrinsically richer).  Marching-squares contours `psi` into ~68
  saddle-shaped wire loops, and the reconstructed `Bz` matches the design
  `Gx*x` to ~0.8% RMS over the DSV.  Uses a numpy Biot-Savart kernel (the
  ObjFlmCur tilted-loop path is unreliable for these non-planar loops);
  single-stroke connection of the nested fingerprint loops is future work.

## References

- Sugahara Lab, "ACA-accelerated stream function method + CMA-ES",
  IEEJ Joint Technical Meeting on Static Apparatus / Rotating Machinery,
  SA-25-020 (manuscript Method 2/3).
- HACApK (ppOpen-HPC, MIT): `src/ext/HACApK/`.
- Validation reference: the f2py-wrapped Fortran `coil_solver.f90`
  (`method_aca_tsvd_1/2`), matched bit-for-bit by
  `tests/test_stream_function.py::test_matches_f90_reference`.
