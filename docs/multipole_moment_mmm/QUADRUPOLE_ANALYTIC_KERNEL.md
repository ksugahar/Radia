# Analytic Closed-Form Moment Kernel (`moment_analytic_kernel`)

Status: **default ON** (2026-07-02). Verified live + correct on method 0/1/2;
the 64-point Gauss kernel remains available as an explicit cross-check path with
`rad.SolverConfig(moment_analytic_kernel=False)`.

## What

The multipole-moment surface-charge kernel `CentroidFieldGradFromFace`
(`src/core/rad_interaction.cpp`) computes, per source face at a target element
centroid, the demag field `H[3]` and its gradient `gH[6]` (the quadrupole
field-gradient = the Hessian of the single-layer potential `I0 = ∫_T 1/|r-r'| dS'`).

- **Default** (`rad.SolverConfig(moment_analytic_kernel=True)`): each face is
  fan-triangulated and integrated with a **closed form**
  (`FieldGradFromChargedTriangleLocal`):
  - `H` = van Oosterom–Strackee (the existing `FieldFromChargedTriangleLocal`),
  - `gH` = the **Mathematica-verified symbolic gradient** of that `H`
    (`gH_ij = ∂H_i/∂obs_j`, the quadrupole field-gradient).
- **Cross-check** (`rad.SolverConfig(moment_analytic_kernel=False)`): a 64-point
  (8×8) Gauss bilinear-quad quadrature for both `H` and `gH`.

The closed-form `gH` is assembled from **only the tangential (log-term) derivatives**
plus tracelessness (`Gzz = -(Gxx+Gyy)`) and symmetry (`Gxz = ∂HH1/∂e3`,
`Gyz = ∂HH2/∂e3`) — so the `atan` (solid-angle) derivative is never needed, which
keeps it well-conditioned near the source plane.

## Why it is correct (derivation + verification)

`gH_ij = ∫_T [δ_ij/r³ − 3 d_i d_j/r⁵] dS' = ∂/∂obs_j ∫_T d_i/r³ dS' = ∂H_i/∂obs_j`,
i.e. the production 64-pt-Gauss `gH` is exactly the Jacobian of the analytic `H`.
Since `H` is already closed form, the closed-form `gH` is its symbolic gradient.

Reproduce / verify (self-testing WolframScript):

- [`quadrupole_hessian_derivation.wls`](quadrupole_hessian_derivation.wls) —
  derives `gH` by symbolic differentiation of the analytic `H` and checks it against
  the production 64-pt Gauss kernel. Result (general-position triangles, off-plane
  observation):
  - `gH` rel vs 64-pt Gauss: **4.6e-6 … 1e-12** (= the Gauss quadrature error; the
    closed form is exact, so it is in fact *more* accurate than the default),
  - tensor **symmetric to 6.2e-16**, **traceless to 8.9e-16** (∇·H = 0, automatic),
  - local→global congruence transform matches the direct global Jacobian to 5.6e-16.
- [`quad_split_validation.wls`](quad_split_validation.wls) — the 2-triangle fan vs
  the bilinear-quad surface: **planar quad exact** (rel 1e-8 … 1e-13 = Gauss error);
  a non-planar (warped) hex face differs by `O(warp)` (~1–4 % at 5 % warp) — the same
  flat-triangulation modeling choice the existing analytic field path already makes.

## Wiring

The single analytic-capable `CentroidFieldGradFromFace` is the one source of truth;
every moment path delegates to it in the analytic branch (the Gauss path keeps its
precomputed 64-sample fast loop). The face corners are cached in
`MomentGeomFaceCache.V4` for the fast paths.

| Method | Path (file `rad_interaction.cpp`) | Wired |
|---|---|---|
| 0 (dense LU) | `BuildCentroidFieldGrad` | yes |
| 1 (matrix-free BiCGSTAB) | `MomentKernelMatVec6x6` | yes |
| 2 (HACApK) | `MomentSystemBlock6x6` / `MomentSystemEntry` | yes |
| mixed hex+wedge/pyramid | `MomentSystemBlockAny` | yes |

Golden lock: [`validation_test/feec/test_moment_analytic_kernel.py`](../../validation_test/feec/test_moment_analytic_kernel.py)
asserts each path is **live** (changes the result vs Gauss — guards against a silent
no-op), reproduces the Gauss physics (ext-B rel < 5e-3), gives the correct cube demag
(`M_z ~ 3·H0`), and that the flag round-trips + defaults on + does not leak.

## Performance (honest)

LAB smoke timing (method-0/2 build, n×n×n hex, best of 3) showed the analytic build
is **~1.2–1.4× faster** than 64-pt Gauss. The mdx knob matrix on 2026-07-02 measured
the dominant method-2 H-matrix build at about **1.5× faster** (`ctype` 28k DOF:
2.97 s → 1.95 s; compact cube 24.6k DOF: 2.76 s → 1.84 s). It is **not** the ~64×
that "64× fewer evaluations/face" would
suggest: the closed form is **transcendental-heavy** (per face ~3 `log` + 3 `atan` +
6 `sqrt` + the rank-2 congruence transform), which offsets most of the fewer-points
advantage, and the method-0 build is further diluted by the O(N³) LU.

So the value is primarily **exactness** (the closed form removes the 64-pt Gauss
quadrature error) plus a **modest** build speedup. The default flip has landed; keep
the Gauss path only as a deliberate cross-check / regression-diagnosis switch.
