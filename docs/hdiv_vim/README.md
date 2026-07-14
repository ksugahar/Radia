# HDiv-Type VIM

HDiv-VIM is Radia's production soft-iron demagnetization route.  It keeps the
magnetization, material state, charge map, and reduced-FEM handoff in NGSolve
mesh/function-space vocabulary, then uses Radia's C++ charge-Gram H-matrix for
the open-boundary integral operator.

## Live API

```python
import radia as rad
import radia.vim as vim

iron = vim.MeshSoftIron(mesh, mu_r=1000)
model = rad.ObjCnt([iron, source])
result = rad.Solve(model, 1e-3, 100, 2, demag_backend="hdiv")
field = rad.Fld(model, "b", [0, 0, 0.02])
```

For direct VIM use, call `radia.vim.Solve(mesh, mu_r=... | bh_table=...,
H_ext=..., image=...)`.

The public three-dimensional VIM contract supports RT1 on pure TET/HEX/WEDGE
meshes and RT2 on flat pure TET meshes.  RT2 uses the same C++ charge-Gram,
mass-Riesz CG, and energy-Newton material paths; HEX/WEDGE, 2D, IMA, and
`rad.Fld` field reconstruction remain RT1-only and fail loudly at RT2.
Curved geometry remains a production route through RT1 on an isoparametric P2
mesh; curved RT2 is gated until its Duffy-build cost is acceptable.
For a geometrically and topologically symmetric reduced/full hex
pair, `rad.Fld` after an image solve must agree with the explicit full solve at
the roundoff contract (`< 10 eps` relative error), not merely within a percent.

The production 3D solve uses symmetric C++ CG on the SPD `W + B^T G B`
system.  The public default is `preconditioner="auto"`: linear solves and
small tet nonlinear solves use the exact mass-Riesz map, while nonlinear
energy-Newton on pure hex/wedge meshes and medium-or-larger tet meshes uses the
exact diagonal of `W + B^T G B`.  The diagonal branch is still the same
symmetric energy-Newton system; it changes only the inner Krylov
preconditioner and avoids one PARDISO phase-33 mass solve per CG iteration.
The result artifact records both `preconditioner_requested` and the resolved
`preconditioner`, plus `preconditioner_policy` so benchmark JSON shows why
`auto` chose that branch.  The current tet switch point is `6000` HDiv DOF
(`RADIA_HDIV_AUTO_JACOBI_TET_NFACE` overrides it for measurement sweeps).
`preconditioner="mass-riesz"` and `preconditioner="jacobi"` remain explicit
diagnostic overrides.

For flat pure-hex `.vol` meshes, the order-1 charge-basis path now follows the
NGSolve reference ordering directly: the Q2 geometry lattice is built from the
linear `.vol` vertices, and the Q1 shape-moment to monomial map is applied as a
cached block-diagonal sparse transform.  Curved `.vol` meshes still use
`GetTrafo` as the geometry source of truth.

The hex charge-Gram build caches the already symmetrized host-pair block
`0.5*(AB + BA^T)`.  For sufficiently far hex host pairs the product quadrature is
orientation-symmetric, so the build uses the one-sided `AB` block directly in the
symmetric H-matrix; near/self pairs still use the explicit `0.5*(AB + BA^T)`
average.  The default far threshold is `1.0*(size_A + size_B)` and can be
disabled with `RADIA_HDIV_HEX_FAR_ONESIDED=0` for diagnostics.  The C++ linear
solve also reports `solve_*` timing fields, and Python passes sparse inputs
through NumPy-array pybind entry points instead of materializing large Python
lists.  The hot hex block-cache hit/miss counters are opt-in via
`RADIA_HDIV_HEX_CACHE_STATS=1`; ordinary timing runs avoid the per-entry atomic
counter overhead.  On the LAB N=20 cube smoke run, the current mass-Riesz path is
about 15.6 s wall time with about 7.9 s in HACApK build, versus about 20.0 s and
12.3 s before the far one-sided symmetric-block optimization.

The first mdx cube timing sweep on 2026-07-08 reached structured-hex N=40
(`1.56M` HDiv DOF) in about `355 s` wall time and `116 s` H-matrix build time.
The matching Netgen tet comparison showed that unstructured tet runs need a
tighter ACA tolerance (`gram_eps=1e-8` to `1e-10`) to keep mass-Riesz CG near
25--30 iterations; too-loose `1e-4` compression can inflate CG to hundreds of
iterations or fail at 4000 iterations.  The hex far one-sided idea also applies
to the flat high-order tet FAR path: `QuadDotFar(a,b)` is a symmetric low-order
double integral, so the direct FAR term is one-sided by default
(`RADIA_HDIV_HO_FAR_ONESIDED=0` restores the diagnostic average).  On mdx this
cut tet H-matrix build by roughly 8--12% with machine-level `M_avg_z`
agreement; the 1.37M-DOF tet cube moved from about `428 s` to `386 s`.  With
the stable tet tolerance and this tet-side optimization, the optimized hex path
is still faster than the current tet path at comparable cube DOF.  The
remaining performance focus is split by regime: charge-Gram build for linear
or small-tet runs, and H-matrix apply count / nonlinear globalization for large
hex-wedge energy-Newton runs.

The 2026-07-09 mdx nonlinear energy-Newton preconditioner sweep showed why the
`auto` policy is element- and size-aware.  For structured hex, the tuned
diagonal branch reduced N=16 from `247.8 s` to `33.7 s` and completed N=20 in
`79.9 s`, with `M_avg_z` matching the mass-Riesz tuned result at the
sub-ppm-to-few-ppm level.  For wedge, N=10 moved from `60.2 s` to `24.5 s`.
For tet, exact mass-Riesz is still faster on small meshes (`maxh=0.20`,
`4314` DOF: `4.22 s` mass-Riesz vs `4.79 s` diagonal), but the diagonal branch
wins by `maxh=0.15` (`8193` DOF: `8.82 s` vs `13.21 s`).  The supporting JSON
and CSV are stored under the manuscript `results/` directory as
`hdiv_energy_newton_preconditioner_mdx_20260709_*`.  The cube benchmark
drivers expose these timing defaults as `--nonlinear-profile mdx-scaling`
(`nl_tol=3e-4`, `newton_continuation=2`, `newton_reuse_tangent_steps=3`,
`preconditioner=auto`); the default `strict` profile keeps solver-regression
settings.

The wedge/prism path uses the same block-serving infrastructure but is less
optimized than the structured-hex path.  The first wedge mdx pass enabled a
matching FAR host-pair one-sided block (`RADIA_HDIV_WEDGE_FAR_ONESIDED=0`
restores the diagnostic average).  At wedge N=8 and N=12 this reduced
H-matrix build by about 6--8% with `M_avg_z` changes below `6e-9` relative.
Wedge now also uses translated host-block reuse by default:
`RADIA_HDIV_WEDGE_TRANS_CACHE=all` / `2` reuses cell-cell, cell-face, and
face-face template blocks; `1` falls back to the older cell-cell-only subset;
`0` disables the translation cache.  The earlier face-bearing drift was traced
to non-deterministic tie breaks in wedge near quadrature.  After stabilizing the
source-site and target-corner choices, a dense N=2 wedge charge-Gram comparison
between `0` and `all` matched to `8.6e-17` relative Frobenius error, and an
N=4 solve matched `M_avg_z` to `3.2e-16` relative while cutting local H-matrix
build from `2.14 s` to `0.315 s`.  The wedge charge-basis build now mirrors the
hex fast path on linear prism meshes: Q2 lattice nodes are interpolated from
mesh vertices, and the small L2/SurfaceL2-to-monomial transforms are cached and
applied as sparse block transforms instead of solving per element.  On mdx this
cut wedge N=12 basis construction from about `59 s` to `0.75 s`.  With both
changes, mdx N=8/N=12 translation-cache A/B at `gram_eps=1e-8` moved from
`12.0/44.2 s` to `2.52/6.59 s` wall time; H-matrix build moved from
`9.31/34.84 s` to `1.04/2.34 s`, with relative `M_avg_z` drift
`4.8e-13/2.5e-8`.  The follow-up apply/solve profile added optional
`hmatvec_*` timing fields (`RADIA_HDIV_HMATVEC_STATS=1`) and limits MKL's
thread count inside HACApK leaf dgemv calls through
`RADIA_HACAPK_MATVEC_MKL_THREADS` (default `1`) so TaskManager remains the
outer parallel layer while PARDISO mass-Riesz solves can still use the process
thread setting.  On mdx, the updated all/default wedge sweep reached N=20
(`293.6k` HDiv DOF) in `25.3 s` wall time with `6.7 s` H-matrix build and
`8.1 s` solve time; the charge-Gram H-matvec portion of that solve was
`0.64 s`, leaving PARDISO mass-Riesz factor/apply as the next solve-side
target.  Wedge also tolerated `gram_eps=1e-4` in earlier cube checks, cutting
build further while moving `M_avg_z` by only `2e-5`--`4e-5` relative to the
`1e-8` reference.  Treat that as a promising wedge timing lane, not yet the
universal default for all wedge geometries.

## Why This Is The Main Route

- `N = B^T G B` is symmetric and loop-free by construction: loops are
  `ker(B)`, so charge-free modes are field-null without hand-built loop bases.
- NGSolve `Mesh`, `GridFunction`, `CoefficientFunction`, `BilinearForm`, and
  `TaskManager` are shared with reduced FEM and motor workflows.
- Curved/high-order geometry uses the same finite-element geometry path instead
  of translating through a separate object discretization.
- TET/HEX/WEDGE and 2D planar validation live in `validation_test/feec/`.

## Validation Expectations

Fast tests belong in `tests/`; heavier numerical checks belong in
`validation_test/feec/` and JSON artifacts labelled with the actual validation
host (`mdx` or `hibino`).  Required HDiv gates:

- analytic demag factors where available;
- nonlinear BH convergence metadata;
- image symmetry checked against an explicitly mirrored full model on truly
  symmetric meshes;
- `rad.Fld` after `rad.Solve(..., demag_backend="hdiv")`;
- RT1/RT2 flat pure-TET accuracy and cost, plus charge-Gram H-matrix build stats and
  memory/timing on idle `mdx` or `hibino`
  for large runs;
- 2D planar motor saliency checks for the motor lane.

## Public Docs

Current user-facing notebooks in this directory should be result-bearing and
paired with synchronized JSON sidecars.  Old migration notes and comparison
archives are not part of the public docs surface; recover them from git history
only if needed.
