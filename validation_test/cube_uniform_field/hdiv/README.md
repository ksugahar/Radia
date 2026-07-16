# HDiv-VIM Cube Benchmark

This directory contains the HDiv-VIM counterpart of the soft-iron cube
benchmark used in the 2026-08-25 static/rotating-machine MMPM manuscript.

The driver is intentionally separate from the older hexahedron/tetrahedron
benchmark scripts because those results came from the retired six-face MMPM
path.  Current Radia soft-iron demag validation should use the BDM1 HDiv-VIM
entry:

```powershell
python validation_test\cube_uniform_field\hdiv\bench_hdiv_cube.py `
  --sizes 2 4 6 `
  --material both `
  --bh-table "W:\02_学会資料\2026年度\2026_08_25_静止器・回転機@八戸\MMMM@菅原・伊田・矢野\原稿\figures\BH.txt" `
  --output "W:\02_学会資料\2026年度\2027_01_静止器・回転機@\HDiv要素MMM@菅原\results\hdiv_cube_initial_lab_20260707.json"
```

Problem definition:

- 1 m cube centered at the origin.
- Structured pure-hex `N x N x N` mesh.
- Applied field `H0 = 200 kA/m` in `+z`.
- Linear `mu_r=1000` or nonlinear BH table.
- Solver: `radia.vim.Solve`, BDM1 HDiv-VIM, analytic HACApK charge Gram.

Recorded fields include `n_el`, `ndof`, `n_charge`, `iters`, `M_avg`,
`demag`, wall time, memory snapshot, and the solver's `hmat_stats`.

The LAB initial run on 2026-07-07 reached `N=6` only.  A later LAB development
smoke on 2026-07-08 reached `N=20` with the mass-Riesz linear path in about
15.6 s wall time (`hmat_stats.build_time` about 7.9 s) after the far one-sided
hex symmetric-block optimization.  These LAB numbers are smoke / setup results,
not publication timing claims.

The first mdx validation sweep on 2026-07-08 used the same linear cube settings
and wrote the consolidated artifacts to
`W:\02_学会資料\2026年度\2027_01_静止器・回転機@\HDiv要素MMM@菅原\results\hdiv_cube_mdx_scaling_20260708.{json,csv}`:

| N | ndof | wall s | H-matrix build s | CG solve s | peak WSet MB |
|---:|---:|---:|---:|---:|---:|
| 20 | 196800 | 34.12 | 9.65 | 17.90 | 4369.9 |
| 25 | 382500 | 61.15 | 23.47 | 25.42 | 8181.5 |
| 30 | 658800 | 130.49 | 58.72 | 41.90 | 13876.3 |
| 35 | 1043700 | 211.57 | 94.58 | 56.48 | 20945.3 |
| 40 | 1555200 | 354.81 | 115.98 | 82.58 | 31796.6 |

The N=20 one-sided optimization check on mdx gave `38.33 s -> 34.12 s` wall
time and `18.60 s -> 9.65 s` H-matrix build time with a relative `M_avg_z`
change of `4.6e-7`.

## Tet Comparison

Use `bench_hdiv_tet_cube.py` for the current BDM1 HDiv-VIM unstructured-tet
comparison.  Do not use the older `validation_test/cube_uniform_field/tetrahedron`
scripts for HDiv timing claims; those belong to the retired Radia object/MMPM
path.

```powershell
python validation_test\cube_uniform_field\hdiv\bench_hdiv_tet_cube.py `
  --maxh-values 0.08 0.07 0.06 `
  --material linear `
  --gram-eps 1e-8 `
  --output "W:\02_学会資料\2026年度\2027_01_静止器・回転機@\HDiv要素MMM@菅原\results\hdiv_tet_cube_mdx_scaling_20260708.json"
```

mdx tet results from 2026-07-08 are stored in
`hdiv_tet_cube_mdx_scaling_20260708.{json,csv}`.  The selected timing lane uses
`gram_eps=1e-8` through `maxh=0.045`, then `gram_eps=1e-10` for larger cases
where `1e-8` started to inflate Krylov iterations.  The deliberately bad
`gram_eps=1e-4` probe is kept in
`hdiv_tet_cube_mdx_eps_sensitivity_20260708.csv`: at `maxh=0.07`, iteration
count was `733` with `1e-4` but `26` with `1e-8`; at `maxh=0.06`, `1e-4`
failed to converge in 4000 CG iterations.

The hex far one-sided lesson was also applied to flat high-order tet FAR
pairs: `QuadDotFar(a,b)` is a symmetric low-order double integral, so the
direct FAR term no longer computes the diagnostic `0.5*(ab+ba)` average by
default.  Set `RADIA_HDIV_HO_FAR_ONESIDED=0` to restore the old A/B path.  The
mdx A/B artifact is `hdiv_tet_ho_oneside_mdx_ab_20260708.{json,csv}`; tested
rows showed machine-level `M_avg_z` agreement and about 8--12% H-matrix build
speedup.

| maxh | ndof | gram eps | wall s | H-matrix build s | CG solve s | iters | peak WSet MB |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.055 | 187929 | 1e-8 | 35.79 | 30.53 | 4.37 | 27 | 1924.0 |
| 0.040 | 435165 | 1e-10 | 123.42 | 105.86 | 11.38 | 25 | 5448.6 |
| 0.035 | 562515 | 1e-10 | 178.67 | 154.17 | 15.33 | 26 | 7351.9 |
| 0.030 | 1372635 | 1e-10 | 386.16 | 314.74 | 42.21 | 30 | 16513.0 |

Nearest-DOF comparison artifacts:

- `hdiv_hex_tet_cube_mdx_comparison_20260708.csv`
- `hdiv_hex_tet_cube_mdx_nearest_pairs_20260708.csv`

Current read: the optimized structured-hex path is already faster than the
current unstructured-tet path at comparable DOF for this cube.  The tet solve
phase is healthy when the Gram ACA tolerance is tight enough; the remaining
gap is mainly H-matrix build cost, not mass-Riesz CG convergence.

## Wedge Comparison

`bench_hdiv_cube.py` also accepts `--mesh-kind wedge` for structured pure-prism
cube meshes:

```powershell
python validation_test\cube_uniform_field\hdiv\bench_hdiv_cube.py `
  --mesh-kind wedge `
  --sizes 8 12 `
  --material linear `
  --gram-eps 1e-8 `
  --output "W:\02_学会資料\2026年度\2027_01_静止器・回転機@\HDiv要素MMM@菅原\results\hdiv_wedge_cube_mdx_scaling_20260708.json"
```

The wedge path shares the hex/wedge block-serving kernel but has less
structured-cache machinery than the pure-hex path.  On 2026-07-08 the same FAR
one-sided idea was enabled for wedge host-pair blocks through
`RADIA_HDIV_WEDGE_FAR_ONESIDED` (default `1`; set `0` for diagnostic
`0.5*(ab+ba)`).  mdx A/B artifacts are
`hdiv_wedge_oneside_mdx_ab_20260708.{json,csv}`.

The next wedge optimization is translation-block reuse.  The default
`RADIA_HDIV_WEDGE_TRANS_CACHE` scope is now `all`/`2`, which caches translated
cell-cell, cell-face, and face-face host-pair blocks; set `1` for the older
cell-cell-only subset or `0` to disable it.  A previous N=4 `all` drift was
traced to near-quadrature tie breaks for translated wedge hosts, not to the
translation key itself.  After stabilizing those tie breaks, the local N=2 dense
charge-Gram `0` vs `all` comparison matches to `8.6e-17` relative Frobenius
error, and the local N=4 solve matches `M_avg_z` to `3.2e-16` relative while
cutting H-matrix build from `2.14 s` to `0.315 s`.  The earlier mdx A/B artifact
for the conservative cell-cell scope is
`hdiv_wedge_trans_cache_mdx_ab_20260708.{json,csv}`.

The wedge charge-basis path was also brought up to the same structural level as
the hex path for linear prism meshes: Q2 lattice nodes come directly from mesh
vertices, and per-orientation monomial transforms are cached and applied through
sparse block transforms.  The mdx artifacts are:

- `hdiv_wedge_trans_cache_fastbasis_mdx_ab_20260708.{json,csv}`
- `hdiv_wedge_all_fastbasis_mdx_scaling_20260708.{json,csv}`
- `hdiv_wedge_all_fastbasis_applysolve_mdx_scaling_20260708.{json,csv}`

| N | ndof | gram eps | wall s | H-matrix build s | CG solve s | iters | peak WSet MB |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 19328 | 1e-8 | 17.18 | 16.31 | 0.59 | 21 | 1529.6 |
| 8 | 19328 | 1e-4 | 13.64 | 12.87 | 0.50 | 22 | 1170.8 |
| 12 | 64224 | 1e-8 | 98.78 | 94.75 | 2.78 | 21 | 5088.6 |
| 12 | 64224 | 1e-4 | 79.00 | 77.37 | 1.32 | 23 | 4525.9 |

With the optimized wedge charge-basis path, translation-cache A/B at
`gram_eps=1e-8` is:

| N | trans off wall s | all/default wall s | trans off build s | all/default build s | build speedup | wall speedup | rel `M_avg_z` diff |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 12.02 | 2.52 | 9.31 | 1.04 | 8.94x | 4.76x | 4.8e-13 |
| 12 | 44.16 | 6.59 | 34.84 | 2.34 | 14.89x | 6.70x | 2.5e-8 |

The all/default wedge scaling lane after apply/solve profiling is:

| N | ndof | wall s | charge-basis s | H-matrix build s | CG solve s | iters | peak WSet MB |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 19328 | 2.11 | 0.24 | 0.97 | 0.62 | 22 | 559.5 |
| 12 | 64224 | 5.09 | 0.73 | 2.09 | 1.92 | 23 | 1521.4 |
| 16 | 151040 | 12.27 | 1.46 | 3.83 | 4.85 | 23 | 3346.3 |
| 20 | 293600 | 25.32 | 2.70 | 6.69 | 8.14 | 24 | 6730.9 |

Set `RADIA_HDIV_HMATVEC_STATS=1` to record `hmatvec_*` fields in the result
JSON.  The current HACApK matvec worker limits MKL to a local single thread by
default (`RADIA_HACAPK_MATVEC_MKL_THREADS=1`) because TaskManager already
parallelizes across H-matrix leaves.  On the N=20 wedge profile this reduced
the charge-Gram H-matvec part of the solve to `0.64 s`.  The next nonlinear
energy-Newton sweep showed that exact PARDISO mass-Riesz is not always the
right large-run default: switching only the inner W-CG preconditioner to the
exact diagonal of `W + N` reduced structured-hex N=16 from `247.8 s` to
`33.7 s`, completed hex N=20 in `79.9 s`, and reduced wedge N=10 from
`60.2 s` to `24.5 s`.  Tet remains size-dependent: mass-Riesz wins at
`maxh=0.20` (`4314` DOF), while the diagonal branch wins at `maxh=0.15`
(`8193` DOF).  The public `preconditioner="auto"` policy follows that split
with a `6000`-DOF tet switch point.  Result JSON records
`preconditioner_policy`; set `RADIA_HDIV_AUTO_JACOBI_TET_NFACE` to sweep the
threshold without editing the benchmark driver.

For publication-style nonlinear timing, run the drivers with
`--nonlinear-profile mdx-scaling`.  This resolves the nonlinear timing defaults
to the measured fast lane (`nl_tol=3e-4`, `newton_continuation=2`,
`newton_reuse_tangent_steps=3`, `preconditioner=auto`) while still allowing any
of those flags to be overridden explicitly.  The default profile remains
`strict` for solver-oriented regression checks.  `--newton-cg-x0` is a separate
diagnostic comparison flag; it is not part of the mdx-scaling profile until a
large-run sweep shows that warm-starting the inner CG reliably reduces the
H-matvec apply count.

For wedge, unlike tet, `gram_eps=1e-4` remained Krylov-stable in these two
checks and moved `M_avg_z` by only `2.2e-5` to `4.2e-5` relative to `1e-8`.
This is a promising timing lane for wedge cube scaling, but the tight `1e-8`
lane remains the accuracy reference until shape/IMA/nonlinear coverage is
broadened.
