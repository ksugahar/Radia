# HDiv-MMM Cube Validation

This directory owns the current mesh-backed uniform-field cube validation for
Radia's BDM1 HDiv-VIM formulation. It replaces the retired six-face Radia
object/MMPM solver comparisons.

## Drivers

- `bench_hdiv_cube.py`: structured HEX or WEDGE, linear or nonlinear material.
- `bench_hdiv_tet_cube.py`: unstructured TET, linear or nonlinear material.
- `bench_hex_nonlattice.py`: structured-cache and non-lattice ablation.

The common problem is a 1 m cube centered at the origin with an applied field
of 200 kA/m in +z. The main drivers write JSON into this directory by default.
Useful sizes are validation workloads: run them on hibino when available, or on
mdx while its CI queue is idle.

```powershell
python validation_test\cube_uniform_field\hdiv\bench_hdiv_cube.py `
  --sizes 2 4 6 --material both

python validation_test\cube_uniform_field\hdiv\bench_hdiv_tet_cube.py `
  --maxh-values 0.4 0.3 --material both
```

For publication runs, pass `--bh-table` with the exact material table and use
`--output` to choose a durable JSON name in this directory. The result records
the mesh, solver parameters, runtime, magnetization, demagnetizing operator,
HACApK statistics, timing, and memory information.

## Committed Evidence

The mdx sweep from 2026-07-08 is preserved as raw JSON:

- `hdiv_cube_mdx_scaling_20260708.json`: structured HEX scaling.
- `hdiv_tet_cube_mdx_scaling_20260708.json`: TET scaling and ACA sensitivity.
- `hdiv_tet_ho_oneside_mdx_ab_20260708.json`: TET one-sided far-entry A/B.
- `hdiv_wedge_cube_mdx_scaling_20260708.json`: WEDGE scaling and sensitivity.
- `hdiv_wedge_oneside_mdx_ab_20260708.json`: WEDGE one-sided far-entry A/B.
- `hdiv_wedge_trans_cache_mdx_ab_20260708.json`: initial translation-cache A/B.
- `hdiv_wedge_trans_cache_all_mdx_ab_20260708.json`: full cache-scope A/B.
- `hdiv_wedge_trans_cache_fastbasis_mdx_ab_20260708.json`: fast-basis A/B.
- `hdiv_wedge_all_fastbasis_mdx_scaling_20260708.json`: fast-basis scaling.
- `hdiv_wedge_all_fastbasis_applysolve_mdx_scaling_20260708.json`: final
  apply/solve profile.

`results_hex_nonlattice.json` records the 2026-07-12 structured-cache ablation.
`evidence_manifest.json` pins every committed record by SHA-256. The imported
raw records intentionally retain their original scratch-path provenance; those
paths are descriptive and are not execution dependencies.

These records support three current observations: structured HEX is faster
than the contemporary unstructured TET lane at comparable DoF for this cube;
TET Krylov behavior requires a tighter charge-Gram tolerance at larger sizes;
and the WEDGE translation/fast-basis caches materially reduce charge-Gram build
cost while preserving the reported magnetization agreement. They are measured
validation evidence, not routine CI thresholds.
