# Hiruma SparseSolv benchmark

This directory is the durable validation lane for the 30 kHz, three-material
Hiruma eddy-current problem used to measure Compact AMS + COCR.  It replaces
the retired `src/ext/sparsesolv/examples/hiruma` location.

The repository tracks the smallest case, `meshes/mesh1_2.5T.vol`, as the
repeatable validation fixture.  The larger 3.5T, 4.5T, 5.5T and 20.5T meshes
are optional scaling inputs because together they add about 78 MB.  Place them
in `meshes/` and use `--all` when a full performance sweep is required.

Run the numerical gate from the repository root:

```powershell
python validation_test/sparsesolv/hiruma/bench_compact_ams.py --verify-baseline
```

The gate fixes mesh hashes and finite-element sizes exactly, while accepting a
small iteration/residual range.  Setup and solve times are recorded but are not
gated because they depend on the machine and concurrent load.  Heavy scaling
runs belong on mdx; the tracked 2.5T case is suitable as a LAB smoke run.

The fixture was converted from the existing laboratory `mesh1_2.5T.msh`
dataset.  `prepare_hiruma_vol.py` preserves its points and volume elements and
adds explicit air/material interface triangles.  The resulting `.vol` passes:

```powershell
check-vol validation_test/sparsesolv/hiruma/meshes/mesh1_2.5T.vol `
  --contract validation_test/sparsesolv/hiruma/hiruma_vol_label_contract.json
```

The assembled complex system before and after adding the interface descriptors
is identical to floating-point precision.  The source and final SHA-256 values,
mesh sizes, golden range, and one observed timing are recorded in
`compact_ams_baseline.json`; the complete mesh check is stored in
`mesh1_2.5T.vol-check.json`.

Additional comparisons remain available as `bench_ams_vs_abmc.py` and
`bench_cocr_vs_gmres.py`.  They write transient result JSON below
`C:\temp\radia-validation`, never into the repository.
