# cubit_panels accel-magnet prune (2026-06-29)

`examples/cubit_panels/accel_magnet` was not a final examples surface.  The
active EM panel sample path is the panel sample/golden-test lane, while the
old examples directory contained earlier development probes:

- `coil_dipole.py` duplicated the shipped panel coil script
  `src/radia/panels/samples/em_sample_coil.py`.
- `experiment_occ_dipole.py` and `experiment_step_fem*.py` were exploratory
  Cubit/OCC/STEP/FEM runs superseded by `calc_accel_magnet.py`,
  `step_mesh_builder.py`, the EM notebook workbench, and the golden panel
  samples.
- `experiment_mmm_ima.py` was a diagnostic for MMM image-method analysis, not
  a live validation gate.
- `BH.txt` is embedded as `STEEL_BH` and also shipped as
  `src/radia/panels/samples/em_sample_bh.txt`.

Panel-only fixtures rescued from the examples tree should live under the
repo-root `panels/` tree.  Reusable calculation or geometry APIs belong in
`src/`; no reusable API was identified in this examples batch beyond code that
already exists in `src`.
