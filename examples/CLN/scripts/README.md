# CLN research scripts (NGSolve / radia consumers)

Python scripts (plus their result JSONs, logs, Mathematica `.wls`, etc.)
that drive the Cauer Ladder Network research line. Migrated from
`W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/` on 2026-05-12 so the
NGSolve / `radia.*` consumers live next to the Radia source and ship with
`pip install radia`.

## Layout

| Subdir | What it is |
|---|---|
| `axifemm/` | Phase B1c/B2 P2 triangle Henrotte FE development — Python reference, NGSolve + C++ benchmarks, edge-convention probes, multi-stage debug reproducers. |
| `ngsolve_validation/` | The bulk of the CLN validation framework — 3D HCurl Hiruma, COMSOL TEAM 28 port, Kameari accumulation, Tanimoto A-T + H-H projection, DD pipeline experiments, sphere/cylinder/cuboid/A1 sweeps, Schöberl-Zaglmayr basis source, and the legacy FP64 references. |
| `tanimoto_canonical/` | M.~Tanimoto 修論 (2025) canonical 4-formulation CLN notebook set (2D H1 scalar / 3D A-T / 3D A-Phi / 3D T-Omega) on a 1 cm Cu cylinder. Reference baseline against which the ngsolve_validation/ scripts cross-check. Mirrored from `S:/NGSolve/谷本/修論/`. See `tanimoto_canonical/README.md` for formula table. |
| `multiconn_loop_method/` | Hiptmair-Ostrowski Loop Method for multiply-connected T-Ω: 1st-cohomology basis construction (`LoopField.py`) + bordered-system solver (`MatrixSolver.AddCoupling`, `SolveCoupled2`) + canonical user-facing notebook (BathPlate with Holes). Required for genus ≥ 1 conductors where the plain T-Ω matrix is singular. Sourced from `S:/NGSolve/EMPY/EMPY_Analysis/`. Algorithm doc in `multiconn_loop_method/README.md`. |

## Why these scripts live here

Per the 2026-05-12 policy decision (memory key `feedback_no_ngsolve_py_in_cln_workdir.md`):
- Every Python script that imports `ngsolve`, `netgen`, or `radia.*` belongs under `examples/`, NOT in the W:/ research working directory.
- Result JSONs, log files, and the Mathematica `.wls` scripts that pair with the Python tools came along in the same bulk move so the validation workflow stays runnable from this checkout.
- LaTeX / research notes / progress PDFs (i.e., the non-code artefacts) stay on W:/.

## Status

Most scripts are *research scripts* in the strict sense: experimental,
sometimes broken, sometimes superseded. Treat anything outside the
"canonical" entry points (`dd_full_pipeline.py`, `dd_sphere_axisym_mp.py`,
`hex_vim_cupy*.py` modulo the DD-canonical policy in memory key
`feedback_dd_canonical_policy.md`) as historical record. See the per-file
memory keys for the verified vs. deprecated lists.

## Working folder for the CLN paper / digest

The IGTE 2026 digest TeX source and Q&A MEMORY.md live in the parent
directory (`examples/CLN/`). New CLN scripts should land directly under
`examples/CLN/scripts/` rather than in W:/.
