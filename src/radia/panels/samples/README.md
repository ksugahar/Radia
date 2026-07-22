# Radia application samples

Canonical sample artifacts shipped with the `radia` wheel.  Each
application exposes one or more analysis modes; each mode requires a
specific combination of input artifacts (mesh `.vol`, coil `.py`
or `.step`, BH curve `.txt`, hysteresis `.hys`, etc.).

This README is the index.  See sub-directories for deeper
documentation (`em/README.md` for the EM-specific corpus).

## Induction Heating block

| Method (UI label) | Calc script (Layer 4) | .jou recipe | Canonical .vol | Coil source |
|---|---|---|---|---|
| PEEC inductance (coil only, STEP) | `calc_inductance.py --coil-solver peec` | `ih_peec_inductance.jou` | (none — coil-only) | `ih_peec_inductance_coil.step` |
| BEM-A inductance (coil only, .vol) | `calc_inductance.py --coil-solver bem-a` | (Cubit-meshed coil .vol) | (none — coil-only) | `<coil>.vol` (source/sink sidesets) |
| PEEC + BEM weak coupling (workpiece) | `calc_inductance.py --coil-solver peec --vol <wp.vol>` | `ih_peec_bem_coarse.jou` | `ih_peec_bem_coarse.vol` | `ih_peec_bem_coarse_coil.step` |
| BEM-A + BEM weak coupling (workpiece) | `calc_inductance.py --coil-solver bem-a --vol <wp.vol>` | (Cubit-meshed coil .vol) | `ih_peec_bem_coarse.vol` | `<coil>.vol` (source/sink sidesets) |
| PEEC coil + FEM wp (SIBC) + Kelvin | `calc_fem_kelvin.py --formulation total --peec-step ...` | (reuses ih_peec_bem_coarse for wp+Kelvin) | `ih_peec_bem_coarse.vol` | `ih_peec_bem_coarse_coil.step` |
| Full simulation (FEM A-V + wp SIBC + Kelvin) | `calc_fem_coilmesh.py` | `ih_fem_kelvin_skin_fine.jou` | `ih_fem_kelvin_skin_fine.vol` | `ih_fem_kelvin_skin_fine_coil.step` |

`calc_peec_inductance.py`, `calc_peec_bem.py`, and
`calc_coil_bem_a_workpiece.py` were unified into a single
`calc_inductance.py` in v4.25.0 (2026-05); see its module docstring
for the migration history.

Golden tests:
- `tests/panels/golden/peec_inductance_torus_50kHz_Cu.json`
- `tests/panels/golden/peec_inductance_3turn_150kHz_Cu.json`
- `tests/panels/golden/peec_bem_coarse_7kHz_Cu.json`
- `tests/panels/golden/fem_coilmesh_gapped_fine_7kHz_Cu.json`

Production paths (per memory `IH panel 最終構成 2026-04-19`): the
**PEEC+BEM weak-coupling** and **FEM Full** paths are the validated
production methods, both verified to give P_wp within 1% of each
other on Cu @ 7 kHz.  The PEEC+FEM Kelvin path solves the same
physics by a different formulation but does not have its own
dedicated sample (it reuses the PEEC+BEM mesh).

Stale / research-only `ih_*` files in this directory (e.g.
`ih_bem_sample.*`, `ih_sample.vol`, `ih_fem_sample.vol`,
`ih_fem_hole.vol`, `ih_fem_kelvin_sample.vol`,
`ih_fem_kelvin_skin.vol`, `ih_closed_torus.*`) are leftovers from
demoted formulations or local test runs.  They are gitignored
(`*.vol`) so they do not enter git, and they do not ship in CI-
built wheels (CI starts from a clean checkout and regenerates only
the canonical .vol files via Cubit batch).  Local LAB checkouts
may retain them for research convenience.

## Electromagnet Simulink block

See `em/README.md` for the full C-yoke corpus (1/1 / 1/2 / 1/4
/ 1/8 reductions, ELF reference, Kelvin Benchmark mode).

Top-level canonical trio (all shipped):

| Slot | File | Role |
|---|---|---|
| Coil | `em_sample_coil.py` (`build_coil() -> CoilBuilder`) | `--coil-script` |
| Mesh | `em_sample.jou` -> `em_sample.vol` | `--vol` |
| BH | `em_sample_bh.txt` (or built-in `STEEL_BH`) | `--bh-file` |

Kelvin Benchmark mode (`Electromagnet` block / `EMDesignSpec` formulation = "Kelvin
Benchmark"): bundled `kelvin_benchmark_sphere_1_2.vol` (1/2
model) + `kelvin_benchmark_sphere_1_4.vol` (1/4 model).  See
`tests/panels/test_kelvin_benchmark_golden.py` for the golden
band (±1.5% at p=2).

## PCB PEEC Simulink block

| Sample | Role |
|---|---|
| `pcb_sample.jou` | Cubit recipe for the canonical PCB geometry |

Calc script: `calc_pcb_peec.py`.  Currently a single sample; the
application mostly targets user-supplied `.jou` files.

## Conventions

- `*.jou` files are tracked in git (canonical Cubit recipes).
- `*.vol`, `*.step`, `*.msh`, `*.sol` are gitignored
  (regenerate via the corresponding `.jou` or build script).
- The wheel ships the .jou + the .vol/.step/.bh that exist in
  this directory at build time; CI regenerates them from .jou
  before each release build.
- `_*.py`, `_*.jou`, `test_*.py`, `compare_*.py` are gitignored
  (private experimentation scripts).

## Deploy verification

The `deploy` skill's golden-range matrix exercises one canonical
sample per application mode end-to-end (see `validation_test/panels/golden/*.json`)
to lock numerical regressions.  Add a sample to `panels/samples/`
only when its golden test passes consistently.
