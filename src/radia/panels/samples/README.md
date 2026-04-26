# Radia panel samples

Canonical sample artifacts shipped with the `radia` wheel.  Each
panel exposes one or more analysis modes; each mode requires a
specific combination of input artifacts (mesh `.vol`, coil `.py`
or `.step`, BH curve `.txt`, hysteresis `.hys`, etc.).

This README is the index.  See sub-directories for deeper
documentation (`em/README.md` for the EM-specific corpus).

## IH panel (`radia_ih.py`) — Induction Heating

| Method (UI label) | Calc script | .jou recipe | Canonical .vol | Coil source |
|---|---|---|---|---|
| PEEC inductance (coil only, STEP) | `calc_peec_inductance.py` | `ih_peec_inductance.jou` | (none — coil-only) | `ih_peec_inductance_coil.step` |
| Fast workpiece heating (PEEC+BEM, 1-way) | `calc_peec_bem.py` | `ih_peec_bem_coarse.jou` | `ih_peec_bem_coarse.vol` | `ih_peec_bem_coarse_coil.step` |
| PEEC coil + FEM wp (SIBC) + Kelvin | `calc_fem_kelvin.py` | (reuses ih_peec_bem_coarse for wp+Kelvin) | `ih_peec_bem_coarse.vol` | `ih_peec_bem_coarse_coil.step` |
| Full simulation (FEM A-V + wp SIBC + Kelvin) | `calc_fem_coilmesh.py` | `ih_fem_kelvin_skin_fine.jou` | `ih_fem_kelvin_skin_fine.vol` | `ih_fem_kelvin_skin_fine_coil.step` |

Golden tests:
- `tests/panels/golden/peec_inductance_torus_50kHz_Cu.json`
- `tests/panels/golden/peec_inductance_3turn_150kHz_Cu.json`
- `tests/panels/golden/peec_bem_coarse_7kHz_Cu.json`
- `tests/panels/golden/fem_coilmesh_gapped_fine_7kHz_Cu.json`

Production paths (per memory `IH panel 最終構成 2026-04-19`): the
**PEEC+BEM 1-way** and **FEM Full** paths are the validated
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

## EM panel (`radia_em.py`) — Accelerator electromagnet

See `em/README.md` for the full C-yoke corpus (1/1 / 1/2 / 1/4
/ 1/8 reductions, ELF reference, Kelvin Benchmark mode).

Top-level canonical trio (all shipped):

| Slot | File | Role |
|---|---|---|
| Coil | `em_sample_coil.py` (`build_coil() -> CoilBuilder`) | `--coil-script` |
| Mesh | `em_sample.jou` -> `em_sample.vol` | `--vol` |
| BH | `em_sample_bh.txt` (or built-in `STEEL_BH`) | `--bh-file` |

Kelvin Benchmark mode (`radia_em.py` formulation = "Kelvin
Benchmark"): bundled `kelvin_benchmark_sphere_1_2.vol` (1/2
model) + `kelvin_benchmark_sphere_1_4.vol` (1/4 model).  See
`tests/panels/test_kelvin_benchmark_golden.py` for the golden
band (±1.5% at p=2).

## PCB panel (`radia_pcb.py`) — PCB inductance

| Sample | Role |
|---|---|
| `pcb_sample.jou` | Cubit recipe for the canonical PCB geometry |

Calc script: `calc_pcb_peec.py`.  Currently a single sample; the
panel mostly targets user-supplied `.jou` files.

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
sample per panel mode end-to-end (see `tests/panels/golden/*.json`)
to lock numerical regressions.  Add a sample to `panels/samples/`
only when its golden test passes consistently.
