# Docs Cruft Prune 2026-06-28

## Lesson

Bulk promotion ledgers and copied legacy source trees do not belong in `docs/`
once their lessons are captured elsewhere. `docs/` should show maintained
theory, result-bearing notebooks, synchronized JSON, and canonical APIs; it
should not become a warehouse for development chronology.

## Pruned

- `docs/examples_consolidation/`: batch-by-batch migration scratch space for
  "next 100" example sweeps. Keep those decisions in memory notes and topic
  notebooks, not as public docs.
- `docs/stream_function/stream_function_examples_archive.*`: source-only topic
  archive that preserved deleted TODO benchmark stubs. The maintained docs now
  live in `docs/stream_function/{theory,regularization,deformation,benchmarks}`.
- `docs/hdiv_vim/vim_examples_archive.*`: full-source inventory of the
  `examples/vim` corpus. The public docs keep `README.md`, productionization
  notes, and result-bearing showcase notebooks; the validation corpus remains
  in `examples/vim` + `validation_test/feec`.
- Topic `*_examples_archive.*` triples for topics that already have maintained
  result-bearing notebooks or validation surfaces. Public docs should point at
  those maintained artifacts, not source-only archive ledgers.
- `docs/electric_machine/electric_machine_validation_archive.*`: validation
  source inventory superseded by `validation_test/electric_machine/` plus the
  result-bearing cogging/skew docs notebook.
- Second-pass archive triples removed for
  `build123d_netgen_gmsh_flow`, `clebsch_hodograph`, `clebsch_legendre`,
  `cubit_mesh_export`, `fem_readable`, `rf_waveguide`, and `visualization`
  after their maintained docs/API/validation surfaces stopped referencing the
  source ledgers.
- `docs/kelvin/legacy_assets/kelvin_transformation/`: old examples mirror and
  debug notes. The one live document, the Kelvin convention, was promoted to
  `docs/kelvin/CONVENTION.md`; old source-level history remains recoverable
  from git if needed.

## Rule

When a docs artifact is only a source archive, migration ledger, failed
attempt log, or old path mirror, prune it before promoting more examples.
Promote the distilled rule/API/result instead.

## IH Simulink Native-Only Contract (2026-07-23)

The production induction-heating interface is the masked Simulink application
block containing `radia_ih_eddy_sfun` and `radia_ih_thermal_sfun`. A heat/power
LUT, lumped thermal plant, or generic discrete state-space block is not an
independent validation route and must not be restored as an IH example, test
fixture, RL environment, or fallback.

The retired MATLAB surfaces were `buildIHControlModel`, `makeIHPlant`, the IH
power/heat LUT constructors and evaluators, the waveform simulators, and the
lumped `makeIHEnvironment` adapter. Generic HCurl reduced state-space support,
the TEAM28 CLN table, and the temperature-dependent BH material block remain
separate canonical capabilities.

The fast regression route uses the actual library block with a synthetic
distributed native operator. It checks current-squared scaling, asymmetric
multi-cell rotation, weighted thermal conservation, temperature feedback,
temperature-dependent operator updates, MAT/JSON loading, fail-fast config
validation, and singular-operator lifecycle recovery. Thermal state advances
only in `mdlUpdate`; advancing it in `mdlOutputs` can apply `dt` more than once
when Simulink reevaluates outputs at one simulation time.
