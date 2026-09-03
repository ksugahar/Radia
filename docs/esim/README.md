# ESIM Documentation Index

**ESIM** = **E**ffective **S**urface **I**mpedance **M**ethod for
nonlinear ferromagnetic eddy-current analysis.  Radia's implementation
combines a 1-D radial cell solver (`esim_cell_problem.py`) with three
outer-solver paths (scalar BIE, FEM-Kelvin, FEM-coilmesh) and an outer
Picard fixed-point ("Karl iteration") on the surface impedance Z_s.

This folder holds the implementation documentation.  Read order depends
on the reader's goal:

## Entry points by goal

| You want to ... | Start with |
|---|---|
| ... use the CLI tools (calc_inductance.py / calc_fem_kelvin.py / calc_fem_coilmesh.py) | [`USAGE.md`](USAGE.md) |
| ... **run an example or reproduce a paper figure** | [`EXAMPLES.md`](EXAMPLES.md) |
| ... understand the math / cite the paper | [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md) |
| ... understand the code architecture / Karl loop internals | [`IMPLEMENTATION.md`](IMPLEMENTATION.md) |
| ... see numerical validation tables | [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) |
| ... understand the honest limits of ESIM validation (and why) | [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) § 1b |
| ... diagnose a Karl loop that didn't converge | [`IMPLEMENTATION.md`](IMPLEMENTATION.md) § 3.4 + [`plot_karl_history.py`](../../validation_test/ih_esim_benchmark/plot_karl_history.py) |
| ... understand WHY this combination of choices (the publication argument) | [`SCALAR_BIE_VS_VECTOR_BEM.md`](SCALAR_BIE_VS_VECTOR_BEM.md) |
| ... diagnose the PEEC vs BEM-A coil-R discrepancy | [`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) |
| ... see the slow STEP-to-filament path analysis | [`PEEC_PERFORMANCE_AND_R_ANALYSIS.md`](PEEC_PERFORMANCE_AND_R_ANALYSIS.md) |

## Document table

| Document | Purpose | Length |
|---|---|---|
| [`README.md`](README.md) | this index | small |
| [`EXAMPLES.md`](EXAMPLES.md) | catalog of runnable ESIM scripts under `examples/` + paper-reproducer | ~ 100 lines |
| [`USAGE.md`](USAGE.md) | CLI guide for the three calc_*.py scripts | ~ 200 lines |
| [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md) | PDE → weak form → quadrature → curvature handling → HCurl weak form → Lorentz reciprocity φ·B form → complex-µ | ~ 970 lines |
| [`IMPLEMENTATION.md`](IMPLEMENTATION.md) | three-solver dispatch + Karl loop + per-DOF Z_s + Lipschitz/Anderson + failure modes + perf chars | ~ 540 lines |
| [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) | Bessel reference + 3-path consistency + curve-order benchmark + **per-element vs scalar headline** + Karl convergence rate | ~ 690 lines |
| [`SCALAR_BIE_VS_VECTOR_BEM.md`](SCALAR_BIE_VS_VECTOR_BEM.md) | **why scalar BIE + curved Tri6 + per-element ESIM is the right combination** — error-rate match, Piola argument, Karl per-iter cost comparison, reviewer Q&A | ~ 200 lines |
| [`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) | focused diagnosis of why PEEC vs BEM-A coil R differs | ~ 130 lines |
| [`PEEC_PERFORMANCE_AND_R_ANALYSIS.md`](PEEC_PERFORMANCE_AND_R_ANALYSIS.md) | PEEC R + STEP→filament slowness deep-dive | ~ 190 lines |

## Quick reference cards

### When does ESIM beat linear SIBC?

| Workpiece | Use |
|---|---|
| Cu / Al / brass (µ_r ≈ 1) | linear SIBC (Dowell) |
| Steel, |H_t| < BH knee | linear SIBC with reasonable µ_r |
| Steel, |H_t| straddles BH knee | **ESIM** + BH curve |
| Steel, strong spatial |H_t| contrast (≥3×) | **ESIM `--esim-per-panel`** |

### Headline numerical result

For steel cylinder at 50 kHz, I_port = 100 A driven through the BH
knee: **per-DOF ESIM reports P_wp = 18.75 W vs scalar uniform ESIM's
30.51 W (-38.5 %)** because local saturation changes the surface
impedance at hot-spot DOFs that the uniform model averages away.

Reproduces from `radia >= 4.67.0` via `calc_inductance.py
--esim-per-panel`.  Full reproducer in [§6b of CROSS_VALIDATION.md](CROSS_VALIDATION.md).

### MCP access

The MCP server `mcp-server-ih` exposes this knowledge via the
`ih_esim(topic="...")` tool.  Topics include `overview`, `bh_file`,
`inductance_cli`, `fem_kelvin_cli`, `fem_coilmesh_cli`, `per_element`,
`convergence`, `json_output`, `troubleshooting`,
`scalar_vs_vector_bem`, `headline_numbers`.  Implementation:
[`packages/radia-mcp/src/radia_mcp/ih/esim_knowledge.py`](../../packages/radia-mcp/src/radia_mcp/ih/esim_knowledge.py).

## Cross-references outside docs/esim/

- [`docs/research/bem_numerics/BEM_SIBC_ESIM_RESEARCH.ipynb`](../research/bem_numerics/BEM_SIBC_ESIM_RESEARCH.ipynb) — research WIP (LAB-only, gitignored).
- [`docs/research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.ipynb`](../research/bem_numerics/NONLOCAL_SIBC_BILICZ_2023.ipynb) — wide-band nonlocal extension roadmap (deferred).
- [`docs/induction_heating/induction_heating_demo_showcase.ipynb`](../induction_heating/induction_heating_demo_showcase.ipynb) — result-saved public ESIM/WPT/RWG demo showcase.
- [`docs/induction_heating/induction_heating_examples_catalog.ipynb`](../induction_heating/induction_heating_examples_catalog.ipynb) — result-saved promotion catalog for the IH docs/API/validation migration.
- [`validation_test/ih_esim_benchmark/`](../../validation_test/ih_esim_benchmark/) — benchmark scripts producing `results.json`.
- [`src/radia/esim_cell_problem.py`](../../src/radia/esim_cell_problem.py) — cell-problem solver source.
- [`src/radia/panels/calc_inductance.py`](../../src/radia/panels/calc_inductance.py) — scalar BIE-SIBC path (path A in [`SCALAR_BIE_VS_VECTOR_BEM.md`](SCALAR_BIE_VS_VECTOR_BEM.md)).
- [`src/radia/panels/calc_fem_kelvin.py`](../../src/radia/panels/calc_fem_kelvin.py) — FEM-Kelvin path (path C).
- [`src/radia/panels/calc_fem_coilmesh.py`](../../src/radia/panels/calc_fem_coilmesh.py) — FEM-coilmesh path (path D).

**Document version**: 2026-05-30 (radia v4.67.0+ dense-sweep baseline).
