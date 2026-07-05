# `radia_mcp.radia_ngsolve` — Radia + NGSolve coupled magnetostatics

**30 MCP tools** — second-largest subpackage. Covers the lab's
production Radia + NGSolve workflow: Kelvin transformation, sparse
solvers (Compact AMS via radia.sparsesolv_ngsolve), Cauer Ladder
Network reduction, PEEC inductance, ngsolve.bem, HDiv-VIM, plus
the closed-form analytical formula bank (Wakao-Igarashi-Fujiwara-
Kameari Part 1-9, cuboid average B, Bessel impedance).

## Quick start

```bash
pip install radia-mcp                # standalone, no radia C++ ext needed
mcp-server-radia-ngsolve             # stdio server
```

```
> radia_ngsolve_status()
> radia_ngsolve_topics()             # 30 topics
> kelvin_transformation(topic="bem-fem-hybrid")
> analytical_formulas(topic="cuboid_average_field")
```

## Tool families

| Family | Examples |
|---|---|
| **Kelvin transformation** | `kelvin_transformation`, `kelvin_identify_pairs`, `kelvin_periodic_bc_setup` |
| **NGSolve usage** | `ngsolve_usage`, `ngsolve_examples`, `ngsolve_recipe_<topic>` |
| **Closed-form analytical** | `analytical_formulas` (11 topics: Wakao Part 1-9, cuboid_average_field, validation_use_cases) |
| **PEEC inductance** | `peec_inductance`, `peec_filament_dispatch` |
| **HDiv-VIM** | `hdiv_vim` (Radia soft-iron demag, charge-Gram H-matrix, Reduced FEM coupling) |
| **sparsesolv** | `compact_ams_preconditioner`, `cocr_solver` (HYPRE-free, TaskManager-native) |
| **CLN SIBC orthogonal** | `cln_sibc_orthogonal` — Hierarchical Cauer SIBC (lab specialty) |
| **Mesh + I/O** | `netgen_workflow`, `vol_format_inspect`, `vol_file_lint` |
| **lint / QA** | `radia_ngsolve_lint` |
| **bibliography** | `radia_ngsolve_bibliography_index` |

Run `radia_ngsolve_status()` for the live list.

## Killer topic: `analytical_formulas` (★)

Closed-form reference layer covering:

- **Wakao-Igarashi-Fujiwara-Kameari Part 1-9** — magnetic shielding,
  plate Joule dissipation, AC thin-shell shielding,
  magnetic-shell interior fields, planar surface impedance,
  Bessel cylindrical-conductor AC impedance, plate eddy current,
  Fabri solenoid, three-phase line, ellipsoid demag/torque
- **Cuboid average field** — sympy-derived G1/G2 antiderivatives
  + 64-corner inclusion-exclusion sum (~40 µs/call, 817× faster
  than Gauss-Legendre baseline) — shipped in radia 4.22.0 C++
- **validation_use_cases** — "given this analysis X, which closed
  form is the trusted reference?" mapping

Use this **before** running any new FE simulation — gives a
ballpark + sanity check.

## Lab specialty: CLN orthogonal SIBC

The `cln_sibc_orthogonal` topic surfaces the lab's Cauer Ladder
Network surface impedance work (Kameari-Ebrahimi-Sugahara-Shindo-
Matsuo 2018 IEEE TMAG canonical paper + Hane-Nakamura 2020 dynamic
hysteresis + Hierarchical Cauer SIBC). This is the only place this
material is curated in MCP-readable form.

## radia-coupled vs radia-free

This server uses `import radia` in 4 places inside knowledge module
function bodies (not at module top). That means:

- `import radia_mcp.radia_ngsolve.server` ✅ works without radia
- `radia_ngsolve_status()` ✅ works without radia
- `radia_ngsolve_usage(topic="kelvin")` ✅ works (text only)
- Anything that exec's example code requires `pip install radia-mcp[radia]`

For pure documentation use, install with `pip install radia-mcp`
(no extras needed).

## Cross-references

- `mcp-server-mathematica` — symbolic verification of analytical
  formulas (paired killer-demo: derive Kelvin → verify in Mathematica)
- `mcp-server-fem` — FEM-formulation theory layer (A-Ω, T-Ω, H,
  Reduced potential, Darwin, MSFEM)
- `mcp-server-bem` — RWG, EFIE/MFIE, Loop-Star, Calderón, Radia
  HDiv-VIM, HACApK, FEM-BEM hybrid
- `mcp-server-matrix-solvers` — Krylov + preconditioner theory
- `mcp-server-mor` — Cauer Ladder Network model-order reduction
- `mcp-server-peec` — PEEC filament/panel, FastHenry, HOIBC
- `mcp-server-differential-forms` — k-form / de Rham math foundations

## Source

- `src/radia_mcp/radia_ngsolve/server.py` — tool registration
- `src/radia_mcp/radia_ngsolve/knowledge/` — per-topic knowledge
- `src/radia_mcp/radia_ngsolve/bibliography_index_knowledge.py` —
  auto-generated bibliography
