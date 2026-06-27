# Equivalence-theorem near-field source -- validation corpus

Schelkunoff / Love equivalence theorem (Stratton-Chu surface integral)
for the Radia + NGSolve stack.  Record EM field on a closed surface
around a source region, then re-evaluate the exterior field at any
external point (or hand off to an external MoM tool via Nastran
format).

Production module: `src/radia/equivalence_source.py` (NearFieldSource
class).  See:
- [`../../docs/equivalence_source/USER_GUIDE.md`](../../docs/equivalence_source/USER_GUIDE.md) — API usage, workflows
- [`../../docs/equivalence_source/CPP_DESIGN.md`](../../docs/equivalence_source/CPP_DESIGN.md) — C++ kernel design
- [`../../docs/equivalence_source/FMM_DESIGN.md`](../../docs/equivalence_source/FMM_DESIGN.md) — Phase D acceleration plan
- [`../../docs/equivalence_source/demos.ipynb`](../../docs/equivalence_source/demos.ipynb) — result-saved rendered showcase
- MCP tool `fem_equivalence_source` — knowledge base

## Phases

### Phase 1 -- `phase1_static_coil.py` (production demonstration)

Static magnetostatic verification:
- Source: circular current loop, radius a = 0.2 m, center at z = 0.6 m,
  current I = 1.0 A.  H field has an analytical closed form
  (Biot-Savart via complete elliptic integrals K(k), E(k)).
- Extraction surface: sphere of radius R = 0.9 m enclosing the loop.
- Workflow: analytic H on sphere -> NearFieldSource -> Stratton-Chu
  reconstruction at exterior points -> compare to analytical Biot-Savart.

Expected result:  max relative error <= 0.5 % at observation points
spanning |r| = 1 m to 5 m, dominated by the sphere triangulation
discretisation error (decreases as 1/N).

Reproduces the Sugahara Lab 2008 axi slide 5-6 verification.

### Phase 3 -- `phase3_e2e_cubit_to_sol.py` (END-TO-END, PASS)

Full production pipeline test: Cubit `.jou` -> `.vol` -> CLI -> NFS
-> CoefficientFunction -> projected GridFunction (`.sol`).

Workflow:
1. Cubit headless: create sphere R=0.5, sideset `nfs_surface`,
   `export netgen` -> `inner_mesh.vol`.
   (Falls back to NGSolve OCC if Cubit 2025.12 plugin isn't
   registered; the plugin needs rebuild against the new Cubit SDK
   but the rest of the pipeline is independent of Cubit.)
2. NGSolve: sample analytic dipole H on VectorH1 inner mesh -> `H_inner.sol`.
3. Subprocess: `calc_equivalence_source.py --vol --sol --surface
   nfs_surface --output nfs.json` -> NFS artifact (~240 KB).
4. Load NFS; build outer mesh (OCC shell R=0.6 to R=3.0);
   `NearFieldSource.project_to_h1_vector(outer_mesh)` -> projects the
   Stratton-Chu reconstruction onto a VectorH1 GridFunction ->
   `H_outer_reconstructed.sol` (~350 KB).
5. Verify: load `H_outer_reconstructed.sol`, evaluate at 8 obs
   points spanning R = 0.7 m to R = 2.5 m, compare to analytic
   dipole.

Result:  max relative error 8.34 %  (threshold 20 %)  --  **PASS**.

Error budget (per stage):
- ~2-3% inner-mesh FEM interpolation of analytic dipole
- ~1-5% Stratton-Chu integral with 1106 surface panels
- ~5-15% outer-mesh order=1 nodal projection (Set is sub-optimal
  for L2; finer mesh or higher order would tighten)

This is the canonical demonstration of the equivalence-theorem
source in the Radia stack: a closed-surface FEM solution (.sol)
becomes a reusable equivalent source artifact that can be replayed
onto an ARBITRARY second mesh as a coefficient function /
GridFunction.

### Phase 2 -- `phase2_wpt_harmonic.py` (time-harmonic)

Time-harmonic Hertzian-dipole demonstration.

Setup:
- Source: small electric dipole at origin, I*l = 1 A*m, omega = 2 pi * 1 MHz.
- Extraction surface: sphere R = 1 m.
- Workflow: analytic (E, H) on sphere -> NearFieldSource ->
  `evaluate()` at exterior obs points.

The harmonic path uses the **full dyadic Green's function**
`(I + grad-grad / k^2) psi`.  The 1 MHz deep-near-field Hertzian
dipole case now passes the 2% band; zero-analytical-H observation
points are checked with an absolute A/m threshold instead of a
singular relative error.

### `null_field_property.py` — equivalence-theorem physics check

The defining property of the equivalence theorem: when the source is
captured on a closed surface, the exterior reconstruction matches
the true field, and the interior reconstruction is ~zero (the
"null-field" property — the equivalent surface sources cancel the
real source contribution everywhere inside).

Setup: magnetic dipole at origin, sphere R = 0.30 m extraction
surface, 6 240 panels.

Result (typical run, see `results_null_field.json`):
- Exterior (|r| from 0.6 m to 1.5 m): max rel err = 0.08 % (threshold 1.5 %)
- Interior (|r| from 0.10 m to 0.20 m): |H_reconstructed| / |H_dipole_true|
  = 0.035 % (threshold 1 %)

The 1 000× contrast between exterior reconstruction and interior
null is the WHOLE POINT of the equivalence theorem.

## Running

```bash
cd S:/Radia/01_GitHub/validation_test/equivalence_source
python phase1_static_coil.py      # ~5 s   (static, analytical coil)
python phase2_wpt_harmonic.py     # ~3 s   (harmonic, Hertzian dipole)
python null_field_property.py     # ~2 s   (interior null / exterior real)
python phase3_e2e_cubit_to_sol.py # ~60 s  (full Cubit -> NGSolve e2e)
python bench_static.py            # ~30 s  (Phase A C++ speedup bench)
```

Each writes `results_<name>.json` next to the script.

Requirements: numpy, scipy (for elliptic integrals in phase 1).
NGSolve is NOT required for these phases (they validate the
Stratton-Chu integrator against analytical sources, no FEM in the
loop).  See `src/radia/panels/calc_equivalence_source.py` for the
NGSolve-coupled version (Stage 2 CLI for the future Cubit panel).

## What this is NOT

- **No IABC / SDI**.  The equivalence theorem doesn't need them.
  Radia's Kelvin transformation handles the FEM open boundary
  (stronger than IABC), and the 2015 Sugahara Lab work showed the
  reconstruction OUTSIDE is insensitive to the inner BC anyway.
- **No competing FEM solver**.  This is a POST-PROCESSING layer on
  top of whatever solver you used (NGSolve, Radia MMM, ngsolve.bem).
