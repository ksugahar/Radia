# Equivalence-theorem near-field source -- examples

CST Microwave Studio "Near-Field Source" equivalent for the Radia +
NGSolve stack.  Record EM field on a closed surface around a source
region, then re-evaluate the exterior field at any external point via
the Stratton-Chu surface integral (or hand off to an external MoM
tool via Nastran format).

Production module: `src/radia/equivalence_source.py` (NearFieldSource
class).  Knowledge reference: MCP tool `fem_equivalence_source`.

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

### Phase 2 -- `phase2_wpt_harmonic.py` (time-harmonic, KNOWN LIMITATION)

Time-harmonic Hertzian-dipole demonstration.  **This test EXPECTEDLY
fails the 2% threshold by ~3x in deep near-field**: the current
`evaluate()` uses the SCALAR Green's-function Stratton-Chu form
which is FAR-FIELD ACCURATE but misses the (1/k^2) grad-grad psi
correction term of the full DYADIC Green's function.  In the deep
near-field of a 1 MHz Hertzian dipole (R_sphere = 1 m << lambda
~ 300 m, R_obs = 10 m = lambda/30), the undershoot is ~factor 3.

Setup:
- Source: small electric dipole at origin, I*l = 1 A*m, omega = 2 pi * 1 MHz.
- Extraction surface: sphere R = 1 m.
- Workflow: analytic (E, H) on sphere -> NearFieldSource -> evaluate()
  at 10 m on z-axis.

**Roadmap**: implement `evaluate_dyadic()` with the full
(I + grad grad / k^2) dyadic Green's function kernel for accurate
near-field harmonic reconstruction.  Until then, the harmonic path
is suitable for FAR-FIELD reconstruction (R_obs >> lambda), not
deep near-field.

The static path (Phase 1, Radia's primary use case) works correctly
out of the box.

## Running

```bash
cd S:/Radia/01_GitHub/examples/equivalence_source
python phase1_static_coil.py      # ~5 s
python phase2_wpt_harmonic.py     # ~3 s
```

Each writes `results_phaseN.json` next to the script.

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
