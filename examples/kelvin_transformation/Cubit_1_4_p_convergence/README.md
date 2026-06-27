# Cubit 1/4 Sphere + Kelvin: p-convergence (VERIFIED)

p-convergence demonstration for the Cubit-meshed `export netgen`
+ NGSolve Kelvin transformation pipeline, on a magnetic sphere
(mu_r = 100) in uniform Hz background.

## Result (verified 2026-04-25)

| order | ne    | ndof   | Hz_origin   | error vs analytical | time |
|------:|------:|-------:|------------:|--------------------:|-----:|
| 1     | 9806  | 2176   | 0.031703    | **+7.79%**          | 0.3s |
| 2     | 9806  | 15131  | 0.029621    | **+0.71%**          | 0.5s |
| 3     | 9806  | 48673  | 0.029571    | **+0.54%**          | 2.3s |

Analytical: `Hz_origin = 3 / (mu_r + 2) * H0 = 0.029412 A/m`

p=2 nails the analytical to <1%, confirming Cubit + Kelvin works
correctly at high order.

## Geometry (1/4 sector: x>=0, y>=0, full z)

- Magnetic sphere at origin, radius `a = 0.05` m, `mu_r = 100`
- Air outer at radius `R = 0.20` m (`kelvin_int` boundary)
- Kelvin sphere at `(0, 0, 0.6)`, radius `R` (`kelvin_ext` boundary)
- BCs: `sym_bn=0_x` (x=0 plane, Natural for Omega),
       `sym_bn=0_y` (y=0 plane, Natural for Omega),
       `GND` (Kelvin sphere centre, Dirichlet Omega = 0)

## Two non-obvious fixes

These are the practical lessons.  The geometry-only debug took several
iterations because each issue was silent (no exception, just a wrong
number that converges to a wrong answer).  Both lessons apply to ANY
Cubit-meshed Kelvin / Omega-Reduced setup.

### 1. `subtract A from B keep` does NOT carve in Cubit

Verified by post-merge probe (`get_relatives("surface", sid, "volume")`)
showing only ONE parent volume per shared surface.  The `keep` flag
preserves A but the minuend B is left unchanged.

**Fix**: drop the `keep` flag, re-create the subtrahend as a fresh
sphere afterwards:

```python
cubit.cmd("create sphere radius R")           # air
cubit.cmd("create sphere radius a")           # mag (subtrahend)
cubit.cmd("subtract volume mag from volume air")   # NO keep -- mag consumed
cubit.cmd("create sphere radius a")           # re-create mag in the cavity
```

After this, imprint+merge correctly identifies the shared mag-air
interface and the netgen `FaceDescriptor` gets `domin = air_index,
domout = mag_index` (two-sided).

### 2. Cubit FD orientation differs from OCC after `export netgen`

Probing both meshes:

|              | OCC `sphere`              | Cubit `sphere`            |
|--------------|---------------------------|---------------------------|
| domin        | 1 (air_inner)             | 2 (air)                   |
| domout       | 2 (magnetic)              | 1 (magnetic)              |
| `+normal`    | OUT of air INTO mag       | OUT of air INTO mag       |

The `+normal` direction is the SAME in both meshes (out of air into
mag).  But `specialcf.normal` returns a normal whose sign depends on
the surface element's local orientation in the mesh, which Cubit and
OCC fix differently.

**Empirically validated rule**: Cubit-meshed `sphere` Neumann
correction needs `-specialcf.normal(mesh.dim)` (same sign as the
OCC reference).  The archived classic reference
`Omega_ReducedOmega/Sphere/3D_sphere_with_Kelvin.py` (preserved in
`docs/kelvin/kelvin_classic_demos_results.json`) already uses
`-specialcf.normal`, so Cubit-meshed solvers should follow the same
sign convention.

## Files

| File | Role |
|------|------|
| `mesh_and_export.py` | Cubit Python: build geometry, mesh, sideset/block, copy-mesh Kelvin (1/4 reduction), export `.vol` at p=1, 2, 3 |
| `solve_p_convergence.py` | NGSolve Python: Omega-Reduced Omega + Kelvin solver, probe Hz at origin, compare to analytical |
| `p_convergence.json` | JSON results (last run) |

## Run

```bash
python mesh_and_export.py --orders "1,2,3"
python solve_p_convergence.py --orders "1,2,3"
```

## Why 1/4 and not 1/8

We attempted 1/8 (x>=0, y>=0, z>=0) first and ran into TWO separate
issues, the first now FIXED:

### Issue 1 (FIXED): copy-mesh anchor was non-deterministic for 1/8

The 1/8 spherical cap has 3 boundary arcs of EQUAL length (all
quarter great-circles of length pi*R/2).  The original
`_add_kelvin_cubit_reduction` used `max(curves, key=length)` to pick
the copy-mesh anchor curve, but the tie-break depended on Cubit's
internal listing order which could pick DIFFERENT arcs on source
(air's cap) vs target (kelvin's cap).  Result: a ~120-deg rotational
mis-projection that left ~5% of vertex pairs ~2 cm off the
translation copy.

**Fix (commit 2026-04-25)**: deterministic geometric predicate.
Pick the arc with the smallest (centroid_z, centroid_y, centroid_x)
tuple.  For the +octant convention this picks the z=0 arc on both
source and target, which is the translation-equivalent curve.  After
the fix, 1/8 vertex pairs match to machine precision (2.6e-16 m,
143/143 paired).  1/4 case is unaffected (was already 5.6e-16 m).

### Issue 2 (deferred): solver gives Hz_origin ~ 0 for 1/8

With the mesh now perfect, the Omega-Reduced Omega solver still
gives Hz_origin ~ -3.8e-4 (analytical 0.029, error -101%) for 1/8
reduction.  Initial diagnosis was the singular reluctivity
`Mu = mu0 * (R/r')^2` near r'=0 lying on the kelvin_far Dirichlet
boundary, but **2026-04-26 sweeps disproved this**: capping `Mu` at
`mu_0` everywhere (i.e. neutralising the Kelvin transformation
entirely) STILL gives -101% on 1/8 -- yet 1/2 (+0.34%) and 1/4
(-0.02%) work fine without the Kelvin reluctivity, because the
Periodic identification + sym BCs alone do most of the open-
boundary work.

So the 1/8 failure is structural to the corner-GND geometry,
not the (R/r')² singularity per se.  The user's observation that
H1 only needs ONE Dirichlet point at r'=0 (since GND already plays
that role) is theoretically correct and was applied to
`calc_kelvin_benchmark.py` (kelvin_far dropped from the Dirichlet
list); this frees ~250 DOFs but does not close the gap on its own.

Open hypotheses for 1/8: corner-discretisation on coarse mesh,
sym_ht=0_z combined with point GND, Periodic identification under
offset_dir="x".  See `memory/feedback_kelvin_1_8_blocker.md`.

### Recommendation

For new Cubit-Kelvin samples that need symmetry reduction, default
to **1/4** with offset_dir along a free axis perpendicular to all
reduction axes -- this avoids the kelvin_far singular plane entirely
and is the production-validated path with p=2 matching analytical
to <1%.

## What this verifies

- **Cubit `export netgen` produces a valid Kelvin-Periodic
  `.vol`** at any order (1, 2, 3 all tested).
- **NGSolve `Periodic` H1 FES correctly slaves high-order DOFs**
  using the Cubit-written point identifications + mesh topology.
- **Omega-Reduced Omega + Kelvin solver gives the analytical answer
  at p=2 to <1% error** on this Cubit mesh.
- **The detect_kelvin_offset / add_periodic_kelvin pipeline in
  `calc_common.py` works end-to-end** for Cubit meshes.
