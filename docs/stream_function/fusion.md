# SF for fusion -- stellarator Stage-2 (NESCOIL / REGCOIL / FOCUS)

The same surface stream function this framework uses for MRI-gradient and
induction-heating coils **is** the object the stellarator-coil community calls
the **winding-surface current potential** -- the unknown of NESCOIL
(Merkel 1987), REGCOIL (Landreman 2017), and the surface part of FOCUS
(Zhu 2018). The "Stage-2" coil problem is mathematically identical to the
gradient/IH design:

```
given a target NORMAL field  B.n  on the PLASMA boundary,
find a current potential psi on a WINDING SURFACE around it whose
Biot-Savart  n.B  reproduces it          (iso-contours of psi = the coils)
```

The design-matrix rows are just the plasma-normal component `n.B` of the
winding-surface Biot-Savart kernel the SF designer already assembles
(`A3` = the 3-component field at the plasma points; `A_n = einsum('mc,mcj->mj',
plasma_normal, A3)`). There is **no** fusion-specific solver code -- the
gradient-coil machinery designs a stellarator coil unchanged.

Two demos cover eight parts:

| demo | parts |
|------|-------|
| [`demo_regcoil_fusion.py`](examples_catalog.ipynb) | forward map  /  REGCOIL L-curve  /  net current  /  VMEC boundary |
| [`demo_regcoil_fusion_advanced.py`](examples_catalog.ipynb) | coil force/stress  /  real li383 wout  /  FOCUS standoff  /  FOCUS winding-shape |

All numbers below are reproduced by the demos and locked by
`tests/panels/test_streamfunction_golden.py::test_regcoil_fusion_*`.

---

## 0. radia's SF solver vs NESCOIL / REGCOIL / FOCUS -- a superset

radia ships the **same object** as these codes (the winding-surface current
potential `psi`, `K = n x grad psi`) through
[`radia.stream_function`](../../src/radia/stream_function.py) +
[`calc_streamfunction.py`](../../src/radia/panels/calc_streamfunction.py). The
honest positioning is **design at parity, deliverable + physics a strict
superset**:

| axis | radia-SF | NESCOIL / REGCOIL / FOCUS |
|------|----------|---------------------------|
| **Design** (`B.n` / current potential) | yes -- `B.n` ~2e-9 on producible targets | yes (the reference) |
| **Scale / method** | dense Biot-Savart **ACA**-compressed (H-matrix) + **ridge-TSVD** pseudo-inverse on **any meshed** surface | dense Fourier least-squares on a parameterised torus |
| **Optimise** | fast (~50 us) folded-TSVD re-solve -> **multi-objective Pareto** over `alpha` **and** winding-surface **shape** (+ sheet-metal deform) | single-`alpha` L-curve |
| **Deliverable** | `psi` -> single-stroke wire -> sheet-metal distort -> **STEP CAD** -> **PEEC** circuit model | stop at `psi` |
| **Physics** | **iron** yoke/shield/core (Kelvin-FEM `M = M_free + M_react`, `--iron-vol`) | free-space (vacuum) only |

So **REGCOIL is a special case** (vacuum, toroidal, design-only) of radia's SF:
design at parity, deliverable + iron a superset, and a complete win for
single-conductor coils (MRI gradient/shim, induction heating -- the same
current-potential math). The design-parity claim is exactly what sections 1-4
and A-D below demonstrate; the four extra axes are:

- **Scale / method -- ACA + ridge-TSVD.** The design matrix is the dense
  Biot-Savart coupling compressed by ACA (an H-matrix) and inverted by a ridge
  (Tikhonov) truncated-SVD pseudo-inverse ([api.md](api.md),
  [regularization.md](regularization.md)). The ridge **is** REGCOIL's `lambda`,
  but the ACA compression + the FE-direct surface make it a *scalable*
  regularised least-squares on an **arbitrary meshed** winding surface, not a
  dense Fourier least-squares on a parameterised torus.
- **Optimise -- fast TSVD -> Pareto + surface deformation.** The TSVD core is
  folded once, so each regularised re-solve is a ~50 us `k x k` core solve
  (`RegularizedTSVD`, `demo_pareto_tikhonov_aca.py`). That cheap re-solve makes a
  whole front affordable, so radia sweeps a multi-objective Pareto over **both**
  the regularisation `alpha` (misfit-vs-energy, plus an L-inf IRLS
  misfit-vs-peak-current front) **and** the winding-surface **shape**
  (`geometry`/geom_scale lever + the sheet-metal `deform`) -- richer than a
  single-`alpha` L-curve. See [regularization.md](regularization.md) +
  [deformation.md](deformation.md).
- **Deliverable -- design-to-manufacture.** NESCOIL/REGCOIL/FOCUS stop at `psi`.
  radia continues: contours -> single-stroke wire (grad-`psi` winding
  orientation, so `l>=2` saddle shims chain without the common series current
  cancelling) -> sheet-metal distort -> STEP CAD (OCC `WriteStep`) -> PEEC. The
  PEEC step is a full circuit-extraction **solver** (`L, R, C, M` + SPICE
  netlist, MMM coupling), not just an inductance number. See
  [single_stroke.md](single_stroke.md).
- **Physics -- iron.** These codes are all free-space (vacuum Biot-Savart).
  radia's material-aware kernel (Kelvin-FEM DtN transfer `M = M_free + M_react`,
  the `--iron-vol` path) designs coils **with** iron yoke/shield/core -- domains
  REGCOIL structurally cannot enter.

**Honest nuance.** The single-stroke win is decisive for single-conductor coils.
Stellarator **modular** coils are intentionally separate coils, so "one wire" is
not their manufacturing step; there radia still adds STEP CAD + PEEC `L` beyond
REGCOIL's `psi`, and its distort is analogous to FOCUS filament-shape
optimisation. Open gaps: no SIMSOPT/STELLOPT Stage-1+2 integration; a
head-to-head at community mode counts is unmeasured.

### Earned by measurement (three goldens, Repository-First)

| golden | what it locks | measured |
|--------|---------------|----------|
| (a) [`test_regcoil_parity_deliverable_golden.py`](../../tests/panels/test_regcoil_parity_deliverable_golden.py) | vacuum **parity** + the deliverable REGCOIL lacks | `B.n` rel **4.9e-9**, STEP **954 kB**, PEEC `L` **3.09 uH** from one run |
| (b) [`test_regcoil_iron_differentiator_golden.py`](../../tests/panels/test_regcoil_iron_differentiator_golden.py) | the **iron** differentiator | free-space **misses >20 %**, material-aware **hits <1e-2** (ratio >10x) |
| (c) [`test_streamfunction_manufacture_e2e_golden.py`](../../tests/panels/test_streamfunction_manufacture_e2e_golden.py) | **manufacture end-to-end** | target `x` -> `psi` **0.042 %** -> single-stroke wire **0.197 %** -> distort **0.175 %** (no regress) -> STEP -> PEEC `L` **81.8 uH** |

The driver for (a) is
[`demo_regcoil_parity_deliverable.py`](examples_catalog.ipynb)
(it reuses the `demo_regcoil_fusion` helpers + `calc_streamfunction`'s
`_write_step_polylines` / `_peec_inductance`); (b) reuses the Kelvin-DtN
`act8_03_general_iron_design` bridge; (c) runs the production
`calc_streamfunction.py --method manufacture ... --distort --step-output --peec`
on the cylinder + DSV fixture.

---

## 1. The forward map is exact

Two **producible** targets -- a uniform vertical field (a PF / equilibrium /
vertical-field coil) and a non-axisymmetric `sin(theta) cos(2 phi)`
(stellarator-like shaping) -- are reproduced to **`B.n` residual ~2e-9**
(machine precision). This is honest: the targets are smooth and within the
winding surface's reach, so the forward map is exact, not hidden behind a
tolerance.

## 2. The REGCOIL trade-off (L-curve)

On a genuinely hard target `sin(3 theta) cos(5 phi)` -- a high mode that decays
across the plasma-coil gap, so it is **not** cheaply producible -- sweeping the
regularisation weight `alpha` traces the classic REGCOIL L-curve:

- large `alpha` -> smooth coil, high `B.n` residual;
- small `alpha` -> low residual, but the peak surface current density
  `|grad psi|` saturates at the surface's representation limit (the vertical leg
  of the L-curve), knee at `alpha_rel ~ 2e-2`.

`(field error, coil complexity)` is exactly the Stage-2 trade-off REGCOIL is
built to expose.

## 3. Net current -- the multivalued / secular term

A single-valued `psi` carries **zero** net current through each hole of the
winding torus. A real coil set carries net current; the full current potential
gains a **secular** term

```
Psi = psi + (G/2pi) zeta + (I/2pi) theta          (zeta toroidal, theta poloidal)
```

whose two extra degrees of freedom are the first cohomology generators of the
winding surface.

- **Their count is topology.** The surface's first Betti number `b1 = 2`
  (genus 1) is **confirmed gmsh-free** via the Euler characteristic of the
  triangulated torus (`b1 = 2 - chi = 2g`).  The *volume* T-Omega cohomology
  CUT engine [`src/radia/cohomology_cut.py`](../../src/radia/cohomology_cut.py)
  is itself gmsh-free (the pure-Python `radia.cohomology` combinatorial Hodge
  Laplacian); on the torus the surface generators are analytic: `grad(zeta)`,
  `grad(theta)` are single-valued vector fields even though the angles are
  multivalued.
- **The TF secular field is verified.** `K_zeta = n x grad(zeta)` is the
  net-**poloidal**-current (TF) sheet; it reproduces the textbook toroidal field
  `B_tor * R = const` inside the tube (Ampere's `1/R`, to **0.2 %**) and `~0`
  outside (no enclosed poloidal current).
- **Key physics -- net current is a PRESCRIBED parameter.** The TF field is
  **tangent** to the plasma boundary, so its `B.n` footprint is **>1000x
  smaller** than the net-toroidal generator's. The net poloidal current is
  therefore *not* fitted from `B.n` -- it is set for the desired on-axis toroidal
  field (1 T at R = 0.30 m needs 1.5 MA). This is exactly why REGCOIL takes
  `net_poloidal_current` as an **input**, not a fitted output.

## 4. A VMEC-shaped plasma boundary

The circular plasma torus is replaced by a non-axisymmetric **rotating-ellipse**
boundary in the VMEC Fourier representation
`R = sum RBC(m,n) cos(m th - n NFP ph)`, `Z = sum ZBS(m,n) sin(...)` (NFP = 3),
with analytic surface normals from the parametric tangents. `--wout wout_*.nc`
reads a real equilibrium's boundary (netCDF4 `rmnc`/`zmns`/`xm`/`xn`/`nfp`);
**stellarator-symmetric only** -- a non-symmetric `lasym=T` file (with
`rmns`/`zmnc`) is **rejected**, not silently truncated.

---

## A. Coil force / stress

A surface current `K` experiences a Lorentz force per unit area

```
f = K x B_avg,        B_avg = 1/2 (B+ + B-)
```

where `B_avg` is the average of the coil's own field on the two sides of the
current sheet (the tangential `B` jumps by `mu0 K` across it; the `+/-eps`
normal average cancels the self-singularity). For the net-poloidal-current (TF)
coil this is the classic magnetic stress: **`|f|` equals the magnetic pressure
`B_tor^2/(2 mu0)`** (ratio ~0.99) and **concentrates on the INBOARD leg**
(~5x the outboard, where `B_tor ~ 1/R` is largest) -- exactly why a tokamak /
stellarator TF coil is inboard-stress-limited.

*Honest scope:* this is the magnetic force per unit area (the stress **driver**,
N/m^2), not a structural hoop-stress model of the conductor.

## B. A real VMEC equilibrium (li383)

The demo designs against the boundary of a genuine free-boundary equilibrium --
the **li383** (NCSX-like, NFP = 3, quasi-axisymmetric) reference `wout` shipped
by [simsopt](https://github.com/hiddenSymmetries/simsopt) (MIT). `--wout PATH`
uses any VMEC output; with no path the demo fetches li383 (a 121 kB netCDF) to a
cache. The Stage-2 forward map reproduces the vertical-field `B.n` on the
genuine 25-mode (m in [0,3], n in [-3,3], R0 = 1.378 m) stellarator boundary to
**~4e-8**. (That a real, strongly-shaped boundary is accepted at all validates
the winding-normal logic: the outward orientation is decided by a majority vote
against the boundary's own `m=0` centre, not a hardcoded major radius.)

## C. FOCUS winding-standoff study

FOCUS / REGCOIL also optimise the **winding surface**. Sweeping the
winding-surface standoff (gap to the plasma) shows that a closer winding surface
couples better, so it reproduces the target with **lower** `B.n` residual **and**
lower peak current density `|grad psi|` (coil complexity) -- the curve is
**monotonic** (~50x over the sweep). The distance optimum is therefore
**constraint-bound**: push the winding surface to the minimum engineering
standoff `d_min` (set by blanket / access / neutron shielding). A bounded
minimisation of coil complexity lands on `d_min`.

## D. FOCUS winding-shape (the core FOCUS contribution)

The deeper FOCUS lever is the winding-surface **shape**, not just its distance.
`_surface_mesh_from_grid` builds an NGSolve surface mesh from an **arbitrary**
`(theta, phi)` point grid (manual netgen `Element2D` + `FaceDescriptor`, two
triangles per quad, periodic in both directions). A manual-mesh SF design
reproduces `B.n` to machine precision (2.6e-10), matching the OCC-revolved
torus -- so the winding surface can be **conformal** to a shaped plasma, not just
a circular torus.

For an elongated (`kappa = 2`) plasma, blending the winding surface from
**circular** (encloses the Z tips -> a large gap at the R sides) to
**conformal** (the plasma offset by the standoff along its normals -> a uniform
gap), at the **same** minimum standoff, cuts coil complexity `peak |grad psi|`
by **~34 %** (circular 2.0e6 -> conformal 1.3e6) at the same `B.n` quality. The
winding **shape**, not only its distance, is a real design lever -- which is
exactly FOCUS's point.

*Honest scope:* the study sweeps a 1-parameter circular->conformal blend; a full
Fourier-mode winding-surface optimiser is the named next step. (The naive
normal-offset can over-concentrate the surface at high-curvature tips, so the
true optimum is *near*-conformal and grid-sensitive -- the committed claim is the
robust "conformal << circular", locked at >15 % reduction.)

---

## Honest scope (summary)

- Parts 1-2 use a single-valued `psi` (correct for PF / RMP / shaping / shim
  fields); part 3 adds the multivalued secular term -- the generator **count** is
  computed gmsh-free (Euler characteristic), the generators are analytic **on the
  torus**, and a general winding surface takes them from the gmsh-free
  `cohomology_cut.py` (`radia.cohomology`).
- Force is the magnetic force per area, not a structural stress.
- We do **not** run VMEC here (no simsopt/desc installed) -- the default
  rotating-ellipse is an analytic model; `--wout` drops in a real equilibrium,
  and the reader is round-trip-verified against the wout schema.
- The winding-shape study is a 1-parameter blend; a full Fourier-mode
  winding-surface optimiser and a structural hoop-stress model are the remaining
  next steps.

## References

- P. Merkel, "Solution of stellarator boundary value problems with external
  currents", *Nucl. Fusion* **27**, 867 (1987) -- NESCOIL.
- M. Landreman, "An improved current potential method for fast computation of
  stellarator coil shapes", *Nucl. Fusion* **57**, 046003 (2017) -- REGCOIL.
- C. Zhu et al., "New method to design stellarator coils without the winding
  surface", *Nucl. Fusion* **58**, 016008 (2018) -- FOCUS.

MCP: `streamfunction("fusion")` (radia-streamfunction server).
