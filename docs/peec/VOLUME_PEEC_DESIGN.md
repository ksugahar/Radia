# Volume PEEC for Pancake Coils — Design Notes

**Status (2026-07-02)**: PROTOTYPED in `C:\temp` (research per the
"run in C:\temp" policy), NOT shipped — because it does NOT close the
gap, and the investigation indicates the gap most likely should NOT
be closed (the PEEC ~4.5 mΩ is probably correct; BEM-A's 15 mΩ is the
outlier).  See the "2026-07-02 outcome" section below.  Prior status
(2026-05-20): DEFERRED, waiting on a 4-terminal Kelvin re-measurement.

## 2026-07-02 outcome — volume PEEC prototyped, does NOT reach 15 mΩ

A cross-section-averaged volume PEEC (graded polar grid, per-cell
sub-area, full Ruehli mutual with cross-section averaging on close
pairs — a real fix over the C++ filamentary mutual, which DIVERGES
with refinement) was built and validated:

- **Isolated straight wire**: reproduces the closed-form Bessel AC
  resistance (the model is sound).  Filamentary (center-only) mutual
  diverges as the mesh refines; cross-section averaging is REQUIRED.
- **3-turn coil**: CONVERGES to ~3.7 mΩ as angular resolution
  increases (n_angular 16/24/32/48 → 3.40/3.53/3.59/3.63 mΩ,
  increments halving).  It does **not** trend toward BEM-A's 15.1 mΩ.

**Why the gap should not be closed** (round-wire proximity is
intrinsically modest):

- Analytic strong-skin proximity `P_prox/P_self = 2 a²/s² ≈ 0.29` per
  neighbour at the coil's `s/2a = 1.32` (2 mm surface gap).  Even the
  numeric 2-wire bundle gives 1.21× at that gap and only 1.30× when
  nearly touching.  A 3-turn coil (middle turn = 2 neighbours, with
  field build-up) caps at ~1.4–1.7× → coil R ≈ 4.5–5 mΩ.
- Three independent methods agree: volume PEEC 3.7, perimeter PEEC
  4.5, analytic ~4.8 mΩ.  **BEM-A's 15.1 mΩ (4.48× self-skin) is the
  3× outlier** — a 4.48× proximity factor is not physical for round
  wires at these gaps (it needs H_ext/H_self > 1.3; here it is ~0.6).
- BEM-A **is** correctly calibrated on self-skin: BEM-A on an
  isolated straight wire = 0.3165 mΩ (maxh 1.5 mm) / 0.3153 (maxh
  1.0 mm) vs Bessel 0.3148 (0.5 % / 0.16 %, mesh-stable).  So the
  coil over-estimate is coil-specific.
- **MECHANISM CONFIRMED (loss map on the R=15.144 mΩ solve):** the
  SIBC loss `Re(Zs)·∫|J|² dS` is pathologically concentrated — **top
  2 % of the surface area carries 71 % of the loss**; peak/mean loss
  density = 5767×.  The near-contact turn-gap tris (21 % of area,
  identified by cross-turn arc-length proximity < 3.5 mm) carry
  **69.7 % of the loss** at ~10⁵× the density elsewhere; the very
  top-density spikes are the source/sink injection edges + faceted-
  STEP edges.  A uniform physical skin would put ~2 % of loss in 2 %
  of area — the 71 % concentration is the **perfect-conductor J
  singularity** at near-touching surfaces + edges.  There the J
  varies on a scale far finer than the skin depth δ, so the
  Leontovich SIBC (which assumes J smooth over δ) **breaks down** and
  over-integrates |J|²; the real finite-σ current spreads over δ and
  does not concentrate this way.  Hence BEM-A's 15.1 mΩ is a
  SIBC-breakdown over-estimate, mesh-divergent (consistent with the
  2.90 → 15.14 refinement jump), NOT the physical AC resistance.
  (A full BEM-A mesh-convergence sweep was impractical — dense EFIE
  on the 634-face coil.stp > 15 min/level; a 2-wire BEM-A hit a
  singular EFIE saddle on disconnected / busbar conductors.)

**BEM-A FIX IMPLEMENTED (2026-07-02): the impedance-EFIE.**  Rather
than compute R post-hoc from the PEC current, put the surface
impedance INTO the saddle system so J is the *finite-impedance*
current and cannot over-concentrate.  The (1,1) block becomes
``j ω μ0 SL + Z_s M`` (complex); Z = Z_s (Jᴴ M J) + j ω μ0 (Jᴴ SL J),
R = Re(Z), L = Im(Z)/ω.  Wired as
``compute_inductance_source_sink(..., impedance_efie=True, omega,
Z_s_complex)`` and the panel flag ``calc_inductance.py
--bema-impedance-efie``.  Validated:

  | geometry | PEC post-hoc R | impedance-EFIE R | truth |
  |---|---|---|---|
  | isolated straight wire (smooth) | 0.3153 mΩ | 0.3153 mΩ | Bessel 0.3148 |
  | kubota 3-turn coil | 15.14 mΩ | **4.63 mΩ** | ~4.5-5 (3 methods) |

  On smooth geometry the two agree (nothing to fix); on the coil the
  impedance-EFIE removes the singular over-concentration and lands
  in the physical band — a **4th independent method** agreeing with
  volume/perimeter PEEC + analytic proximity.  Why the small
  resistive term (~1/Q of the reactance at Q≈27-90) moves R 3×: the
  singular concentration contributes almost nothing to the magnetic
  energy (JᴴSL J) but enormously to the loss (Jᴴ M J), so even a
  weak resistive penalty makes J avoid it.  Golden:
  ``validation_test/bem/test_coil_bem_a_impedance_efie.py`` (locks
  R_imp ≤ R_pec + L unchanged on the gapped-torus fixture).  It is
  **opt-in** (default off) pending a decision to flip the BEM-A R
  default (the PEC path is a known ~3× over-estimate for tightly-wound
  coils, so flipping is recommended).

**Recommendation**: keep perimeter PEEC (~4.5 mΩ) as the screening
R; do NOT trust BEM-A R for tightly-wound multi-turn coils without a
mesh-convergence check.

**The 2-terminal LCR ~15 mΩ is ALSO not a trustworthy arbiter** — and
not only because of contact/lead resistance.  At 150 kHz this coil is
reactance-dominated: `ωL = 2π·150k·430n ≈ 0.405 Ω`, so R is the small
real part of an impedance ~30–90× larger in the imaginary part
(**Q = ωL/R ≈ 27 at R=15 mΩ … 90 at R=4.5 mΩ**).  Then
`dR ≈ ωL·dδ` → a **1° loss-angle (phase) error = ~7 mΩ error in R**,
and distinguishing 4.5 vs 15 mΩ needs ~1.5° phase accuracy (or a
dissipation factor D = R/ωL of 0.011 vs 0.037 resolved to ±0.002),
which is at/below a bench LCR meter's floor at 150 kHz.  So **"BEM-A
matches the measurement" is most likely a coincidence of two
independent upward biases** (BEM SIBC-breakdown + LCR phase/contact
error) landing near 15 mΩ for different reasons — NOT mutual
validation.  The best current estimate is the first-principles
~4.5–5 mΩ; there is no solid experimental anchor yet.

**Phase-independent ground-truth routes** (in preference order):
(1) **DC 4-terminal R_dc** (robust, ~0.35 mΩ) × the simulated
`R_ac/R_dc ≈ 10–13`; (2) **calorimetric** loss at a known current
(no phase dependence — the gold standard for high-Q loss);
(3) **resonance/Q-bandwidth** method; all with open/short/load
compensation + a known ~5 mΩ reference-resistor check.  Volume
PEEC via the constant-current parallel-filament bundle also mildly
UNDER-counts the coil's along-length-varying proximity (conical
helix → nearest-neighbour direction rotates), so its 3.7 mΩ is a
lower bound; a proper volume PEEC would need a full 3-D PEEC MNA
(nodes at each station, radial current exchange) — a much larger
effort than the ~1 week estimated below, and not worth it while
BEM-A already brackets from above and perimeter PEEC is near the
true value.  Prototype: `C:\temp\volpeec_proto\`
(wire_bessel_clean.py, coil_volume_vec.py, two_wire_prox.py).

What actually shipped in v4.57.0: **perimeter PEEC + proximity-iterative
solver** ([`peec_proximity.py`](../../src/radia/peec_proximity.py)),
which lifts R from 3.68 mΩ (self-skin only) to 4.48 mΩ (proximity
factor 1.218×) on the 3-turn 150 kHz case.  Sweep shows the factor
is constant ~1.20–1.22× across 1 kHz – 1 MHz — the **structural
ceiling of perimeter PEEC**, which captures surface-Leontovich
proximity but not transverse eddy loops in the wire interior.

If the Kelvin re-measurement confirms the coil R is truly ~15 mΩ
(i.e. lead R is small) then this Volume PEEC design needs to ship.
If the re-measurement returns ~4–5 mΩ (lead R was the bulk of the
discrepancy), proximity-iterative perimeter PEEC is sufficient and
this design stays archived.

## Motivation

The current PEEC pipeline places `n_peri` filaments on the cross-section
PERIMETER only.  At 150 kHz on the 3-turn Cu pancake (a/δ = 18.5) this
gives `R_coil = 4.48 mΩ` with the v4.57.0 proximity-iterative solver
(or 3.68 mΩ with `--no-peec-proximity`).  The LCR hi-tester read
~15 mΩ; the gap is unresolved — either proximity-driven transverse
eddy currents the perimeter bundle cannot see, or lead/contact
resistance in the 2-terminal probe setup.

For a 3-turn pancake **if proximity were ~3–4×**, the "correct" R
would be ~10–15 mΩ.  To capture interior eddy loops the PEEC bundle
needs filaments that span the **interior** of the cross-section, not
just the perimeter -- because the proximity-induced eddy current is a
TANGENTIAL surface current with non-uniform azimuthal distribution
that can only be represented if the bundle's mutual L matrix sees
the asymmetric coupling between INTERIOR positions of adjacent turns.

This design is preserved against the case the Kelvin measurement
shows lead R is small.

## Why the existing `_filaments_via_coil_builder` (walker) path
does NOT work

`filaments_from_step(step, n_peri=None, nwinc=N, nhinc=M, ...)`
dispatches to `_filaments_via_coil_builder` which uses the **walker**.
The walker steps along the spine, sections the solid at each step,
and identifies cross-section corners.

CLAUDE.md / `coil_from_cad.py` doc comments:
> the walker hangs or natively crashes on multi-turn loft STEPs
> (Kubota's 3turncoil.stp: walker hangs netgen.occ > 5 min;
>  on 100号機 the subprocess exits with an unhandleable native
>  error code).

The 3-turn pancake is exactly the geometry where the walker fails.
Smoke-tested 2026-05-19: `nwinc=3, nhinc=3` on `3turnCoil_work_coil.step`
fails with `RuntimeError("seed plane has no usable cross-section")`
in `coil_from_step.py:428`.

## Proposed design: long-edge spine + volume-grid expansion

Reuse the **longest-edge UV path** (which works on 3turncoil) for
centerline extraction.  Then expand each `(station, perimeter_angle)`
sample to a `nwinc × nhinc` grid INSIDE the cross-section.

### Pipeline

```
STEP solid
  -> extract_centerline_from_step (longest-edge)        # works on 3turncoil
  -> n_stations × n_peri (perimeter) samples            # existing
  -> per-station local basis (tangent, u_hat, v_hat)    # existing
  -> for each (station, perim) point:
       expand to nwinc × nhinc grid in (u, v) plane     # NEW
  -> K = nwinc * nhinc * n_peri filaments
     each filament traces n_stations points along the conductor
  -> PEECBuilder.add_connected_segment for each grid filament
  -> build_loop_bundle_impedance → R_f (K×K diag, includes per-filament
     R_DC weighted by cross-section share), L_f (K×K mutual)
  -> solve_loop_bundle → complex bundle solve, redistribution
     captures proximity
```

### Cross-section grid placement

For a **round wire** (cross_section_kind == "circle", radius a_eq):

```
Polar grid (n_radial × n_angular):
  r_i = (i + 0.5) * a_eq / n_radial     for i in [0, n_radial)
  θ_j = 2π * j / n_angular              for j in [0, n_angular)
  offset_ij = r_i * (cos(θ_j) u_hat + sin(θ_j) v_hat)
  sub_area_ij = π a_eq² / (n_radial * n_angular)  (uniform)

K = n_radial × n_angular total filaments per cross-section.
```

For a **rectangular wire** (cross_section_kind == "rect", w × h):

```
Cartesian grid (nwinc × nhinc):
  u_i = (i + 0.5 - nwinc/2) * w / nwinc       for i in [0, nwinc)
  v_j = (j + 0.5 - nhinc/2) * h / nhinc       for j in [0, nhinc)
  offset_ij = u_i * u_hat + v_j * v_hat
  sub_area_ij = w * h / (nwinc * nhinc)

K = nwinc × nhinc total filaments per cross-section.
```

### Per-filament cross-section share for R_DC

`build_bundle_solver(filament_paths, ..., cell_wh=cell_wh)` already
handles per-filament cross-section.  Each grid filament gets:

```
cell_wh[k][i] = (sub_w, sub_h)
where:
  for cartesian: sub_w = w/nwinc, sub_h = h/nhinc
  for polar:     sub_w = sub_h = sqrt(A_sub / 1)  (equivalent square)
```

PEECBuilder uses `R = ρ * L / (sub_w * sub_h)` per segment.  Total
parallel R recovers the wire's bulk R_DC.

### Filament chain (port-to-port)

The K filaments are CONNECTED IN PARALLEL between the port nodes,
identical to the existing perimeter case.  Each filament k traces:

```
filament_paths[k] = [
    (offset_k_at_station_0, offset_k_at_station_1),
    (offset_k_at_station_1, offset_k_at_station_2),
    ...
]
```

where `offset_k_at_station_i = centerline_pos[i] + Δu_k(i) · u_hat[i]
+ Δv_k(i) · v_hat[i]`.

The local `u_hat`, `v_hat` come from the parallel-transport frame
(existing `_compute_rmf_frame`).

### API surface

New args to `filaments_from_step`:

```python
def filaments_from_step(
    step_path,
    sigma=5.8e7,
    n_peri=None,             # legacy perimeter-only count
    n_radial=1, n_angular=None,   # NEW polar grid for round
    nwinc=None, nhinc=None,  # NEW cartesian grid for rect
    ...
):
```

Dispatch:
* If `n_peri` only -> perimeter (current behaviour, preserves
  v4.55.4 self-skin Bessel result).
* If `n_radial × n_angular` given AND cross_section is round ->
  polar grid, K = n_radial × n_angular.
* If `nwinc × nhinc` given AND cross_section is rect -> cartesian
  grid, K = nwinc × nhinc.
* Mutually exclusive: per CLAUDE.md "No Fallbacks", raise if the
  user passes both n_peri and (n_radial, n_angular).

### Panel UI (calc_inductance.py)

Add to argparse:
```
--peec-n-radial INT  (round wire interior, default 1 = perimeter-only)
--peec-n-angular INT (round wire angular, default = same as n_peri)
--peec-nwinc INT     (rect wire interior width, default 1)
--peec-nhinc INT     (rect wire interior height, default 1)
```

Panel adds 4 spin boxes; existing `n_peri` becomes either the
"perimeter samples around the wire" knob (legacy) or the "angular
divisions" of the polar grid (new).  Default values reproduce the
current self-skin Bessel R for backward compatibility.

### Expected results on 3-turn pancake

Convergence study target (n_radial × n_angular):

| n_radial × n_angular | filaments | R [mΩ] | t_solve | comment |
|---|---|---|---|---|
| 1 × 16  (current)     | 16  | 3.7  | <1 s  | self-skin only |
| 2 × 16                | 32  | 6-9  | ~1 s  | partial proximity |
| 3 × 16                | 48  | 9-12 | ~2 s  | most proximity |
| 4 × 16                | 64  | ~13  | ~3 s  | near-converged |
| 5 × 16                | 80  | ~14  | ~5 s  | over-resolved? |

(Numbers are best-guess from a/δ = 18.5; will refine with actual
runs.)  Target: match measured 15 mΩ to within 20%.

### Existing infrastructure reuse

* `extract_centerline_from_step` (longest-edge UV path) -- works
  on 3turncoil, returns centerline + per-station basis.
* `_compute_rmf_frame` -- existing parallel-transport frame.
* `build_bundle_solver(filament_paths, ..., cell_wh=cell_wh)` --
  existing PEEC bundle with per-segment cross-section.
* `solve_loop_bundle(..., Zs_fil=Zs_fil)` -- existing bundle solve
  with optional per-filament SIBC.  For volume PEEC we may NOT
  need Zs_fil because the mutual L_f matrix captures skin
  redistribution naturally IF the discretization resolves the
  skin depth.  At a/δ = 18.5 we need ~37 × 37 ≈ 1400 grid points
  for true skin resolution, which is computationally expensive.
  Hybrid (coarse grid + Zs_fil) is the practical recipe.

## Implementation phases

1. **R1**: Modify `filaments_from_polyline_uv` (or equivalent) to
   accept `n_radial`, `n_angular` and emit interior grid filaments.
2. **R2**: Update `cell_wh` computation for interior cells.
3. **R3**: Add panel argparse + UI spinboxes.
4. **R4**: Golden test for proximity convergence on 3turncoil.

Effort: 2-3 days of focused work.  No new C++.

## Related

* Proximity effect physics: Dowell 1966, Ferreira 1994.
* FastHenry's volume-filament discretization (this is the same
  approach but driven by build_bundle_solver instead of nodal MNA).
* PEEC-MCS (Multi-Conductor System): same idea, more degrees of
  freedom (per-segment current variation along axis).  Not needed
  for IH coils where each turn is short relative to wavelength.

---
**Document version**: 2026-05-19 (radia 4.56.x).
