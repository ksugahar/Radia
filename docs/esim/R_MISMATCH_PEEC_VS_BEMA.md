# R Discrepancy: PEEC inductance vs BEM-A inductance

**Status (2026-05-20, proximity-aware)**: PEEC now (a) injects per-filament
Bessel self-skin `Zs_fil` (v4.55.4, 2026-05-19), and (b) augments it
with the Leontovich surface dissipation evaluated from the actual
Biot-Savart H field at each filament's wire-surface position via
`solve_proximity_iterative` (v4.57.0, 2026-05-20, default ON).

On the 3-turn pancake (Cu, n_peri=16, 150 kHz):

| Path | R_coil | L_coil | Notes |
|---|---|---|---|
| pre-2026-05-19 (R_DC only) | 0.3945 mΩ | 426.30 nH | wrong — DC at all f |
| 2026-05-19 (Bessel self-skin) | 3.6752 mΩ | 430.14 nH | round-wire AC asymptote |
| **2026-05-20 (+ proximity)** | **4.4793 mΩ** | **431.32 nH** | **production default** |
| LCR hi-tester measurement | ~15 mΩ | — | see §6 |

**Update 2026-07-02 (3-turn measurement case reopen)**: on a
converged BEM-A run of the SAME 3-turn coil geometry, BEM-A returns
**R = 15.14 mΩ** matching the LCR hi-tester ~15 mΩ.  So the
4.48 → 15 mΩ gap is NOT lead/contact resistance -- BEM-A reproduces
the measurement without leads.  The gap is a real PEEC vs BEM-A
physics discrepancy for tightly-packed multi-turn coils.

The structural ceiling of perimeter PEEC + proximity iteration is
~1.2× the self-skin value and is a **formulation ceiling, not an
algorithm-parameter ceiling**.  Investigated 2026-07-02:

- Increased ``n_peri`` from 16 to 32 to 64 (more parallel perimeter
  filaments): R stays 4.3--4.5 mΩ, no trend toward 15 mΩ.
- Increased eval-point density per filament from 1 (centroid) to 32,
  128, 256 axial samples: same 4.4 mΩ saturation.
- Switched from Zs_fil-embedded iteration to direct dissipation
  summation (``R = R_bessel + ∫|H_prox|² / |I_port|²``): same 4.4 mΩ.
- Varied the near-source exclusion (source segments within
  ``R_exclude_factor × wire_radius`` of an eval point are dropped
  to avoid double-counting the self-skin baked into Bessel):
  factor 1.0 gives 4.3 mΩ, factor 8.0 gives 3.7 mΩ.  Peak is 4.3 mΩ.
- Shifted the eval point radially outward from the wire surface
  by ``shift_factor × wire_radius``.  Motivated by an
  independently-verified observation: at ``r = a`` exactly on the
  discrete filament shell, Biot-Savart returns ``I/(4πa)`` --
  precisely half of Ampere's ``I/(2πa)`` due to the surface-current
  jump condition.  Naively this suggested a factor-2 miss in the
  self-skin contribution.  Investigation:

    * Calibrated on an isolated straight wire: ``shift_factor = 0.13``
      gives ``R_from_H / R_bessel_analytic = 1.028`` (correct).
    * At the calibrated ``shift_factor = 0.13`` on the 3-turn coil:
      R = 4.55 mΩ, ratio 1.24× Bessel (same as the current
      centroid-based iter, no improvement).
    * At ``shift_factor = 0.05`` on the 3-turn coil: R = 11.98 mΩ
      (looks close to 15) -- but the same shift on an isolated wire
      gives ``R_from_H / R_bessel = 3.0`` (3x over-count).  The
      "improvement" is an artefact of near-source singularity in the
      Biot-Savart when eval is too close to the surface filaments;
      it does not reflect real proximity physics.

Root cause: at the wire surface the discrete-filament Biot-Savart
sum from N line currents distributed around the perimeter yields
tangential H on the order of 1300 A/m at N=16, versus Ampere's
``I/(2πa) = 5052 A/m`` for the isolated wire.  Ampere's limit is
approached only as N → ∞ (log-divergently).  For N in the practical
range (16--64), the surface H that Biot-Savart produces is
fundamentally under the surface H that BEM-A's EFIE saddle J
distribution produces on the same surface -- the two formulations
are not equivalent for surface impedance dissipation.

**Correct tool for this problem class**: BEM-A (surface EFIE
saddle) captures the physics that perimeter-filament PEEC
structurally cannot.  Volume PEEC (deferred, see
[`VOLUME_PEEC_DESIGN.md`](../peec/VOLUME_PEEC_DESIGN.md), ~1-week
effort) would add radial filaments to represent the transverse eddy
loops inside the wire cross-section; combined with the existing
perimeter filaments this may close the gap.

To get the self-only Bessel R/L for cross-checks: pass
`--no-peec-proximity` to `calc_inductance.py`.

---

## TL;DR (current behaviour, 2026-05-19+)

| Solver | What R actually is | Where computed |
|---|---|---|
| **PEEC** | **Full Bessel `Z_cyl(ω)`** for a round-wire bundle.  Per-filament `Zs_fil_k = n_peri · (Z_cyl(ω) − R_DC_per_m) · L_k` is added to the loop-bundle diagonal so the bundle impedance reaches `Z_cyl(ω) · L_filament` at all frequencies (R_DC at ω=0, SIBC asymptote at high ω). | [`calc_inductance.py:_solve_coil_peec`](../../src/radia/panels/calc_inductance.py#L93) → [`peec_bundle.solve_loop_bundle(..., Zs_fil=Zs_fil)`](../../src/radia/peec_bundle.py#L240) |
| **BEM-A** | **AC SIBC**: `R = Re(Z_s) · J^T M J = (1/(σδ)) · ∫_S \|J\|^2 dS` where J is the perfect-conductor surface current from the real-valued EFIE saddle. | [`coil_inductance_ngsolve.py:218-226`](../../src/radia/bem/coil_inductance_ngsolve.py#L218-L226) |

For a round wire both solvers should now agree to within a few %
across the full frequency range.  The remaining gap is mesh / quadrature.

## TL;DR (historical, pre-2026-05-19, before the fix)

| Solver | What R returned | Issue |
|---|---|---|
| **PEEC** | **`R_DC = ρL/A`** at every frequency | `solve_loop_bundle` called without `Zs_fil`; the Bessel / Dowell skin term was never injected. |
| **BEM-A** | AC SIBC (correct) | OK |

The earlier (pre-2026-05-19) version of this document had the
interpretation **reversed**: it claimed PEEC included Dowell and that
BEM-A under-estimated.  That was wrong.  See the "What changed"
section at the bottom.

---

## 1. PEEC's R: Bessel skin-corrected loop-bundle (current behaviour)

PEEC discretises each filament as a series chain of segments.  Each
segment carries a DC resistance from the C++ `PEECBuilder`:

```python
# build_bundle_solver -- per segment
builder.add_connected_segment(node_a, node_b, w_i, h_i, sigma=sigma)
# C++ builder fills topology['R'] with rho * L_i / A_i (DC).
```

`build_loop_bundle_impedance` collapses to filament-level diagonal R:

```python
# peec_bundle.py:216
R_f[k, k] = sum_{i in seg_of_filament[k]} R_dc[i, i]
```

`calc_inductance.py:_solve_coil_peec` then **injects the Bessel skin
contribution** as `Zs_fil`:

```python
# calc_inductance.py:_peec_skin_impedance_per_filament
a_eq    = _equivalent_wire_radius_m(topo, n_peri)         # round-wire radius
Z_ac    = cylinder_ac_impedance(a_eq, sigma, omega)        # Bessel, per unit length
R_dc_m  = cylinder_dc_resistance(a_eq, sigma)              # ρ/(πa²)
dZ_per_m = Z_ac - R_dc_m                                   # skin contribution only
for k:
    Zs_fil[k] = n_peri * dZ_per_m * L_filament_k

# calc_inductance.py:_solve_coil_peec → solve_loop_bundle:
I_fil, V_port = solve_loop_bundle(R_f, L_f, args.frequency,
                                  I_port=args.current,
                                  Zs_fil=Zs_fil)
```

`solve_loop_bundle` then builds `Z_fil[k,k] = R_f[k,k] + jωL_self_k +
Zs_fil[k]`, and the parallel-of-n_peri equivalence makes the bundle
self-impedance reach `Z_cyl(ω) · L_filament` exactly — both real and
imaginary parts.

**Frequency response on the in-repo 3-turn coil** (n_peri=16, Cu,
3turnCoil_work_coil.step):

| frequency | `a/δ` | `R_coil` [mΩ] | `L_coil` [nH] |
|---:|---:|---:|---:|
| 100 Hz   | 0.48  | 0.355 | 463.2 |
| 1 kHz    | 1.51  | 0.398 | 459.8 |
| 150 kHz  | 18.5  | 3.675 | 430.1 |

The 150 kHz value matches the high-skin analytic `R_AC ≈ R_DC·a/(2δ) ≈
3.98 mΩ` (Bessel asymptote) within ~8 %; the L drop from 463 → 430 nH
is the internal-inductance shrinking (`μ_0/(8π)` at DC → 0 at high
ω, integrated over 0.78 m wire length ≈ 39 nH).

## 2. BEM-A's R: AC SIBC on perfect-conductor J

BEM-A solves the **real-valued** EFIE saddle system on the coil
surface ([`coil_inductance_ngsolve.py:202-211`](../../src/radia/bem/coil_inductance_ngsolve.py#L202-L211)):

```
[SL     D^T] [J]   [0]
[D      0  ] [p] = [g]
```

`Z_s` does NOT appear in this matrix.  The recovered J is the
**perfect-conductor surface current**.  The AC resistance is then
computed post-hoc
([line 218-226](../../src/radia/bem/coil_inductance_ngsolve.py#L218-L226)):

```python
R_coil = Re(Z_s) · J^T M J
       = (1 / (σ δ)) · ∫_S |J|^2 dS
```

This is the standard Leontovich SIBC surface integral: dissipation
density `½ Re(Z_s) |H_t|^2` (with `H_t = n × J`) integrated over the
conductor surface.

In the **strong-skin limit** (`a / δ >> 1`) BEM-A recovers the round-
wire Bessel asymptote to within mesh quadrature error.  In the
**weak-skin limit** (`a / δ << 1`) BEM-A → 0 because `Re(Z_s) → 0` as
`δ → ∞`; it does NOT pick up R_DC.

## 3. Why PEEC and BEM-A now agree

PEEC's per-filament `Zs_fil_k = n_peri (Z_cyl(ω) − R_DC/m) L_k` gives
a parallel-bundle self-impedance of `Z_cyl(ω) · L`.  The real part is
the full Bessel AC resistance.

BEM-A's `Re(Z_s) ∫|J|² dS = ρL/(δ·P)` is the high-skin asymptote of
the same Bessel formula.

In the high-skin limit (the IH regime, `a/δ ≥ 3`) the two formulas
give the same R within mesh and `n_peri` discretisation error.  In
the weak-skin / DC limit PEEC stays at R_DC (correct), BEM-A drops to
0 (an artefact of SIBC ≠ DC).

## 4. Empirical numbers on the 3-turn pancake coil

Test geometry: 3-turn pancake, R_avg ≈ 35 mm, r_wire = 3.15 mm, two
lead bars 60 mm each at y = ±12.5 mm, Cu (σ = 5.8 × 10⁷ S/m), 150 kHz.

Analytic anchors at 150 kHz:

```
δ              ≈ 0.171 mm
a / δ          ≈ 18.5      (strong skin)
R_DC = ρL/A    ≈ 0.355 mΩ  (using the PEEC-resolved path length ~ 0.643 m;
                            the 0.43 mΩ figure in older docs assumed a
                            shorter "naïve" 0.78 m chain.)
R_AC ≈ R_DC · a/(2δ)  ≈ 3.28 mΩ  (round-wire high-skin asymptote)
```

Measured by each solver:

| Solver | `L_coil` [nH] | `R_coil` [mΩ] | comment |
|---|---|---|---|
| Analytic R_DC | — | 0.355 | reference |
| Analytic R_AC (Bessel a/(2δ)) | — | 3.28 | reference |
| **PEEC** (`n_peri=16`, after 2026-05-19 fix) | 430.1 | **3.675** | full Bessel; ~12 % above the simple a/(2δ) asymptote because Bessel `J_0/J_1` has higher-order terms |
| **BEM-A** (.vol sample) | 244.0\* | 2.90 | AC SIBC, surface mesh under-converged |
| FEM A-V (`calc_fem_coilmesh.py`) | TBD | TBD | true 3-D reference |

\* The PEEC L (430 nH) and BEM-A L (244 nH) differ because the BEM-A
path consumes a pre-meshed `_peec_bem_bem.vol` whose surface mesh
does not match the PEEC STEP-derived filament geometry.  That is a
separate mesh / topology issue.

**Pre-2026-05-19 historical PEEC** (for comparison only):

| Solver (historical) | `L_coil` [nH] | `R_coil` [mΩ] | what was wrong |
|---|---|---|---|
| PEEC pre-2026-05-19 | 421.8 | 0.379 | R was R_DC at every frequency — `Zs_fil` not passed |

## 5. Roadmap / open items

- **Non-round cross-sections** (rectangular bars, ribbons): the
  Bessel `cylinder_ac_impedance` falls back to "equivalent circle"
  via mean cross-section area.  This is OK for square-like
  cross-sections (within ~10 %) but errs for high-aspect-ratio
  ribbons.  A proper Dowell formula for rectangular bars is in
  [`dielectric_solver._apply_dowell_correction`](../../src/radia/dielectric_solver.py)
  but not wired into `_solve_coil_peec` yet.  Once wired, choose by
  `topo['cross_section_kind']` (round → Bessel; rect → Dowell).
- **ESIM nonlinear steel filaments**: when filaments are themselves
  ferromagnetic with `μ_r(H)` dependence, the per-filament Z_s should
  come from the ESIM cell solver, not the linear Bessel formula.
  Branch on a `--filament-impedance esim` flag (future).
- **(Future) Complex-Z_s EFIE for BEM-A**: add `(Z_s / jω) · M` to
  the saddle SL block so J redistributes under the full Leontovich
  impedance.  Recovers the J-redistribution physics inside the BEM-A
  formulation.  Effort: ~1 week.

## What changed 2026-05-19

The earlier version of this document (2026-05-18) had the physics
interpretation backwards:
- "PEEC's filament-level skin formula increases R_DC by the Dowell
  factor" → wrong, no skin formula was applied; `Zs_fil` was not
  passed by `calc_inductance.py`.
- "BEM-A under-estimates R by ~60 %" → wrong, BEM-A was correctly
  applying the SIBC integral; PEEC was sitting on R_DC and was the
  one off by ~10 ×.
- "Prefer PEEC's R for engineering screening" → wrong, the opposite
  was true.

Same day (2026-05-19), the PEEC code was fixed: `_solve_coil_peec`
now computes `Zs_fil_k = n_peri (Z_cyl(ω) − R_DC_per_m) L_k` and
passes it to `solve_loop_bundle`.  PEEC and BEM-A now agree on R to
within mesh / `n_peri` discretisation error in the high-skin regime;
PEEC additionally recovers R_DC correctly at low frequency where
BEM-A goes to zero.

---

**Document version**: 2026-05-19 fixed (radia v4.55.4+).
