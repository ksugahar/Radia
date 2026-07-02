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

> **FOLLOW-UP 2026-07-02 (volume-PEEC + proximity-physics probe --
> the discrepancy now points the OTHER way; see
> [`VOLUME_PEEC_DESIGN.md`](../peec/VOLUME_PEEC_DESIGN.md)
> "2026-07-02 outcome").**  A cross-section-averaged **volume PEEC**
> (validated: reproduces Bessel on an isolated wire) CONVERGES to
> **~3.7 mΩ** on this coil (n_angular 16→48: 3.40/3.53/3.59/3.63),
> NOT toward 15.  Round-wire proximity is intrinsically modest:
> analytic `2a²/s² ≈ 0.29`/neighbour and a numeric 2-wire bundle
> 1.21× at the coil's 2 mm gap → a 3-turn coil caps at ~1.4–1.7× →
> **~4.5–5 mΩ**.  Three methods agree (volume 3.7, perimeter 4.5,
> analytic 4.8); **BEM-A 15.1 (4.48×) is the 3× OUTLIER** and 4.48×
> is unphysical for round wires at these gaps.  BEM-A self-skin IS
> calibrated (isolated wire = Bessel to 0.5 %, mesh-stable).
> **MECHANISM CONFIRMED by a loss map on the 15.144 mΩ solve:** the
> SIBC loss `Re(Zs)·∫|J|²dS` is pathologically concentrated — **top
> 2 % of the surface area = 71 % of the loss** (peak/mean density
> 5767×); the near-contact turn-gap tris (21 % area) carry 69.7 % of
> the loss at ~10⁵× the density elsewhere, plus source/sink injection
> + faceted-edge spikes.  That is the perfect-conductor **J
> singularity** at near-touching surfaces/edges, varying far below
> the skin depth δ where the Leontovich SIBC breaks down — so BEM-A
> over-integrates |J|² and the 15.1 mΩ is a SIBC-breakdown,
> mesh-divergent over-estimate (matches the 2.90→15.14 jump), NOT the
> physical R.  The 2-terminal LCR ~15 mΩ is **also not a trustworthy
> arbiter**: at 150 kHz the coil is reactance-dominated (Q = ωL/R ≈
> 27–90), so R is a tiny real part and a **1° phase error ≈ 7 mΩ** —
> distinguishing 4.5 vs 15 mΩ is at/below a bench LCR's phase floor,
> on top of contact/lead R.  So "BEM-A matches the measurement" is
> most likely a coincidence of two independent UPWARD biases (BEM
> SIBC-breakdown + LCR phase/contact error), NOT mutual validation.
> **Ground truth needs a phase-independent measurement** — DC
> 4-terminal R_dc × simulated R_ac/R_dc, calorimetric loss, or
> resonance/Q (see VOLUME_PEEC_DESIGN.md).  Treat the "BEM-A is the
> correct tool / 15 mΩ is right" conclusion below as SUPERSEDED for
> tightly-wound round-wire coils.
>
> **BEM-A FIX SHIPPED (default ON): the impedance-EFIE.**  Putting Z_s into
> the saddle system ((1,1) = jωμ0 SL + Z_s M, complex) makes J the
> finite-impedance current so it cannot over-concentrate; on the coil
> this gives **4.63 mΩ** (vs the PEC path's 15.14), a 4th independent
> method in the ~4.5-5 mΩ band, and reproduces Bessel on a smooth
> wire.  The panel/CLI uses this path unconditionally for AC BEM-A;
> the legacy PEC post-hoc comparison path was removed rather than kept
> as a compatibility mode.
> (`compute_inductance_source_sink(..., omega, Z_s_complex)`).  Golden:
> `validation_test/bem/test_coil_bem_a_impedance_efie.py`.  See
> VOLUME_PEEC_DESIGN.md "2026-07-02 outcome".

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

## TL;DR (current behaviour, 2026-07-02+)

| Solver | What R actually is | Where computed |
|---|---|---|
| **PEEC** | **Full Bessel `Z_cyl(ω)`** for a round-wire bundle.  Per-filament `Zs_fil_k = n_peri · (Z_cyl(ω) − R_DC_per_m) · L_k` is added to the loop-bundle diagonal so the bundle impedance reaches `Z_cyl(ω) · L_filament` at all frequencies (R_DC at ω=0, SIBC asymptote at high ω). | [`calc_inductance.py:_solve_coil_peec`](../../src/radia/panels/calc_inductance.py#L93) → [`peec_bundle.solve_loop_bundle(..., Zs_fil=Zs_fil)`](../../src/radia/peec_bundle.py#L240) |
| **BEM-A** | **Impedance-EFIE SIBC** (sole formulation): `Z_s = (1+j)/(σδ)` sits inside the complex saddle `jωμ0·SL + Z_s·M`; `R = Re(Z_s)·(Jᴴ M J)` on the finite-impedance J, `L = μ0·(Jᴴ SL J)` (external).  The legacy PEC post-hoc `R = Re(Z_s)·JᵀMJ` was REMOVED (over-estimated ~3× on tightly-wound coils; see §2). | [`coil_inductance_ngsolve.py`](../../src/radia/bem/coil_inductance_ngsolve.py) `compute_inductance_source_sink` |

For a round wire both solvers agree to within a few % across the full
frequency range, and on tightly-wound multi-turn coils they now agree
too (kubota 3-turn: BEM-A 4.63 vs PEEC 4.48 mΩ).  The remaining gap
is mesh / quadrature.

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

## 2. BEM-A's R: the impedance-EFIE (Z_s inside the saddle)

Since 2026-07-02 BEM-A solves the **complex impedance-EFIE** saddle
on the coil surface (`coil_inductance_ngsolve.py`, sole formulation):

```
[jω μ0 SL + Z_s M   D^T] [J]   [0]
[D                   0 ] [p] = [g]        (complex, ω > 0)
```

with the complex Leontovich `Z_s = (1+j)/(σ δ)` in the (1,1) block.
The recovered J is the **finite-impedance surface current** -- it does
NOT over-concentrate at near-contact gaps or edges.  Then

```python
R_coil = Re(Z_s) · (Jᴴ M J)          # Leontovich dissipation
L_coil = μ0 · (Jᴴ SL J)              # EXTERNAL inductance (geometry)
```

(the internal surface reactance `Im(Z_s)·(Jᴴ M J)/ω` regularizes J but
is deliberately not folded into the reported L).  At `ω = 0` the solve
reduces to the real vacuum saddle `[SL, D^T; D, 0]` with R = 0.

**Historical (removed 2026-07-02)**: the original implementation
solved the real perfect-conductor saddle (no `Z_s` in the matrix) and
evaluated `R = Re(Z_s)·JᵀMJ` post-hoc from the PEC current.  On smooth
geometry that agreed with the impedance-EFIE and with Bessel; on
tightly-wound coils the PEC J concentrates singularly at near-contact
turn gaps and edges (varying below δ, where the Leontovich integral
breaks down) and over-estimated R ~3× (kubota 3-turn: 15.14 vs
4.63 mΩ; a loss map put 71 % of the loss on 2 % of the area).  The
path was deleted per the Discard-the-PoC / No-Fallbacks policies.

In the **strong-skin limit** (`a / δ >> 1`) BEM-A recovers the round-
wire Bessel asymptote to within mesh quadrature error.  In the
**weak-skin limit** (`a / δ << 1`) the Leontovich SIBC itself loses
validity and R → 0 as `δ → ∞`; BEM-A does NOT pick up R_DC (use PEEC
or frequency=0 + analytic R_DC there).

## 3. Why PEEC and BEM-A agree

PEEC's per-filament `Zs_fil_k = n_peri (Z_cyl(ω) − R_DC/m) L_k` gives
a parallel-bundle self-impedance of `Z_cyl(ω) · L`.  The real part is
the full Bessel AC resistance.

BEM-A's `Re(Z_s) (Jᴴ M J)` on the impedance-EFIE current reproduces
the same Bessel asymptote on smooth wires (0.16-0.5 % measured), and
-- unlike the removed PEC path -- stays physical on tightly-wound
coils: kubota 3-turn 4.63 mΩ vs perimeter PEEC 4.5 / volume PEEC 3.7 /
analytic proximity 4.8 mΩ.

In the weak-skin / DC limit PEEC stays at R_DC (correct), BEM-A drops
to 0 (Leontovich SIBC ≠ DC).

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
- **Complex-Z_s EFIE for BEM-A — SHIPPED 2026-07-02 and now the SOLE
  formulation** (see §2): the shipped (1,1) block is
  `jω μ0 SL + Z_s M` (equivalent to the `(Z_s/jω)·M` proposal up to
  the overall `jω μ0` scaling), J redistributes under the full
  Leontovich impedance, R = Re(Z_s)(Jᴴ M J), and L keeps the
  external-inductance convention `μ0 (Jᴴ SL J)`.  The legacy PEC
  post-hoc path was removed in the same change.

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
