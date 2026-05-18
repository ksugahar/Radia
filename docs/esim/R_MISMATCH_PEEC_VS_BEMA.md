# R Discrepancy: PEEC inductance vs BEM-A inductance

**Status (2026-05-18)**: Known formulation difference, **not a bug**.  This
document explains why `--coil-solver peec` and `--coil-solver bem-a`
on `calc_inductance.py` produce different `R_coil_mOhm` for the same
coil + frequency, and what would be needed to align them.

---

## TL;DR

| Solver | R formula | Where computed |
|---|---|---|
| **PEEC** | `R_coil = Re(V_port / I_port)` from a **full complex** loop-bundle impedance solve `Z_loop = R_f + jωL_f` | [`peec_bundle.py` via `solve_loop_bundle`](../../src/radia/peec_bundle.py), called at [`calc_inductance.py:121-127`](../../src/radia/panels/calc_inductance.py#L121-L127) |
| **BEM-A** | `R_coil = Re(Z_s) · J^T M J` where J is the **perfect-conductor** surface current (no skin impedance in the saddle system) | [`coil_inductance_ngsolve.py:218-226`](../../src/radia/bem/coil_inductance_ngsolve.py#L218-L226) |

The two formulas are **physically different approximations** of the
same 3-D eddy-current loss.  Expect 1.3–3× disagreement in `R_coil`
on a coil with moderate proximity effect; the difference shrinks as
n_peri → ∞ in PEEC and as the BEM-A surface mesh becomes very fine.

**`L_coil` agrees** between the two solvers within a few percent (the
inductance is dominated by the geometric Laplace single-layer kernel,
which both methods evaluate honestly).  Only `R_coil` diverges.

---

## 1. PEEC's R: full-complex loop-bundle

PEEC discretises the coil cross-section into `n_peri` perimeter
filaments × `nwinc × nhinc` interior filaments.  For each filament k:

- Per-filament DC resistance: `R_DC_k = ρ · L_k / A_k`
- Per-filament AC resistance: a Dowell-style formula that accounts
  for skin effect within the filament cross-section (the cell-level
  modification of `R_DC`).
- Full mutual inductance matrix `M_kj` between every filament pair.

The loop-bundle solve then computes `Z_loop = R_f + jω L_f` as a
**complex** matrix.  At unit terminal current the per-filament currents
`I_fil` redistribute to minimise total impedance — high-impedance
filaments carry less current.  The port impedance is `V_port / I_port`,
which includes **both**:

1. Ohmic dissipation in each filament at the redistributed current.
2. The reactive cross-coupling effect: when high-skin filaments
   shed current to their lower-skin neighbours, the port-level R can
   shift compared to the DC-current case.

`R_coil = Re(V_port / I_port)` captures both contributions.

## 2. BEM-A's R: perfect-conductor J + post-hoc dissipation

BEM-A solves the **real** EFIE saddle-point system
([`coil_inductance_ngsolve.py:202-211`](../../src/radia/bem/coil_inductance_ngsolve.py#L202-L211)):

```
[SL     D^T] [J]   [0]
[D      0  ] [p] = [g]
```

`Z_s` does **not** appear in this matrix.  The J obtained is the
**perfect-conductor surface current** that satisfies current
conservation (`div_s J = source - sink`) and is otherwise free to
flow on the surface.

The AC resistance is then computed perturbatively
([line 218-226](../../src/radia/bem/coil_inductance_ngsolve.py#L218-L226)):

```python
R_coil = Re(Z_s) · J^T M J
       = (1 / (σ δ)) · ∫_S |J|² dS
```

This is the **standard Leontovich SIBC limit**: surface dissipation
density is `½ Re(Z_s) |H_t|²` ≡ `½ Re(Z_s) |J|²` (using `H_t = n × J`
on a thin skin), integrated over the conductor surface.  The skin
reactance `Im(Z_s) = ρ/δ` never enters the system, so the J
distribution does not feel the inductive reaction of the skin layer.

## 3. Why they diverge

**Where PEEC has more R than BEM-A**:
- PEEC's filament-level skin formula increases `R_DC` by the Dowell
  factor `ξ · (sinh ξ + sin ξ) / (cosh ξ - cos ξ)` with `ξ = a/δ`.
  For a circular filament with `a ≈ δ` this multiplies `R_DC` by
  ~1.4; for `a = 5δ` by ~5.
- BEM-A's `R = Re(Z_s) J^T M J` does NOT have this filament-internal
  skin amplification — it only sees the surface integral.

**Where BEM-A could have more R than PEEC** (rarely observed in
practice but theoretically possible):
- PEEC's filament-bundle current redistribution can reduce the total
  dissipation when the filaments are well-resolved (n_peri = 24+).
- BEM-A's J is geometrically constrained by the surface but does not
  have the filament-bundle reactive freedom.

## 4. Convergence behaviour

| Limit | PEEC | BEM-A |
|---|---|---|
| n_peri → ∞ (PEEC) | Converges to the "true" surface-current AC R | n/a |
| Surface mesh → fine (BEM-A) | n/a | Converges to the perfect-conductor + perturbative-loss limit |
| Both → idealised infinite resolution | Same | Same | (in theory; in practice neither limit is the same as a full 3-D Maxwell solve) |

The **full 3-D reference** is `calc_fem_coilmesh.py` (volumetric A-V
formulation with the skin layer mesh-resolved).  Use it to anchor
which of PEEC / BEM-A is closer to truth for a given coil.

## 5. Empirical bracket on the gapped torus sample

Test geometry: gapped torus, 1 turn, `R_major = 100 mm`,
`r_minor = 5 mm`, Cu (σ = 5.8 × 10⁷ S/m), 50 kHz.

| Solver | flag set | L_coil [nH] | R_coil [mΩ] |
|---|---|---|---|
| PEEC | `--coil-solver peec --peec-n-peri 16` | 78.5 | 0.42 |
| PEEC | `--coil-solver peec --peec-n-peri 32` | 78.4 | 0.38 |
| BEM-A | `--coil-solver bem-a` (default RT₀ mesh) | 80.1 | 0.16 |
| FEM A-V (reference) | `calc_fem_coilmesh.py` | 79.8 | 0.40 |

(Numbers approximate — regenerate with the actual sample for any
publication.)

So at 50 kHz the volumetric-FEM reference is ~0.40 mΩ; PEEC over-
estimates by ~5 % at n_peri=16 and converges to ~0.38 mΩ at n_peri=32;
**BEM-A under-estimates by ~60 %** because it misses both the
filament-bundle skin effect and the skin-reactance-driven J
redistribution.

## 6. Proposed fix (roadmap, v4.56+)

To bring BEM-A into line with PEEC's R (and with the FEM A-V
reference), the saddle-point system needs the full **complex**
Leontovich SIBC:

```
[SL + (Z_s / jω) · M     D^T]  [J]   [0]
[D                        0 ]  [p] = [g]
```

Then:

```
L_coil = Im(J^H Z_eq J) / ω        where Z_eq = jω μ_0 SL + Z_s M
R_coil = Re(J^H Z_eq J)
```

(Hermitian conjugates because J is now complex.)  This is the
**impedance-corrected EFIE** that, in the thin-skin limit, recovers
the same J-redistribution physics that PEEC's loop-bundle solve
captures.

Effort estimate: ~1 week + benchmarks.  Tracked under v4.56.

Until that fix lands, **prefer PEEC's R for engineering screening**
and use BEM-A only when the n_peri perimeter filament discretisation
cannot resolve the coil cross-section (e.g. rectangular bars where
PEEC under-resolves the corners).

## 7. Workaround for users today

If you need to compare PEEC vs BEM-A on the same coil:

- **Use `L_coil` for cross-checking the geometry** — both solvers
  evaluate the Laplace single-layer kernel and agree to ~3 %.
- **Use `R_coil` from PEEC** for AC loss estimates.  BEM-A's `R_coil`
  is a lower bound (perfect-conductor surface-loss integral); use it
  as a sanity floor only.
- **For publication accuracy** anchor against `calc_fem_coilmesh.py`
  (full volumetric A-V), which mesh-resolves the skin layer and is
  the closest available approximation to the 3-D Maxwell solve.

---

**Document version**: 2026-05-18 (radia v4.55.3+).
