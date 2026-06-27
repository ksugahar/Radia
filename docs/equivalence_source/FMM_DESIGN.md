# Acceleration design — Equivalence-theorem near-field source via NGSolve.bem

**Status**: design only.  Implementation deferred until a real
use case demands `N_face × N_obs > 1e9` direct pairs.

**Key architectural commitment**: any acceleration of the
equivalence-source one-shot evaluator goes through **NGSolve.bem**
existing Multilevel Expansion (FMM-equivalent) machinery, NOT a
Radia-vendored FMM library.  Consistent with `CLAUDE.md`:

- "Complement NGSolve, don't compete with it"
- "Do NOT Implement... custom H-matrix algorithms, mesh generation
  wrappers, full-wave BEM, ..."

## 1. What NGSolve.bem already provides (audited 2026-05-26)

`from ngsolve.bem import ...` exposes the full BEM stack on
NGSolve 6.2.2603:

### 1.1 Multilevel (FMM-equivalent) kernel summation

| Symbol | Purpose | Kernel |
|---|---|---|
| `BiotSavartCF` | Direct Biot-Savart from filaments / surface current | `(1/(4π)) ∇(1/R) × J` |
| `BiotSavartRegularMLCF` | **Multilevel-accelerated** Biot-Savart far-field | regular-expansion FMM |
| `BiotSavartSingularMLCF` | Multilevel-accelerated Biot-Savart near-field | singular-expansion FMM |
| `RegularExpansionCF` / `SingularExpansionCF` | Single-level FMM expansion of any kernel | generic |
| `RegularMLExpansion(CF)` / `SingularMLExpansion(CF)` | Multilevel hierarchical FMM | generic |

### 1.2 Layer-potential operators (kernel summation)

| Symbol | Kernel | Use |
|---|---|---|
| `LaplaceSL` / `LaplaceDL` | `1/(4πR)` and its normal-derivative | Phase A static H, the `(n·H) ∇(1/R)` term |
| `HelmholtzSL` / `HelmholtzDL` | `e^{-jkR}/(4πR)` scalar | Phase B scalar Stratton-Chu term |
| `HelmholtzSingleLayer/DoubleLayerPotentialOperator` | scalar Helmholtz layer | Phase B  |
| `HelmholtzHypersingularOperator` | normal-derivative of DL | Phase B (n·E grad-psi term) |
| `MaxwellSingleLayer/DoubleLayerPotentialOperator` | vector Helmholtz dyadic | **Phase B dyadic correction** (closed-form, no `∇∇` issues) |
| `MaxwellSingleLayerPotentialOperatorCurl` | `curl(SL)` of Maxwell SL | E/H interchange |

The Maxwell SL/DL operators internally use the full dyadic Green's
function `Ḡ_e = (I + ∇∇/k²)ψ` — exactly what we hand-coded in
`rad_equivalence_source.cpp::EvaluateHarmonic`.  NGSolve.bem already
solves the FP-precision issues at low `k` via the standard EFIE/MFIE
weak-form trick.

### 1.3 What NGSolve.bem does NOT provide

- `HMatrix` / H-matrix compression — **not in 6.2.2603**.  Acceleration
  is via Multilevel Expansions (FMM), not via algebraic H-matrix
  compression.  Consistent with the user's "FMM not H-matrix" call:
  one-shot evaluation needs kernel-summation acceleration, not
  matrix-recompression.

### 1.4 Lucy Weggler's stabilized EFIE — out of scope here, in scope elsewhere

NGSolve.bem ships **Lucy Weggler's product-space EFIE stabilization**
(`HDivSurface × SurfaceL2` formulation with `κ² V_κ` re-scaling) for
low-frequency BEM **solves** — see the lab notebook
`to_developers/ngsolve/low_freq_efie_ngbem_applications.ipynb`
(Sugahara 2026-02) and the upstream demo at
[Weggler/docu-ngsbem `Maxwell_DtN_Stabilized.ipynb`](https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb).

Weggler's stabilization fixes the **operator conditioning** that
breaks down as `κ → 0`.  It is the right tool for PEEC-style
circuit extraction, scattering, and any BEM **solve** in the low-
frequency regime — i.e. the world of `ngsolve.bem` BilinearForm +
GMRES.

For the equivalence-source one-shot evaluator the Weggler
stabilization is **out of direct scope** (we don't solve a system);
its conceptual sibling — Laplace-ML routing of the kernel summation
— addresses the FMM-side LF pathology instead.  See §3.4 for the
distinction.

## 2. Cost model (unchanged from previous draft)

Break-even (FMM > direct C++ kernel): **`N_face × N_obs > 10⁹`**
(≈5 s direct walltime on 8 cores).

| N_face | N_obs | Direct (Phase A C++) | NGSolve.bem ML |
|---|---|---|---|
| 10⁴ | 10⁴ | 0.5 s | 0.01 s |
| 10⁵ | 10⁵ | 50 s | 0.5 s |
| 10⁶ | 10⁶ | 1.4 h | 50 s |

Today's examples all sit at `N_face × N_obs ≤ 10⁶`.  The direct C++
kernel is faster than spinning up a NGSolve `BilinearForm` for that
size.  ML wins only at larger problems.

## 3. Mapping equivalence_source → NGSolve.bem

The Phase A kernel:

    H_A(r) = (1/(4π)) ∮ [∇(1/R) × J_s - (n·H) ∇(1/R)] dS'

decomposes into:

- **Term 1** `∇(1/R) × J_s`: a Biot-Savart of the surface current `J_s = n × H_s`.
  → `ngsolve.bem.BiotSavartRegularMLCF(j_s_grid_function, mesh)`
  evaluated at `obs_points` — gives ML-accelerated H from surface
  current, with the `1/(4π)` and the curl built in.
- **Term 2** `(n·H_s) ∇(1/R)` (scalar source, gradient kernel): a
  **gradient of a Laplace single layer** with scalar density `(n·H_s)`.
  → `ngsolve.bem.LaplaceSL` applied to the GridFunction, then take
  the spatial gradient via `Grad(...)` CF.

Both pieces are CoefficientFunctions, so the result is itself a CF
suitable for `gf_H.Set(...)` into a volumetric mesh — automatic
integration with NGSolve post-processing.

The Phase B production kernel:

    J_s = n × H,     M_s = n × E
    E(r) = ∮ [-jωμ_0 Ḡ_e · J_s + ∇ψ × M_s] dS'
    H(r) = ∮ [+jωε_0 Ḡ_e · M_s + ∇ψ × J_s] dS'

decomposes into:

- `Ḡ_e · J_s`, `Ḡ_e · M_s`: the Maxwell SL operator gives this
  exactly (it IS the Stratton-Chu dyadic).
  → `MaxwellSingleLayerPotentialOperator(space, kappa)`
- `∇ψ × M_s`, `∇ψ × J_s`: curl of Helmholtz SL.
  → `MaxwellSingleLayerPotentialOperatorCurl` (provides curl of the
  Helmholtz SL).
The scalar-charge terms are not added separately in this dyadic form;
the longitudinal contribution is contained in the Maxwell SL operator.

The dyadic `(1/k²) ∇∇` correction is handled internally by Maxwell SL
— no special low-k cancellation handling needed at the *operator*
level.  However, the **multilevel expansion** behind the operator
still depends on the frequency regime — see §3.3.

### 3.3 Frequency-dependent FMM applicability

The break-even analysis (§2) silently assumed the FMM is *numerically
stable* at the operating frequency.  That is **not always true**.
Standard Helmholtz FMM exhibits the well-known **low-frequency
breakdown** (Greengard, Huang, Rokhlin, Wandzura 2002; Cheng et al.):
when `kR → 0`, `J_n(kR) → 0` too rapidly with `n`, the multipole
expansion becomes numerically ill-conditioned, and the Helmholtz FMM
gives garbage well before the static limit is reached.

For the Radia / equivalence_source target regimes:

| Regime | Example | `kR_max` (R≈1 m) | FMM path |
|---|---|---|---|
| **Static** | DC, permanent magnets | 0 | Laplace ML (`BiotSavartRegular/SingularMLCF`, `LaplaceSL/DL`).  No breakdown. |
| **Low frequency** | IH 10 kHz–500 kHz | 2·10⁻⁵ – 1·10⁻² | **LF breakdown territory** for Helmholtz FMM.  Route through Laplace ML; ω-dependent imaginary part is small (verified: Phase B vs Phase A static-limit agree to 5e-15 at ω = 1 Hz). |
| **Low / mid** | WPT 100 kHz | 2·10⁻³ | Borderline.  Use Laplace ML + direct imaginary correction. |
| **Mid** | WPT 13.56 MHz | 0.28 | **Standard Helmholtz FMM works.**  Use `MaxwellSL/DL` ML directly. |
| **High** | Antenna GHz | 21 (at R=1 m) | HF-FMM (plane-wave expansion).  Use `MaxwellSL/DL` ML; NGSolve.bem internally selects HF routine when kappa is large. |

**Operational decision rule for Phase D autoswitch**:

```python
if omega == 0:
    # Static reduction — Laplace ML always safe.
    use_path = "laplace_ml"
elif k * R_obs_max < 0.01:
    # Low-frequency Radia regime (IH).  Helmholtz FMM breaks down here;
    # the imaginary correction is < 1% of the real part, so route the
    # static reduction through Laplace ML and ADD the imaginary part
    # via the direct C++ kernel (cheap when ω is small).
    use_path = "laplace_ml + direct_im_correction"
elif k * R_obs_max < 1e3:
    # Standard Helmholtz FMM regime (mid-frequency WPT through low-GHz).
    use_path = "helmholtz_ml"
else:
    # High-frequency regime; NGSolve.bem internally picks the
    # plane-wave-based HF-FMM.  Same operator API.
    use_path = "helmholtz_ml_hf"
```

The `kR < 0.01` threshold is conservative; the actual LF breakdown
of Helmholtz FMM kicks in around `kR ~ 10⁻³ – 10⁻⁵` depending on the
multipole order `p` used.  Below the cutoff, the static-Laplace path
plus a direct linear-in-ω imaginary correction is *both* faster
*and* more accurate than trying to push Helmholtz FMM out of its
regime.

**NGSolve.bem visible symbol gap** (audited 2026-05-26): the API
exposes `BiotSavartRegular/SingularMLCF` (Laplace, well-conditioned
at all frequencies) and Maxwell/Helmholtz `SL/DL` operators with a
`kappa` argument (Helmholtz, subject to LF breakdown at small kappa).
A dedicated **LF-FMM** symbol that handles `kR ≪ 1` Helmholtz cleanly
is **not visible**.  This is consistent with the operational rule
above: in NGSolve.bem the low-frequency regime is handled by routing
through the Laplace ML path, not by a separate LF-FMM kernel.

The Phase D plan (§4) therefore assumes the Laplace-ML / Helmholtz-ML
split per the rule above, NOT an "everything through Helmholtz FMM"
approach.

### 3.4 Two distinct LF breakdowns — Weggler vs Greengard-Huang

Two **different** numerical pathologies hide under the name
"low-frequency BEM breakdown".  They affect different parts of the
pipeline and need different fixes:

| Pathology | Where | Fix in NGSolve.bem |
|---|---|---|
| **Operator-level**: standard EFIE block matrix scales `V_κ / κ²` and becomes ill-conditioned (`O(κ⁻²)`) as `κ → 0`; Krylov solve fails. | BEM **solve** stage (BilinearForm + GMRES) | **Lucy Weggler's product-space formulation** `HDivSurface × SurfaceL2` with `κ² V_κ` re-scaling.  `O(1)` conditioning at all frequencies. Lab reference: `to_developers/ngsolve/low_freq_efie_ngbem_applications.ipynb` (Sugahara 2026-02). |
| **FMM-level**: Helmholtz multipole expansion `J_n(kR) → 0` underflow as `kR → 0`; multipole-to-local translation loses precision. | Kernel **summation** stage (multilevel expansion machinery) | **Greengard-Huang-Rokhlin LF-FMM** (2002).  Or, simpler: route the low-frequency case through Laplace ML and add a small imaginary correction directly. |

For the equivalence-source **one-shot evaluator**:

- It does NOT solve a system → operator-level breakdown is irrelevant.
  Weggler's product-space stabilization is for **EFIE/MFIE solves**
  (used by `radia.peec` / `ngsolve.bem` users who actually invert
  the BEM matrix), not for our Stratton-Chu evaluation.
- It DOES sum kernels at obs points → FMM-level breakdown is the
  one that bites.  Mitigation: §3.3 rule above (Laplace ML for
  `kR < 0.01`).

### 3.5 Conceptual unity (and why the equivalence-source LF rule is right)

Although Weggler's product space and our Laplace-ML routing solve
different problems (solve vs evaluation), they share the same
underlying observation:

> **At low frequency, the static (Laplace) physics dominates and the
> wave correction is a small perturbation.**

Weggler exploits this by reformulating the EFIE matrix so the
solenoidal (`HDivSurface`, `A_κ`) block reduces cleanly to the
Laplace inductance operator as `κ → 0` — see Section 2 of the
notebook: `LaplaceSL(HDivSurface)` projected onto the harmonic
subspace gives the analytic air-core inductance to 0.3 %.

The equivalence-source LF rule (§3.3) exploits the same observation
on the evaluation side: at `kR ≪ 1`, the real (Laplace-like) part
of the Stratton-Chu integrand is captured by Laplace ML, and the
imaginary `O(kR)` correction is small enough to add via direct
kernel summation.

In other words: **Weggler's `κ → 0` limit of the `A_κ` block is the
SOLVE-side analogue of our Laplace-ML routing of the EVALUATE-side
kernel.**  Same physics, different machinery, both correct.

There is one situation where Weggler's path matters directly for
equivalence_source: if a future user wants to **compute** (E, H) on
the equivalence surface from a BEM solve (instead of an FEM solve),
the `ngsolve.bem` BEM solve should use Weggler's product-space
formulation at low frequency.  At the moment the design assumes the
(E, H) on the surface comes from an FEM solve (e.g.
`calc_fem_kelvin.py`), so the Weggler stabilization is in the FEM
side and we never touch it.

## 4. Phase D delivery plan (when actually built)

| Step | What | Acceptance |
|---|---|---|
| D1 | Surface mesh extraction: take `NearFieldSource` panels and construct a Netgen surface mesh + `Surface FESpace` from them | `Mesh.nedge`, `mesh.nv` populate; integration `Integrate(CF(1), mesh, BND)` matches `sum(areas)` to 1e-12 |
| D2 | `evaluate_static_H_ml(obs)`: route Phase A (ω=0) through `BiotSavartRegularMLCF` + `LaplaceSL + Grad`.  ‖result - direct‖∞ / ‖direct‖∞ < 1e-6 | unit test on N=10⁴ panel sphere |
| D3 | Autoswitch in `NearFieldSource.evaluate_static_H`: ML when `N_face*N_obs > 1e9`, else direct C++ | benchmark: N=10⁵×10⁴ ML ≈ 0.05 s vs direct ≈ 5 s |
| D4a | **Low-frequency** branch (`kR_max < 0.01`, IH typical): route Phase B through the **Laplace** ML path (D2) for the real part + direct C++ kernel for the linear-in-ω imaginary correction.  Bypass Helmholtz FMM entirely (LF breakdown). | unit test at ω = 2π × 50 kHz, R = 0.5 m: ‖ML - direct‖∞ / ‖direct‖∞ < 1e-5 |
| D4b | **Mid / high-frequency** branch (`kR_max ≥ 0.01`, WPT MHz and above): route Phase B through `MaxwellSL + MaxwellSLCurl + HelmholtzSL`.  Machine-precision match vs direct on a small case. | unit test at ω = 2π × 13.56 MHz, R = 1 m: same tolerance |
| D4c | `evaluate(...)` autoswitch picks D4a / D4b based on `kR_max` computed from `omega * c_inv * max(‖obs‖)`. | unit test sweeping ω: smooth transition through `kR = 0.01` with < 0.5% discontinuity |
| D5 | Docs update: mark Phase D delivered, drop the in-tree FMM speculation | docs PR |

**Zero new C++ in Radia.**  All acceleration code is glue Python that
calls existing NGSolve.bem operators (Laplace ML for low frequency,
Helmholtz ML for mid / high).  Matches CLAUDE.md "Complement NGSolve"
+ "Do NOT Implement custom H-matrix algorithms".

## 5. Out of scope (do NOT do this)

- ~~Re-vendor ExaFMM-t into `src/ext/`~~ — that was the WRONG idea
  in this doc's previous draft.  NGSolve.bem already has Multilevel
  Expansion (FMM-equivalent) for all needed kernels.
- HACApK for equivalence_source — there is no matrix to compress.
- Custom Multilevel Expansion in Radia C++ — duplicates NGSolve.bem.

## 6. Open questions

1. NGSolve.bem ML primitives work on **NGSolve meshes + GridFunctions**.
   The `NearFieldSource` carries raw `(centroids, normals, areas)` +
   per-face `(E, H)`.  D1 needs a clean bridge: either
   (a) reconstruct a `netgen.meshing.Mesh` from the panel triangulation
   (panels were originally from a Netgen surface mesh anyway), OR
   (b) write a CF that wraps the raw arrays and call ML primitives
   point-by-point (kills the FMM advantage — likely the wrong choice).
2. Real concrete user case: do we ever build a 200³ obs grid on a
   N_face=10⁴ extraction surface (`N_face × N_obs = 8e10`)?  This is
   the trigger.
3. NGSolve.bem ML is on the **NGSolve roadmap** — version 6.2.2603
   ships it.  Upstream evolution may add more multipole/local
   translation options; track via NGSolve release notes.
4. **LF-FMM**: does NGSolve.bem internally use a stabilised low-
   frequency Helmholtz routine when `kappa → 0`, or does the user
   need to detect the breakdown regime and route through Laplace ML
   manually?  §3.3 + §4 (D4a vs D4b) assumes the *manual* split.
   If NGSolve.bem ships a stabilised `MaxwellSL(kappa→0)` upstream,
   D4a / D4b can collapse into a single call with autoselection
   inside the operator.  Probe before implementing D4a.
5. **D4a imaginary correction**: at `kR < 0.01`, the imaginary part
   of the harmonic kernel is `O(kR)` smaller than the real part.
   Does it actually need full kernel summation, or does a coarse
   numerical estimate (e.g. dipole-moment approximation of the
   surface sources) suffice for the imaginary part?  This affects
   D4a cost: full direct C++ at `O(N_face × N_obs)` vs dipole
   approximation at `O(N_obs)`.  Decide based on a user case where
   the imaginary part matters (eddy current loss in IH workpiece
   reconstruction far from the source).

## 7. Recommendation

**Defer Phase D until a concrete user case crosses the break-even.**
Until then:
- The direct C++ kernel (`EvaluateStaticH` / `EvaluateHarmonic`) is
  the reference path.
- 50-120× speedup over Python at N≤10⁵ is already enough for all
  current examples.
- This document is the anchor for the future implementation; no code
  is shipped now.

When the case arrives, do D1-D3 (~few days work), validate, ship as
an internal optimisation behind the existing `evaluate_static_H` /
`evaluate` API — no API surface change.

---

Previous draft of this doc proposed re-vendoring ExaFMM-t into
`src/ext/`.  That violated CLAUDE.md "Complement NGSolve" architecture
and is **rejected** in favour of using NGSolve.bem.  Kept the cost
model and decision rule from the previous draft.
