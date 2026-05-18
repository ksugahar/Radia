# IGTE Symposium 2026 — Paper Outline

**Working title.**  "A Curved-Element, Per-Element Nonlinear Surface
Impedance Method for Boundary Element Eddy-Current Analysis"

**Audience.**  IGTE Symposium 2026 (Aug 2026, Graz / virtual) →
target journal IEEE Trans. Magn. selected papers.

**Author.**  K. Sugahara (Kindai University) + collaborators.

**Code release.**  All numerical results in this paper are
reproducible from radia ≥ 4.55.3 (PyPI):
[github.com/ksugahar/Radia](https://github.com/ksugahar/Radia).

---

## 1. Contribution Statement (1 paragraph)

Hollaus, Kaltenbacher & Schöberl (IEEE TMAG, 2025,
DOI [10.1109/TMAG.2025.3613932](https://doi.org/10.1109/TMAG.2025.3613932))
introduced the nonlinear Effective Surface Impedance Method (ESIM) in
a **magnetic scalar potential FEM formulation**: a 1-D cell problem
through the conductor depth returns `Z_s(|H_t|)`, used as a surface
boundary condition in the outer FEM volume solve, coupled via a
Picard fixed-point iteration ("Karl iteration").

This paper extends the ESIM idea to **Boundary Element Method (BEM)
analysis** of induction-heating workpieces, with three new
contributions:

1. **BEM formulation** with the Hollaus cell problem as the surface
   impedance source — eliminates the volume mesh of air and the
   Kelvin / PML truncation needed by the FEM path.
2. **Curved isoparametric Tri6 elements** for the boundary integral
   equation, capturing surface curvature exactly at second order
   without remeshing.
3. **Element-by-element (per-DOF) nonlinear surface impedance**:
   each surface degree of freedom carries its own `Z_s` computed from
   the locally extracted `|H_t|`, capturing spatial saturation
   patterns that a single mesh-RMS scalar `Z_s` cannot resolve.

The combination of these three with the Hollaus-type damped Picard
iteration enables high-fidelity nonlinear eddy-current analysis on
unbounded-domain induction-heating problems at a fraction of the
FEM-path cost.

---

## 2. Background and State of the Art

### 2.1 Linear SIBC (Leontovich; Krähenbühl-Muller 1993)

Surface impedance boundary conditions reduce 3-D eddy-current problems
with thin skin layer to surface PDEs.  Linear SIBC: `Z_s = (1+j)ρ/δ`
(Dowell-style closed form) for non-magnetic conductors.  Linear-μ
extension via `δ = √(2ρ/(ω μ_0 μ_r))`.

### 2.2 Nonlinear cell problem (Hollaus-Kaltenbacher-Schöberl 2025)

For ferromagnetic workpieces (silicon steel, electrical steel) the
saturation pattern `μ(|H|)` must be resolved.  Hollaus et al. propose:

- 1-D radial PDE through the conductor: `(ρ/r) ∂_r[r ∂_r H] + jω μ(|H|) H = 0`
- Driven by `H(R) = H_t` on the surface
- Returns `Z_s(H_t) = 2(P' + jQ')/|H_t|²`
- **Outer FEM scalar potential** consumes Z_s as Robin BC
- **Picard fixed-point** ("Karl iteration") couples cell ↔ outer

### 2.3 What's missing in Hollaus 2025

(a) **Volume FEM mesh** of the entire air domain → expensive for
unbounded problems; needs Kelvin transformation or PML truncation.
(b) **Single scalar `Z_s`** per Karl outer iteration → saturation
pattern not spatially resolved.
(c) **Linear (flat) elements** typical → curvature captured by mesh
refinement rather than basis order.

This paper addresses all three.

---

## 3. Method (BEM + curved + per-element + nonlinear)

### 3.1 Boundary integral formulation

Scalar magnetic potential `φ` on the workpiece surface `Γ` satisfies
the BIE

$$
\bigl(\tfrac{1}{2} M + \mathrm{DL} + \tfrac{j\omega}{Z_s}\,\mathrm{SL}\cdot M^{-1}\cdot K\bigr)
\varphi
=
\varphi_{\mathrm{inc}}
$$

where SL, DL, M, K are the single-layer, double-layer, mass, and
tangential-gradient stiffness operators on `Γ`.  Solution proceeds via
either dense LU (for small Γ) or HACApK-ACA-compressed GMRES (for
large Γ).

### 3.2 Curved isoparametric Tri6 elements

Each surface element carries 6 nodes (3 vertex + 3 edge midpoint).
The Lagrange-P2 basis `{φ_i^P2}_i^6` is integrated against the
isoparametric mapping `x(u, v) = Σ φ_i^P2 x_i`, giving quadratic
geometric representation.  Cubit's `radia_export netgen "f.vol" order 2`
embeds the curved nodes; NGSolve `Mesh()` automatically loads them.

The in-tree assembler at `bem_sibc_solver.py` builds SL, DL, M, K
on Tri6 with Sauter-Schwab Duffy quadrature for the singular pair
(`intree_singular_n_q = 6`) and Gauss-Legendre at degree 7 for
regular pairs.  Convergence rate: **O(h³)** in `||φ - φ_h||_2` for
smooth Γ + smooth Z_s, vs O(h²) for flat Tri3.

### 3.3 Per-element (per-DOF) ESIM

Instead of a single scalar Z_s applied uniformly to the BIE, each
DOF `i ∈ {1, ..., n_DOF}` carries an independent `Z_s[i] ∈ ℂ`,
computed from the local tangential field amplitude.  The BIE row `i`
becomes:

$$
\Bigl(\tfrac{1}{2} M_{i,:} + \mathrm{DL}_{i,:} + \frac{j\omega}{Z_s[i]}\,(\mathrm{SL} M^{-1} K)_{i,:}\Bigr)\varphi = (\varphi_{\mathrm{inc}})_i.
$$

Per-DOF `|H_t|` is extracted variationally:

$$
|H_t|^2_i = \frac{|\varphi_i \cdot (K\varphi)_i|}{(M\mathbf{1})_i}
$$

(the i-th summand of `φᵀ K φ` localised + divided by the i-th
lumped-mass / row-sum-of-M factor).  This formula is consistent with
the global mesh-RMS `H_t_rms² = (φᵀKφ)/A_wp` in the limit of uniform
basis-function gradient.

### 3.4 Karl iteration on the per-element Z_s vector

Outer fixed-point on the **ndof-dimensional vector** `Z_s ∈ ℂ^n`:

1. Seed: `Z_s^{(0)}[i] = E(H_t = 5 A/m)` for all i (conservative low-H).
2. Solve BIE at current `Z_s` → extract `H_t[i]` per DOF.
3. Cell-solver call per DOF: `Z_s^{(k+1)}_new[i] = E(H_t[i])`.
4. Damped update: `Z_s^{(k+1)}[i] = α Z_s_new[i] + (1-α) Z_s^{(k)}[i]`.
5. Stop when `max_i |dZ[i]| / |Z_s[i]| < tol`.

The ndim-version converges with the same geometric rate as the scalar
case (empirical Lipschitz `L ≈ 1` for steel workpieces; `α = 0.5`
gives `αL = 0.5 < 1`, strict contraction).  Cost: `N_DOF · t_cell`
extra per Karl iteration; on production meshes (166-500 DOFs) this
overhead is ~5-30 s, negligible compared to the BIE solve.

---

## 4. Numerical Results

### 4.1 Cell-problem validation: Bessel reference (§ 2 of [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md))

| ξ = R/δ | Frequency | Re(Z_s)_anal [Ω] | Re(Z_s)_num [Ω] | Rel. err |
|---|---|---|---|---|
| 4   | 1 kHz   | 7.50e-3 | 7.50e-3 | <1e-5 |
| 14  | 10 kHz  | 2.49e-2 | 2.49e-2 | <1e-4 |
| 45  | 100 kHz | 7.95e-2 | 7.95e-2 | 4e-4 |
| 140 | 1 MHz   | 2.50e-1 | 2.50e-1 | 1.3e-3 |

The cell solver matches the closed-form Bessel-`I_0/I_1` reference
(Wakao-Igarashi-Fujiwara-Kameari Part 5) to <0.13 % across `ξ ∈ [4, 140]`.

### 4.2 Internal consistency: 3 outer-solver paths (§ 3)

PEEC filament + BEM-SIBC ↔ PEEC filament + FEM-Kelvin ↔ FEM A-V
volumetric coil on a canonical steel-cylinder + gapped-torus test
case.  At convergence:

- **PEEC-BEM ↔ FEM-Kelvin** agree on `|Z_s|` to <2 % at all four
  test frequencies (10, 50, 100, 500 kHz).
- **FEM-coilmesh** tracks at 10 kHz (<1 %) but diverges at higher
  frequencies due to coil-mesh under-resolution (`δ_coil → h_coil`),
  not a method bug.

### 4.3 Curved-element benefit (h vs p refinement)

The BEM curved element / per-element ESIM is implemented at
`--h1-order 2` with `intree_geom_order=2`.  On the gapped-torus
benchmark at 50 kHz:

- Flat Tri3 (`--h1-order 1`, curve_order=1): ~166 BIE DOFs → P_wp = X
- Curved Tri6 (`--h1-order 2`, curve_order=2): ~580 BIE DOFs → P_wp = Y

(Actual values TBD — see [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md)
§ 6 placeholder.)

### 4.4 Per-element vs scalar Z_s comparison

On the same gapped-torus benchmark, compare:

- Scalar Karl: single Z_s_rms from mesh-averaged H_t → P_wp_scalar
- Per-element Karl: 166 Z_s values from per-DOF |H_t| → P_wp_per_element

Per-element captures the saturation pattern on the workpiece — the
edge regions (where H_t is highest) get a saturated, low-|Z_s|; the
center stays in the high-μ regime.  Expected effect: 5-15 % difference
in P_wp when the workpiece's |H_t| range spans the BH knee.

(Actual numerical comparison TBD; needs a workpiece geometry that
intentionally drives the spatial saturation contrast.  The
`em_sample_bh.txt` curve at I_port = 1 A places the entire workpiece
in the linear-μ regime, so the per-element advantage is small for the
default benchmark.)

### 4.5 Karl convergence rate (§ 3 history table)

`dZ/Z` decays geometrically with ratio ~0.4 per iteration at `α=0.5`,
converging in 5-7 iterations to `tol=1e-3`.  Same rate for scalar and
per-element variants.

---

## 5. What this paper does NOT do

For honest scope-limiting:

- **No measurement vs simulation** validation (lab-built steel
  cylinder + B-H meter).  All validation is internal-consistency or
  against closed-form references.
- **No saturation-regime benchmark** (Stoll 1974 envelope; Lavers-
  Biringer 1985 plate).  These are roadmap.
- **No transient analysis** — only steady-state phasor.
- **No coupling to NGSolve.bem upstream FMM** — the BIE matrix is
  assembled either dense or HACApK in-tree.  See
  [`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) for
  the per-element + HACApK gap.

These are signalled as future work in § 6.

---

## 6. Conclusion + Future Work

(Standard journal closing — ~150 words.  Key points:)

- Curved-element BEM with per-element nonlinear Z_s is a viable
  alternative to FEM A-V volumetric eddy-current analysis for
  induction-heating workpieces.
- The Hollaus-Kaltenbacher-Schöberl 2025 cell problem can be cleanly
  decoupled from its FEM context and re-applied to BEM with no loss
  of accuracy in the regime where SIBC is valid.
- Per-element Z_s captures saturation patterns that scalar mesh-RMS
  cannot resolve, important for workpieces with high `|H_t|`
  contrast.
- Karl iteration extends naturally from scalar to vector-valued
  fixed-point with the same geometric convergence rate.

Future work:
- Anderson acceleration (m=2-3) on the Karl loop for deep-saturation
  cases.
- HACApK backend support for per-element `Z_s` array (currently
  dense-only).
- Wide-band extension via nonlocal SIBC (Bilicz-Badics-Pávó 2023).
- Measurement validation on physical IH workpieces.

---

## 7. Reproducibility

All figures and tables are reproducible from:

```bash
pip install 'radia[cubit,gui]==4.55.3' radia-mcp==0.55.3

# Bessel reference (Table § 4.1)
python examples/ih_esim_benchmark/analytical_bessel_baseline.py

# Three-path consistency (Table § 4.2)
python examples/ih_esim_benchmark/benchmark.py --frequencies "1e4,5e4,1e5,5e5"

# Curved-element comparison (Table § 4.3) — TBD: workflow at
# tests/panels/test_h1_p2_curved.py is planned

# Per-element vs scalar (Table § 4.4)
python src/radia/panels/calc_inductance.py \
    --coil-solver peec --coil-step <coil.step> --vol <wp.vol> \
    --impedance-model esim --bh-file <bh.txt> \
    [--esim-per-panel] \         # toggle for per-element vs scalar
    --frequency 50000 ...
```

The full benchmark JSON outputs are committed to
[`examples/ih_esim_benchmark/results.json`](../../examples/ih_esim_benchmark/results.json)
and `tests/panels/golden/`.

---

## 8. References

(Subset of [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md) § 8;
the IGTE digest will use 8-10 refs.)

1. **K. Hollaus, M. Kaltenbacher, J. Schöberl**, "A Nonlinear
   Effective Surface Impedance in a Magnetic Scalar Potential
   Formulation," IEEE Trans. Magn., 2025.
   DOI: 10.1109/TMAG.2025.3613932.
2. **L. Krähenbühl, D. Muller**, "Thin layers in electrical
   engineering — Example of shell models in analysing eddy-currents
   by boundary and finite element methods," IEEE Trans. Magn. 29(2),
   1993.
3. **J. D. Lavers, P. P. Biringer**, "An efficient calculation of
   effective surface impedance for nonlinear ferromagnetic
   materials," IEEE Trans. Magn. 21(5), 1985.
4. **R. L. Stoll**, *The Analysis of Eddy Currents*, Oxford
   University Press, 1974.
5. **R. F. Harrington**, *Time-Harmonic Electromagnetic Fields*,
   McGraw-Hill, 1961. (Reaction integral / Lorentz reciprocity.)
6. **E. Dlala, A. Belahcen, A. Arkkio**, "Optimal Convergence of the
   Fixed-Point Method for Nonlinear Eddy-Current Problems," IEEE
   Trans. Magn. 44(6), 2008.
7. **S. Yuferev, N. Ida**, *Surface Impedance Boundary Conditions:
   A Comprehensive Approach*, CRC Press, 2009.
8. **S. Bilicz, Z. Badics, J. Pávó**, "Nonlocal surface impedance
   boundary condition for wide-band eddy-current problems," ISEM 2023.
9. **S. Wakao, H. Igarashi, K. Fujiwara, A. Kameari**, "Various
   Verifications of Eddy Current Analysis (Parts 1-9)," T.IEE Japan,
   2008-2018.  (Part 5 used for the Bessel cylinder baseline.)

---

## Cross-references to supporting docs

- [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md) — full
  formulation (PDE derivation, weak forms, Lorentz reciprocity,
  HCurl weak form, complex-μ, limitations).
- [`IMPLEMENTATION.md`](IMPLEMENTATION.md) — three-solver dispatch,
  Karl loop details, per-DOF Z_s mechanics, Lipschitz / Anderson.
- [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) — Bessel parity, three-
  path consistency, 2-D axisym, mesh convergence (placeholder),
  worked example.
- [`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) —
  why PEEC vs BEM-A produce different coil R (relevant to the
  paper's "future work" discussion of vector vs scalar BEM).
- [`USAGE.md`](USAGE.md) — CLI guide.

---

**Document version**: 2026-05-18 (radia v4.55.3+).
**Next steps**: regenerate § 4.3 (curved vs flat element comparison)
and § 4.4 (per-element vs scalar) with real benchmark data before
IGTE submission deadline.
