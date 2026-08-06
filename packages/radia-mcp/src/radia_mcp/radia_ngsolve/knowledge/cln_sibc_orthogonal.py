"""CLN expansion-point + SIBC orthogonal-residual theory.

POLICY NOTE (2026-06-12): The Warburg-Schur termination described in
this module was SUPERSEDED by the Mixed Galerkin (CLN bulk + HOIBC
surface envelope, Schur-coupled, no `d` parameter) framework.  The
Warburg-Schur code was hard-deleted from the repo on 2026-06-12; see
`memory/project_warburg_schur_deprecated_2026_06_12.md` for the
history.  This knowledge module is kept for context but its
"Warburg-block closure" approach is deprecated.  New work should
reference `validation_test/mixed_galerkin/` and use the s-dependent HOIBC
envelope, not the rational K √s/(s+d) Warburg block.  A runnable,
VERIFIED results showcase (cylinder Senior tower 0.04% -> 1e-5%; sphere
gamma_1 0.11% -> 0.001%, 102x) is at
`docs/mixed_galerkin/mixed_galerkin_results.ipynb` (it imports the kept
corpus's `_references/`).

Knowledge module for the new CLN theory established in IGTE 2026 work
(Sugahara session 2026-05-16): N stages of Cauer ladder + SIBC analytical
residual form an L^2-orthogonal complement that spans the full frequency
response.  Avoids the precision frontier of high-stage Cauer extraction
(Nagamine 2026 verified affine arithmetic) by replacing the high-stage
tail with a closed-form Warburg / SIBC asymptote anchored by physical
K_SIBC = sqrt(sigma/mu) times a geometric factor.

Verified analytically with Mathematica on:
  - 1D slab (c = 1 mm)
  - 2D rectangular prism (a x b = 17.72 x 2 mm, Dirichlet A_z = 0,
    Rong-Cangellaris near-field BC)

Scripts:
  repo:/validation_test/maglev/research_cln/A1_SIBC_CLN_orthogonal_slab.wls
  repo:/validation_test/maglev/research_cln/A1_SIBC_CLN_orthogonal_rect2D.wls
  repo:/validation_test/maglev/research_cln/A1_SIBC_CLN_orthogonal_*.py (plotters)

Companion to the rigorous Nagamine et al. (JJIAM 2026 submitted) verified
N-stage Cauer extraction for a 2D square prism with affine arithmetic;
together they form a precision-frontier-aware framework: Nagamine = verified
high-stage extraction with intervals, this theory = small-N + asymptotic
SIBC tail to skip the high-stage frontier entirely.

Also includes the "xfem_vs_sibc" decision framework (2026-05-25 β-1 Phase
1-3 benchmark): when to use Hiruma 2023 XFEM at the FE level vs the
augmented CLN at the ROM level vs classical SIBC at the boundary, and
the port-driven scope condition that makes augmented CLN inapplicable
to volume-source problems without XFEM stacking.
"""


CLN_SIBC_ORTHOGONAL_OVERVIEW = """
# CLN expansion-point with SIBC-orthogonal residual

## One-line claim

For a linear eddy-current body Ω with Foster eigenmodes {φ_n}, the
L^2-orthogonal decomposition
    L^2(Ω) = span{φ_0, ..., φ_{N-1}} ⊕ span{φ_N, φ_{N+1}, ...}
induces an admittance decomposition
    Y_exact(jω) = Y_CLN^N(jω) + Y_SIBC^⊥(jω)
in which Y_SIBC^⊥(jω) → K_SIBC · (jω)^{-1/2} as ω → ∞ (Warburg / SIBC).
A small N (e.g. 3) plus the closed-form asymptote covers the entire
frequency band without needing to extract high-stage Cauer rungs.

## Why this matters

High-stage Cauer extraction is precision-bound: even Nagamine et al. (2026)
with 192-bit MPFR + affine arithmetic show ~60×/stage interval growth, and
ordinary 30-dps mpmath grows ~300×/stage.  In the IGTE 2026 work the FEM
Tanimoto-AT + QD-Wynn pipeline becomes non-monotonic past stage 3-4 even at
mesh+p-convergence.  Replacing the high-stage tail with an analytical SIBC
asymptote removes the precision-frontier dependence:

    accuracy limit moves from   "N → ∞ extraction"   to
                                "N small + SIBC residual fidelity"

## Construction

1. Solve the eddy-current eigenvalue problem in Ω:
       -∇² φ_n = λ_n φ_n,    φ_n|_{∂Ω} = 0  (Dirichlet / Rong-Cangellaris)
   yielding orthonormal eigenfunctions and τ_n = μσ/λ_n.

2. Compute coupling g_n^2 = |⟨A_ext, φ_n⟩|^2 / ⟨φ_n, φ_n⟩.

3. Y_exact(jω) = Σ_{n} g_n^2 / (1 + jωτ_n)  (Mittag–Leffler / Foster).

4. Choose N (typically 2–4); define
       Y_CLN^N(jω)    = Σ_{n=0}^{N-1} g_n^2 / (1 + jωτ_n)
       Y_SIBC^⊥(jω)   = Y_exact - Y_CLN^N = Σ_{n≥N} g_n^2 / (1 + jωτ_n)
   The two summands live in orthogonal L^2 subspaces by Foster orthogonality.

5. At high ω, Y_SIBC^⊥(jω) → K_SIBC (jω)^{-1/2}, with K_SIBC determined by
   geometry / physics (NOT continuity), e.g.
       K_SIBC^{slab}    = √(σ/μ)
       K_SIBC^{2D prism} = P_perim · √(σ/μ)   (P_perim = 2(a+b) for rect)
       K_SIBC^{3D body}  = surface integral involving Z_s = √(jωμ/σ).

   This anchors the unified admittance without inheriting the truncation
   bias of the empirical "continuity seam" Y_CLN^N(jω_N) = K_eff·(jω_N)^α
   used in the original Matsuo expansion-point construction.

## Verification signature

arg Y_SIBC^⊥(jω) → -45° (Warburg phase) at high ω, monotonically.

|Y_SIBC^⊥(jω)| · |jω|^{1/2} = "effective K_eff(ω)" must converge to the
physically anchored K_SIBC as ω → ∞.  Convergence rate increases with N:
in the 2D rectangle 17.72 × 2 mm, N=3 reaches arg → -44° at f ≈ 10^7 Hz,
N=20 reaches the same at f ≈ 10^5 Hz.
"""


CLN_SIBC_ORTHOGONAL_RELATION_TO_MATSUO = """
# Relation to Matsuo's expansion-point theory (SA-26-014)

The original Matsuo construction connects the truncated Cauer at ω_N to
a CPE termination Z = K(jω)^{-α} by **continuity**:
    Y_CLN^N(jω_N) = K_eff · (jω_N)^{-α},  α = (k_R+1)/(k_L+k_R+2)
where ω_N = √(2^(α/2)-1) · R_N/L_N and (k_L, k_R) come from the asymptotic
power law of Cauer elements.

This construction is **biased** for small N because Y_CLN^N(jω_N) is
already in the truncation-error regime at ω_N (the truncation error in
Y_CLN^N is the missing higher modes, which the seam is supposed to capture).
Therefore K_eff inherits a 10–30% bias relative to the physical K_SIBC.

The CLN+SIBC orthogonal-residual reformulation:

  - replaces the empirical "K from continuity" with the physical
    K_SIBC = √(σ/μ) × geometric factor (no continuity needed),
  - exhibits the Cauer truncation and the SIBC tail as L^2-orthogonal
    components of Y_exact in the Foster eigenmode basis,
  - reduces the seam jump to a quantitative measure of the missing
    high-stage contribution Σ_{n≥N} g_n^2/(1+jωτ_n) at ω_N.

This is the **basis-family orthogonalization** at the expansion point
that connects different basis spaces (Cauer polynomial / FE on Ω vs SIBC
boundary trace on ∂Ω) implicit in Matsuo's framework.
"""


CLN_SIBC_ORTHOGONAL_MATH = """
# Mathematical derivation

## Foster decomposition completeness

For a Dirichlet eigenvalue problem on a Lipschitz bounded Ω, the Foster
eigenfunctions {φ_n} are complete in L^2(Ω) and form an orthonormal basis:

    Σ_n φ_n(x) φ_n(y) = δ(x - y)  (completeness),
    ⟨φ_n, φ_m⟩ = δ_{nm}            (orthonormality, L^2(Ω)).

The Foster admittance series Σ_n g_n^2/(1+sτ_n) is the Mittag-Leffler
representation of the slab/prism/body admittance, convergent for every s
in the right-half plane.  Partial sums up to N preserve the Foster
orthogonality of captured modes.

## Y_CLN^N vs Y_SIBC^⊥ orthogonality

The truncated and residual admittances are **functionals** of different
mode subspaces:

    Y_CLN^N(s) = Σ_{n<N} g_n^2/(1+sτ_n)    lives in span{φ_n}_{n<N}
    Y_SIBC^⊥(s) = Σ_{n≥N} g_n^2/(1+sτ_n)   lives in span{φ_n}_{n≥N}

Since ⟨φ_a, φ_b⟩ = 0 for any a < N ≤ b, the two contributions are
L^2-orthogonal in the field representation:

    A_red(x; jω) = jω Σ_n c_n(jω) φ_n(x)
    A_red^CLN     = jω Σ_{n<N} c_n φ_n  ⟂  jω Σ_{n≥N} c_n φ_n = A_red^SIBC

This automatic orthogonality is the underlying reason the seam is well-
defined: it is the Riesz projection onto complementary closed subspaces,
not a heuristic continuity.

## SIBC asymptote of the residual

For large ωτ_n,
    g_n^2/(1+jωτ_n) ≈ g_n^2/(jωτ_n) = (g_n^2/τ_n) (jω)^{-1}.

Summing the residual with the Foster density g_n^2/τ_n ~ const · 1/(2n+1)
(slab, n→∞) gives a logarithmic divergence in (1/jω) coefficient, which
must be regularized; the correct asymptotic representation is the
Warburg form:

    Σ_{n≥N} g_n^2/(1+jωτ_n)
      ~ ∫_{N}^∞  (2σc/((n+1/2)^2 π^2)) / (1 + jω · μσc^2/((n+1/2)^2 π^2)) dn
      = K_SIBC (jω)^{-1/2}   as  ω → ∞,
   with  K_SIBC = √(σ/μ)·(geometric factor).

For 1D slab: factor = 1 → K_SIBC = √(σ/μ) ≈ 6.79×10^6 (Cu).
For 2D rect: factor = perimeter 2(a+b) → K_SIBC = 2(a+b)·√(σ/μ).
For 3D body: factor = total wetted surface S, with edge/corner corrections.

## K_eff(ω) convergence to K_SIBC

Define
    K_eff(ω) := Y_SIBC^⊥(jω) · (jω)^{+1/2}.
By construction K_eff(ω) → K_SIBC and arg K_eff → 0 as ω → ∞.
Convergence rate is set by the largest retained τ_N: |K_eff/K_SIBC - 1|
≲ O(1/(ωτ_N)) for ω > 1/τ_N.

Practical regime: with N = 3 on Cu cuboid 17.72×2 mm cross-section,
arg K_eff = -45° (Warburg phase) is reached at f ≈ 10^7 Hz; with N = 20
at f ≈ 10^5 Hz.
"""


CLN_SIBC_ORTHOGONAL_VERIFICATION = """
# Mathematica verification (2026-05-16)

## c_basis analytical for sphere (Stoll closed form)

Sphere with uniform B_0 applied (radius a, Cu, mu = mu_0):
  - Stoll susceptibility: 1+chi(s) = (3/alpha^2)(alpha coth(alpha) - 1),
    alpha = a Sqrt[s mu sigma]
  - High-omega asymptote: 1+chi -> 3/alpha -> K_SIBC^{Stoll} = 3/(a Sqrt[mu sigma])
    Numerically (a=10mm Cu):  K_SIBC^{Stoll} = 35.14
    Numerical verify at omega=10^12: 35.140 (0.001% agreement)
  - Tanimoto-AT R_0 = 15 / (2 pi sigma B_0^2 a^5)
    Numerically: R_0 = 411.6 (Tanimoto units)
  - c_basis^{sphere} = Y_Tanimoto(0)/Y_Stoll(0) = (1/R_0)/1 = 0.00243
  - K_SIBC^{Tanimoto, sphere} = c_basis * K_phys = 0.0854

Universal vs basis-specific structure:
  Universal:      Foster L^2 orthogonality, K_SIBC^{phys} = sqrt(sigma/mu) * G
                  (G = geometric factor: slab->c, sphere->3a, cuboid->4ac)
  Basis-specific: c_basis depends on observable normalisation + geometry
                  Sphere: 0.00243; Cuboid AC fit: ~22000 (7 orders apart)

  Script: sphere_unified_basis_derivation.{wls, json}

## Foster N=3..100 visualisation (SIBC = high-mode collective)

For 1D slab, K_eff(N, omega) = (Y_exact - Y_Foster^N) * Sqrt[j omega]:

  omega [rad/s]   N=3     N=10    N=50    N=100
  10^6  (160 kHz) 0.53    0.17    0.035   0.017
  10^8  (16 MHz)  0.95    0.84    0.34    0.17
  10^10           0.995   0.98    0.92    0.84

For each N, K_eff approaches K_SIBC as omega rises, but the transition
frequency itself shifts with N (N=3: ~10^7 Hz, N=100: ~10^11 Hz).  This
shows the SIBC asymptote is the *collective* effect of modes n >= N, not
a single-mode behaviour.  The 20%-wall of 1-parameter Cole-Davidson
(Section above) arises specifically because, at omega ~ 10^6 rad/s with
N=3, the residual is in the middle of its asymptotic transition and a
1-parameter form cannot capture the non-trivial transition shape.

  Script: slab_foster_N_visualization.{wls, py, png, pdf}

## 1D slab (c = 1 mm, Cu)

  - Foster spectrum first 10 modes: τ_n = c^2 μσ/((n+1/2)^2 π^2)
    τ_0 = 29.54 μs, τ_1 = 3.28 μs, τ_2 = 1.18 μs, ...
  - g_n^2 = 2σc/((n+1/2)^2 π^2) sums to σc = Y(0) (Mittag-Leffler ✓)
  - L^2 orthogonality of φ_n on (-c, c): δ_{nm} verified symbolically

  Frequency sweep |Y/Y(0)|:
    f       |Y_exact|  |Y_CLN^3|  |Y_SIBC^⊥|  |residual|/|Warburg|
    1 Hz     1.000      0.933      0.067       1.4e-3
    100 kHz  0.148      0.103      0.065       0.44
    1 MHz    0.047      0.013      0.038       0.81
    10 MHz   0.015      0.001      0.014       0.94 ← residual ≈ Warburg

  K_eff convergence (target K_SIBC = √(σ/μ) = 6.79×10^6):
    f = 10 MHz:  |K_eff| = 6.38e6, arg = -3.8°
    f = 1 GHz:   |K_eff| = 6.75e6, arg = -0.4°

  Script: A1_SIBC_CLN_orthogonal_slab.wls
  Output: A1_SIBC_CLN_orthogonal_slab.{json,png,pdf}

## 2D rectangular prism (a = 17.72 mm, b = 2 mm, Cu, Dirichlet A_z = 0)

  Aspect a/b = 8.86; eigenvalue spectrum non-degenerate.
  Top 12 Foster modes by g^2 (B applied in x-hat):
    rank 1: (m=1, n=1)  τ = 29.17 μs   g^2 ranked by descending value
    rank 2: (1, 2)      τ = 7.36 μs
    rank 3: (3, 1)      τ = 26.50 μs
    rank 4: (1, 3)      τ = 3.28 μs
    ...

  Y/Y(0) sweep + phase: see Mathematica output. Residual phase reaches
  -45° at f ~ 10^7 Hz with N=3, at f ~ 10^5 Hz with N=20.

  Script: A1_SIBC_CLN_orthogonal_rect2D.wls
  Output: A1_SIBC_CLN_orthogonal_rect2D.{json,png,pdf}

## Open implementations

- **2D square (Nagamine geometry, a = b)**: eigenvalue degeneracy
  (1,3)/(3,1) of identical τ → bundle the 2D eigenspace and project
  applied B onto its physical representative.  Same K_SIBC theory holds.

- **3D cuboid (A1 17.72 x 17.72 x 2 mm)**: requires HDivDivFreeHex basis
  surface traces + Green's first identity to express the K_SIBC inner
  product as a pure surface integral on ∂Ω (no interior fields, since
  SIBC has no internal solution).  Volume residual term vanishes if
  basis functions extend harmonically into Ω (Trefftz-type), otherwise
  remains as a correction.
"""


CLN_SIBC_ORTHOGONAL_NAGAMINE_LINK = """
# Relationship to Nagamine et al. (JJIAM 2026 submission)

Nagamine, Yamaguchi, Sugahara, Hiruma, Mifune, Matsuo provide verified
numerical computations of the Cauer network representation for a 2D
infinite-length square prism conductor (Dirichlet A_z = 0 BC,
Rong-Cangellaris formulation):

  - 192-bit MPFR + interval arithmetic + affine arithmetic + INTLAB
  - 16 inductor stages (L_1 ... L_31) and 15 resistor stages (R_2 ... R_30)
    with verified interval bounds
  - ~60× interval width growth per stage from QD-Wynn + affine arithmetic
  - Final stage L_31 interval width ≈ 1.90×10^{-5}

This is the **rigorous reference** for high-stage Cauer extraction in 2D.
Our orthogonal-residual framework is **complementary**:

   Nagamine:  rigorous N = 15 Cauer rung extraction
   Ours:      small N = 2-4 Cauer rungs + analytic SIBC tail

Together they form the "verified small-N + analytical SIBC" path, where:
  1. Nagamine's verified rungs are the trusted ground for the Cauer head.
  2. Our orthogonal residual provides the matching high-frequency tail
     via the Warburg asymptote.
  3. The seam consistency Y_Cauer^N(jω_N) ≈ K_SIBC(jω_N)^{-1/2} can be
     used to check the missing-modes contribution at any N.

For IGTE 2026 narrative the two-paper coupling is:
  - Nagamine = "verification of Cauer extraction algorithm"
  - This work = "operational framework that does not need high N"

Both rest on the same physical CLN representation; together they cover
the precision frontier (Nagamine direct, this work via asymptote skip).
"""


CLN_SIBC_ORTHOGONAL_KURIYAMA = """
# Kuriyama 2019 multi-expansion-point framework (the canonical method)

The expansion-point orthogonalisation discussed above was rigorously
formalised by Kuriyama, Kameari, Ebrahimi, Fujiwara, Sugahara, Shindo
and Matsuo in IEEE Trans. Magn. 55(6), paper WO1-2 (2019):

  "Cauer Ladder Network With Multiple Expansion Points for Efficient
   Model Order Reduction of Eddy-Current Field"
  Handle: http://hdl.handle.net/2433/244835
  Open access: https://core.ac.uk/download/pdf/237666689.pdf

## Setup (A-formulation)

For expansion point s_i = 2 pi f_i, define the coefficient matrix
    K_i = C^T nu C + sigma s_i
where C is the FEM edge-face incidence matrix, nu the reluctivity
matrix, sigma the conductivity matrix.

## Two expansion points (s_0, s_1) — the core procedure

Stages 0 .. M-2 use the standard single-point CLN with s_0.  At
stage M-1, the expansion point switches to s_1 via eqs (9)-(13):

    K_1 a_{2M-1} - K_0 a_{2M-3} = (1/mu_{2M-2}) sigma e_{2M-2}     (9)

    e_{2M} - e_{2M-2} = -(1/nu_{2M-1}) (a_{2M-1} + sum p_{2i-1} a_{2i-1})  (10)

The orthogonalisation coefficient is the KEY equation:

    p_{2i-1} = -(1/mu_{2i-1}) a^T_{2M-1} K_0 a_{2i-1}    (1 <= i <= M-1)  (13)

This is Gram-Schmidt orthogonalisation of the new (s_1-derived) basis
function a_{2M-1} against the captured (s_0-derived) basis functions
{a_1, a_3, ..., a_{2M-3}} using the K_0 (= s_0-evaluated) inner product.

The corrected basis a'_{2M-1} = a_{2M-1} + sum p_{2i-1} a_{2i-1} is then
used in the network expansion (eqs 15, 16).

## Additional inductor at the transition (eq 17)

When a^T_{2M-3} K_1 a'_{2M-1} != 0 an extra L element is required:

    L = -nu_{2M-1} / [(s_1 - s_0) sigma_{2M-2} (nu_{2M-1} + p_{2M-3} nu_{2M-3})]

## Numerical results

  - Iron bar rectangular cross-section: single expansion s_0 = 0 gives
    inaccurate impedance even at 100 stages; with (s_0, s_1) = (0, 2 pi f_1)
    just 4 stages cover several decades.
  - Wireless power transfer at 85 kHz: (s_0, s_1) = 2 pi (0, 85 kHz) with
    only 6 stages reproduces the impedance accurately from 0 to 1 MHz.

## Relation to SIBC connection

Kuriyama 2019 uses FINITE expansion points s_i = 2 pi f_i.  The SIBC limit
corresponds to s = infinity, which the paper does not treat explicitly.
Extension paths:

  - Shift-and-invert at s = infinity: replace K_infinity (singular)
    with the matrix-vector operation A v = (1/sigma) C^T nu C v.  The
    Krylov subspace span{v, A v, A^2 v, ...} captures the high-frequency
    asymptote without ever forming K_infinity.
  - Adaptive s_1 schedule: take s_1 large enough that ω s_1 dominates,
    use Kuriyama's procedure unchanged.
  - Circuit-level interpretation: the s_infinity stage corresponds to a
    Warburg element (CPE alpha = 1/2) appended to the Cauer ladder.

The Foster L^2 orthogonality + SIBC asymptote argument given in the
Math section above is the field-theoretic counterpart of the Kuriyama
algebraic procedure.  Both achieve the same goal of bridging different
expansion-point basis families via Gram-Schmidt orthogonalisation.

## Why this is the canonical method

  - Kuriyama's procedure is exact at each expansion point (matches all
    moments at s_0 and at s_1 simultaneously).
  - The "20% mid-band wall" of single-parameter closed-form residual
    approximations does not arise — by construction the multi-point
    expansion matches local Taylor moments at every chosen s_i.
  - The Wireless Power Transfer numerical example shows accuracy at
    every frequency in 0-1 MHz with only 6 stages, demonstrating that
    the framework operates well past the 1-MHz scale where our slab
    Cole-Davidson hit the 20% wall.

For our SIBC-CLN unification programme this paper is the *operational*
standard; the Foster L^2 orthogonality view (this knowledge module's
Math section) is the *theoretical* counterpart.  Implementing Kuriyama
2019 directly is the recommended path to break the 20% mid-band wall.

## s = infinity extension (planned, 2026-05-16 directive)

Kuriyama 2019 treats only finite s_i = 2 pi f_i.  The natural extension
to s = infinity yields the SIBC connection:

  - K_infinity = C^T nu C + sigma * infinity is singular; do not form
    it directly.
  - Shift-and-invert: build the Krylov subspace span{v, A v, A^2 v, ...}
    with A = sigma^{-1} C^T nu C.  This captures the high-frequency
    asymptote without ever inverting K_infinity.
  - The first basis vector A v corresponds to the SIBC (Warburg)
    contribution; subsequent vectors generate higher-order corrections.

Open question (user 2026-05-16): *Is a single SIBC element sufficient
at s = infinity, or are corner/edge corrections required?*

  - Sphere: a single SIBC mode suffices in the asymptotic limit
    (K_SIBC = 3/(a sqrt(mu sigma)), Stoll closed form).
  - Cuboid: corner/edge singularities may demand higher-order SIBC
    (HOSIBC) or multiple "near-infinity" expansion points selected
    adaptively.
  - This is the central research question for the journal extension
    of the IGTE 2026 draft (currently held).

Concrete next-session task:
  1. Implement Kuriyama eq (13) on the 1D slab with s_0 = 0 and
     finite s_1 (e.g. 2 pi * 1 MHz); demonstrate the 20% mid-band
     wall is broken.
  2. Generalise to s_1 -> infinity via shift-and-invert.
  3. Test on Nagamine-style 2D square (verified reference) and on
     A1 cuboid (AC ground truth at three meshes available).
"""


CLN_SIBC_ORTHOGONAL_OUTLOOK = """
# Outlook: 2D Nagamine + 3D cuboid completion

## Path 1 — 2D Nagamine square (a = b)

Geometry: square prism a × a, z infinite, Dirichlet A_z = 0 on ∂Ω.

Eigenvalues λ_{m,n} = (π/a)^2 (m^2 + n^2), DEGENERATE under (m,n)↔(n,m)
when m ≠ n.  Degenerate pairs (1,3)/(3,1), (1,5)/(5,1), ... carry the
same τ_{m,n} and must be bundled into a single Foster "stage" with
combined g^2 = g_{m,n}^2 + g_{n,m}^2.

With applied B in x-hat (A_ext_z = B_0 y), the symmetric combination
[(m,n) + (n,m)] / √2 couples; the antisymmetric combination is in the
null space of the projection.  Thus a "stage" in the symmetry-adapted
basis is a single mode.

K_SIBC^{2D square} = 4a · √(σ/μ) (perimeter formula).

Verification: extend A1_SIBC_CLN_orthogonal_rect2D.wls with a=b and a
symmetry-adapted projection.  Compare phase residual convergence to the
non-degenerate rectangle case; convergence rate should be similar.

## Path 2 — 3D cuboid (A1 17.72 × 17.72 × 2 mm)

Eigenvalue problem on a 3D box with Dirichlet A·n = 0 (or Rong-Cangellaris
on conductor surface):

    -∇² A_red - jωμσ A_red = jωμσ A_ext  in Ω_c,
    A_red·n = 0  on ∂Ω_c.

Foster eigenmodes are HDiv-div-free vector fields with eigenvalues
λ_{m,n,p}; in our HDivDivFreeHex basis they are not the trigonometric
separated solutions but polynomial combinations of order ≤ p.

K_SIBC^{3D body} requires the surface admittance bilinear form
    K^{SIBC}(u, v) = ∮_∂Ω Z_s(jω) (u·τ)·(v·τ) dS,
with Z_s = √(jωμ/σ).  By Green's first identity for vector fields,

    ∫_Ω curl u · curl v dV = ∫_Ω u · curl curl v dV + ∮_∂Ω (n × u) · curl v dS,

and for SIBC the contract is "no interior fields", so the volume term must
be replaced by either:
  (a) a Trefftz-type basis (free-space harmonic extension, curl curl v = 0
      inside ∂Ω) → only the surface integral remains;
  (b) a residual ∫_Ω u·∇²v dV that quantifies the deviation from Trefftz
      for the actual basis (e.g., HDivDivFreeHex).

The natural decomposition for 3D cuboid:
  - 2 large faces (z = ±t/2, area 17.72^2 mm^2)
  - 4 small faces (lateral, area 17.72 × 2 mm^2 each)
  - 12 edges + 8 corners (1-D / 0-D singular contributions)

For B applied in z-hat, only lateral faces carry surface currents (the
top/bottom faces have n parallel to H so n × H = 0).  K_SIBC^{3D, B∥z}
≈ (4 · 17.72 mm · 2 mm) · √(σ/μ) / (lateral perimeter factor) with
edge/corner corrections.

Verification target: |Y_residual_N(jω) → K_SIBC^{3D}(jω)^{-1/2}| at high ω,
agreement with NGSolve AC sweep ground truth at intermediate ω.

## Implementation order

  1. ✅ 1D slab (sanity, done)
  2. ✅ 2D rectangle 17.72 × 2 mm (done)
  3. ⏳ 2D Nagamine square a × a (next, needs symmetry-adapted basis)
  4. ⏳ 3D A1 cuboid 17.72 × 17.72 × 2 mm (Green + HDivDivFreeHex)

Done = "if results agree at 2D Nagamine and 3D cuboid the theory is
complete enough for IGTE 2026 + the eventual journal paper."
"""


CLN_SIBC_ORTHOGONAL_XFEM_COMPARISON = """
# XFEM vs SIBC vs augmented CLN: three abstraction layers of the same
# SIBC asymptote — when to use which

Three methods produce the high-frequency K_SIBC/√s asymptote of an
eddy-current admittance Y(s); they are **complementary, not competing**,
because each operates at a different abstraction layer of the same
physics.  The decision is governed by the *use case* (single-frequency
3D field solve vs frequency-sweep / time-domain / SPICE export) and
the *source structure* (port-driven vs volume-source).

## The three methods

| Method | Abstraction layer | What it modifies |
|--------|-------------------|------------------|
| Classical SIBC | FE boundary condition | Replaces interior with surface Robin BC Z_s = √(jωμ/σ) |
| XFEM (Hiruma 2023) | FE trial space | Enriches V_h with ψ(x) = exp(-γξ), γ = √(jωμσ) |
| Augmented CLN (Sugahara 2026) | Reduced matrix (ROM) | Adds 1×1 Schur block z(s) = (s+d)/(K_SIBC √s) to K_r(s) |

## Axis-by-axis comparison

| Axis | classical SIBC | XFEM | Augmented CLN |
|------|----------------|------|---------------|
| Conductor interior mesh | NO (boundary only) | YES coarse (δ not resolved) | NO (after ROM extraction) |
| DC limit | BROKEN (Z_s → 0) | exact (ψ → 1) | exact (z(0) = ∞) |
| K_SIBC/√s asymptote | exact | exact | exact (structural) |
| Wall band r/δ ~ 1 | leading-order, O(δ/L) | excellent | depends on base ROM basis |
| Sharp edges / corners | BROKEN — needs HOIBC | OK (Hiruma demonstrated) | absorbed in Foster modes |
| Curvature corrections | Senior/Mitzner HOIBC (c1=-H, c2=(K_g-H^2)/2, ...) | local ξ absorbs curvature | scalar K_SIBC = S √(σ/μ) |
| Volume-source (uniform J_0 = σ) | OK | OK | **BROKEN** — FE residual σ(|Ω| - b_r^T M_r^{-1} b_r) does not vanish |
| Port-driven (boundary excitation) | OK | OK | OK (Y_CLN(∞) → 0 naturally) |
| Frequency sweep | Z_s(ω) refactor per ω | γ(ω) re-enrich + refactor per ω | single ROM, s → jω is one eval |
| Time-domain / SPICE export | hard (Fourier round-trip) | hard (φ depends on ω) | natural (Schur 1-element → Foster-of-Cauers) |
| DOF count | smallest (∂Ω only) | mid (88 DOF reached 0.14% at r/δ=15 on Hiruma cyl.) | smallest ROM (4-5 DOF) |
| Implementation | trivial — every FE library | compound H1×H1 space + product-rule grad + bonus_intorder | one Schur row appended |
| Nonlinear σ(B), μ(B) | needs ESIM extension | XFEM can carry it inside | σ(T) bilinear coupling breaks (ω·τ_thermal ≫ 1) — use table |

## Decision tree

```
Q1.  Are you solving for spatial field distribution at one frequency?
     YES → XFEM (best 3D field, handles volume-source AND corners)
     NO  → continue.

Q2.  Do you need a frequency sweep, time-domain, or SPICE export?
     YES → Augmented CLN (single ROM, all-s validity)
     NO  → continue.

Q3.  Is the conductor body driven from a boundary port, or by an
     interior volume source?
     port-driven → Augmented CLN works (Y_CLN(∞) → 0).
     volume-source (e.g. uniform J_0 = σ, Hiruma 2023 cylinder, IH
                    workpiece in uniform inducing field) →
                    XFEM at FE level, OR re-formulate as port-driven.

Q4.  Pure industrial high-f only (DC accuracy not needed)?
     YES → classical SIBC is the simplest, cheapest option.
```

## Use-case routing (matches Paper 1 §VIII table)

| Application | Recommended |
|-------------|-------------|
| Single-frequency 3D, volume source, corner detail (Hiruma cyl., IH workpiece) | **XFEM** at FE level |
| Single-frequency 3D, port-driven, corner detail | **XFEM** |
| Frequency sweep + port-driven (R, L, port impedance vs f) | **Augmented CLN** (compact, fast) |
| Time-domain transient + port-driven (PWM, dead-time, pulse) | **Augmented CLN → SPICE** export |
| Multi-conductor BEM (5-turn coil cross-section, multi-strand) | **Augmented CLN** per element + BEM coupling |
| Pure industrial high-f, no DC | **classical SIBC** (lowest cost) |
| Volume-source workpiece that needs SPICE | **XFEM at FE → CLN extraction** (stack) OR re-formulate port-driven |

## Stacking: XFEM ∘ Augmented CLN (Paper 1 open direction iii)

For volume-source problems where SPICE export IS needed, the natural
hybrid is

```
                     XFEM      base CLN      Schur block
   continuous FEM → enriched → reduced  →  + (s+d)/(K_SIBC √s)
   ────────────    ────────   ────────     ───────────────
   physics solve   FE level   ROM         ROM tail
                   cures      reduction   cures rational
                   volume-                tail at high f
                   source FE
                   residual
```

Phase 3 of the β-1 benchmark (`docs/hiruma_xfem_comparison/`)
verified that:
  - Hiruma XFEM at 88 DOF gives 0.14% error at r/δ=15 on the
    volume-source cylinder — confirms FE-level cure of the high-f
    asymptote in volume-source.
  - Augmented CLN alone gives Y_CLN(∞) = 53 (44 DOF mesh) → 13 (3558
    DOF mesh), algebraic decay only — fails the C-axis ratio
    r(f) → 1 in volume-source.
Stacking (XFEM-enriched FE → ROM extraction → +Schur block) is the
correct route when the problem mixes volume-source physics AND
SPICE-export ambitions.

## Reproduction status

| Method | Implementation | Status |
|--------|----------------|--------|
| Hiruma XFEM ψ = exp(-γξ) on cylinder | the phase2 section of `docs/hiruma_xfem_comparison/hiruma_xfem_comparison.ipynb` (NGSolve compound H1×H1 + manual product-rule grad + bonus_intorder=10) | ✓ 0.14% at r/δ=15, 88 DOF |
| Augmented CLN Galerkin-Krylov | the phase3b section of `docs/hiruma_xfem_comparison/hiruma_xfem_comparison.ipynb` | ✓ d = 4.78×10⁴ rad/s (4% from analytical wall band); ✗ r(f) → 1 fails for volume-source |
| Classical SIBC Robin BC | `src/radia/panels/calc_fem_kelvin.py` | production-grade for port-driven |

## Key lesson (port-driven scope)

Theorem 1 of Paper 1 (single-DOF asymptote-preserving augmentation)
provides Y_R = Y_CLN + K_SIBC √s / (s+d) which preserves the SIBC
asymptote ONLY when Y_CLN(s → ∞) → 0.  This holds automatically when
the source vector b is supported on a boundary port (b_r^T M_r^{-1} b_r
= |∂Ω_port| in the discrete limit) but NOT when the source is a
volume body (J_0 = σ over Ω) — then Y_CLN(∞) = σ(|Ω| - b_r^T M_r^{-1} b_r)
is a FE-discretization residual that decays algebraically with mesh,
not exponentially.  Paper 1 §III now states this as a Remark
after Theorem 1; Paper 1 §VIII routes volume-source problems to XFEM.

## References

  - Hiruma, IEEE TMag 59(5), 2023 — XFEM ψ = exp(-γξ) on 2D cylinder.
  - Hiruma, 2024 — XFEM extension / corner handling.
  - Senior, IRE Trans. AP, 1962 — original HOIBC (c1 curvature term).
  - Mitzner, 1967 — Mitzner SIBC (c2, c3 terms).
  - Stoll, "Surface Impedance Concept" — classical SIBC reference.
  - Sugahara 2026 Paper I, Theorem 1 — uniqueness of single-DOF Schur
    augmentation, port-driven scope.
  - β-1 Phase 1-3 benchmark (2026-05-25):
      docs/hiruma_xfem_comparison/hiruma_xfem_comparison.ipynb
      memory/project_augmented_cln_scope_port_driven.md
"""


CLN_SIBC_EVRS_ESIM = """
# Production HCurl EVRS + local ESIM-SIBC mixed Galerkin

This is the current Radia production interpretation of the SIBC bridge.  It is
separate from the deprecated scalar Warburg closure described historically in
other sections of this document.

## Topology-aware space

Start from a high-order NGSolve `HCurl(p)` parent space and classify every
conductor face by its neighboring materials:

- conductor-air/exterior: eligible for surface-Omega/SIBC modes;
- conductor-conductor: retain conductive graph-cycle bridge modes;
- conductor-insulator: retain a non-SIBC blocking trace;
- ordinary conductor interior: compress by the Eddy-Visible Response Space
  (EVRS), then quotient directions carrying no independent `curl(T)` current.

Thus `p` controls parent-space admissibility near corners and edges, while the
retained EVRS/CLN rank controls reduced-model size.  An SIBC face is never
inferred merely from being a boundary face.

## Local ESIM surface operator

For a nonlinear magnetic conductor, solve the one-dimensional normal ESIM cell
at the local tangential-field magnitude and tabulate
`Z_s^ESIM(|H_t|,s)`.  Do not replace those values by their surface average.
For surface-current modes `K_m`, assemble

    [G_Gamma(Z_s)]_mn = int_Gamma Z_s(x,s) K_m(x) . K_n(x) dGamma.

Radia exposes this contract as `SurfaceImpedanceGram` and
`AssembleSurfaceImpedanceGram`.  The latter projects per-quadrature-sample
impedances and checks the surface basis and weights.  If
`Re Z_s(x,s) >= 0` pointwise, the Hermitian dissipative part of the Gram is
positive semidefinite.  Passing an unnamed array as if it were a modal diagonal
is rejected because it would hide the required surface projection.

`LocalESIMSurfaceModel` owns the validated BH curve and SciPy one-dimensional
cell controls. `SolveLocalESIMSurfaceVIM`, or
`TopologyAwareHybridVIM.solve_local_esim`, runs the fail-loud HCurl outer loop
for exactly one physical excitation.  For a fixed bulk HDiv magnetic operator,
`CoupledHDivHybridVIMSystem.solve_frequency_local_esim` runs the same update
around the complete magnetic/current mixed solve and returns
`CoupledHDivHCurlLocalESIMSolution`. The returned fields and
`SurfaceImpedanceGram` always belong to the same linearization point.

For design sweeps, `BuildLocalESIMSurfaceLUT` solves the cell grid offline and
`LocalESIMSurfaceLUT` persists `Z_s(f,|H_t|)` in a pickle-free NPZ archive.
Attach a loaded table with `LocalESIMSurfaceModel.with_lut`. The online Karl
iteration then reports `cell_solve_count = 0`. Interpolation is bilinear in
log-frequency and log-field for `log(|Z_s|)` and continuous impedance phase;
extrapolation is forbidden. A SHA-256 material/cell signature rejects reuse with
a different BH curve, conductivity, or numerical cell setup. Establish LUT
adequacy by field/frequency node refinement. Add state
dimensions, rather than forcing a two-dimensional LUT, for hysteresis,
temperature-dependent data not already represented by BH/sigma, rotational
magnetization, or multidimensional corner cells.
`ValidateLocalESIMSurfaceLUT` automates this check with direct cell solves at
field midpoints and at both stored and midpoint frequencies; choose the table
density from its reported maximum relative interpolation error.

The nonlinear solve uses an outer Karl iteration:

1. solve the reduced HCurl system for the current local surface Gram;
2. recover the solution-shaped `|H_t|` profile on conductor-air faces;
3. interpolate the converged ESIM cell table;
4. rebuild `G_Gamma(Z_s)` and repeat to the impedance-change tolerance.

The converged Gram is solved directly inside the fixed-operator
HDiv-MMM/HCurl coupling. Do not overclaim the scope:
`NgsolveBDMEddyBubbleVIM` currently builds the BDM-MMM response reduction from
scalar `mu_r`. A simultaneous constitutive iteration that updates an ordinary
nonlinear HDiv magnetic operator together with the local ESIM Gram is still an integration task.

## Exact mixed-Galerkin elimination

Partition the reduced operator into eliminated ordinary bulk coordinates `b`
and retained HDiv/bridge/SIBC coordinates `k`:

    A = [[A_bb, A_bk], [A_kb, A_kk]].

The frequency-specific trial and test maps are

    P_trial = [-A_bb^-1 A_bk; I],
    P_test  = [-A_bb^-T A_kb^T; I],
    P_test^T A P_trial = A_kk - A_kb A_bb^-1 A_bk.

For a nonzero eliminated-block right-hand side, reconstruction is affine:

    x = P_trial c + [A_bb^-1 f_b; 0],
    rhs_reduced = P_test^T f.

Omitting the particular lift loses the directly excited bulk field even though
the Schur matrix itself is correct.  `MixedGalerkinOrthogonalization.solve`
implements the complete affine reconstruction.

## Recorded p=6 broadband gate

`validation_test/cln/evrs_esim_sibc_mixed_notched_p6.json` records an mdx run
on a re-entrant notched conductor from 1 kHz to 1 MHz:

- 3557 active parent-HCurl DoFs;
- Krylov depth 8: 10 parent-T responses -> 8 current-Gram bulk modes;
- 14 conductor-cycle modes + 3 exterior-SIBC modes;
- 25 retained eddy coordinates, 0.703% of the parent DoFs;
- local-ESIM port difference from depth 22: at most 4.95e-8;
- mixed-Galerkin port difference from the direct reduced solve: at most 6.44e-16;
- uniform ESIM vs local ESIM: 0.229%--0.655%;
- linear SIBC vs local ESIM: 6.98%--38.31%;
- the 4-by-96 LUT has maximum direct-midpoint relative error 8.17e-4;
- all 16 outer solves converged in 10--17 updates with passive local Grams;
- every online outer update used the LUT and performed zero cell solves.

The linear-to-ESIM difference is a constitutive-model comparison, not an error
estimate.  Use local ESIM-SIBC when the high-frequency layer is thin and locally
one-dimensional.  Retain volume EVRS/VIM for low-frequency penetration,
genuinely three-dimensional edge/corner current, or rotational/hysteretic state
outside the ESIM cell model.
"""


def get_cln_sibc_orthogonal_documentation() -> str:
    """Return full markdown documentation for the CLN+SIBC orthogonal theory."""
    return "\n\n".join([
        CLN_SIBC_ORTHOGONAL_OVERVIEW,
        CLN_SIBC_ORTHOGONAL_RELATION_TO_MATSUO,
        CLN_SIBC_ORTHOGONAL_KURIYAMA,
        CLN_SIBC_ORTHOGONAL_MATH,
        CLN_SIBC_ORTHOGONAL_VERIFICATION,
        CLN_SIBC_ORTHOGONAL_NAGAMINE_LINK,
        CLN_SIBC_ORTHOGONAL_OUTLOOK,
        CLN_SIBC_ORTHOGONAL_XFEM_COMPARISON,
        CLN_SIBC_EVRS_ESIM,
    ])


def get_cln_sibc_orthogonal_section(name: str = "list") -> str:
    """Retrieve a specific section of the CLN+SIBC orthogonal documentation.

    Sections:
        "list"           - list available section names
        "overview"       - one-page summary of the construction
        "matsuo"         - relation to Matsuo expansion-point theory
        "kuriyama"       - Kuriyama 2019 multi-expansion canonical method
        "math"           - derivation, orthogonality, asymptote
        "verification"   - Mathematica results for 1D / 2D + sphere c_basis
        "nagamine"       - link to Nagamine 2026 verified extraction
        "outlook"        - 2D Nagamine square + 3D cuboid completion path
        "xfem_vs_sibc"   - XFEM / classical SIBC / augmented CLN decision
                           framework, port-driven scope, stacking strategy
        "evrs_esim"      - production HCurl EVRS, local ESIM surface Gram,
                           and exact affine mixed-Galerkin reconstruction
        "all"            - full documentation
    """
    table = {
        "overview":     CLN_SIBC_ORTHOGONAL_OVERVIEW,
        "matsuo":       CLN_SIBC_ORTHOGONAL_RELATION_TO_MATSUO,
        "kuriyama":     CLN_SIBC_ORTHOGONAL_KURIYAMA,
        "math":         CLN_SIBC_ORTHOGONAL_MATH,
        "verification": CLN_SIBC_ORTHOGONAL_VERIFICATION,
        "nagamine":     CLN_SIBC_ORTHOGONAL_NAGAMINE_LINK,
        "outlook":      CLN_SIBC_ORTHOGONAL_OUTLOOK,
        "xfem_vs_sibc": CLN_SIBC_ORTHOGONAL_XFEM_COMPARISON,
        "evrs_esim":    CLN_SIBC_EVRS_ESIM,
        "all":          get_cln_sibc_orthogonal_documentation(),
    }
    if name == "list":
        return "Available sections: " + ", ".join(
            ["overview", "matsuo", "kuriyama", "math", "verification",
             "nagamine", "outlook", "xfem_vs_sibc", "evrs_esim", "all"])
    if name not in table:
        return (f"Unknown section '{name}'. "
                f"Use 'list' to see available sections.")
    return table[name]
