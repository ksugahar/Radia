"""
ESIM IGTE Symposium 2026 paper knowledge.

Source of truth: s:/Radia/01_GitHub/docs/esim/*.md (7 files, ~3300 lines).
This module condenses the paper-relevant content into topic-keyed
strings that MCP tools serve.  When the canonical docs change, this
file should be re-curated -- it is NOT auto-generated.

Author note: this is the IGTE 2026 Symposium paper "Curved-Element
Per-Element Nonlinear Surface Impedance Method for Boundary Element
Eddy-Current Analysis" by K. Sugahara (Kindai University) +
collaborators.  Target journal: IEEE Trans. Magn. selected papers.

Canonical source docs:
  docs/esim/IGTE_PAPER.md           — paper outline (abstract -> refs)
  docs/esim/MATHEMATICAL_ANALYSIS.md — formulation
  docs/esim/IMPLEMENTATION.md        — architecture
  docs/esim/CROSS_VALIDATION.md      — numerical results
  docs/esim/R_MISMATCH_PEEC_VS_BEMA.md — focused diagnosis
  docs/esim/PEEC_PERFORMANCE_AND_R_ANALYSIS.md — deep dive
  docs/esim/USAGE.md                 — CLI guide
"""


ESIM_PAPER_OVERVIEW = """
# IGTE 2026 Paper: BEM Curved-Element Per-Element Nonlinear ESIM

**Target**: IGTE Symposium 2026 (Aug 2026, Graz) -> IEEE Trans. Magn.
**Authors**: K. Sugahara (Kindai University) + collaborators
**Status**: Draft (2026-05-18). LaTeX file at
  W:/02_学会資料/2026年度/2026_09_IGTE_Symposium/ESIM@久保田/digest/igte_symposium_2026.tex

## What the paper proposes

Hollaus, Kaltenbacher & Schoberl (IEEE TMAG 2025,
DOI 10.1109/TMAG.2025.3613932) introduced the nonlinear Effective
Surface Impedance Method (ESIM) for the FEM scalar magnetic potential
formulation.  This paper extends ESIM to BEM with three additions:

  1. BEM formulation (no air volume mesh; no Kelvin/PML truncation)
  2. Curved isoparametric Tri6 (O(h^3) vs flat O(h^2))
  3. Element-by-element (per-DOF) nonlinear Z_s
  4. Karl-style damped Picard on the ndof-dimensional Z_s vector

## Headline numerical result

At I_port = 100 A on a steel cylinder workpiece (sigma=2e6, BH curve,
50 kHz):

  Scalar Karl (mesh-RMS):  P_wp = 30.6 W   (converged, 15 iter)
  Per-element Karl:        P_wp = 45.4 W   (+48 % vs scalar)
  Z_s spatial ratio (max/min): 3.30x
  H_t range:                250-2951 A/m

The per-element formulation captures the spatial saturation pattern
that the scalar mesh-RMS averaging eliminates.  This is the
paper's headline contribution.

## When to use ESIM (induction heating context)

ESIM applies when:
  * workpiece is simply-connected (cylinder, plate, shaft, bar)
  * conductive enough that delta << workpiece thickness (thin-skin)
  * material has saturation behaviour (steel, electrical steel)

Out of scope:
  * DC / very low frequency where delta >= workpiece thickness
  * Hard-magnetic / branching-hysteresis materials
  * Anisotropic / laminated cores
  * Multiply-connected workpieces (rings, gear teeth as separate ring)
  * >100 MHz wave regime
"""


ESIM_PAPER_CONTRIBUTION = """
# Contribution Statement (paper abstract material)

## Background: Hollaus 2025 FEM ESIM

Hollaus-Kaltenbacher-Schoberl 2025 (DOI 10.1109/TMAG.2025.3613932)
introduced ESIM in a magnetic scalar potential FEM formulation:

  1. 1-D cell problem in radial direction returns Z_s(|H_t|) for the
     local surface tangential field.
  2. Outer FEM scalar potential solve uses Z_s as a Robin BC.
  3. Karl iteration (damped Picard) couples cell <-> outer.

  Hollaus uses a SCALAR Z_s (single value per Karl outer iteration,
  derived from the mesh-RMS H_t).

## Our extension to BEM with per-element Z_s

This paper extends to BEM with three new contributions:

### 1. BEM scalar magnetic potential formulation
   The workpiece is represented by a scalar potential phi on its
   surface (BIE): (1/2 M + DL + (jw/Z_s) SL M^-1 K) phi = phi_inc.
   No air volume mesh; no Kelvin/PML truncation.  Suitable for
   unbounded domains.

### 2. Curved isoparametric Tri6 elements
   Lagrange-P2 basis with isoparametric mapping x(u,v) =
   sum_i phi_i^P2 * x_i (6 nodes per triangle).  O(h^3) convergence
   for smooth Gamma + smooth Z_s, vs O(h^2) for flat Tri3.
   Cubit's `radia_export netgen "f.vol" order 2` embeds the curved
   nodes; the in-tree assembler at `src/radia/bem_sibc_solver.py`
   builds SL, DL, M, K on Tri6 with Sauter-Schwab Duffy quadrature.

### 3. Per-element (per-DOF) nonlinear Z_s
   Each DOF i carries an independent Z_s[i] = E(|H_t|[i]) computed
   from the LOCAL H_t.  Replaces Hollaus's single scalar Z_s.

   The BIE row i becomes
     (1/2 M[i,:] + DL[i,:] + (jw/Z_s[i]) (SL M^-1 K)[i,:]) phi = phi_inc[i]
   so Z_s enters the system matrix by direct row-scaling.

   Per-DOF |H_t| extracted variationally:
     |H_t|^2_i = |phi_i (K phi)_i| / (M 1)_i
   (i-th summand of phi^T K phi divided by i-th lumped mass).
   Consistent with the global mesh-RMS in the uniform-gradient limit.

### 4. Karl iteration on the ndof-dimensional Z_s vector
   Same damped Picard structure as Hollaus, but on the vector Z_s.
   Same geometric convergence rate; same Lipschitz analysis applies
   per-component.

## Why these four pillars matter

For induction heating workpieces, the BH curve introduces strong
spatial saturation contrasts: regions near the coil's current
concentration see |H_t| >> BH knee while distant regions stay in
the linear-mu regime.  Scalar mesh-RMS Z_s averages these away;
per-element captures them.

The benchmark in CROSS_VALIDATION.md § 6b shows 48 % under-estimate
of P_wp by the scalar method vs per-element on a steel cylinder at
I = 100 A.
"""


ESIM_PAPER_FORMULATION = """
# Mathematical Formulation (paper § II material)

## Cell problem (1-D radial PDE for the workpiece skin layer)

Inside the workpiece (homogeneous conductor, MQS limit, azimuthal
ansatz H = H_phi(r) phi_hat):

  Strong form:
    (rho/r) d/dr[r dH/dr] + jw mu(|H|) H = 0   in r in [0, R]
    H'(0) = 0       (axis regularity)
    H(R) = H_t      (driven by outer surface tangential field)

  Weak form (with radial measure r dr):
    int_0^R rho r (dH/dr)(dv*/dr) dr
    - jw int_0^R mu(|H|) H v* r dr  =  0
    forall v in H^1_0([0, R])

Implementation: src/radia/esim_cell_problem.py (FD on uniform mesh
recovers exact Galerkin P1 stencil with midpoint quadrature on
stiffness + vertex quadrature on mass).

Output: Z_s(H_t) = 2(P' + j Q') / |H_t|^2 (Ohm).
Coupled to outer via inner Picard iteration (Picard on mu(|H|)).

## Outer scalar BIE (workpiece surface)

  (1/2 M + DL + (jw/Z_s) SL M^-1 K) phi = phi_inc on Gamma_wp

where SL, DL, M, K are the single-layer, double-layer, mass,
tangential-gradient stiffness operators on Gamma_wp.  phi_inc is
the incident scalar potential from the coil (PEEC filaments via
Biot-Savart, or BEM-A surface current).

Implementation: src/radia/bem_sibc_solver.py:ScalarBIESIBCSolver.

Solver: dense LU (small Gamma, <~ 25k DOFs) or HACApK-ACA-compressed
GMRES (large Gamma).

## Per-element Z_s as Z_s in R^{ndof}

Replace scalar Z_s with Z_s[i] for each DOF i.  System matrix row i:

  A[i, :] = 1/2 M[i, :] - DL[i, :] + (jw/Z_s[i]) (SL M^-1 K)[i, :]

This is implemented by element-wise multiplication on the diagonal
of (jw/Z_s) before forming the SL contribution.

## Per-DOF |H_t| extraction

The H1 BIE solves for phi at H1 P1 DOFs.  Local |H_t|^2 at DOF i:

  |H_t|^2_i = |phi_i * (K phi)_i| / (M @ 1)_i

This is the i-th summand of the global phi^T K phi (which equals
A_wp * |H_t|^2_rms in the scalar formulation) divided by the
lumped-mass i-th entry (which equals the effective area of basis
function i).

abs() needed because each per-DOF summand can be signed; only the
sum is provably non-negative.  In the uniform-gradient limit the
per-DOF formula recovers the scalar mesh-RMS exactly.

## Lorentz reciprocity for Delta L (paper "Telegen" terminology)

The induced workpiece surface current J_s gives a port impedance
perturbation Delta Z = jw Lambda, where

  Lambda = (1/I_port^2) * integral_Gamma phi (n . B_inc) dS
         = (1/I_port^2) * integral_Gamma J_s . A_inc dS   (continuum identical)

This is the **Lorentz reciprocity reaction integral** (Harrington
1961 sec.3.8; lab nicknames it "Telegen" but the literature term
is Lorentz).  The first form (phi . n.B) is gauge-invariant; the
second form (J_s . A) is gauge-invariant in continuum but
gauge-dependent in a P1 FE discretisation.  Radia uses the first.

Decomposition:
  Delta L  = Re(Lambda)
  Delta R  = -w * Im(Lambda)
  Delta Z  = jw * Lambda  =  Delta R + jw * Delta L

Implementation:
  src/radia/workpiece_surface.py:delta_L_telegen_phiB
  Each workpiece triangle T contributes
    Lambda|_T = phi_avg_T * (n_T . B_inc(c_T)) * A_T / I_port^2
  with single-point centroid quadrature on B_inc (smooth on the
  workpiece scale) and exact integration of phi_h (P1 -> phi_avg
  is the corner mean).

## Karl iteration on the Z_s vector

For k = 0, 1, ...:
  1. Solve BIE at current Z_s^(k)  ->  phi^(k)
  2. Extract |H_t|^(k)[i] per DOF (variational formula above)
  3. Cell-solver call per DOF:
       Z_s_new[i] = E(|H_t|^(k)[i])     (a 1-D PDE solve per DOF)
  4. Damped update:
       Z_s^(k+1)[i] = alpha * Z_s_new[i] + (1 - alpha) * Z_s^(k)[i]
  5. Stop when max_i |Z_s^(k+1)[i] - Z_s^(k)[i]| / |Z_s^(k)[i]| < tol

Default alpha = 0.5, tol = 1e-3.  Empirical Lipschitz L ~ 1 in
linear regime; alpha * L = 0.5 gives strict contraction.

In deep saturation (H_t straddles BH knee), L can rise above 1 at
some DOFs and Karl fails to converge at alpha = 0.5.  Workarounds
under investigation: lower alpha (0.3 or 0.2), Anderson acceleration
on Z_s vector (m = 2-3 history), per-DOF adaptive alpha_i.
"""


ESIM_PAPER_VALIDATION = """
# Numerical Validation (paper § IV material)

All numbers below are REAL benchmark data from radia v4.55.3.

## Cross-check A: Linear-mu Bessel reference

Cell solver matches the Wakao-Igarashi-Fujiwara-Kameari Part 5
analytical Bessel I_0/I_1 reference to < 0.13 % on a cylindrical
steel workpiece (mu_r = 100) across 4 frequencies (1 kHz - 1 MHz,
xi = R/delta from 4 to 140).

  xi    f     Z_s_anal     Z_s_num    err
   4    1k    7.50e-3      7.50e-3   <1e-5
  14   10k    2.49e-2      2.49e-2   <1e-4
  45  100k    7.95e-2      7.95e-2   4e-4
 140    1M    2.50e-1      2.50e-1   1.3e-3

Source: examples/ih_esim_benchmark/{benchmark.py, analytical_bessel_baseline.py}
Reproducible: python examples/ih_esim_benchmark/analytical_bessel_baseline.py

## Cross-check B: Three independent paths (PEEC-BEM / FEM-K / FEM-CM)

Same geometry, 4 frequencies; data from results.json:

  f [kHz] | Path        | |Z_s| [Omega] | L_total [nH] | P_wp [uW] | iter
  -------+-------------+--------------+--------------+-----------+-----
   10    | PEEC-BEM    |   1.103e-2   |    97.77     |    21.1   |  7
   10    | FEM-Kelvin  |   1.101e-2   |   163.27     |    19.0   |  7
   10    | FEM-coilmesh|   1.097e-2   |    93.21     |    18.9   |  7
   50    | PEEC-BEM    |   2.572e-2   |    97.48     |   162.5   |  6
   50    | FEM-Kelvin  |   2.543e-2   |   161.74     |   148.2   |  6
   50    | FEM-coilmesh|   2.321e-2   |    35.83     |    66.6   |  7
  100    | PEEC-BEM    |   3.770e-2   |    97.31     |   305.8   |  5
  100    | FEM-Kelvin  |   3.722e-2   |   160.76     |   282.3   |  5
  100    | FEM-coilmesh|   3.200e-2   |    26.34     |   101.7   |  7
  500    | PEEC-BEM    |   1.314e-1   |    97.08     |   561.6   |  4
  500    | FEM-Kelvin  |   1.293e-1   |   159.41     |   527.9   |  4
  500    | FEM-coilmesh|   8.55e-2    |     9.35     |    95.3   |  7

PEEC-BEM <-> FEM-Kelvin agree on |Z_s| to < 2 % at all 4 freq.
FEM-coilmesh tracks at 10 kHz then diverges (53 % at 500 kHz) due
to coil mesh under-resolution (coil_h_max ~ 5.4 mm vs delta_coil
~ 0.09 mm @ 500 kHz).  Not a Karl-loop bug.

Karl convergence at 50 kHz PEEC-BEM: dZ/Z decays geometrically
with ratio ~0.4, converges in 6 iter (default alpha = 0.5).

## Cross-check C: External 2-D axisymmetric reference

Cu cylinder + gapped torus @ 7 kHz, I = 100 A.  Reference value
P_wp^ref ~ 6.63e-5 W from external 2-D axisym Mathematica solve.

  Path         | P_wp [W]    | Delta vs ref
  -------------+-------------+--------------
  PEEC-BEM     | 6.48e-5     |  -2.3 %
  FEM-coilmesh | 6.543e-5    |  -1.3 %

Golden tests: tests/panels/golden/peec_bem_coarse_7kHz_Cu.json,
fem_coilmesh_gapped_fine_7kHz_Cu.json (both PASS).

## HEADLINE: Per-element vs scalar Z_s

Steel cylinder workpiece (sigma=2e6, BH curve em_sample_bh.txt) at
50 kHz, I_port = 100 A (chosen to push H_t through the BH knee
at ~1000 A/m).

  Quantity          | Scalar Karl    | Per-element Karl
  ------------------+----------------+------------------
  Convergence       | 15 iter, dZ=1e-3 | 15 iter, dZ_max=0.34 (NOT converged)
  P_wp [W]          | 30.60          | 45.36   (+48 %)
  L_total [nH]      | 87.96          | 84.76   (-3.6 %)
  Delta L [nH]      | -10.01         | -13.21  (+32 % more negative)
  H_t_rms [A/m]     | 680.7 (single) | 1193 mean, 2951 max
  |Z_s|             | 2.25e-2        | 1.13e-2 to 3.74e-2 (3.30x ratio)

Physical interpretation: surface H_t varies 12x across the workpiece
(250 A/m far face -> 2951 A/m near face).  Steel's BH knee at
~1000 A/m causes mu_r to drop from ~100 (linear regime) to ~5
(saturated).  Z_s = (1+j) rho/delta sqrt(mu_r) thus varies by 3.3x
in magnitude.  Scalar mesh-RMS Z_s averages this away; per-element
samples it correctly.

Why P_wp differs by 48 %: P_wp ~ Re(Z_s) * |H_t|^2 is nonlinear in
H_t; the integral over Gamma of Re(Z_s(|H_t|)) * |H_t|^2 does not
equal Re(Z_s(<|H_t|>_rms)) * <|H_t|^2>_rms * A_wp.

Open observation: per-element Karl with alpha=0.5 does NOT
converge in 15 iter at this drive level (dZ_max stalls at 0.34).
Root cause: BH-knee-straddling DOFs have local Lipschitz L ~ 2-3,
exceeding the contraction condition alpha*L < 1.  Workarounds
(alpha=0.3 + max_iter=60, Anderson acceleration) under investigation;
this is itself a publishable observation.

## Curve-order benchmark

Same Cu cylinder mesh (2150 BND tris), 3 (curve_order, basis_order)
combinations, linear SIBC at 50 kHz, I=1 A:

  Case    | curve  | basis | P_wp [uW] | t_asm | t_solve
  --------+--------+-------+-----------+-------+--------
  p1_h1   | flat   | P1    | 167.13    | 1.27s | 0.19s
  p1_h2   | flat   | P2    | 168.24    | 4.96s | 7.19s   (+0.7 %)
  p2_h2   | curved | P2    | 168.81    | 4.96s | 7.21s   (+0.3 %)

Smooth-cylinder workpiece -> curvature effect < 1 % because the
mesh is already in the asymptotic regime.  Larger effect expected
on highly-curved geometries (gear roots, fillets) -- future
benchmark.

## Reproducibility

  # Cross-check A (Bessel parity)
  python examples/ih_esim_benchmark/analytical_bessel_baseline.py

  # Cross-check B (three-path consistency)
  python examples/ih_esim_benchmark/benchmark.py --frequencies "1e4,5e4,1e5,5e5"

  # Cross-check C (2-D axisym)
  pytest tests/panels/test_peec_bem_golden.py tests/panels/test_fem_coilmesh_golden.py

  # HEADLINE per-element vs scalar
  for flag in "" "--esim-per-panel"; do
    python src/radia/panels/calc_inductance.py \\
      --coil-step src/radia/panels/samples/ih_fem_kelvin_demo_coil.step \\
      --coil-solver peec \\
      --vol src/radia/panels/samples/ih_bem_sample_p1.vol --wp-label sibc \\
      --sigma 2e6 --mu-r 100 --half-thickness 0.005 \\
      --frequency 50000 --current 100.0 --coil-sigma 5.8e7 \\
      --impedance-model esim --bh-file src/radia/panels/samples/em_sample_bh.txt \\
      --esim-max-iter 15 --esim-tol 1e-3 --esim-relax 0.5 \\
      --h1-order 1 --wp-bem-backend intree-dense $flag
  done
"""


ESIM_PAPER_REFERENCES = """
# Paper References (Hollaus headlined)

1. K. Hollaus, M. Kaltenbacher, J. Schöberl, "A Nonlinear Effective
   Surface Impedance in a Magnetic Scalar Potential Formulation,"
   IEEE Trans. Magn., 2025.  DOI: 10.1109/TMAG.2025.3613932.
   ** Canonical reference for the ESIM cell-problem method and Karl
   iteration.  The "Karl" in "Karl iteration" is Karl Hollaus's
   first name (lab shorthand; literature term: "Hollaus-type damped
   Picard relaxation"). **

2. L. Krähenbühl, D. Muller, "Thin layers in electrical engineering
   — Example of shell models in analysing eddy-currents by boundary
   and finite element methods," IEEE Trans. Magn. 29(2), 1993.
   Foundational thin-shell SIBC reference.

3. J. D. Lavers, P. P. Biringer, "An efficient calculation of
   effective surface impedance for nonlinear ferromagnetic
   materials," IEEE Trans. Magn. 21(5), 1985.  Closed-form
   envelope reference for cylinder + slab nonlinear ESIM.

4. R. L. Stoll, The Analysis of Eddy Currents, Oxford University
   Press, 1974.  Reference for analytical Bessel function
   comparisons.

5. R. F. Harrington, Time-Harmonic Electromagnetic Fields,
   McGraw-Hill, 1961.  Section 3.8 — Lorentz reciprocity /
   reaction integral.  NOT "Telegen" (Telegen is a network
   theorem); the integral form is correctly called Lorentz
   reciprocity.

6. E. Dlala, A. Belahcen, A. Arkkio, "Optimal Convergence of the
   Fixed-Point Method for Nonlinear Eddy-Current Problems,"
   IEEE Trans. Magn. 44(6), 2008, pp. 1318-1321.  Optimal
   contraction-factor analysis for the same Picard scheme; informs
   the alpha = 0.5 default.

7. S. Yuferev, N. Ida, Surface Impedance Boundary Conditions: A
   Comprehensive Approach, CRC Press, 2009.  Textbook reference
   for nonlinear SIBC iteration schemes.

8. S. Bilicz, Z. Badics, J. Pávó, "Nonlocal surface impedance
   boundary condition for wide-band eddy-current problems,"
   ISEM 2023.  Wide-band nonlocal extension (deferred to roadmap).

9. S. Wakao, H. Igarashi, K. Fujiwara, A. Kameari, "Various
   Verifications of Eddy Current Analysis (Parts 1-9)," T.IEE
   Japan, 2008-2018.  Part 5 used for the Bessel cylinder
   baseline cross-check.
"""


ESIM_PAPER_IMPLEMENTATION = """
# Implementation (paper § III material)

## Three-solver architecture

All three production paths use the same Hollaus cell solver but
differ in outer-solver formulation:

  calc_inductance.py        — BEM scalar BIE workpiece (paper's focus)
  calc_fem_kelvin.py        — HCurl FEM workpiece + Kelvin (FEM analogue)
  calc_fem_coilmesh.py      — full FEM A-V volumetric coil (anchor reference)

The paper's headline contribution lives in calc_inductance.py.
calc_fem_kelvin.py and calc_fem_coilmesh.py serve as cross-validation
references and as "the per-element extension also works in FEM"
supporting evidence.

## Cell solver internals (src/radia/esim_cell_problem.py)

Two classes:
  ESIMFiniteSlabSolver (line 339) — used by all 3 production paths
  ESIMCellProblemSolver (line 816) — legacy infinite-domain solver

Discretisation: finite difference on uniform radial mesh
(n_nodes=100 default), tridiagonal system solved with
scipy.sparse.linalg.spsolve.

Axis singularity (cylinder mode, r=0): L'Hopital limit gives
diag_main[0] = -4 rho/h^2 + jw mu[0], diag_upper[0] = 4 rho/h^2.

Initial guess: analytical Bessel I_0 (cylinder) or cosh (slab) for
the constant-mu linear case, then Picard iterates on mu(|H|).

Inner Picard convergence: tol=1e-6, max_iter=50, relaxation=0.5.

Output dict: Z (complex impedance), P_prime (active power density),
Q_prime (reactive), P_magnetic (hysteretic loss for complex mu),
R_ratio (R_ac/R_dc; used by peec_bundle), H_solution, converged,
iterations.

## BEM-SIBC scalar solver (src/radia/bem_sibc_solver.py)

class ScalarBIESIBCSolver.

Assembly:
  - SL, DL, M, K matrices on H1 P1 (basis_order=1) or H1 P2
    (basis_order=2) on the workpiece surface
  - Curved isoparametric mapping when intree_geom_order=2 (Tri6)
  - Sauter-Schwab Duffy quadrature for the singular (i,j) pair
    (intree_singular_n_q=6)
  - Gauss-Legendre degree 7 for regular pairs

solve(phi_inc, Z_s, omega): assembles
  A_sys = 0.5 M - DL + (jw/Z_s) SL M^-1 K

For per-DOF Z_s, Z_s is an ndarray[ndof]:
  A_sys[i, :] = 0.5 M[i, :] - DL[i, :] + (jw/Z_s[i]) (SL M^-1 K)[i, :]
implemented via diag-multiplication.

Dense LU at small ndof (< 25000); HACApK-ACA + GMRES at large ndof.

Note: HACApK does NOT yet support per-DOF Z_s ndarray (only scalar);
this is a roadmap limitation (--esim-per-panel requires
--wp-bem-backend intree-dense).

## Karl iteration on Z_s vector

Implemented at calc_inductance.py:_solve_workpiece_weak_coupled
(lines 417-840), specifically lines 598-712 for the Karl loop body.

Seed: Z_s_seed = E(H_t=5 A/m); seed every DOF with the same scalar
to start.

Each iter:
  1. solve_bem(phi_inc, Z_s) returns phi_vec
  2. Extract per-DOF |H_t| via variational formula
     (lines 652-670: phi_i * (K phi)_i / M_lump_i)
  3. Cell solver per DOF (lines 660-665): N_DOF calls to esim.solve
  4. Damped update (line 671-672)
  5. Convergence check (line 686-687): max_i |dZ[i]| / |Z_s[i]| < tol

Cost: N_DOF * ~5ms per Karl iter (cell solver is fast).  For
N_DOF ~ 166 (typical IGTE benchmark) the overhead is < 1 s per
iter; dominates only above ~5000 DOFs.

## Wall-time breakdown (LAB, Windows, MKL)

  STEP -> filaments (PEEC topology):  ~3 s  (constant across freq)
  BIE assembly (P1, 166 DOFs):         ~0.2 s
  BIE assembly (P2 Tri6, 580 DOFs):    ~5 s
  BIE solve (dense LU, 166 DOFs):      ~0.05 s
  BIE solve (dense LU, 580 DOFs):      ~7 s
  Cell solver (single call):            ~5 ms
  Karl loop (5-7 iter):                 ~1 s + assembly cost

Total wall time: ~5 s for 1 PEEC-BEM evaluation (P1, n_DOF~166).
Frequency sweep over N points = N * 5 s (no caching today).
Topology caching for sweeps (roadmap): ~3 s savings per freq beyond
the first.
"""


ESIM_PAPER_CAVEATS = """
# Limitations and Scope (paper § V material)

## In-scope regimes

### Frequency
  10 Hz - 1 MHz : production IH range; SIBC valid

### Material (workpiece)
  Non-magnetic conductors (Cu, Al, brass): use --impedance-model sibc
  Linear ferromagnetic (constant mu_r > 1): use --impedance-model sibc --mu-r M
  Soft-magnetic with single-valued BH curve (silicon steel,
    electrical steel, soft iron): use --impedance-model esim --bh-file BH.txt

### Geometry (workpiece)
  Idealised cylindrical / planar workpieces: scalar BIE perfect
  Complex shapes with locally-varying curvature: scalar BIE OK; per-element
    captures spatial variation

### Coupling
  Weak coupling (coil current fixed, workpiece reacts): production mode

## Out of scope

### Frequency
  DC / < 1 Hz: delta >> workpiece size; SIBC limit broken
  > 100 MHz: wave effects in air (displacement current);
    use ngsolve.bem with Helmholtz kernel

### Material
  Hard-magnetic / hysteretic with branching BH loops (permanent magnets,
    detailed Preisach steel modelling): scalar mu(|H|) cannot capture
  Anisotropic / laminated cores: scalar cell solver is isotropic;
    use Hollaus-Hannukainen-Hiptmair homogenisation method
  Lossy ferrite (Mn-Zn, Ni-Zn) with complex mu = mu' - j mu":
    cell solver SUPPORTS this (complex_mu kwarg in
    ESIMCellProblemSolver) but CLI does NOT yet expose it

### Geometry
  Sharp corners (gear roots, fillets): per-element BEM should help;
    benchmark not yet run -- future work
  Very thin workpieces (a < delta): thin-shell SIBC formulation needed
    (Krähenbühl-Muller 1993)
  Multiply-connected workpieces (rings, gear teeth as separate rings):
    scalar phi becomes multi-valued; use cut-surface scheme OR
    vector RWG workpiece formulation -- out of scope for this paper

### Coupling
  Strong coupling (workpiece reaction modifies coil current
    distribution): use calc_fem_coilmesh.py full A-V

## Known open items (paper future work)

  1. Per-element Karl convergence at deep saturation:
     alpha=0.5 fails at high drive (BH-knee straddling DOFs).
     Workarounds: Anderson acceleration on Z_s vector,
     per-DOF adaptive alpha_i.

  2. HACApK + per-DOF Z_s: HACApK backend currently scalar-only;
     extending to per-DOF requires assembler modification.
     Effort: ~1 week.

  3. PEEC vs BEM-A coil R discrepancy: full complex Leontovich SIBC
     in BEM-A saddle matrix.  Effort: ~1 week.
     (see R_MISMATCH_PEEC_VS_BEMA.md for full diagnosis)

  4. Wide-band extension: nonlocal SIBC (Bilicz-Badics-Pavo 2023)
     for cases where SIBC limit breaks (very thin / very thick skin).

  5. Measurement validation: lab-built steel cylinder + B-H meter.
"""


TOPICS = {
    "all":            None,  # special: concatenate all
    "overview":       ESIM_PAPER_OVERVIEW,
    "contribution":   ESIM_PAPER_CONTRIBUTION,
    "formulation":    ESIM_PAPER_FORMULATION,
    "implementation": ESIM_PAPER_IMPLEMENTATION,
    "validation":     ESIM_PAPER_VALIDATION,
    "caveats":        ESIM_PAPER_CAVEATS,
    "references":     ESIM_PAPER_REFERENCES,
}


def get_esim_paper_documentation(topic: str = "all") -> str:
    """Return paper content for the requested topic."""
    topic = (topic or "all").strip().lower()
    if topic == "all":
        return "\n\n".join(
            f"## {key}\n{value}" if value else ""
            for key, value in TOPICS.items() if key != "all"
        )
    if topic in TOPICS:
        return TOPICS[topic]
    # Unknown topic: list available
    return (
        f"Unknown topic: {topic!r}.\n\n"
        f"Available topics: {', '.join(k for k in TOPICS if k != 'all')}.\n"
        f"Pass topic='all' for the complete paper outline."
    )
