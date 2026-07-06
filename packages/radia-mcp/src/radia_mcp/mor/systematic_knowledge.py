"""Systematic MOR (Model Order Reduction) knowledge for radia_mcp.mor.

Distilled from the 3-volume MOR Handbook in the Sugahara lab library
(Volumes 1+2+3 totaling 1221 pages, deGruyter 2020-2021, edited by
Benner-Grivet-Talocia-Quarteroni-Rozza-Schilders-Silveira), plus the
Kiss-Orosz 2024 Energies review specifically for rotating machines,
plus the Sugahara-Wakao CLN papers already covered in
`radia_mcp.mor.cln_knowledge`.

Coverage:
  mor_taxonomy          — the full MOR landscape, three families
  projection_pod_rb     — POD + Reduced Basis (Vol 2 Ch 2, 4)
  projection_krylov     — Krylov / moment-matching (Vol 1 Ch 3)
  pgd                   — Proper Generalized Decomposition (Vol 2 Ch 3)
  hyperreduction        — DEIM / EIM / hyper-reduction (Vol 2 Ch 5)
  system_theoretic_bt   — Balanced Truncation, Hankel norm (Vol 1 Ch 2)
  data_driven_dmd_oi    — DMD, Operator Inference, Loewner (Vol 1 Ch 6, Vol 2 Ch 7)
  parametric_pmor       — Parametric MOR (Vol 2 Ch 1)
  cln_sugahara          — CLN cross-reference + lab heritage
  em_specific_ioan      — Ioan Ch 5 EM-specific MOR (Vol 3)
  rotating_machines_kiss_orosz_2024 — recent review specifically
                          for motors (POD/PGD/OIM/CLN comparison)
  software_lab          — MOR software (Vol 3 Ch 13) + lab tools
  lab_recommendation    — decision guide for radia + NGSolve use
"""

from __future__ import annotations


MOR_TAXONOMY = """\
## MOR taxonomy — three families

MOR techniques fall into three structurally distinct families.  Each
makes different assumptions about what is known about the system
being reduced.

### Family 1: Projection-based MOR (intrusive, knows the equations)

Assumes you have the **full-order system equations** (the FE matrices
K, M, C and right-hand side b).  Build a reduced basis V of size r
much less than ndof; project:

  K_r = V^T K V,    M_r = V^T M V,    b_r = V^T b

Solve the reduced system, recover full solution as u = V u_r.

| Method | Basis V construction | Reference |
|--------|------------------------|-----------|
| POD | SVD of snapshot matrix | Vol 2 Ch 2 |
| Reduced Basis (RB) | Greedy parameter selection + Galerkin | Vol 2 Ch 4 |
| Krylov / Arnoldi | K^{-1}M iteration around expansion point | Vol 1 Ch 3 |
| Lanczos (symmetric) | Three-term recurrence | Vol 1 Ch 3 |
| PGD | Separable greedy expansion | Vol 2 Ch 3 |
| **CLN** (Sugahara) | Lanczos on PEEC + magnetic coupling | `cln_knowledge` |

### Family 2: System-theoretic MOR (intrusive, knows transfer function)

Treats the full system as an LTI state-space model `(A, B, C, D)`.
Reduce to keep input-output behavior (transfer function `H(s) =
C(sI-A)^{-1}B + D`) close.

| Method | Idea |
|--------|------|
| Balanced Truncation (BT) | Equalize observability and controllability Gramians, truncate small Hankel singular values |
| Optimal Hankel norm | Best approximation in Hankel norm |
| H_inf optimal | Minimize sup_omega ||H(jomega) - H_r(jomega)|| |
| H_2 optimal (IRKA) | Minimize integral_omega ||H - H_r||^2 domega |

Strength: a priori guaranteed error bounds (especially BT).
Weakness: assumes LTI; doesn't extend to nonlinear directly.

### Family 3: Data-driven MOR (non-intrusive, knows snapshots only)

Assumes you have **snapshots / simulation data** of the full system
but not the equations themselves.  Extracts a low-dimensional model
purely from data.  Bridge to ML.

| Method | Idea |
|--------|------|
| DMD (Dynamic Mode Decomposition) | Best-fit linear operator on snapshot pairs |
| Operator Inference (OI) | Least-squares fit of polynomial operator |
| Loewner framework | Construct rational interpolant from frequency-response data |
| Manifold interpolation | Geodesic in matrix manifold |
| Vector fitting | Pole-residue rational fit |
| Kernel / Kriging surrogates | Gaussian process regression |
| Sparse manifold (PINN MOR) | Autoencoder + ODE for latent dynamics |

The boundary with Family 1 blurs (e.g. POD computed from snapshots
is technically data-driven but yields a projection basis).

### When to use which family

- Have access to FE matrices, parametric study → Family 1 (POD/RB)
- Have only black-box simulator, large parameter sweep → Family 3 (DMD/OI)
- Need certified error bounds → Family 1 RB (a posteriori bounds) or
  Family 2 BT (a priori bounds)
- Need real-time inference for control → Family 1 (if affine
  parameter dependence) or Family 3 (otherwise)
"""

PROJECTION_POD_RB = """\
## POD and Reduced Basis methods (Vol 2 Ch 2, 4)

### Proper Orthogonal Decomposition (POD)

Given M snapshots `{u_i ∈ R^N}_{i=1}^M` of the full-order system at
M parameter values (or time instances), assemble the snapshot
matrix X = [u_1, u_2, ..., u_M] ∈ R^{N x M}.

Compute SVD: X = U Σ V^T.  Keep first r columns of U as the
**POD basis**: V_POD = [u_1, ..., u_r].

Optimality (Karhunen-Loève): the r-dim POD basis MINIMIZES the
total snapshot reconstruction error in Frobenius norm among ALL
r-dim subspaces.

### Choosing r

Standard criterion: keep r such that the **energy ratio**

  E(r) = Σ_{i=1}^r σ_i² / Σ_{i=1}^M σ_i²  > 0.999

(99.9% of snapshot energy captured).  For smooth problems r ~ 10-50
typically suffices.  For convection-dominated or shock problems,
POD struggles (large Kolmogorov n-width).

### Reduced Basis (RB) method

Greedy alternative to POD: instead of computing snapshots at a dense
grid of parameters, build the basis **incrementally**.  At each
step, pick the parameter mu_k+1 that maximizes a cheap error
estimator on the current basis, solve full-order problem at mu_k+1,
orthogonalize the new snapshot, add to basis.

Advantages over POD:
1. **Adaptive parameter selection** — focus on hard parts of
   parameter space, skip easy ones.
2. **A posteriori error bounds** — rigorous bounds via residual
   norm + stability constant.
3. **Offline-online decomposition** — all expensive computation
   moved to offline phase; online phase is O(r^3).

Requires:
- Affine parameter dependence: A(mu) = Σ_q theta_q(mu) A_q
- Stability constant alpha_LB (lower bound on coercivity constant
  — computed by Successive Constraint Method)

### Application to EM / motor analysis

POD/RB is the natural choice for:
- Parametric magnetostatic studies (vary coil current, mu_iron, gap)
- Reduced-fidelity inverse design (compute many candidate designs)
- Frequency sweeps in linear regime

For nonlinear ν(B), POD/RB require **hyper-reduction** (DEIM/EIM,
see `hyperreduction` topic) because the projection `V^T f(V u_r)`
still costs O(N) per nonlinear evaluation without it.

Sugahara lab usage:
- Hollaus-Schöberl-Schöbinger 2020 MSFEM+MOR uses POD on MSFEM
  snapshots (cf. `radia_mcp.motor.hollaus_eddy_knowledge.msfem_plus_mor`)
"""

PROJECTION_KRYLOV = """\
## Krylov / Arnoldi / Lanczos moment-matching (Vol 1 Ch 3)

### Moment-matching idea

The transfer function H(s) of an LTI system can be Taylor-expanded
around an expansion point s_0:

  H(s) = M_0 + M_1 (s - s_0) + M_2 (s - s_0)^2 + ...

where M_k are the moments `M_k = (-1)^k C (s_0 I - A)^{-k-1} B`.

A reduced model H_r(s) **matches the first r moments** of H(s) if
the reduced basis V is constructed from the Krylov subspace

  K_r(M_inv, M_inv b) = span{M_inv b, M_inv^2 b, ..., M_inv^r b}

where M_inv = (s_0 I - A)^{-1}.

### Algorithm

**Arnoldi** for general (non-symmetric) A:
- Iteratively build orthonormal basis V via modified Gram-Schmidt
- Produces an upper-Hessenberg reduced operator H_r
- Cost: r matrix-inverse applications + O(N r^2) orthogonalization

**Lanczos** for symmetric A:
- Three-term recurrence (more efficient than Arnoldi)
- Produces a tridiagonal reduced operator T_r
- This is the algorithm used in **CLN** (Sugahara/Wakao)

**Two-sided Lanczos** for asymmetric A with both inputs and outputs:
- Builds V (for state) and W (for output), with W^T V = I_r
- Produces tridiagonal T_r
- Used in **PRIMA** (Padé via Lanczos, Odabasioglu et al. 1998)

### Multipoint expansion (rational Krylov)

Single expansion point matches moments well near s_0 but poorly
elsewhere.  Pick MULTIPLE expansion points {s_0, s_1, ..., s_K},
build Krylov subspaces at each, concatenate the bases.

Sweet spot: expansion points spread across the frequency range of
interest.  IRKA (Iterative Rational Krylov Algorithm) iteratively
optimizes the expansion points for H_2 optimality.

### Why moment-matching is fast

- No SVD of dense snapshot matrix (unlike POD)
- No greedy parameter selection (unlike RB)
- Cost dominated by r linear solves (with shared factorization if
  the matrix is fixed)

For LTI systems with O(10^6) DOFs, moment-matching reaches r ~ 30
with ~30 PARDISO back-substitutions — minutes, not hours.

### Connection to CLN (Sugahara)

CLN is essentially **Lanczos applied to the PEEC LR-Hamiltonian-coupled
system** with specific physical interpretation of the reduced
tridiagonal matrix as a **Cauer ladder circuit** (R-L cells in
series + parallel).  This makes the reduced model directly usable
in SPICE-style circuit simulators.

See `radia_mcp.mor.cln_knowledge` for the lab's CLN implementation.
"""

PGD = """\
## Proper Generalized Decomposition (Vol 2 Ch 3)

### The idea — separable expansion

Approximate a multi-variable function (space x parameter) by a
**sum of rank-1 separable products**:

  u(x, mu) ≈ Σ_{i=1}^r  F_i(x) · G_i(mu)

This is structurally similar to a POD modal expansion, but:
- POD: G_i is fixed at "snapshot value at mu_i"
- PGD: G_i is itself unknown, computed by **greedy alternating
  minimization** in (F, G) space

### Algorithm sketch

Greedy: for k = 1, 2, ..., r:
  Initialize F_k^0, G_k^0 arbitrarily.
  Repeat:
    Fix G_k.  Solve for F_k by reducing the residual (1D-in-x problem).
    Fix F_k.  Solve for G_k by reducing the residual (1D-in-mu problem).
  Until convergence.
  Subtract rank-1 term from residual.

### Why this matters

PGD is **a priori** — no offline snapshot generation needed.  The
reduced model is built incrementally as part of the solution
process.

For high-dimensional parameter spaces (UQ with 20+ parameters), PGD
is **the** practical choice (POD/RB suffer the curse of
dimensionality).

### EM/motor applications

- Parametric motor design with 5-20 geometric parameters (slot
  width, tooth height, magnet thickness, ...) — PGD handles this
  natively
- Stochastic ν(B) (manufacturing-tolerance UQ on iron grade)
- Multi-physics coupled (E&M + thermal + mechanical) — separable in
  the physics dimension

Sugahara lab status: ✗ not implemented; potential future direction.
"""

HYPERREDUCTION = """\
## Hyper-reduction: DEIM, EIM, gappy POD (Vol 2 Ch 5)

### Why hyper-reduction is needed

Linear projection-based MOR gives O(r^2) reduced solve cost, BUT
the **nonlinear** term `V^T f(V u_r)` still costs O(N) per
evaluation because:
1. Reconstruct u = V u_r (cost O(Nr))
2. Evaluate f(u) at all N degrees of freedom (cost O(N))
3. Project f back: V^T f (cost O(Nr))

For nonlinear problems (ν(B), hysteresis), this kills the
performance gain entirely.

### DEIM (Discrete Empirical Interpolation Method)

Chaturantabut-Sorensen 2010.  Pick a small set of **interpolation
indices** {p_1, ..., p_m} (where m is small) and approximate

  f(u)  ≈  U_DEIM (P_DEIM^T U_DEIM)^{-1} P_DEIM^T f(u)

where U_DEIM is a POD basis of f-snapshots and P_DEIM is the
selector matrix for the interpolation indices.

The key savings: `P_DEIM^T f(u)` only evaluates f at the m chosen
indices, NOT all N.  Combined with `V_DEIM = V^T U_DEIM` precomputed,
the entire nonlinear evaluation costs O(m^2) per step instead of O(N).

### Choosing interpolation indices

Greedy: pick p_1 to maximize |U_DEIM[:, 1]|; pick p_2 to maximize
the residual after fitting at p_1; and so on.  Yields m points with
provable a priori interpolation error bounds.

### EM/motor application

The HOLLAUS 2023 MSFEM+MOR+DEIM paper (IEEE TMAG 2023) is exactly this:
- MSFEM reduces 3D-laminated to 2D-1D (geometry hyper-reduction)
- POD on MSFEM snapshots reduces to r ~ 30
- DEIM hyper-reduces the nonlinear ν(B) evaluation
- Result: real-time-feasible motor simulation with nonlinear iron

This is the **state of the art** as of 2024 for laminated machine
transient analysis.

### Variants

- EIM (Empirical Interpolation Method, Maday-Nguyen-Patera): the
  continuous analog of DEIM, used in Reduced Basis context.
- Gappy POD (Everson-Sirovich 1995): the original idea, less
  rigorous than DEIM but simpler.
"""

SYSTEM_THEORETIC_BT = """\
## Balanced Truncation + H_inf / H_2 optimal (Vol 1 Ch 2)

### Balanced Truncation (BT, Moore 1981)

For an LTI state-space (A, B, C, D), compute:
- **Controllability Gramian** P satisfying A P + P A^T + B B^T = 0
- **Observability Gramian** Q satisfying A^T Q + Q A + C^T C = 0

In the balanced realization, P = Q = Σ = diag(σ_1, ..., σ_n) with
σ_1 ≥ σ_2 ≥ ... ≥ σ_n ≥ 0 = **Hankel singular values**.

Truncate to first r dimensions → reduced (A_r, B_r, C_r, D).

**A priori error bound**:

  ||H - H_r||_{H_inf}  ≤  2 (σ_{r+1} + σ_{r+2} + ... + σ_n)

This is a STRICT a priori bound, available before doing the
reduction.  Unique to BT in the MOR universe.

### Computational cost

Standard BT requires SOLVING the two Lyapunov equations, each costing
O(N^3) in dense arithmetic.  For N > 10000 this is prohibitive.

Modern alternatives use **low-rank Lyapunov solvers** (ADI, RKSM)
that exploit the rank structure of B B^T (B is N × p_inputs with
p_inputs << N).  Cost reduces to O(N · r^2).

### H_2-optimal MOR (IRKA)

Iterative Rational Krylov Algorithm (Gugercin-Antoulas-Beattie 2008):
- Pick initial expansion points
- Build Krylov bases
- Compute reduced model
- Move expansion points to mirror images of reduced eigenvalues
- Iterate to convergence

Yields a reduced model that is H_2-optimal: minimizes
`integral_omega ||H(jomega) - H_r(jomega)||^2 domega`.

Practical advantages over BT:
- No Lyapunov solve required
- Cost: O(r) linear solves with the full operator
- Works for problems where Lyapunov is infeasible

### When system-theoretic MOR fits EM/motor

EM/motor problems are usually NOT presented as LTI state-space —
they're parametric PDEs.  But after spatial discretization + linear
coupling (no nonlinear iron), the resulting frequency-domain solver
IS LTI in the input-output sense.

Use BT/IRKA for:
- Linear eddy current modes (resonant winding coupling)
- Linear inductor characterization with multiple ports
- Cable + motor system EMI (LTI is the right framework)

For nonlinear iron, fall back to Family 1 (POD with hyper-reduction).
"""

DATA_DRIVEN_DMD_OI = """\
## DMD, Operator Inference, Loewner (Vol 1 Ch 6, Vol 2 Ch 7)

### Dynamic Mode Decomposition (DMD, Schmid 2010)

Given two snapshot matrices X = [u_0, u_1, ..., u_{M-1}] and
Y = [u_1, u_2, ..., u_M] (consecutive in time), DMD computes the
best-fit **linear operator** A_DMD such that

  Y ≈ A_DMD X

via least squares: A_DMD = Y X^+ (pseudoinverse).  Then eigen-
decomposition of A_DMD gives **DMD modes + frequencies**.

For nonlinear systems, DMD doesn't capture nonlinearity, but its
modes still reveal dominant **coherent structures** in the data
(eddies in fluid flow, current ripple modes in motor PWM, ...).

### Operator Inference (Peherstorfer-Willcox 2016)

Generalize DMD to **polynomial operators**: fit

  Y ≈ A u + B u^2 + C u^3 + ... + D control

by least squares.  Adds nonlinearity in a data-driven way without
needing the underlying PDE.

Particular advantage: **physics-informed regularization** can be
added (energy conservation, dissipation inequality, ...) to make the
fitted operator obey known constraints.

### Loewner framework (Mayo-Antoulas 2007)

Given transfer-function samples `{H(jomega_k)}_{k=1}^K` from a
black-box simulator, construct a rational interpolant of minimum
degree.  Uses the **Loewner matrix** L_ij = (H(omega_i) -
H(omega_j)) / (omega_i - omega_j).

The SVD of L gives the McMillan degree (minimum required model
order).  No equations needed — pure frequency-response data.

### EM/motor applications

- **DMD** on PWM motor transient: extracts dominant time-harmonic
  modes of current/torque ripple — useful for controller design.
- **Operator Inference** on temperature-dependent ν(B,T) coupling:
  fit a polynomial T-dependent operator from snapshot data.
- **Loewner** for vendor-supplied transformer / motor frequency
  responses: build an LTI model from S-parameter measurements.

For radia/NGSolve specifically, DMD is the most directly useful —
the snapshot data from a transient calc_motor_transient.py run
naturally feeds into a DMD post-processor.
"""

PARAMETRIC_PMOR = """\
## Parametric MOR (Vol 2 Ch 1)

### The parametric extension

The standard RB / POD framework handles one fixed PDE.  Parametric
MOR (PMOR) handles a family of PDEs

  L(mu) u(mu) = b(mu)

where mu varies in a parameter space P ⊆ R^d.  The goal: a reduced
basis V valid across ALL of P.

### Affine parameter dependence (the lucky case)

If A(mu) = Σ_q theta_q(mu) A_q (sum of parameter-functions times
mu-independent matrices), then the reduced system

  A_r(mu) = Σ_q theta_q(mu) (V^T A_q V)

is assembled cheaply at each new mu (just sum r×r matrices).  Online
cost: O(r^3) PER mu, independent of N.

### Non-affine (the common case)

If A(mu) is non-affine (e.g. nonlinear material), use EIM (Empirical
Interpolation) to approximate it as affine plus error:

  A(mu) ≈ Σ_q theta_q(mu) A_q  + O(eps_EIM)

EIM picks the A_q from a greedy snapshot procedure.  This is the
"hyperreduction" sibling of DEIM (cf. `hyperreduction` topic).

### Greedy parameter sampling

To build V efficiently, use a greedy loop:
1. Start with one parameter mu_1; solve full → u_1; V = [u_1].
2. For each candidate mu in a training set:
   - Cheaply estimate error eta(mu) from residual + stability constant
3. Pick mu_k+1 = argmax eta; solve full; orthogonalize → enlarge V
4. Stop when max eta < tolerance

### EM/motor use

PMOR is the natural framework for design exploration:
- Vary slot dimensions, magnet thickness, iron grade
- Each design takes O(r^3) instead of O(N^3)
- Train basis on ~50-200 design candidates, then evaluate thousands

The reduced model can then drive an outer optimization loop
(genetic algorithm, gradient descent) at real-time speed.
"""

CLN_SUGAHARA = """\
## CLN (Cauer Ladder Network) — Sugahara/Wakao lab specialty

Cross-reference to `radia_mcp.mor.cln_knowledge` which has the
full lab-specific treatment.

### Quick summary in MOR taxonomy

CLN sits in **Family 1 (projection-based)** as a SPECIALIZED
Lanczos variant:

- Two-sided Lanczos applied to the PEEC + magnetic-coupling system
- Reduced model is **structurally a tridiagonal R-L state-space**
  with physical R, L coefficients
- Tridiagonality is the structural fingerprint of Lanczos; the
  physical interpretation as **R-L ladder** is the lab insight

### Why CLN is special

Generic Lanczos gives an abstract reduced (A_r, B_r, C_r) with no
particular structure beyond tridiagonality.  CLN's contribution is
showing that for **PEEC-discretized magnetic-conductor systems**:

1. The Lanczos reduction PRESERVES energy / passivity (provable)
2. The reduced state variables have INTERPRETABLE units (currents
   through ladder branches)
3. The reduced system is directly EXPRESSIBLE as a SPICE netlist
   without further translation

This makes CLN the natural pre-processing step for **circuit-
simulator-coupled** motor / power-electronics analysis.

### Sugahara contributions (multiple papers)

The Sugahara lab has co-authored at least 5+ CLN papers:
- Original CLN for PEEC matrices
- Extension to magnetic-material coupled systems
- Combined CLN+SIBC for IH workpieces
- Verified-arithmetic CLN (DD precision, double-double)
- 3D CLN extensions

Cross-reference `radia_mcp.mor.cln_knowledge` for detailed
formulations and CLAUDE.md "PRIMA Model Order Reduction" for the
production implementation in `radia/lanczos_reduction.py`.

### Position in the MOR taxonomy (CLN as a research-niche)

| Aspect | Generic Lanczos | CLN (Sugahara) |
|--------|------------------|----------------|
| Reduced structure | Tridiagonal | Tridiagonal |
| Physical interpretation | None | R-L ladder |
| Passivity guarantee | No (without effort) | Yes (built in) |
| SPICE export | Manual | Automatic |
| Verification arithmetic | Standard FP | DD/Interval supported |
| Lab application | None | IH, transformer, motor coil |

CLN is **not a competitor** to POD/RB for general MOR — it's a
specialized tool for the PEEC-discretized circuit-coupled regime.
"""

EM_SPECIFIC_IOAN = """\
## EM-specific MOR — Ioan et al. (Vol 3 Ch 5, 56 pages)

Reference: D. Ioan et al., "Complexity reduction of electromagnetic
systems", in: Model Order Reduction Vol 3 Applications, deGruyter,
2021, pp. 145-209 (56 pages).

This is the canonical EM-specific MOR reference.  Author Daniel
Ioan is a long-time computational EM person at Politehnica
Bucharest.

### Structure

1. **EM modeling pipeline** (steps 1-7: conceptual → AAM → numerical
   → simulation → MOR → V&V)
2. **EM fundamental quantities** + Tonti diagram / De Rham
   sequences / Maxwell's House
3. **Discretization methods** (FEM, FIT, BEM, ...) and what each
   does to the MOR space
4. **EM-specific MOR challenges**: passivity, structure preservation,
   gauge invariance, low-frequency stability
5. **Survey of EM MOR methods** with application examples
6. **EM-specific MOR techniques**: Krylov for FE, AWE/PRIMA for PEEC,
   POD for parametric design

### Why Tonti diagram matters for MOR

The Tonti diagram (also discussed in
`radia_mcp.differential_forms.maxwell_stress`) classifies EM
quantities into a 2x2 grid of (configuration, source) x (primal,
dual) cells linked by exact differential operators.  This structure
SHOULD be preserved by any good MOR scheme:

- Reduced inputs and outputs should sit in the right cells
- Reduced operators should preserve the kernel/range relationships
- Passivity (energy positivity) must survive reduction

Generic MOR (POD applied blindly) does NOT preserve this structure
and can produce reduced models that violate physical conservation
laws.  Ioan emphasizes **structure-preserving MOR** specifically
for EM.

### Key EM-MOR results from this chapter

- **Frequency-domain Maxwell + PEEC** → directly amenable to PRIMA
  (Padé-via-Lanczos): preserves passivity automatically
- **Time-domain transient FEM** → POD with hyperreduction for the
  nonlinear ν(B) coupling
- **Parametric magnetostatic** → RB with stability constants
- **Coupled multi-physics** → PGD for separability across physics

### Lab implication

For radia + NGSolve MOR work, follow Ioan's structure-preservation
guidance:
- For PEEC reduction → use the Sugahara CLN (which IS PRIMA + ladder
  realization)
- For FE eddy-current reduction → use POD + DEIM (Hollaus 2023)
- For parametric design → use RB / PGD

This module + `cln_sugahara` + `hyperreduction` topics give the
combined picture.
"""

ROTATING_MACHINES_KISS_OROSZ_2024 = """\
## Kiss-Orosz 2024 — MOR review for rotating machines (Energies)

Reference: K.L. Kiss, T. Orosz, "Model Order Reduction Methods for
Rotating Electrical Machines: A Review", Energies 17:5145, 2024.
DOI: 10.3390/en17205145.
Open access (Creative Commons CC BY).

### Why a dedicated review now

E-mobility is driving demand for **real-time control** of high-
performance machines.  Real-time control requires REDUCED motor
models (full FE is too slow at 10-100 kHz inverter rates).
Industrial state-of-the-art (lookup tables) is limited; ML / MOR are
the alternatives.

### Methods surveyed

The paper systematically reviews:

| Method | Mentioned in keywords | Coverage |
|--------|------------------------|----------|
| POD | ✓ | extensive |
| PGD | ✓ | extensive |
| OIM (Orthogonal Interpolation Method) | ✓ | moderate |
| **CLN (Cauer Ladder Network)** | ✓ | dedicated section |
| Krylov / moment-matching | (implied) | brief |
| RB | (implied) | brief |
| Balanced truncation | not in keywords | minimal |
| DMD / OI | not in keywords | minimal |

CLN being in the abstract keywords is significant: it confirms CLN
has reached canonical-method status in the rotating-machine MOR
literature.

### Practical conclusions (paper's perspective)

1. **No published industrial real-time MOR implementation** (as of
   2024).  Gap between academic methods and production deployment.
2. **POD + DEIM is mature** for offline parametric design.
3. **CLN is the right tool** for control-loop-coupled simulation
   where the reduced model must talk to a circuit simulator.
4. **PGD is under-utilized** despite theoretical advantages — has
   not crossed the academic-industrial gap.

### Lab implication

The Sugahara lab is **already practicing what Kiss-Orosz 2024
recommend** for the control-loop-coupled case (CLN).  The
gap-of-2024 (no published industrial implementation) is an
**opportunity**:

- Take CLN out of the lab and into a production motor controller
- Pair with Lange-Henrotte-Hameyer transient framework
  (`calc_motor_transient.py`) for the offline FE side
- Result: end-to-end open-source motor control simulation toolchain

This would be a high-impact paper / open-source release direction.
"""

SOFTWARE_LAB = """\
## MOR software landscape (Vol 3 Ch 13)

### Open-source MOR toolboxes

| Software | Language | Strengths |
|----------|----------|-----------|
| pyMOR | Python | RB, POD, EIM; well-documented; active community |
| MORLab | MATLAB | TBR, Krylov, multiple input/output; mature |
| MOR Wiki | Online | Benchmark problems repository |
| RBniCS | Python (FEniCS) | RB integrated with FEniCS FEM |
| EMORE | C++ | EM-specific; legacy at TU Berlin |

### Lab-specific tooling

The Sugahara lab has in-house CLN implementation at
`radia/lanczos_reduction.py` (CLAUDE.md "PRIMA Model Order
Reduction" section).  Key classes:

```
SPICEExtractionConfig     - configuration object
PRIMASchurExtractor       - PRIMA + Schur complement for coupling
LoopStarMagneticCoupled  - magnetic-magnetic-material coupling layer
```

### Recommended workflow

For new MOR work in radia + NGSolve:

1. **Prototype** in pyMOR (Python, fast iteration)
2. **Production** in `radia.lanczos_reduction` if CLN-specific, else
   write a custom NGSolve + numpy implementation
3. **Benchmark** against MOR Wiki test cases when applicable

### Integration points with radia / NGSolve

NGSolve `BilinearForm.mat` is a sparse matrix that pyMOR can
consume directly.  The pattern:

```python
# In a calc_*.py script:
from ngsolve import BilinearForm, ...
from pymor.algorithms.pod import pod
from pymor.bindings.ngsolve import NGSolveVectorSpace, NGSolveVisualizer

# Solve at training parameters
snapshots = NGSolveVectorSpace(fes).empty()
for mu in training_set:
    gfu = solve_full(mu)
    snapshots.append(gfu.vec.FV().NumPy().copy())

# Build POD basis
pod_basis, sigmas = pod(snapshots, atol=1e-6)

# Use pod_basis as V for projection
```
"""

LAB_RECOMMENDATION = """\
## Sugahara lab MOR recommendation

### Decision flowchart for a new MOR project

```
Is the system LINEAR (no nu(B))?
  YES → Is it PEEC-discretized (coil-only or coil-conductor)?
    YES → CLN (radia.lanczos_reduction)
    NO  → POD or PRIMA depending on input/output count
  NO  → Is it MOR for a CONTROL loop (real-time)?
    YES → Combine FE + CLN with a lookup table for nonlinear iron
    NO  → POD + DEIM (offline parametric study)
```

### Method-by-application matrix

| Application | Recommended MOR | Notes |
|-------------|------------------|-------|
| IH coil + workpiece (linear) | CLN | Sugahara specialty |
| IH coil + workpiece (nonlinear via ESIM) | CLN + ESIM lookup | Frequency-domain |
| Transformer characterization | CLN or RB | depends on # ports |
| Motor parametric design | POD + DEIM | Hollaus 2023 reference |
| Motor real-time control | FE → CLN → SPICE | Kiss-Orosz 2024 gap-closer |
| Accelerator magnet field map | Direct evaluation | no MOR needed (small N) |
| EMC / EMI cable + motor | CLN + BEM (ngsolve.bem) | Helmholtz at high f |

### What is NOT in scope

These would be over-engineering for the lab:
- Generic black-box LTI MOR (use system identification tools instead)
- Quantum / many-body MOR (different mathematical universe)
- Network MOR (network/graph MOR Vol 3 Ch 11 — different from EM)

### What IS in scope but unimplemented

Worth considering for future investment:
- PGD for parametric motor design (5-20 geometric params, currently
  done via RB which struggles past 5)
- DMD for transient post-processing (extracts dominant modes from
  simulation history, useful diagnostic)
- Loewner from measurement (combine measured transformer S-params
  with FE-derived MOR for vendor-supplied components)
"""

SECTIONS = {
    "mor_taxonomy": MOR_TAXONOMY,
    "projection_pod_rb": PROJECTION_POD_RB,
    "projection_krylov": PROJECTION_KRYLOV,
    "pgd": PGD,
    "hyperreduction": HYPERREDUCTION,
    "system_theoretic_bt": SYSTEM_THEORETIC_BT,
    "data_driven_dmd_oi": DATA_DRIVEN_DMD_OI,
    "parametric_pmor": PARAMETRIC_PMOR,
    "cln_sugahara": CLN_SUGAHARA,
    "em_specific_ioan": EM_SPECIFIC_IOAN,
    "rotating_machines_kiss_orosz_2024": ROTATING_MACHINES_KISS_OROSZ_2024,
    "software_lab": SOFTWARE_LAB,
    "lab_recommendation": LAB_RECOMMENDATION,
}


def get_systematic_mor_knowledge(topic: str = "mor_taxonomy") -> str:
    """Return systematic MOR knowledge.

    Args:
        topic: section name; "all" returns everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
