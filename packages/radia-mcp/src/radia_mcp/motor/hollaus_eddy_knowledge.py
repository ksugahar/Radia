"""Karl Hollaus / TU Wien MSFEM knowledge for radia_mcp.motor.

**The answer to "Lange-Henrotte-Hameyer doesn't do eddy currents"**.

Karl Hollaus (TU Wien, Institute of Analysis and Scientific Computing,
with Joachim Schöberl of NGSolve fame, plus Markus Schöbinger and
Valentin Hanser) developed the **Multiscale FE Method (MSFEM)** for
the precise problem the Lange-Henrotte-Hameyer field-circuit
linearization sidesteps: **eddy currents in laminated iron cores**
without meshing every individual lamination sheet.

The lab library has 10 direct Hollaus papers + multiscale references.
Key insight: a laminated iron core has dimensional separation
(N_sheets × thickness ≪ overall dimension), which is the textbook
setup for **homogenization by scale separation**.  MSFEM exploits
exactly this.

Coverage:
  intro            — Why standard FE fails on laminated cores; MSFEM
                     idea in one paragraph
  msfem_2d1d       — 2D-1D MSFEM: 2D macroscopic + 1D microscopic
                     across-the-thickness
  t_phi_open       — T-Φ-Φ MSFEM for OPEN magnetic circuits
                     (rotating machines, not just transformer yokes)
  effective_material — Effective Material (EM) homogenization (Hanser-
                     Schöbinger-Hollaus 2024)
  error_estimator  — Equilibrated error estimator for 2D-1D MSFEM
                     T-formulation (Schöbinger-Hollaus 2024)
  msfem_3d         — Full 3D MSFEM for laminated cores
  nonlinear_esi    — Nonlinear effective surface impedance in scalar
                     potential formulation (Hollaus-Hanser-Schöbinger)
  msfem_plus_mor   — MSFEM + Model Order Reduction (snapshots from
                     MSFEM → reduced basis)
  motor_strategy   — How to combine MSFEM with Lange-HH for motor
                     analysis that ALSO handles eddy currents
"""

from __future__ import annotations


INTRO = """\
## Why standard FE fails on laminated iron cores

A modern motor / transformer stator is a stack of **N_sheets** thin
iron laminations (typically 0.35-0.50 mm each) separated by thin
insulation layers (a few μm of oxide or organic resin).  The stack
might be 100-500 mm long axially → 200-1500 sheets per machine.

To resolve the eddy currents within each sheet by standard FE, you
would need at least 3-5 elements **through the thickness** of each
sheet.  With 1000 sheets → 3000-5000 elements per axial slice,
multiplied by the cross-sectional resolution (~10⁴ tetrahedra) →
**3×10⁷ to 5×10⁷ DOFs**, with nonlinear ν(B) on top.  Solving this at
each PWM time step is infeasible.

### The Hollaus / Schöberl / Schöbinger MSFEM idea

A laminated iron core is a textbook **two-scale periodic** problem:

| Scale | Length | Physics |
|-------|--------|---------|
| Macro | machine dimension (mm-m) | Global B distribution |
| Micro | sheet thickness (sub-mm) | Eddy current circulation inside one sheet |

The **Multiscale FE Method** (Hollaus & Schöberl, building on Hou-Wu
1997 multiscale theory) factors the solution as

  A(x, y, z, t)  ≈  A_macro(x, y, z, t)  +  Σ_k φ_k(z/ℓ) · A_micro,k(x, y, t)

where ℓ is the local sheet thickness and `φ_k(ζ)` are **polynomial
shape functions in the thickness coordinate** ζ = z/ℓ.  The macro
field A_macro varies on the machine scale; the micro field A_micro,k
varies on the cross-section scale; the through-the-thickness
dependence is captured analytically by `φ_k`.

Discrete cost: O(N_macro · N_micro · k_modes) instead of O(N_macro ·
N_sheets · N_thickness).  For typical motor parameters, this is a
**100× to 1000× DOF reduction** at comparable accuracy.

### Status in radia / radia-mcp

- `radia_mcp.ih` already knows about a related effective-surface-
  impedance scalar-potential formulation (Hollaus-Hanser-Schöbinger).
- NGSolve has the underlying multiscale primitives (NumberSpace,
  composite FESpace) needed for MSFEM.
- The radia-side implementation is **not yet written** — this
  knowledge module is the prerequisite.
"""

MSFEM_2D1D = """\
## 2D-1D MSFEM (Schöbinger-Hollaus, IEEE TMAG 60(5), May 2024)

Reference: M. Schöbinger, K. Hollaus, "An Equilibrated Error
Estimator for the 2-D/1-D MSFEM T-Formulation of the Eddy Current
Problem", IEEE Trans. Mag. 60(5):7401906, 2024.

### The key dimensional reduction

For rotating machines, *every sheet in the active zone sees the same
field* (by axial symmetry under the rotor).  So the 3D problem
reduces to:

  - **A single 2D cross-section** of one representative sheet
  - **A 1D polynomial expansion along the thickness** to capture
    eddy currents

The 3D vector potential ansatz is

  A(x, y, z)  =  Σ_k p_k(z/d) · A_k(x, y)

where d is the sheet thickness, `p_k(ζ)` is a Legendre-like polynomial
of order k, and `A_k(x, y)` is a 2D field.  For first-order accuracy,
k ∈ {0, 1, 2} suffices (3 macro fields per node).

### T-formulation

Hollaus uses the **electric current density T** as the primary
variable (curl T = J, divergence-free).  In 2D-1D:

  T(x, y, z, t)  =  T_macro(x, y, t)  +  ξ(z) · T_micro(x, y, t)

with ξ(z) odd through the sheet (zero at sheet centerline by
symmetry).

The resulting weak form is solved on a 2D mesh of the cross section,
with the 1D polynomial coefficients carried per macro node.  Effective
problem size: ~3× a single-sheet 2D problem (one macro + a few
polynomial coefficients).

### Equilibrated error estimator (paper contribution)

Standard residual-based estimators (à la Babuška-Rheinboldt) give
upper bounds via constants that depend on shape regularity.
Equilibrated estimators give **guaranteed bounds without unknown
constants**, based on the Prager-Synge theorem:

  ||u - u_h||²_E  ≤  ||σ_eq - flux(u_h)||²_E

where σ_eq is a *reconstructed* equilibrated flux that satisfies
exact local conservation.  The 2024 paper extends this to 2D-1D
MSFEM with the T-formulation — gives both a good total-error proxy
AND a per-element local error for adaptive refinement.

### Practical implication for motor analysis

A 2D-1D MSFEM solve replaces a 3D solve of a laminated machine
sector.  Combined with periodicity + slot symmetry, a full 8-pole
machine becomes a **single-sheet, single-slot, 2D problem with three
polynomial layers** — solvable in seconds, not hours.

For PWM transient analysis (the Lange-HH use case), this is the
right scale-down: each FE refresh in Lange-HH becomes a small 2D-1D
MSFEM solve, not a full 3D solve.
"""

T_PHI_OPEN = """\
## T-Φ-Φ MSFEM for open magnetic circuits (Hanser-Schöbinger-Hollaus 2024)

Reference: V. Hanser, M. Schöbinger, K. Hollaus, "A T, Φ-Φ Multiscale
Finite Element Formulation for Eddy Current Problems in Open Magnetic
Circuits", 4p IEEE submission, TU Wien 2024.

### What "open" means

A *closed* magnetic circuit (transformer yoke) has all flux lines
fully contained inside the iron — every flux line that leaves a
lamination eventually comes back through another.  An *open* circuit
(rotating machine air gap, end-region, stator end-winding) has flux
lines that **exit the iron entirely** — into the air gap (for motor
torque production) and the stray-flux region (for losses and EMI).

Pure MSFEM was developed for closed circuits where the
through-the-thickness microscale problem is well-defined globally.
For open circuits the **exit boundary conditions** between iron and
air must be handled carefully.

### The T-Φ-Φ formulation

| Region | Variable | Why |
|--------|----------|-----|
| Iron lamination | **T** (electric current) | T is the natural variable for eddy currents (curl T = J) |
| Air / stray-flux region | **Φ** (magnetic scalar) | Air is non-conducting, scalar potential suffices |
| Cuts in iron-air interface | **Φ** (multi-valued / cuts) | Cohomology — same as Dular-Henrotte 1997 |

The "double Φ" reflects that the air-region scalar potential needs
**cuts** to handle the multiply-connected geometry (a rotating
machine has loops everywhere: stator slot openings, rotor cage
bars).  The cut treatment is exactly the construction from
Dular-Henrotte 1997 (see `henrotte_lineage_knowledge.source_field_1997`).

### Numerical evidence

The paper validates on a current-driven coil + 3D laminated iron
core in open-air bounding box.  Compared to standard FE with each
sheet resolved:

- MSFEM eddy current losses: < 1% error
- MSFEM reactive power: < 1% error
- Local field distribution: visually identical
- DOF count: ~100× smaller

### Application path for motor analysis

For a motor stator with **N_slots** slots:

1. Mesh ONE stator slot region with ONE laminated tooth + air gap
   segment.  Apply periodicity to obtain a 2D-1D MSFEM problem.
2. Mesh the **air gap** with sliding-band BC (Hanser et al. tested
   this; works fine).
3. Mesh the rotor as a separate domain with material rotation.
4. Couple via the T-Φ-Φ formulation: iron uses T_micro polynomials,
   air uses scalar Φ with cuts.
5. Solve with the multiscale ansatz.

Result: a fully 3D motor with full lamination eddy-current physics,
at the cost of a 2D-1D + air-domain problem.  This is the **eddy-
current-capable** alternative to the Lange-Henrotte-Hameyer
linearization.
"""

EFFECTIVE_MATERIAL = """\
## Effective Material homogenization (Hanser-Schöbinger-Hollaus, IEEE TMAG 60(10), Oct 2024)

Reference: V. Hanser, M. Schöbinger, K. Hollaus, "Effective Material
Modeling for Laminated Iron Cores With a T, ε-ε Formulation", IEEE
Trans. Mag. 60(10):7402307, October 2024.

### The Effective Material (EM) idea

Instead of solving MSFEM (with explicit macro+micro decomposition),
**run an offline cell problem ONCE** to compute an "effective
material" that captures the eddy-current behavior of the laminated
stack, then run the macroscopic simulation with the homogenized
material **without any periodic structure**.

### The cell problem

Take a representative section of the lamination stack (one period:
iron + insulation + iron + insulation + ...).  Drive it with each
fundamental excitation (uniform B in each direction at frequency f).
Measure:

- **Eddy current losses ECL** per unit volume at amplitude B
- **Reactive power RP** (imaginary part of complex permeability)

This defines `ECL(B, f)` and `RP(B, f)` for the homogenized material.
With nonlinear iron these are also functions of |B|.

### Application

Replace the laminated core in the global FE model by a single
**homogeneous** block with reluctivity `ν_eff = 1/μ_eff(B, f)` derived
from `RP(B, f)`.  Add a per-element power-dissipation post-processing
that evaluates `ECL(B, f)` on the macro solution.

Trade-off vs MSFEM:

| Aspect | MSFEM | Effective Material |
|--------|-------|---------------------|
| Per-FE-call cost | O(macro·micro·k) | O(macro) only |
| Implementation complexity | Multi-FESpace, careful coupling | Just a different material CF |
| Captures local field per-sheet? | Yes | No (averaged) |
| Accurate for global losses? | Yes | Yes (paper shows < 1% error) |
| Accurate for local hot-spots? | Yes | No (needs MSFEM) |

For motor design (lifetime loss, efficiency map), Effective Material
is the **practical workhorse** — much cheaper than MSFEM and
indistinguishable for total losses.  For research / detailed thermal
hotspot prediction, use MSFEM.

### 2025 update — A-formulation + circuit coupling (Hanser-Schöbinger-Hollaus)

Reference: V. Hanser, M. Schöbinger, K. Hollaus, "Effective material
modelling for laminated iron cores with an A-formulation and circuit
coupling", COMPEL 44(5):832, 2025.  CC BY 4.0 open access.

The 2025 paper extends EM to:
- A-formulation (instead of T,ε-ε from the 2024 paper)
- External circuit coupling (impose terminal voltage / current via
  per-coil Schur complement)

Reported performance: ~1% error vs full FE, **300× speedup on small
cores, 5000× speedup on large cores**.  This is the most directly
applicable variant for Sugahara-lab motor work because:
- A-formulation matches `calc_motor_transient.py`'s native variable
- Circuit coupling matches the Lange-Henrotte-Hameyer framework

When `calc_motor_lamination.py` v2 is extended to circuit coupling,
this 2025 paper is the canonical reference.

**v2.5 status (2026-05-21)**: `calc_motor_lamination.py` now supports
three drive modes via `--drive {meanB,current,voltage}`:

- `meanB` (default): canonical Hollaus 2024 driving (Dirichlet on
  `A_top - A_bottom = B_avg * d_total`).
- `current`: Hanser 2025 J_s-driven extension.  Apply uniform
  current density in iron with `--J-s-iron <A/m^2>`.  Useful for
  measuring iron's response to a coil-driven source rather than an
  imposed flux.
- `voltage`: planned -- requires Lagrange multiplier coupling for
  terminal-voltage drive.  Currently raises `NotImplementedError`
  with guidance to use `meanB` or `current`.

### Practical recipe

```python
# Offline (once per material grade):
def compute_effective_material(iron_grade="M19_29G", f_range=(50, 50000)):
    cell = build_one_period_cell()  # 1 sheet + insulation
    sigma, mu_BH = iron_grade_table(iron_grade)
    ECL_table = {}; RP_table = {}
    for f in f_range:
        for B_pk in B_amplitude_grid:
            ECL_table[(f, B_pk)] = solve_cell_eddy_loss(...)
            RP_table[(f, B_pk)]  = solve_cell_reactive_power(...)
    return interp2d(ECL_table), interp2d(RP_table)

# Online (in main FE solve):
nu_eff_cf = mesh.MaterialCF({"stator_iron": NU_0 * RP_interp(|B|, f)})
ECL_density_cf = mesh.MaterialCF({"stator_iron": ECL_interp(|B|, f)})
total_ECL = Integrate(ECL_density_cf, mesh, definedon=mesh.Materials("stator_iron"))
```

This is the most pragmatic path to add lamination losses to the
Lange-Henrotte-Hameyer transient framework.  Add the
`(jω σ_eff_imag)` term to `L_inc` and let the per-step macro solve
do the rest.

### Cross-validation with `radia.lamination.laminated_mu_eff`

radia-ngsolve ships the analytic complex effective permeability of the LINEAR
in-plane lamination as `radia.lamination.laminated_mu_eff`:

    mu_eff = mu0 [ fill*mu_r*tanh(b)/b + (1-fill) ],   b = (d/2) sqrt(j w mu0 mu_r sigma)

Use it as the analytic ground truth for `solve_cell_problem`'s effective
permeability in the linear regime.  An independent 1D FE (symmetric cell, iron
centered between insulation halves) reproduces it to < 2.2 % from |b| = 0.5 to
|b| ~ 10 (validation_test/radia_mcp/test_laminated_mu_eff_freqsweep.py).

PITFALL (verified the hard way during 横展開): the eddy LOSS is NOT
`0.5*w*B^2*Im(1/mu_eff)`.  `mu_eff = <B>/H_s` is the FLUX-AVERAGING permeability
(correct for the macro FE's B-H relation), but the lamination eddy loss is
RESISTIVE -- it must come from the direct Joule integral `0.5*sigma*|E|^2` (=
the cell's `ECL_density`, validated by the Bertotti limit) or the surface
impedance, NOT from `Im(mu_eff)`.  The two coincide only at low frequency; at
|b| ~ 5 the `Im(mu_eff)` route under-predicts the loss by ~70 %.  This is exactly
why the EM method keeps `RP` (-> nu_eff, the field part) and `ECL` (the loss
part) as SEPARATE cell outputs -- do not derive one from the other.
"""

ERROR_ESTIMATOR = """\
## Equilibrated error estimator (Schöbinger-Hollaus 2024)

Reference: M. Schöbinger, K. Hollaus, IEEE Trans. Mag. 60(5):
7401906, May 2024.

### What it provides

An error estimator η_h that gives:
1. **Guaranteed upper bound**: ||u - u_h||²_E ≤ η_h²
2. **Local distribution**: η_h = sqrt(Σ_K η_K²) where η_K is per-element
3. **No unknown constants** (unlike residual-based estimators)

### Theoretical foundation — Prager-Synge

For elliptic problems with bilinear form a(u,v) = ⟨A∇u, ∇v⟩, if you
can construct an **equilibrated flux** σ_eq satisfying:
  - `div σ_eq = -f` (local equilibrium)
  - `σ_eq · n = 0` on Neumann boundary

Then the error in the H¹-energy norm is **bounded above** by the
"distance" between σ_eq and the discrete flux A∇u_h:

  ||u - u_h||²_E  ≤  ||A^{-1/2}(σ_eq - A∇u_h)||²_{L²}

The estimator η_h = ||A^{-1/2}(σ_eq - A∇u_h)||_{L²} is computable.

### Extension to 2D-1D MSFEM T-formulation

The Schöbinger-Hollaus 2024 paper:
1. Constructs σ_eq for the T-formulation eddy-current problem.
2. Adapts to the 2D-1D MSFEM ansatz (the equilibrium condition is
   modified by the polynomial-thickness expansion).
3. Shows local error captures both **discretization error** (macro
   mesh too coarse) and **truncation error** (too few polynomial
   modes in the thickness expansion).

### Application to motor analysis

Adaptive mesh refinement driven by η_K:
- If macro-mesh error dominates → refine the 2D cross-section mesh
- If polynomial-truncation error dominates → add one more polynomial
  mode k → k+1

This is critical for production-grade motor simulation: the user
can't be expected to know a priori how fine the mesh should be or
how many polynomial modes are needed.  The estimator decides.
"""

MSFEM_3D = """\
## Full 3D MSFEM for laminated iron cores (Hollaus 2020s monograph-style paper)

Reference: K. Hollaus, "A MSFEM to simulate the eddy current problem
in laminated iron cores in 3D", COMPEL (16-page comprehensive paper),
sole author: Karl Hollaus, TU Wien.

### Where this fits

The 2D-1D version covers machines with axisymmetric end effects (the
"main active zone" of a rotating machine).  But end-regions, leakage
flux into end-shields, and 3D-asymmetric structures (e.g. an
electromagnet with vertical lamination plus a top-yoke perpendicular
lamination) need **true 3D MSFEM**.

### MSFEM ansatz in 3D

  A(x, y, z, t)  =  A_macro(x, y, z, t)  +  Σ_k φ_k(ζ) A_micro,k(x, y, z, t)

where ζ = z/d is the local through-thickness coordinate.  Now BOTH
A_macro and A_micro are 3D — the macro is on a coarse mesh that
ignores the lamination, the micro is on a finer mesh oriented with
the lamination.

### Implementation variants studied

1. **Frequency domain** with Biot-Savart source (linear materials).
2. **Higher-order MSFEM** — Σ_k with k up to 4 polynomial modes.
3. **Time stepping** (for nonlinear ν(B), DC offset).
4. **Harmonic balance** (for nonlinear ν(B) with single-frequency
   excitation; gives stationary periodic solution directly).

### Key engineering result

For a U-core transformer with 0.3 mm sheets at 50 Hz:
- Standard 3D FE (every sheet meshed): ~10 million DOF, impractical
- MSFEM: ~50,000 DOF, < 1% error vs sheet-by-sheet reference
- Including edge effects (where eddy currents close at sheet edges)
- Faithful reproduction of distributed losses

### Application path for radia / NGSolve

NGSolve has all the primitives needed:
- `H1(..., order=p)` for both macro and micro spaces
- `FESpace([fes_macro, fes_micro_1, fes_micro_2, ...])` composite
- `CoefficientFunction` for the polynomial φ_k(ζ) projections
- `BlockMatrix` for the coupled system

The 16-page paper has the full discrete weak-form derivation; this
is a 1-2 month implementation project starting from there.
"""

NONLINEAR_ESI = """\
## Nonlinear effective surface impedance in scalar potential (Hollaus-Hanser-Schöbinger)

Reference: K. Hollaus, V. Hanser, M. Schöbinger, "A Nonlinear
Effective Surface Impedance in a Magnetic Scalar Potential
Formulation", 7p (in IH application folder).

### When to use this

Scenario: very large conductive (and possibly ferromagnetic) component
where the penetration depth δ is << overall dimensions, AND the
material has a nonlinear B-H curve.

Example: a thick solid-iron magnetic shield at 1 kHz.  δ ~ 5 mm in
mild steel, shield thickness 50 mm.  Eddy currents fill only the
outer ~10% of the shield volume.  Meshing this thin layer with FE
costs 10⁶+ DOFs; with effective-surface-impedance, you replace it by a
boundary condition on the shield surface and solve a much cheaper
magnetic scalar potential everywhere.

### Cell problem for nonlinear effective SI

For LINEAR material, the effective surface impedance is given
analytically by `Z_s = (1+j) / (σ δ)`.

For NONLINEAR material with `μ(|H|)`, no analytical form exists.  The
Hollaus et al. recipe:

1. **Cell problem**: 1D semi-infinite domain perpendicular to the
   surface, drive with sinusoidal H_t on the surface, solve the
   1D nonlinear PDE for the in-plane current density.
2. **Apparent power matching**: compute the *complex* power
   `S = E_t · H_t* / 2` averaged over one period.
3. **Extract Z_s(|H_t|, f)**: define the surface impedance that gives
   the same complex power in the equivalent linear problem.

This is structurally identical to **ESIM** (Effective Surface
Impedance Method, in `radia_mcp.ih.esim_cell_problem`) which the
Sugahara lab already implements in Radia for induction-heating
workpieces!

### Equivalence with Radia's existing ESIM

| Aspect | Radia ESIM (existing) | Hollaus Nonlinear ESI |
|--------|------------------------|----------------------|
| Cell problem | 1D PDE, B-input drive | 1D PDE, H-input drive |
| Output | Z_s(|H|, f) table | Z_s(|H|, f) table |
| Coupling to global | SIBC Robin on workpiece | SIBC on shield surface |
| Formulation | A-formulation | Magnetic scalar potential |
| Nonlinearity | ✓ | ✓ |
| Complex permeability | ✓ | (implicit) |

The two are **equivalent up to formulation choice**.  Radia ESIM
+ A-formulation handles IH workpieces; Hollaus ESI + scalar potential
handles magnetic shields and similar large-conductive-piece
problems.  Same physics, different application.

### Cross-reference

- `radia_mcp.ih.esim_cell_problem` — Radia's implementation
- `radia_mcp.motor.hollaus_eddy_knowledge.nonlinear_esi` — this section
- The Sugahara lab is already producing this kind of work.  This is a
  natural place to cite Hollaus 2024 as the related TU Wien
  implementation in the scalar-potential setting.
"""

MSFEM_PLUS_MOR = """\
## MSFEM + MOR (Hollaus-Schöberl-Schöbinger, IEEE TMAG 56(2), Feb 2020)

Reference: K. Hollaus, J. Schöberl, M. Schöbinger, "MSFEM and MOR
to Minimize the Computational Costs of Nonlinear Eddy-Current
Problems in Laminated Iron Cores", IEEE Trans. Mag. 56(2):7508104,
2020.

### Cascading two reductions

MSFEM already reduces 3D-laminated to 2D-1D-polynomial (~100×
reduction).  But for nonlinear time-stepping with many time steps
(transformer harmonics, motor PWM), the MSFEM cost per step is still
significant.

**Idea**: use **MSFEM to cheaply generate snapshots** for a reduced
basis (POD / Krylov), then **MOR substitutes a small reduced model
for the full MSFEM** at most time steps.

### Workflow

1. Run MSFEM on the full nonlinear time history at coarse Δt → get
   snapshots `{A(t_1), A(t_2), ..., A(t_M)}` of the MSFEM solution.
2. SVD the snapshot matrix → keep N_r dominant modes (typically
   N_r << ndof_MSFEM).
3. Project the full MSFEM system onto the N_r-dim reduced basis.
4. Run the **reduced** model with fine Δt for the production
   simulation — each step is just an N_r × N_r solve.
5. Recover full-MSFEM field via the basis expansion when needed.

### Quantitative result (paper)

For a small transformer with 50 Hz nonlinear excitation:
- Full MSFEM: 75 s per simulation
- MSFEM + MOR (with N_r = 30 modes from 5-cycle snapshots): 2.4 s
- Loss accuracy: < 0.5%
- Flux-density distribution: visually identical

### Connection to radia_mcp.mor and Sugahara lab CLN work

The Hollaus-Schöberl-Schöbinger snapshot/POD approach is **one of**
the canonical MOR techniques.  The Sugahara lab's specialty is
**Cauer Ladder Network (CLN)** MOR (Sugahara as co-author on 5+ CLN
papers; see `radia_mcp.mor.cln_knowledge`) — a **structurally
different** Lanczos-based MOR.

For motor analysis specifically:
- POD-snapshot from MSFEM → Hollaus 2020 approach.  Works for any
  problem.
- CLN from PEEC-coupled FE → Sugahara approach.  Requires loop-star
  decomposition, but gives directly-interpretable RL ladder network
  (good for circuit simulator hand-off).

Both approaches can be combined: use MSFEM for the FE side, hand off
to CLN for the circuit side.
"""

MOTOR_STRATEGY = """\
## Eddy-current-capable motor analysis — synthesis

The user's original concern: "Lange-Henrotte-Hameyer doesn't handle
eddy currents."  This is **correct as stated** — Lange-HH explicitly
drops `σ ∂A/∂t = 0` in the FE subproblem.  But eddy currents matter:
- Cage rotor bars of induction motors → core operating principle
- Lamination iron loss → 30-50% of total motor loss
- PM eddy currents in PMSM → demagnetization, hot-spots
- Stator coil AC copper losses (skin, proximity) → high-frequency
- Shaft voltage / common-mode currents → bearing failures

### Solution map for motor analysis with eddy currents

| Eddy-current source | Solution | radia_mcp module |
|---------------------|----------|------------------|
| Lamination losses (homogeneous) | Effective Material (Hanser-Schöbinger-Hollaus 2024) | `motor.hollaus_eddy_knowledge.effective_material` |
| Lamination losses (local hot-spots) | MSFEM 2D-1D | `motor.hollaus_eddy_knowledge.msfem_2d1d` |
| Cage rotor bars (induction motor) | Per-bar conducting region with σ ∂A/∂t = jωσA (FD) or BDF (TD) | Standard NGSolve, ONELAB im_3kW.pro pattern |
| Stator coil AC loss | Carstensen analytical Dowell + post-process | `motor.henrotte_lineage_knowledge.carstensen_2007` |
| PM eddy currents | Solid-PM mesh with σ_PM, FD/TD | Standard NGSolve, requires PM region with finite σ |
| Shaft voltage / common-mode | Darwin TD with parasitic C | `motor.darwin_model_knowledge.formulation` |
| Magnetic shield (thick solid) | Nonlinear ESI in scalar potential | `motor.hollaus_eddy_knowledge.nonlinear_esi` |

### Recommended hybrid stack for production motor analysis

```
┌─────────────────────────────────────────────────────────────────┐
│  Outer loop: scipy ODE on circuit + mechanical state            │
│    L_inc · di/dt + R i = v - e_bemf,   J dω/dt = T_em - T_load  │
│    (Lange-Henrotte-Hameyer 2009 framework)                      │
└──────────────────────────┬──────────────────────────────────────┘
                           │ at each Δt_FE
┌──────────────────────────▼──────────────────────────────────────┐
│  FE refresh layer                                               │
│  ───────────────────────────────────────────────────────────    │
│  1. Stator iron: Effective Material with ECL(|B|, f) homog.     │
│     (Hanser-Schöbinger-Hollaus 2024)                            │
│  2. Rotor iron: Effective Material similarly                    │
│  3. Cage bars (if IM): explicit σ ∂A/∂t with jω at slip freq    │
│  4. PM (if PMSM): rotated Br + per-magnet σ ∂A/∂t for high freq │
│  5. Air gap: sliding-band BC (or moving-band with remesh)       │
│  6. Stator winding: analytical Carstensen Dowell post-process   │
└─────────────────────────────────────────────────────────────────┘
```

### Implementation roadmap (next Stage-2 steps)

In `calc_motor_transient.py`:
1. **Add lamination loss option** `--lamination-effective` that uses a
   table-based `ν_eff(|B|, f)` and post-processes ECL.
2. **Add rotor cage option** `--rotor-cage` that injects σ on the
   `rotor_bar_*` regions, requiring jω-based complex FE (or BDF TD).
3. **Add PM eddy current option** `--pm-sigma 1e6` for solid PM σ.
4. **Add shaft voltage option** `--darwin-shaft-cap-pF` that adds the
   Darwin scalar potential block.

Each is incremental and can be golden-tested independently against
ONELAB cage-rotor IM templates (im_3kW.pro etc.).

This roadmap is what makes radia_mcp.motor genuinely competitive
with JMAG / ANSYS Maxwell for eddy-current-aware motor analysis.
"""

NONASYMPTOTIC_HOMOGENIZATION = """\
## Schöbinger-Hollaus-Tsukerman 2020 — Nonasymptotic Homogenization

Reference: M. Schöbinger, K. Hollaus, I. Tsukerman, "Nonasymptotic
Homogenization of Laminated Magnetic Cores", IEEE Trans. Mag. 56(2):
7509504, Feb 2020.  DOI: 10.1109/TMAG.2019.2943463.

### Why "nonasymptotic" is in the title

Classical homogenization theory (Bensoussan-Lions-Papanicolaou 1978)
assumes the **period of the microstructure goes to zero relative to
the wavelength / penetration depth** — the so-called asymptotic limit.
For laminated iron cores at high frequency this assumption FAILS:

- Lamination period: ~0.4 mm (iron + insulation)
- Penetration depth δ at 1 kHz in iron (μ_r=5000, σ=4e6 S/m): ~0.5 mm
- Period / δ ≈ 0.8, NOT << 1

So classical homogenization gives **wrong effective permeability** in
this regime.  Schöbinger-Hollaus-Tsukerman 2020 provides a
**NONasymptotic** analytical homogenization that stays correct even
when period ~ δ.

### Two competing definitions of effective permeability

The paper makes a foundational philosophical point: **there are TWO
non-equivalent natural definitions** of μ_eff for laminations:

1. **Reaction-field-optimized μ_eff^(R)**: defined so that the
   B field OUTSIDE the homogenized core matches the actual
   laminated core (an "outside observer" can't distinguish them).
   This is the right definition for **interaction with the stator
   winding** — the externally measurable reaction.

2. **Loss-optimized μ_eff^(L)**: defined so that the eddy-current
   losses in the homogenized core equal the actual losses.

These are NOT equal in general!  At low frequency (δ >> period) they
coincide, but in the nonasymptotic regime they diverge — and the
paper shows the imaginary part of μ_eff^(L) develops a **peculiar
resonance** that does NOT appear in μ_eff^(R).

### Practical implication for motor analysis

When using the Effective Material approach
(`effective_material` topic), the practitioner must **choose** which
definition to use based on the question being asked:

| Question | Use |
|----------|-----|
| Stator winding inductance, terminal voltage | μ_eff^(R) (reaction-field) |
| Iron-loss prediction, thermal hotspots | μ_eff^(L) (loss-optimized) |
| Both? | Solve cell problem twice with different optimization |

The Hollaus-Hanser-Schöbinger 2024 "Effective Material" paper
(see `effective_material`) uses the loss-optimized definition by
default; the 2020 paper explains WHY this matters.

### Cell problem setup (paper §II)

- Cylindrical or planar lamination cell
- Periodic BC at z = ±a/2 (a = lamination period)
- Tangential H field specified at r = r_min (inner radius)
- Homogeneous Dirichlet on tangential H at r = r_max (outer)
- Net B_z = 0 (no through-thickness flux, typical for active zone)

The paper notes the last restriction ("no through-thickness flux")
is lifted in **subsequent publications** — including the 2026
Frljić-Hanser-Hollaus-Schöbinger paper on perpendicular-flux
homogenization (just found via Crossref, not yet in lab library).

### Why this paper matters for radia_mcp.motor

The current `calc_motor_lamination.py` cell problem uses the
A-formulation with imposed mean B (v2.1 + v2.5).  This corresponds
to the **reaction-field** definition of μ_eff (because the imposed
mean B fixes the externally-measurable flux).  For loss prediction
specifically — which is what most users actually want from the cell
problem — there's a subtle systematic bias documented by this paper.

**TODO for v3**: implement BOTH cell-problem variants (reaction +
loss optimized), expose via `--effective-target {reaction,loss}`
CLI flag.
"""

MSFEM_MOR_DEIM_DETAIL = """\
## Hollaus-Schöbinger 2024 — MSFEM + MOR + DEIM detailed formulation

Reference: K. Hollaus, M. Schöbinger, "MSFEM With MOR and DEIM to
Solve Nonlinear Eddy Current Problems in Laminated Iron Cores",
IEEE Trans. Mag. 60(3):7401004, March 2024.
DOI: 10.1109/TMAG.2023.3314835.

### Why MSFEM alone is insufficient (paper's central admission)

> "It turned out that the reduction of the computational costs with
> the aid of the multiscale finite element method (MSFEM) is not
> sufficient."

Even MSFEM (100× DOF reduction vs full sheet-by-sheet FE) is too
expensive for routine design iteration.  Need ANOTHER 100× reduction.

### The DEIM approach (paper §II.C)

Start with the nonlinear ODE system from MSFEM discretization:

  A(ν(u)) u + M du/dt = b(i)             (Eq. 3)

where A(ν) is the n×n nonlinear stiffness matrix and M is the
linear mass matrix.

Apply the fixed-point method (FPM, Hantila-style polarization):

  A(ν_FP) u^(k,l) + (1/Δt) M u^(k,l)
    = b(i^(k)) + A(ν_FP - ν(u^(k,l-1))) u^(k,l-1) + (1/Δt) M u^(k-1)
                                                    (Eq. 4)

The LHS is now LINEAR (ν_FP is constant), so it can be LU-factored
once.  The nonlinearity moves to the RHS — but evaluating it still
costs O(n) per iteration, which is the bottleneck.

### MOR via snapshots (paper §II.C, Eq. 5-8)

Collect m snapshot solutions {u^(1), ..., u^(m)} from a "training"
time-stepped run.  Form snapshot matrix S ∈ R^{n × m}.

Project (4) onto S:
  S^T (A + (1/Δt) M) S y = S^T f                    (Eq. 7)
  K y = g                                              (Eq. 8)

with K ∈ R^{m × m}, m << n.  Reduced solve is now O(m^3).

### The DEIM trick for the RHS (paper §II.D)

But `S^T f` still requires evaluating f at ALL n points, then
projecting.  **DEIM** (Chaturantabut-Sorensen 2010) approximates f
by interpolating at p << n SELECTED "magic points":

  f ≈ U_DEIM (P^T U_DEIM)^{-1} P^T f                   (DEIM eq)

where U_DEIM is a POD basis of f-snapshots and P selects the p
magic points.  Cost: evaluate f only at p points, NOT n.

### The structure-preserving twist (paper's contribution)

> "We propose structure preserving MOR to avoid the known large
> errors due to applying the DEIM to update the nonlinear right-hand
> side of the fixed point method."

Naive DEIM applied to f breaks the **symmetry / passivity** of the
underlying problem.  The 2024 paper introduces a structure-
preserving DEIM that:

1. Splits f into f_lin (linear part) + f_nl (nonlinear part)
2. Applies DEIM ONLY to f_nl
3. Preserves the energy structure of the original FPM

This avoids the "known large errors" without sacrificing the speed
benefits of DEIM.

### Speedup figures (paper §III experiment)

For a 3-phase transformer with 50 Hz nonlinear excitation:
- Full MSFEM nonlinear solve: reference (slow)
- MSFEM + POD (no DEIM, projection only): 5× faster
- MSFEM + POD + naive DEIM: 50× faster but ~5% error
- **MSFEM + POD + structure-preserving DEIM: 50× faster + <1% error**

This is the state-of-the-art for routine industrial nonlinear motor /
transformer iron simulation as of 2024.

### Practical implication for `calc_motor_lamination.py`

The current implementation uses linear ν (mu_r_iron constant).
Extending to nonlinear ν(B²) WITHOUT DEIM would make it slow for
practical use.  The Hollaus-Schöbinger 2024 method is the
canonical path forward:

**Stage-4 roadmap** (not yet implemented):
1. Add `--bh-curve FILE` for nonlinear ν(B²)
2. Implement fixed-point method (FPM) iteration on the cell problem
3. POD on solution snapshots from a training transient
4. Structure-preserving DEIM on the RHS
5. Validate against full nonlinear MSFEM on a small test problem

Estimated effort: 1-2 weeks for someone familiar with FPM + DEIM.
"""

HANSER_2025_CIRCUIT_COUPLING = """\
## Hanser-Schöbinger-Hollaus 2025 — A-formulation + circuit coupling

Reference: V. Hanser, M. Schöbinger, K. Hollaus, "Effective material
modelling for laminated iron cores with an A-formulation and
circuit coupling", COMPEL 44(5):832, 2025.
DOI: 10.1108/COMPEL-12-2024-0520 (open access, CC BY 4.0).
Funded by Austrian Science Fund grant 10.55776/P36395.

### What this paper adds beyond 2024 effective material

The 2024 paper (`effective_material` topic) developed the EM
homogenization with the **T,ε-ε formulation** and **imposed B**
driving.  The 2025 paper extends to:

1. **A-formulation** — uses magnetic vector potential A as primary
   variable (matches the radia_mcp.motor toolchain convention).
2. **External circuit coupling** — terminal voltage or current
   driven via an additional Schur-complement block, matching the
   Lange-Henrotte-Hameyer 2009 framework directly.

### Speedup figures (paper abstract)

- ~1% error vs full sheet-by-sheet FE reference
- **300× speedup on small cores** (transformer test)
- **5000× speedup on large cores** (machine test)

These numbers are dramatic and motivate full integration into the
radia_mcp.motor production stack.

### Schur-complement formulation for circuit coupling

For a coil winding with N_coil terminals, the coupled system is:

  [ K_FE    G_coil ] [ A_FE   ]   [ 0        ]
  [ G_coil^T   Z_c ] [ I_coil ] = [ V_term   ]

where:
- K_FE ∈ R^{n×n}: the global FE stiffness (with EM homogenized iron)
- G_coil ∈ R^{n×N_coil}: coupling between FE DOFs and coil currents
- Z_c ∈ R^{N_coil×N_coil}: per-coil lumped impedance
- V_term: terminal voltage source

Eliminate A_FE via Schur complement:
  (Z_c - G_coil^T K_FE^{-1} G_coil) I_coil = V_term
       ↑ effective coil impedance (N_coil×N_coil) — small!

This is identical in structure to the Lange-Henrotte-Hameyer 2009
field-circuit coupling (cf. `femm_transient_knowledge`), with the
key difference that K_FE here uses the EFFECTIVE (homogenized) iron
material — NOT the full sheet-by-sheet discretization.

### Why the speedup is so dramatic

The bottleneck in standard motor FE is:
- 50,000+ DOF for stator iron (with sheets resolved)
- per-time-step solve: O((50000)^2) for direct, faster for iterative
- combined with 1000+ time steps in a transient: hours per simulation

With EM homogenization + A-formulation + Schur:
- 500 DOF for stator iron (homogenized, no sheet structure)
- per-time-step solve: O((500)^2) = 10000× cheaper
- Schur reduction to N_coil DOF for the actual transient: ~3 DOF
  for a 3-phase machine
- 5000× speedup is consistent with this back-of-envelope

### Direct applicability to `calc_motor_lamination.py` v3

The paper IS the implementation roadmap for the next phase of
`calc_motor_lamination.py`:

1. **Cell problem** (already done in v2.5): solve once for
   μ_eff, σ_eff per (|B|, f).
2. **Global FE with EM** (already done in v2): replace stator iron
   ν with ν_eff in the global solve.
3. **Add Schur block for coil terminals** (NEW): expose
   `--drive voltage --terminal-V <val>` mapping to the Schur
   complement.  Compute effective coil R, L, plus iron losses.
4. **JSON output**: terminal current waveform, total iron loss,
   coil-equivalent impedance for the controller designer.

Estimated effort: 3-5 days starting from the existing v2.5 codebase.
The cell-problem v2.5 `--drive voltage` stub already raises
NotImplementedError pointing here — the 2025 paper IS the citation
for the implementation.

### Free / open-access status

CC BY 4.0 (Emerald Publishing).  Full PDF lives at
`public-safe curated corpus
material modelling for laminated iron cores with an A-formulation
and circuit coupling.pdf` (16 pages, 3.9 MB).
"""

HOLLAUS_2014_NONLINEAR_TWO_SCALE = """\
## Hollaus-Hannukainen-Schöberl 2014 — Nonlinear two-scale FE

Reference: K. Hollaus, A. Hannukainen, J. Schöberl, "Two-Scale
Homogenization of the Nonlinear Eddy Current Problem with FEM",
IEEE Trans. Mag. 50(2):7010104, Feb 2014.
DOI: 10.1109/TMAG.2013.2282334.
(Hannukainen is from Aalto University, Finland — Schöberl is
NGSolve's author.)

### The foundational MSFEM-for-nonlinear paper

This is the chronologically EARLIEST entry in the lab's Hollaus
collection (2014), predating the 2D-1D MSFEM development.  It
establishes the **two-scale ansatz** that all subsequent Hollaus
papers build on:

  ~A(x)  =  A_0(x)  +  φ(x/ε) · (0, A_1, 0)  +  ∇(φ · w)         (Eq. 4)

where:
- A_0(x): the **macro-scale** vector potential (slow variable)
- A_1(x), w(x): **micro-scale** correction fields
- φ(z) = micro-shape function, periodic in z, encoding the
  through-the-thickness lamination structure (Fig. 2 of paper)
- ε = lamination period

Substituting (Eq. 4) into the weak form of the nonlinear curl-curl
equation gives a coupled system in (A_0, A_1, w).

### Key technical insights

1. **Neglecting derivatives of A_1**: numerical experiments showed
   that keeping ∂A_1/∂x_macro made the solution LESS accurate.  So
   A_1 is treated as a constant inside each macro cell.
2. **Polynomial-order assignment per de Rham**: A_0 uses H(curl,Ω),
   A_1 uses L²(Ω_m), w uses H¹(Ω_m).  These are chosen per the
   de Rham complex (FEEC), giving an inf-sup stable system.
3. **Adapted integration rules**: the nonlinear ν(B²) coupling
   makes assembly the bottleneck.  The paper introduces special
   quadrature rules that accelerate FE assembly while preserving
   the macro / micro coupling.

### How this paper sits in the lineage

| Year | Paper | Builds on this |
|------|-------|----------------|
| 2014 | **Hollaus-Hannukainen-Schöberl** (this) | — (foundation) |
| 2019 | Hollaus 3D MSFEM (`msfem_3d` topic) | 2014 ansatz extended to 3D |
| 2020 | Hollaus-Schöberl-Schöbinger MSFEM+MOR (`msfem_plus_mor`) | adds snapshot/POD |
| 2024 | Hollaus-Schöbinger MSFEM+MOR+DEIM (`msfem_mor_deim_detail`) | adds nonlinear hyperreduction |

So this paper IS the seed.  Anyone wanting to learn MSFEM from
scratch should start here, not with the more recent papers.

### Practical use

Already cited in the existing `msfem_2d1d` and `msfem_3d` topics
as the original two-scale theory.  This `hollaus_2014_*` topic
captures the formulation directly from the source paper.
"""

HIERARCHICAL_ERROR_ESTIMATOR = """\
## Schöbinger-Hollaus 2021 — Hierarchical error estimator (3D MSFEM)

Reference: M. Schöbinger, K. Hollaus, "A Hierarchical Error
Estimator for the MSFEM for the Eddy Current Problem in 3-D",
IEEE Trans. Mag. 57(5):7401005, May 2021.
DOI: 10.1109/TMAG.2021.3062041.

### Why this is needed

The 2024 paper (`error_estimator` topic) provided an **equilibrated**
error estimator for the 2D-1D MSFEM (Prager-Synge based).  This 2021
paper is the earlier 3D version using a different mathematical
machinery: **hierarchical** error estimation via p-refinement.

### Hierarchical idea

1. Solve the MSFEM problem with polynomial order p (your chosen
   discretization).
2. Build a **richer space** with order p+1 (the "hierarchical
   enrichment").
3. The difference between the p+1 solution and the p solution is
   the estimated local error.

Mathematically: if the discrete bilinear form `a_h` satisfies the
**saturation assumption** `||u-u_{p+1}||_a ≤ γ ||u-u_p||_a` for some
γ < 1, then

  η_h := ||u_{p+1} - u_p||_a   satisfies   (1-γ²)^{1/2} η_h ≤ ||u-u_p||_a ≤ η_h

So η_h is a **two-sided bound** on the actual error.

### Adapted to MSFEM specifically (paper §III)

The enrichment p+1 must be done **separately for the macro nodal
elements and the micro edge elements**.  The paper adapts:

- **Nodal hierarchical elements** for A_0 (macro vector potential
  components in H¹-style spaces)
- **Edge hierarchical elements** for A_1, w (micro-scale fields
  in H(curl))

### Numerical evidence (paper §IV)

For a 3-laminate stack test problem:
- Uniform p-refinement: convergence rate O(h^p)
- Adaptive p-refinement driven by η_h: convergence rate ~O(h^(p+0.5))
- DOF count: ~3-5x reduction for the same accuracy

Paper concludes: "It is the first error estimator presented for the
3-D MSFEM for the eddy current problem."

### Relationship to 2024 equilibrated estimator

| Aspect | 2021 hierarchical | 2024 equilibrated |
|--------|-------------------|---------------------|
| MSFEM type | 3D | 2D-1D |
| Math basis | p-refinement saturation | Prager-Synge flux equilibration |
| Bound type | Two-sided (γ-dependent) | Guaranteed upper |
| Constants | Saturation γ (problem-dependent) | None |
| Implementation effort | Easier (just solve p+1 too) | Harder (reconstruct equilibrated flux) |

For PRODUCTION code, prefer the 2024 equilibrated approach (no
problem-dependent γ).  For research / quick prototyping, 2021
hierarchical is easier to add.

### Lab implication

If/when `calc_motor_lamination.py` v3 implements MSFEM proper
(not just Effective Material), 2021 + 2024 papers BOTH provide
canonical adaptive-refinement strategies.  Pick one based on
priorities (rigor vs. ease).
"""

FRLJIC_2026_PERPENDICULAR_FLUX = """\
## Frljić-Hanser-Hollaus-Schöbinger 2026 — Perpendicular flux in open cores

Reference: S. Frljić, V. Hanser, K. Hollaus, M. Schöbinger,
"Homogenization of the Eddy Current Problem in Laminated Open-Type
Cores Accounting for Perpendicular Flux", IEEE Trans. Mag. 62(5):
6300408, May 2026.  DOI: 10.1109/TMAG.2026.3680721.
First author S. Frljić from University of Zagreb, Croatia
(international collaboration extension of the TU Wien program).

### The "open-type core" problem

Most prior MSFEM/EM work assumed **closed-type cores** (e.g.
toroidal transformer cores) where flux circulates entirely inside
the iron — `<B_perpendicular_to_lamination> ≈ 0`.  But for
**open-type cores** used in:

- Power voltage transformers (substation service)
- Rural electrification distribution transformers
- E-core single-phase machines

a significant fraction of the flux is **perpendicular to the
laminations** (i.e. through-the-thickness B_⊥).  Standard
homogenization either ignores this (wrong) or requires resolving
each lamination (expensive).

### The 2026 contribution

Schöbinger-Hollaus-Tsukerman 2020 explicitly noted (cf.
`nonasymptotic_homogenization` topic):

> "We assume that the net magnetic flux in the z-direction is zero...
> This restriction will be lifted in subsequent publications."

This 2026 paper is THE subsequent publication.  Key extensions:

1. **Two-vector-potential formulation** — instead of single A,
   use two potentials so that B_parallel and B_perpendicular
   can be approximated independently.
2. **Anisotropic effective material** — the homogenized μ_eff is
   no longer a scalar; it's a 2x2 tensor with different
   components for parallel and perpendicular field directions:

   μ_eff = [μ_parallel  0          ]
           [0            μ_perpendicular]

   where μ_parallel is the standard sheet-conducting homogenization,
   and μ_perpendicular incorporates the through-thickness flux
   correction.
3. **Coarse mesh in stacking direction** — because the perpendicular
   flux is now properly handled by the anisotropic material, the
   global mesh can be relatively coarse in the stacking direction
   (where it previously had to resolve each lamination).
4. **Validation**: numerical example with discretized-lamination
   reference shows "very good agreement".

### Scope: linear, frequency-domain only

Paper §I: "Only the linear case in the frequency domain using the
phasor representation is considered."

So **nonlinear extension is future work** — likely combining with
Hollaus-Schöbinger 2024 DEIM (cf. `msfem_mor_deim_detail`) for
nonlinear ν(B).

### Relevance to motor analysis

For **closed slot motors** (standard PMSM, IM with conventional
stator winding), `<B_⊥> ≈ 0` and the 2020 / 2024 EM approach is
fine.  For:

- **Open-slot stators** (preferred for high-speed motors due to
  airflow / cooling)
- **E-core single-phase motors** (washing-machine, fan motors)
- **Stator end-region** (where the flux path turns out-of-plane)

the 2026 anisotropic μ_eff is necessary for accurate iron loss.

### Implementation in radia_mcp

`calc_motor_lamination.py` v5 (future) should expose:

- `--anisotropic-mu` flag to use μ_⊥ ≠ μ_∥
- Cell problem variant that solves for perpendicular response
- Tensor-valued MaterialCF for the global FE

Estimated effort: 1 week for the linear case, more for nonlinear.
"""

SECTIONS = {
    "intro": INTRO,
    "msfem_2d1d": MSFEM_2D1D,
    "t_phi_open": T_PHI_OPEN,
    "effective_material": EFFECTIVE_MATERIAL,
    "error_estimator": ERROR_ESTIMATOR,
    "msfem_3d": MSFEM_3D,
    "nonlinear_esi": NONLINEAR_ESI,
    "msfem_plus_mor": MSFEM_PLUS_MOR,
    "nonasymptotic_homogenization": NONASYMPTOTIC_HOMOGENIZATION,
    "msfem_mor_deim_detail": MSFEM_MOR_DEIM_DETAIL,
    "hanser_2025_circuit_coupling": HANSER_2025_CIRCUIT_COUPLING,
    "hollaus_2014_nonlinear_two_scale": HOLLAUS_2014_NONLINEAR_TWO_SCALE,
    "hierarchical_error_estimator": HIERARCHICAL_ERROR_ESTIMATOR,
    "frljic_2026_perpendicular_flux": FRLJIC_2026_PERPENDICULAR_FLUX,
    "motor_strategy": MOTOR_STRATEGY,
}


def get_hollaus_eddy_knowledge(topic: str = "intro") -> str:
    """Return Karl Hollaus MSFEM eddy-current knowledge.

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
