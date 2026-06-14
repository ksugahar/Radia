"""Henrotte–Hameyer–RWTH research-arc knowledge for radia_mcp.motor.

This module consolidates the lineage of energy-based, Whitney-form-friendly
electromagnetic machine-FE methods developed at Université de Liège
(Henrotte's origin), Université catholique de Louvain (his later
affiliation), and RWTH Aachen (Hameyer's IEM where he ended up).  All
papers are present in the Sugahara-lab library.

Coverage:
  axisym_1993      — Henrotte et al., A New Method for Axisymmetrical
                     Linear and Nonlinear Problems (foundation of
                     `radia.axifem` Henrotte basis)
  source_field_1997 — Dular, Henrotte et al., A Generalized Source
                      Magnetic Field Calculation Method (h-φ formulation
                      source-field with cuts)
  energy_hys_2006  — Henrotte-Nicolet-Hameyer, Energy-based vector
                     hysteresis (friction-like pinning force)
  variational_2013 — François-Lavet-Henrotte-Stainier-Noels-Geuzaine,
                     Variational FE-friendly extension
  jacques_2018     — K. Jacques (ULB/ULg PhD), 254p comprehensive
                     monograph
  carstensen_2007  — C. Carstensen (RWTH PhD, Hameyer Korreferat),
                     SRM winding eddy currents (187p)
  research_arc     — Synthesis: how these papers connect, and how to
                     use them together for motor analysis

For the 2009 field-circuit-coupling capstone of this arc, see
`femm_transient_knowledge` (Lange-Henrotte-Hameyer 2009).
For the 2004 EM force-density paper that bridges to forces, see
`radia_mcp.differential_forms.forces_knowledge`.
"""

from __future__ import annotations


AXISYM_1993 = """\
## Henrotte 1993 — Axisymmetrical Linear and Nonlinear Problems

Reference: F. Henrotte, H. Hedia, N. Bamps, A. Genon, A. Nicolet,
W. Legros, "A New Method for Axisymmetrical Linear and Nonlinear
Problems", IEEE Trans. Mag. 29(2):1352-1355, March 1993.

### The problem (slot-1)

The classical axisymmetric A-formulation with auxiliary potential
`V = A_φ / r` (i.e. λ = -1 in their notation) and first-order finite
elements suffers from **large saw-tooth induction errors** near the
symmetry axis whenever a magnetic material touches the axis.  The
error is **not** simply discretization error — refining the mesh does
not help in any practical sense.

The root cause: with `V = A_φ / r`, the radial derivative produces a
`1/r` term in the induction:

  B_z = (1/r) ∂(r A_φ)/∂r = (1/r) ∂(r²V)/∂r = 2V + r ∂V/∂r

The `2V` term has no physical meaning (only enforces continuity), but
*does* dominate the gradient near r = 0 when V is approximated by
*linear* shape functions.  Result: spurious oscillations.

### Henrotte's solution — the `s = r²` coordinate transformation

Change variable to `V = r·A_φ` (i.e. λ = +1) AND change coordinate to
`s = r²/2`:

  V(s, z) = a + b·s + c·z + d·s·z       (bilinear in (s, z))

This is **first-order in (s, z)** but **bilinear in (r, z)**.
Crucially, with this special quadrilateral element on axis, the
spurious `1/r` term **vanishes identically**:

  B_z(s, z) = 2(b + d·z)                (free of singularity!)

So linear shape functions on `(s, z)` give *exact* piecewise-constant
induction reproduction for axisymmetric problems with piecewise
constant materials — the same property linear elements have for plane
2D problems.

### Practical implication

This `{1, r², z}` basis (also called the **Henrotte basis** in the
modern literature) is the basis used in `radia.axifem` for the
magnetic A_φ curl-curl problem.  See
`docs/axifemm/FORMULATION.md` sections 5-6.

The lab CLAUDE.md policy "Axisymmetric FE: Henrotte for Magnetic,
Standard H1 for Scalar (FEMM-Canonical)" rests on this 1993 paper:
- Magnetic A_φ → Henrotte basis (this paper)
- Scalar Laplacian → standard H1 with 2πr Jacobian (FEMM convention)

### When this matters for motor analysis

- Single-phase shaded-pole motor analysis (axisymmetric simplification)
- End-winding leakage inductance via axisymmetric reduction
- Solenoidal accelerator-grade actuator analysis
- Any motor cross-section reducible to axisymmetric (some BLDC
  outrunner geometries)

For mainstream 2D radial-cut motor analysis (PMSM, IM, SRM), the
plane 2D A_z formulation suffices and Henrotte basis is unnecessary.
"""

SOURCE_FIELD_1997 = """\
## Dular-Henrotte 1997 — Generalized source-field construction

Reference: P. Dular, F. Henrotte, F. Robert, A. Genon, W. Legros,
"A Generalized Source Magnetic Field Calculation Method for Inductors
of any Shape", IEEE Trans. Mag. 33(2):1398-1401, March 1997.

### What problem this solves

In an h-φ (magnetic-scalar potential + h-field) formulation, you need
a *source field* `h_s` satisfying `curl h_s = j_s` (the imposed
inductor current).  Three classical paths exist:

1. **Biot-Savart law**: `h_s(x) = (1/4π) ∫ j_s(y) × (x-y)/|x-y|³ dy`.
   Exact but expensive (O(N_obs × N_src)).  Naïvely O(N²); requires
   acceleration to be tractable.
2. **Magnetic-vector-potential post-processing**: solve `curl A = B`,
   take `h_s = (1/μ) (curl A)`.  Fast but requires a second FE solve.
3. **Edge-element source field**: build `h_s` directly in the
   Whitney 1-form FE space.  Cheapest, but requires careful
   topological treatment.

This paper develops path (3) and extends it to **multiply-connected
conductors** (e.g., a coil loop is topologically a donut).

### Key insight: cuts handle the cohomology

For a simply-connected non-conducting region Ωcc adjacent to a
conductor Ωc, the curl-free h-field can be derived from a scalar
potential φ: `h = -grad φ`.  The challenge is that for a **loop**
(multiply-connected) conductor, the scalar φ would be multivalued —
you must introduce **cuts** Γ_i to make the domain simply connected.

In the discrete edge-element setting:

- Define edge basis functions `S_α` for each edge α
- Identify edges on cuts → contribute jump-DOFs `α_i` (the cohomology
  classes, one per generator of `H¹(Ωcc)`)
- The source field decomposes as:

      h = Σ_k h_k S_k  (interior edges of conductor)
        + Σ_n φ_n ψ_n  (boundary scalar potential)
        + Σ_i α_i c_i  (cut contributions — one per topological cycle)

  where each `c_i` is supported on a thin transition layer adjacent
  to the cut.

The `α_i` is **exactly** the coil current through that loop.  This is
how a circuit current `I_coil` enters the FE system as a discrete-
cohomology basis coefficient, not as a volumetric source term.

### Practical implication for radia_mcp.radia_ngsolve / motor

This construction is what ONELAB / GetDP do internally when you
declare a `Stator_Ind_*` region with a phase current `I_A`.  The
machinery is hidden in `machine_magstadyn_a.pro` — the user just
writes `Stator_Ind_Ap = Region[...]` and the framework constructs the
edge-element source field with cuts.

In NGSolve, the equivalent machinery requires:

- A `Periodic` or `Trace`/`Mortar` to label coil-loop master/slave
  faces (the "cut" surface)
- A `NumberSpace` per loop carrying the discrete-cohomology DOF `α`
- The source field `h_s` assembled as a `LinearForm` over coil
  volume

This is **NOT yet** implemented as a generic motor primitive in
radia_mcp.radia_ngsolve — currently the lab uses explicit
volumetric `J_s = I_coil · N_turns / A_coil · ê_coil` on the coil
region (the "simple" approach valid for stranded conductors without
skin effect).  Implementing the Dular-Henrotte 1997 source-field
construction would be a research-grade contribution.

### When to care

- Stranded-coil rated-frequency operation: simple volumetric J_s
  suffices.
- Solid-bar conductors with skin effect: Dular-Henrotte gauge matters
  (otherwise the discrete `div h_s ≠ 0` pollutes the AC solution).
- Multi-conductor cable / litz-wire: cut-handling is essential for
  per-strand current accounting.

See also: Bossavit-Henrotte 2004 Handbook, Chapter 5 ("Discrete
Magnetostatics") — same construction in textbook form.
"""

ENERGY_HYS_2006 = """\
## Henrotte-Nicolet-Hameyer 2006 — Energy-based vector hysteresis

Reference: F. Henrotte, A. Nicolet, K. Hameyer, "An energy-based
vector hysteresis model for ferromagnetic materials", COMPEL 25(1):
71-80, 2006.  DOI: 10.1108/03321640610634344.

### Why a new hysteresis model

Existing models (1990s state of the art):

- **Preisach**: phenomenologically powerful, but **no energy
  interpretation** — requires extra assumptions to use in coupled
  thermal/mechanical/circuit problems.
- **Jiles-Atherton**: starts from an energy interpretation but
  loses it during the differential development; does **not extend
  naturally to vector** (only scalar).
- **Bergqvist 1997**: shares some ideas with this paper's approach.

The Henrotte 2006 model is **intrinsically vector** (not a
vectorization of a scalar model) and remains energy-consistent
throughout.

### Core idea — friction-like pinning force

Decompose the applied field h into a *reversible* part `h_r` (derived
from a Helmholtz free-energy `U(J)`) and an *irreversible* part `h_i`
(the pinning force):

  h = h_r + h_i,        h_r = ∂_J U,        |h_i| ≤ χ (constant)

The pinning force `h_i` derives from a **non-differentiable convex
functional** (a friction term):

  h_i = χ · (˙J / |˙J|)        (when ˙J ≠ 0)
       ∈ {v : |v| ≤ χ}         (otherwise — the sub-differential)

Graphically: the tip of the applied field h is allowed to move freely
inside a sphere of radius χ centered at `h_r(J)`.  Magnetization
changes only when h leaves the sphere (the pinning threshold) — this
is the source of the hysteresis.

### Multi-cell extension

A single cell gives a too-rectangular B-H loop.  Stack `N` cells with
different pinning thresholds `χ_1 < χ_2 < ... < χ_N` and free-energy
parameters `μ_1, ..., μ_N`.  Each cell contributes proportionally to
the macroscopic magnetization.  This is structurally identical to the
Preisach decomposition but with energy-consistent cell-blocks.

### Implementation simplicity

The non-differentiability becomes a single `if` statement in code:

```
if |h - h_r(J)| < χ:  # inside the pinning sphere → reversible
    dJ/dt = 0
else:                  # on the boundary → friction-like update
    dJ/dt = (h - h_r(J) - χ · (h - h_r(J))/|h - h_r(J)|) / damping
```

This is *much* simpler than Preisach's history-dependent staircase
operator or Jiles-Atherton's nonlinear ODE.

### Connection to radia material models

This is the same family as Radia's `MatPlayHysteresis` and
`MatEnergyHysteresis` — the **Play model** (commonly attributed to
Mayergoyz / Bobbio) and the **Energy model** (Egger / Henrotte
school) are both implementations of the energy-based friction
framework.  CLAUDE.md "Hysteresis Materials (Play and Energy Models)"
documents the Radia C++ implementations.

For motor analysis:
- Stator iron lamination: use Play (faster) or Energy model
- Anisotropic grain-oriented iron: use the vector form of either
- Rotor PM demagnetization curve: Play model is sufficient (low
  hysteresis)
"""

VARIATIONAL_2013 = """\
## François-Lavet et al. 2013 — Variational FE-friendly extension

Reference: V. François-Lavet, F. Henrotte, L. Stainier, L. Noels,
C. Geuzaine, "An energy-based variational model of ferromagnetic
hysteresis for finite element computations", J. Comput. Appl. Math.
246:243-250, 2013.

### What this adds to Henrotte 2006

The 2006 model is **differential** (an ODE for J̇) — perfectly fine
for time-stepped solvers but awkward to embed in implicit FE.  The
2013 paper recasts it as a **variational** problem:

> At each time step, the new state (J, internal variables) is the
> **minimizer of a thermodynamic potential** that accounts for both
> the stored magnetic energy and the dissipated energy.

This has three advantages:

1. **Variational consistency** — the discrete Euler-Lagrange equations
   inherit the continuous symmetries (energy conservation, dissipation
   inequality).
2. **FE-friendly** — the minimization fits naturally into a Newton-
   Raphson nonlinear solve at every Gauss point in every FE
   element.  No separate update routine needed.
3. **Both energies always available** — stored magnetic energy and
   dissipated energy are tracked at *all* times, not just after a
   closed loop (which is the limitation of differential models).

### Thermodynamic foundation (paper Section 2)

First and second laws of thermodynamics for the magnetic system:

  u̇ = ḣ · ḃ - div q       (1st law, energy conservation)
  ṡ ≥ - div (q/T)          (2nd law, entropy)

Define internal energy `u(b, χ)` where `χ` is a set of internal
variables (the per-cell pinning states).  At each time step, with
applied field `h_n+1` known, find the new internal state by minimizing
the incremental energy + dissipation functional:

  (J_{n+1}, χ_{n+1}) = argmin {u(J, χ) - h_{n+1} · J + D(χ - χ_n)}

The dissipation `D` is the friction-like term from Henrotte 2006.

### Vector graphical interpretation

The applied field h is decomposed as `h = h_r + h_i` with:
- `h_r = ∂_J u`  (reversible part)
- `h_i = χ · ˙J/|˙J|`  (irreversible, |h_i| ≤ χ)

Geometrically, the tip of h lives inside a sphere of radius χ centered
at the tip of h_r.  Multi-sphere = multi-cell hysteresis.

### Practical implication for motor analysis

For NGSolve implementation, this variational form is the cleanest
path because:
- NGSolve's Newton iteration already does Gauss-point function
  evaluation
- The per-Gauss-point hysteresis update becomes one line: `J_new =
  minimize(...)` via SciPy `optimize.minimize` (or a custom Newton
  using the analytical Jacobian)
- Energy bookkeeping is automatic — useful for iron-loss accounting
  and torque-via-virtual-work consistency checks

This is also the approach behind the **Egger Schur-complement Newton**
inside Radia's `MatEnergyHysteresis` C++ implementation
(CLAUDE.md, Hysteresis Materials section).
"""

JACQUES_2018 = """\
## Jacques 2018 — Energy-Based Magnetic Hysteresis Models (PhD, 254p)

Reference: Kevin Jacques, "Energy-Based Magnetic Hysteresis Models —
Theoretical Development and Finite Element Formulations", PhD
dissertation, Université Libre de Bruxelles / Université de Liège,
November 2018.

This is the **canonical comprehensive reference** for the
energy-based hysteresis school.  Supervised by Christophe Geuzaine
(co-author of GetDP).

### Outline

| Ch | Title | What you'll find |
|----|-------|-------------------|
| 1 | Introduction | Motivation, dissertation goals |
| 2 | Magnetic Losses Modelling — State of the Art | Preisach / Play-Stop / Jiles-Atherton / Bergqvist comparison |
| 3 | **Energy-Based Hysteresis Model** | Fundamentals, reversible / irreversible split, single-cell, multi-cells, homogenization |
| 4 | **Discrete Energy-Based Implementations** | Vector Play model, variational functional minimization, full differential / angle-searching, numerical comparison |
| 5 | Anhysteretic curves and identification | How to fit the model parameters to measurements |
| 6 | Anisotropy and rotational hysteresis | Multi-axial extension, rotational loss |
| 7 | FE coupling | Insertion into a 2D / 3D FE solver, time-stepping, convergence |
| 8 | Numerical examples | Stator yoke section, transformer core, etc. |
| 9 | Conclusion + perspectives | |

### Why this dissertation is the right reference

- Side-by-side comparison of all major hysteresis-model families
- Detailed FE implementation (not just theory)
- Open-source code released alongside (GetDP/ONELAB)
- Validation against measurement data for several commercial iron
  grades

### When to consult

- Implementing a new hysteresis model in NGSolve from scratch →
  start here for the canonical reference structure.
- Choosing between Play / Vector-Play / Variational / Full-
  Differential implementations → Chapter 4 has the numerical
  comparison.
- Calibrating model parameters from measured B-H loops → Chapter 5.
- Lab library location:
  `W:/03_文献・論文/00_電磁界解析/磁気特性/ヒステリシス/12_Energy_Based/`
"""

CARSTENSEN_2007 = """\
## Carstensen 2007 — SRM Winding Eddy Currents (PhD, 187p)

Reference: Christian Carstensen, "Eddy Currents in Windings of
Switched Reluctance Machines", PhD dissertation, RWTH Aachen, ISEA
(Institut für Stromrichtertechnik und Elektrische Antriebe).
Advisor: Rik W. De Doncker.  Second reviewer: **Kay Hameyer**.
November 2007.

This monograph is the canonical reference for **AC copper losses in
slot windings**, including skin effect, proximity effect, and their
analytical decomposition.

### Outline

| Ch | Title | Key content |
|----|-------|-------------|
| 1 | Introduction | Motivation: SRM coil losses dominate at high frequency |
| 2 | **SRM Design and Verification** | Common design approaches, current density, thermal, end-effects, lamination, acoustics, winding design (number of turns, mechanical design), measurements (efficiency, loss distribution, performance) |
| 3 | **Fundamentals of Analytic Eddy Current Loss Calculation** | Skin effect in round conductor, current displacement in slot windings (skin + proximity), analytic loss calculation, definition of winding parameters, arbitrary current waveforms |
| 4 | **Finite Element Model** | Simulation domain, winding modeling, search coils, mesh, electrical circuit, materials, automated model generation |
| 5+ | Validation, parameter studies, etc. | (full TOC available via PyMuPDF) |

### Why this matters for motor analysis (Sugahara lab)

SRM is the harshest case because:
- Square-wave / chopped current → very high di/dt
- High-frequency harmonic content reaches into the kHz range
- Slot fill factor + bundled coils → severe proximity effect
- Single-tooth concentrated winding → strong inter-turn coupling

But the **methods generalize**: skin/proximity decomposition,
analytic loss models, FE search-coil flux determination — all are
directly usable for PMSM / SynRM / IM with inverter-fed PWM
operation.

### Recipe for AC copper loss in radia_mcp.motor

1. **DC base resistance** R_DC from `(σ · A_coil)^-1 · l_turn · N_turns`.
2. **Skin factor** k_skin from analytical Bessel formula (Carstensen Ch 3.1):
   - Round conductor in air: `R_AC / R_DC = ξ · ber·ber' + bei·bei') / (ber'² + bei'²)`
     with `ξ = d_w / (δ √2)`, δ = skin depth.
3. **Proximity factor** k_prox from slot-position Dowell formula
   (Carstensen Ch 3.2): a layer of N conductors in a slot has
   per-layer factor that depends on the slot-magnetic-field strength
   at that layer's position.  Use the FE solution's `|H_slot|` map.
4. **Total AC copper loss** = R_DC · (k_skin + k_prox · |H_slot|²/I²)
   per harmonic, summed over harmonics.

This avoids needing to mesh individual strands (which would explode
the DOF count to millions for a 100-turn coil).

### Connection to PEEC

The radia_mcp.peec subpackage handles this differently — it uses
*surface* meshes on solid conductors + SIBC (Surface Impedance
Boundary Condition).  For bundled stranded coils (where each strand
is too thin to mesh individually), the Carstensen analytical-Dowell
approach is the right tool; PEEC + SIBC is for solid conductors with
resolvable skin depth.

### Library location

`W:/03_文献・論文/00_電磁界解析/MOR_モデル縮約/Eddy_Currents/Eddy
Currents in Windings of Switched Reluctance Machines.pdf` (187 pages)
"""

RESEARCH_ARC = """\
## The Henrotte–Hameyer–RWTH research arc — synthesis

The 11 Hameyer/Henrotte papers in the lab library are not random
selections — they are **a single coherent research arc** from the
mid-1990s ULg/RWTH school developing **energy-consistent E&M FE
methods**.  Read in this order:

```
1993  Axisymmetric s = r²/2 basis           [Henrotte et al.]
       ↓
1997  Source-field construction with cuts   [Dular, Henrotte et al.]
       ↓
2004  EM Force Density (Maxwell stress)    [Henrotte, Hameyer]
       ↓                                    (in differential_forms.forces_knowledge)
2006  Energy-based vector hysteresis        [Henrotte, Nicolet, Hameyer]
       ↓
2007  SRM winding eddy currents (PhD)       [Carstensen, Hameyer]
       ↓
2009  Field-circuit coupling (linearization)[Lange, Henrotte, Hameyer]
       ↓                                    (in femm_transient_knowledge)
2013  Variational hysteresis (FE-ready)    [François-Lavet et al.]
       ↓
2018  Energy-based hysteresis PhD (254p)    [Jacques]
```

Common philosophical thread:

1. **Energy is the primary quantity**.  Constitutive relations
   derive from energy / co-energy functionals.  No ad hoc additions
   that break energy bookkeeping.
2. **Discrete cohomology is structural**, not a workaround.  Cuts,
   loop currents, scalar potentials are organized around
   `H^k(Ω)` — Whitney-form / FEEC compatible.
3. **Vector formulations are intrinsic**, not vectorized scalars.
   This is what distinguishes the school from Preisach / Jiles-
   Atherton.
4. **FE-readiness is non-negotiable**.  Every paper provides explicit
   discrete formulations, not just continuous theory.

### How to use this knowledge for motor analysis

For a new motor analysis project, the lineage gives you a **default
methodology stack** built entirely on lab-library papers:

| Layer | Method | Lineage paper |
|-------|--------|---------------|
| Geometry  | 2D radial cut (default) or axisym | Henrotte 1993 if axisym |
| Source    | Edge-FE source-field with cuts    | Dular-Henrotte 1997 |
| Material  | Energy-based hysteresis (vector)  | Henrotte 2006 / 2013 / Jacques 2018 |
| Iron loss | Bertotti + variational hysteresis | Jacques 2018 |
| Copper loss | Carstensen analytical Dowell    | Carstensen 2007 |
| Field-circuit | Temporary linearization      | Lange-Henrotte-Hameyer 2009 |
| Force/torque | Maxwell stress + energy method | Bossavit-Henrotte 2004 (differential_forms.forces_knowledge) |

### Sugahara-lab implementation status

| Method | Status in radia / radia-mcp |
|--------|------------------------------|
| Henrotte basis (axisym magnetic) | ✓ implemented in `radia.axifem` |
| Energy-based hysteresis | ✓ Play and Energy C++ models in radia (`MatPlayHysteresis`, `MatEnergyHysteresis`) |
| Dular-Henrotte source-field | ✗ not yet — currently use simple volumetric J_s |
| Lange-Henrotte-Hameyer field-circuit | ✗ Stage-2 implementation planned (`calc_motor_transient.py`) |
| Carstensen analytical Dowell | ✗ not yet — PEEC SIBC handles solid conductors only |
| Variational hysteresis (FE) | partial — Energy model is differential, not variational |
| Maxwell stress + energy force | ✓ in differential_forms.forces_knowledge |

This list is the natural **roadmap for radia_mcp.motor maturation**.
"""

SECTIONS = {
    "axisym_1993": AXISYM_1993,
    "source_field_1997": SOURCE_FIELD_1997,
    "energy_hys_2006": ENERGY_HYS_2006,
    "variational_2013": VARIATIONAL_2013,
    "jacques_2018": JACQUES_2018,
    "carstensen_2007": CARSTENSEN_2007,
    "research_arc": RESEARCH_ARC,
}


def get_henrotte_lineage_knowledge(topic: str = "research_arc") -> str:
    """Return Henrotte-Hameyer-RWTH research-arc knowledge.

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
