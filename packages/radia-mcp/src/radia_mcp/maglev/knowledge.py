"""Magnetic levitation knowledge (incl. the lab's Radia-IEM + CLN maglev research).

Scope: MAGNETIC LEVITATION, unified -- both the maglev SYSTEMS (EMS/EDS
trains, SCMaglev, Halbach/Inductrack, magnetic wheels, PM/SC bearings)
AND the levitation FORCE physics (induction/eddy-current lift, EML
melting, active magnetic bearings, superconducting Meissner/pinning,
diamagnetic, Earnshaw + loopholes, and how to compute the force).  The
former linear-drive (LIM/LSM, end-effect) material was removed.  (The
force-physics half was previously a separate radia_mcp.levitation server;
it was consolidated into this one server.)

The lab's own maglev research line (CAE-AI Lab, Yano Takaaki + Sugahara)
is the headline content -- two topics, `radia_iem_fem` and
`cln_mor_control`, distilled from:
  - Yano & Sugahara, conf. digest E-3-1 "Construction of an IEM-FEM
    hybrid computation method toward high-accuracy magnetic levitation
    analysis" (CAE-AI Lab) -- Radia IEM (MMM/MSC) <-> reduced-potential
    FEM weak coupling for moving-magnet eddy-current levitation force.
  - Yano, master poster "Implementation of a 3D Multiport Cauer Ladder
    Network method for accelerating control-coupled analysis" -- CLN
    model-order reduction for real-time control-coupled maglev (TEAM 28).

Distilled from the lab archive (W:/.../99_applications/07_magnetic-levitation/;
folder/file names below are English translations of the Japanese originals)
plus open literature for standard maglev topics that the
lab does not have direct PDFs for (those sections are clearly marked
"open literature, not lab PDF" so users know to cross-check).

Primary lab PDFs read for this revision:

  - 07_magnetic-levitation/research/magnetic-levitation-by-eddy-currents.pdf
      (Murata, Yamagata HS) - 1-page Arago-disk demo, NdFeB array
      over rotating Al disk, eddy currents lift small NdFeB pendulum.

  - 07_magnetic-levitation/04_Kansai-Univ/PM-based-magnetic-levitation_original.pdf
      Saiki, Ino, Inada (Kansai Univ., 2021): "Optimization on Pole
      Pitch of Magnetic Wheels and Thickness of Metal Plate for
      Floating and Propulsion System using Permanent Magnets". 2D
      analytical EDS model with skin depth delta_s, secondary B/J
      fields, F_L = J x B, F_L/P scaling proportional to sqrt(pole
      pitch tau). Builds on Fujii (1995, 1997) IEEJ "Revolving
      Permanent Magnet Type Magnet Wheel" lineage.

  - 07_magnetic-levitation/01_patents/maglev_Sumitomo-Heavy-patent-1.pdf
      JP Pub. Patent 7-327337 (Sumitomo Heavy Industries): PM-repulsion
      axial magnetic bearing with eddy-current brake ring (Al
      annulus) for overspeed prevention. Glass-wool thermal
      insulator between magnet and ring.

  - 07_magnetic-levitation/01_patents/maglev_Sumitomo-Heavy-patent-2.pdf
      JP Pub. Patent 2007-215264 (Sumitomo Heavy Industries): planar
      mover with rotating-PM "drive head" (9 reluctance paths),
      air-bearing levitation, eddy currents in conductive base
      plate provide thrust. 4-head omnidirectional version with
      wireless DC-motor control. Stepper / LSM variants.

NOT read on this pass (size budget):
  - 09_linear-drive/*.pdf  (each 10-44 MB)
  - 03_textbooks/magnetic-levitation-and-magnetic-bearings.pdf  (155 MB textbook)
  - research/RocketScience_Dec2014.pdf  and other > 2 MB items

Sections marked "(open literature)" are textbook-level summaries
included so the knowledge file is self-contained; cross-check
against the unread lab PDFs when those become accessible.
"""


# Authoritative topic enum for the dispatcher tool (wired into
# `maglev_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "Magnetic levitation -- systems (EMS/EDS/PM/SC/Halbach) + force physics (induction/EML/AMB/SC/diamagnetic/Earnshaw) + the lab's Radia-IEM / CLN research line",
    # -- maglev systems + the lab's Radia/CLN research --
    "radia_iem_fem": "Radia IEM (MMM/MSC) <-> reduced-potential FEM weak coupling for moving-magnet eddy-current levitation force; A-phi/T-Omega; no re-mesh on magnet motion (Yano bachelor, lab research)",
    "cln_mor_control": "Cauer Ladder Network (CLN) model-order reduction for real-time control-coupled maglev: ~1/500 speedup, multiport matrix-CLN, 3D gauge A-phi/T-Omega/A-T, TEAM 28 (Yano master, lab research)",
    "physical_tensor_rom": "Physical (exterior-matched) polarizability tensor alpha(s) as a passive, stable LTI: AAA discovers the Stoll poles + NNLS passive residues, fitting the per-frequency 3D HCurl tensor; Kameari+Kelvin accumulation BREAKS DOWN for the general 3D body (rom_fit.py, lab research)",
    "pm_maglev_zero_power": "Zero-power passive PM levitation: Maxwell-Earnshaw constraint, axial PM bearings, halbach diamagnetism",
    "eddy_current_maglev": "Eddy-current EDS: Arago-disk physics, magnetic wheels (Fujii/Kansai 2D model), Inductrack",
    "sumitomo_heavy_industrial": "Sumitomo Heavy patents: PM axial bearing with brake ring (JP 7-327337); planar eddy-current mover (JP 2007-215264)",
    "kansai_research": "Kansai University Saiki/Fujii magnetic-wheel lineage: skin depth, F_L/P proportional to sqrt(pole pitch)",
    "scmaglev_eds": "Superconducting EDS (SCMaglev / Chuo Shinkansen): null-flux figure-8 levitation + guidance",
    "halbach_arrays": "Halbach array maglev / Inductrack: passive track coils, threshold-velocity self-stability",
    # -- levitation FORCE physics (force that lifts + stationary/industrial uses) --
    "induction_levitation": "Eddy-current (AC) lift force: jumping/Thomson ring, mutual-inductance circuit model, force vs frequency, perfect-conductor image limit",
    "eml_melting": "Electromagnetic levitation MELTING (EML): containerless molten metal, lift-vs-heat coupling, conical coils, microgravity TEMPUS/ISS-EML (ties to radia_mcp.ih)",
    "magnetic_bearings": "Active magnetic bearings (AMB): F = k_i*i + k_x*x linearization, negative position stiffness, bias current, radial/axial, flywheels; passive PM + SC bearings",
    "superconducting": "Superconducting levitation: Meissner (unstable) vs flux pinning (stable 3D), HTS bulk YBCO, field-cooled vs ZFC, frozen-image force model, hysteresis",
    "diamagnetic": "Diamagnetic levitation: f = (chi/2mu0)*grad(B^2), pyrolytic graphite over PM array, Geim water/frog, the passive+static Earnshaw loophole",
    "earnshaw_stability": "Earnshaw's theorem + its 5 loopholes (diamagnet/eddy/feedback/spin/SC-pinning); Braunbek limit; stability stiffness matrix",
    "force_computation": "Compute levitation force in Radia/PEEC/NGSolve: Maxwell stress tensor, virtual work dW/dz, time-average <J x B>",
    "benchmarks": "Validation references: TEAM 28 electrodynamic levitation, jumping-ring analytic, sphere-over-coil, EML lift coefficient",
    "all": "Concatenation of every topic above",
}


OVERVIEW = r"""
# Magnetic levitation landscape

> **POLICY (analysis method)**: solve maglev with **Radia (IEM) +
> NGSolve (FEM) weak coupling**, NOT standalone FEM.  Pure FEM must mesh
> the large PM<->conductor air gap and RE-MESH on every magnet move;
> Radia supplies the open-boundary external field analytically (no air
> mesh) and NGSolve reduced-potential FEM solves only the eddy reaction
> field.  See topics `radia_iem_fem` (the weak coupling) and
> `cln_mor_control` (CLN MOR for real-time control).  Mirrors the project
> CLAUDE.md policy "Maglev Analysis: Radia + NGSolve, Not FEM Alone".

## Levitation principles

| Method | Force source | Examples | Stability |
|--------|--------------|----------|-----------|
| **EMS** (Electromagnetic Suspension) | Attractive (controlled current) | Transrapid, JR maglev (low-speed) | Active (unstable, needs control) |
| **EDS** (Electrodynamic Suspension)  | Repulsive (induced eddy in track) | SCMaglev (Chuo Shinkansen)        | Self-stable above threshold velocity |
| **Superconducting EDS**              | Repulsive (SC coil on vehicle)    | SCMaglev                           | Self-stable |
| **PM repulsion (Halbach)**           | Repulsive (Halbach array on car)  | Inductrack (LLNL)                 | Self-stable above threshold |
| **PM repulsion (rotating wheel)**    | Repulsive eddy (PM dipole rotates over conductor) | Kansai-Univ magnetic wheel | Self-stable above threshold |
| **Diamagnetic**                      | Repulsive (bismuth, graphite)     | physics demos, SC pinning         | Stable but small force |
| **Hybrid (PM + EM)**                 | Bias by PM, control by EM         | Modern industrial bearings        | Active near zero current |
| **Passive PM (axial)**               | Repulsive (face-to-face PM)       | Sumitomo bearings (JP 7-327337)   | Axially stable, radially unstable (Earnshaw) |

## The lab's maglev research line (headline content)

The Sugahara / CAE-AI lab attacks the hard part of maglev design --
**fast, accurate eddy-current force between a moving magnet and a
conductor** -- with two complementary methods built on Radia + NGSolve:

| Topic | What | Why it matters for maglev |
|-------|------|----------------------------|
| `radia_iem_fem`   | Radia IEM (MMM/MSC) computes the open-boundary external field; reduced-potential FEM computes only the eddy reaction field; weak-coupled, fed back to demagnetisation | The magnet MOVES -> only the external field updates, **no re-mesh** of the air gap; IEM removes the air-region discretisation error |
| `cln_mor_control` | Cauer Ladder Network (CLN) model-order reduction turns the 3D eddy-current FEM into a compact equivalent circuit | Real-time control-in-the-loop maglev design at **~1/500** of full-FEM time; TEAM 28 benchmark |

These are validated on the standard eddy-current benchmarks: the
eddy-current FORCE/loss solver against TEAM Problem 7 (the canonical
"Asymmetrical Conductor with a Hole", 50 Hz harmonic, NOT a levitation
problem), and the levitation itself against TEAM Problem 28 (the
electrodynamic levitation device).

## Lab focus

- **Radia-IEM + FEM weak coupling** for moving-magnet eddy-current force
  (Yano & Sugahara) -- topic `radia_iem_fem`
- **CLN model-order reduction** for control-coupled maglev (Yano) --
  topic `cln_mor_control`
- **Magnetic-wheel EDS** -- Kansai-Univ collaboration (Saiki et al. 2021)
- **PM axial bearings** -- Sumitomo Heavy lineage (heavy-industrial)
- **Bearingless motor + WPT** (cross-link `radia_mcp.motor`, `radia_mcp.pcb`)

## Levitation FORCE physics (also in this server)

Beyond the maglev systems above, this server also covers the levitation
FORCE physics itself -- the force that lifts and its stationary /
industrial / laboratory uses:

- `induction_levitation` -- eddy-current (AC) lift; jumping/Thomson ring
- `eml_melting`          -- electromagnetic levitation melting (ties to IH)
- `magnetic_bearings`    -- active magnetic bearings (AMB), flywheels
- `superconducting`      -- Meissner vs flux pinning, HTS bulk
- `diamagnetic`          -- grad(B^2) levitation (graphite / water-frog)
- `earnshaw_stability`   -- Earnshaw's theorem + its 5 loopholes
- `force_computation`    -- Maxwell stress / virtual work / time-avg <J x B>
- `benchmarks`           -- TEAM 28, jumping-ring analytic, EML lift coeff

## Cross-references

- `radia_mcp.mor` (mor_cln) -- the CLN / Cauer ladder MOR theory
- `radia_mcp.fem` (potential_formulations) -- A-phi / T-Omega / A-T gauges
- `radia_mcp.team_benchmark.force_motion.problem_28` -- TEAM 28 Electrodynamic Levitation benchmark
- `radia_mcp.motor` -- analogous rotary motor analysis
- `radia_mcp.electromagnet` -- DC magnet design (EMS pole face B)
"""


RADIA_IEM_FEM = r"""
# Radia IEM <-> reduced-potential FEM weak coupling for moving-magnet maglev
*(Lab research: Yano & Sugahara, CAE-AI Lab, conf. digest E-3-1.)*

Maglev design hinges on the eddy-current ELECTROMAGNETIC FORCE between a
permanent magnet and a moving conductor across a large air gap.  Neither
standard method handles this well alone:

| Method | Strength | Weakness for maglev |
|--------|----------|---------------------|
| **FEM** | eddy currents, nonlinear materials | large PM<->conductor air gap must be meshed; **magnet motion forces re-meshing**; air-region discretisation error |
| **IEM** (integral element method = Radia MMM/MSC) | NO air mesh, exact open boundary [Chadebec 2006] | does not solve eddy currents efficiently |

## The weak coupling (the lab's answer)

Run IEM and FEM SEQUENTIALLY, exchanging fields (weak coupling):

```
  IEM (Radia MMM/MSC)                 reduced-potential FEM
  magnet region                       conductor region Omega_c
  --------------------                --------------------------
  analytic external field      ---->  source term
    A_ext, H_ext                        (reduced-potential RHS)
                                      solve ONLY the reaction field
                                        H_r, J  (eddy currents)
  demagnetisation  <----  feedback  ----'
```

- IEM computes the external field (A_ext, H_ext) ANALYTICALLY from the
  magnets -- no air mesh, open boundary handled exactly.
- That external field becomes the FEM source term via the **reduced
  (total-minus-external) potential** formulation [Biro 2000]; the FEM
  solves only the unknown REACTION field.
- The FEM reaction field (eddy current J, reaction H_r) feeds back to the
  demagnetisation analysis.
- **Because only the external field updates when the magnet moves, there
  is NO re-meshing** -- the decisive advantage for the moving-magnet
  maglev problem.

## The two reduced-potential formulations (both validated)

**A-phi method** [Biro 2000] -- uses the IEM A_ext.  Split the magnetic
vector potential A = A_ext + A_r; in the conductor Omega_c:

```
  sigma * d(A_r)/dt + curl( (1/mu) curl A_r ) = -sigma * d(A_ext)/dt
                                                  \_____ source _____/
```

**T-Omega method** [Biro 2000] -- uses the IEM H_ext.  Split H = H_ext +
H_r, with the reaction field via an electric vector potential T and a
magnetic scalar potential Omega:  J = curl T,  H_r = T - grad(Omega):

```
  curl( rho * curl T ) + mu * d/dt( T - grad Omega ) = -mu * d(H_ext)/dt
                                                         \_____ source ____/
```

## Validation (Yano & Sugahara digest E-3-1)

Test problem: a permanent magnet ROTATING + TRANSLATING above a copper
plate (the magnet-wheel / Arago class), dt = 0.0111 s, 181 steps.  The
two formulations are physically consistent but feed DIFFERENT physical
quantities (A_ext vs H_ext) to the FEM -- so cross-checking them is an
independent validation:

| Quantity     | Mean rel. error | Max rel. error |
|--------------|-----------------|----------------|
| Joule heat   | 3.38 %          | 10.65 %        |
| Lorentz force| 4.81 %          | 13.20 %        |

All quantities agree within ~5 % mean -> the coupled method is sound.
The eddy-current force/loss solver was also validated against TEAM
Problem 7 (eddy-current "Asymmetrical Conductor with a Hole") and the
levitation against TEAM Problem 28 (electrodynamic levitation device).

## Mapping to the Radia / NGSolve stack

| Role | Tool |
|------|------|
| IEM external field A_ext, H_ext | **Radia MMM/MSC** (`rad.Fld(obj,'a'|'h',pts)`; ObjHexahedron/ObjTetrahedron magnets) -- exact analytic, open boundary |
| reduced-potential FEM reaction field | **NGSolve** A-phi / T-Omega eddy-current solve on the conductor mesh |
| coupling of Radia field into FEM | `radia_mcp.fem.equivalence_source` (NearFieldSource), `radia_mcp.radia_ngsolve` RadiaField CoefficientFunction |

## Cross-references

- `radia_mcp.fem` (potential_formulations) -- A-phi / T-Omega / A-T gauges
- `radia_mcp.radia_ngsolve` -- Radia field -> NGSolve coupling
- `radia_mcp.team_benchmark.force_motion.problem_28` -- TEAM 28 validation
- topic `cln_mor_control` -- speeding the SAME eddy-current FEM up for control
- topic `eddy_current_maglev` -- the moving-magnet-over-plate physics
- Refs: Chadebec et al. 2006 (IEM open boundary); Biro 2000 (reduced
  potential A-phi / T-Omega); Yano & Sugahara digest E-3-1.
"""


CLN_MOR_CONTROL = r"""
# Cauer Ladder Network (CLN) model-order reduction for control-coupled maglev
*(Lab research: Yano, CAE-AI Lab master thesis, "Implementation of a 3D
Multiport Cauer Ladder Network method for accelerating control-coupled
analysis".)*

Magnetic levitation is OPEN-LOOP UNSTABLE, so it needs real-time feedback
control.  Designing and tuning that controller requires an
electromagnetic model running INSIDE the control loop -- and a full 3D
eddy-current FEM is far too slow for that.  The lab's answer is
model-order reduction (MOR) via the **Cauer Ladder Network**.

## The Cauer ladder equivalent circuit

CLN replaces the 3D eddy-current FEM model with a continued-fraction
(Cauer ladder) equivalent circuit whose port impedance is:

```
  Z(s) = V(s)/I(s)
       = R0 + 1 / ( 1/(s*L1) + 1 / ( R2 + 1 / ( 1/(s*L3) + ... ) ) )
```

The CLN basis is generated by an alternating recurrence (Kameari et al.
2018, IEEE TMag 54(3):7201804) that produces orthogonal current/field modes directly from the FEM
operators; truncating at N ladder stages gives an N-th-order reduced
model, and the error decreases MONOTONICALLY with the number of stages.
CLN is closely tied to matrix continued fractions.

## Result: TEAM 28, ~1/500 the time

The 1-axis (vertical, Z) control-coupled simulation of TEAM Workshop
Problem 28 (the electrodynamic levitation device) was completed with CLN:
dt = 0.1 ms, 20000 steps, at about **1/500** of the wall-clock time of a
conventional full eddy-current analysis -- fast enough for desktop
control design.  (Sugahara et al. 2023 extended CLN to problems with
CONDUCTOR MOVEMENT using constant basis functions, which is what makes
the moving-magnet maglev tractable.)

## Verified in-repo: CLN-reduced levitation FORCE vs height (2026-06-04)

A worked, verified example lives in
`examples/CLN/scripts/team28_levitation/`.  The coil-driven axisymmetric
eddy problem is `(K + s*N) X = F` (K = s-independent magnetostatic
operator, N = conductivity term, F = coil source); the CLN/Cauer
reduction is the Krylov subspace built from the COIL SOURCE
(`V0 = K^-1 F`, `V_{k+1} = orthonormalise(K^-1 (N V_k))`).  A **6-stage
CLN reproduces the full-FEM levitation force vs height to < 0.1%** and
recovers the levitation equilibrium **dZ = +4.1 mm** (lift == disk weight
~1.055 N; lab full-FEM ~+4 mm); the force converges in ~5 stages
(stage 3 = 0.14%, stage 5 = 0.000%).  The repo full-FEM baseline matches
the lab axisymmetric ground truth to 0.01%.  This is the first CLN
reduction carried through to the actual TEAM 28 levitation force (prior
CLN-on-TEAM28 work was decay-spectrum-only).  See `radia_mcp.mor`
mor_cln (applications) for the same example from the CLN-theory side.

## Multiport CLN (matrix Cauer ladder) -- for multi-axis maglev

A real maglev needs more than the vertical axis: lateral (XY) guidance
must be co-simulated.  CLN extends to a MATRIX continued fraction
(multiport), with admittance:

```
  Y(s) = ( G0^-1 + ( (s*L1)^-1 + ( G2^-1 + (s*L3)^-1 + ... )^-1 )^-1 )^-1
```

The multiport Cauer ladder is realised in Simulink; the extracted
Y-matrix converges monotonically to the truth value as ladder stages are
added (8-stage CLN taken as reference).  Two enabling steps:
- **CLN 3-dimensionalisation** (3D formulation)
- **CLN multiport-isation** (matrix CLN)

## 3D formulation: the gauge choice matters

In 3D, B = curl A leaves A non-unique (gauge freedom); a Coulomb gauge
(div A = 0) is imposed to reduce DOF.  Yano compared three Coulomb-gauged
formulations -- **A-phi, T-Omega, and A-T** -- for the CLN circuit
constants of a square-prism conductor against an accuracy-guaranteed
analytic solution.  The **A-T method matched the analytic solution best**
(Tanimoto, Yano, Sugahara & Nagamine 2025).  Picking the right gauge is
what makes the 3D multiport CLN accurate.

## Mapping to the Radia / NGSolve stack

| Role | Tool |
|------|------|
| CLN / Cauer ladder MOR theory | `radia_mcp.mor` (mor_cln, mor_cln_multiport) |
| Radia transient / reduced-order core | radia `cln_core`, `lanczos_reduction` (PRIMA) |
| eddy-current FEM the CLN reduces | A-phi / T-Omega / A-T -- topic `radia_iem_fem` |
| gauge formulations | `radia_mcp.fem` (potential_formulations) |
| benchmark | TEAM 28 (`radia_mcp.team_benchmark.force_motion.problem_28`) |

## Cross-references

- `radia_mcp.mor` (mor_cln / mor_cln_multiport) -- the CLN MOR theory
- topic `radia_iem_fem` -- the eddy-current FEM that CLN compresses
- `radia_mcp.fem` (potential_formulations) -- A-phi / T-Omega / A-T
- `radia_mcp.team_benchmark.force_motion.problem_28` -- TEAM 28
- Refs: Kameari, Ebrahimi, Sugahara, Shindo & Matsuo 2018 (CLN
  representation of eddy-current fields, IEEE TMag 54(3):7201804);
  Sugahara, Tanimoto, Takahashi & Matsuo 2023 (CLN with constant basis
  for conductor movement); Matsuo, Fujiwara, Kuriyama & Shindo 2019
  (multiport matrix CLN); Tanimoto, Yano, Sugahara & Nagamine 2025 (3D
  CLN gauge comparison).  (See `radia_mcp.mor` mor_cln for the canonical
  CLN paper list.)
"""


PHYSICAL_TENSOR_ROM = r"""
# Physical polarizability tensor alpha(s) as a passive, stable LTI
*(Lab research: radia.levitation, the "physical Stoll spectrum -> CLN/LTI"
route; src/radia/levitation/mixed_galerkin/rom_fit.py, 2026-06-20.)*

The maglev FORCE on a moving conductor is F ~ Re[alpha(s)] grad(B^2); the
conductor is fully described, per direction, by its eddy-current
polarizability alpha_i(s).  For control-in-the-loop simulation we want
alpha(s) as a PASSIVE, STABLE state-space LTI (drops into Simulink /
SPICE), NOT a per-frequency solve.  This is the maglev counterpart of the
IH/port CLN ladder (topic `cln_mor_control`): there the object is a port
admittance Z(s); here it is the open-boundary polarizability TENSOR
alpha_ij(s) of an isolated body in free space.

## The honest obstacle: Kameari + Kelvin accumulation BREAKS DOWN here

The natural idea -- get the physical (exterior-matched, free-decay /
Stoll) eddy spectrum by a Kameari A-T accumulation with a Kelvin open
boundary, then read off a Cauer ladder -- DOES NOT WORK for a general 3D
isolated-conductor-in-vacuum body.  The lab's own canonical script
`examples/levitation/research_cln/ngsolve_validation/
cuboid_521_kameari_kelvin_v15_canonical.py` is a documented BREAKDOWN demo:
even with the two known bug fixes applied, the iteration sign-flips L_1 at
stage 1 and the Schmidt energy norm grows x15.  This is structural for the
HCurl + vacuum-coupled BC class (the A_ext gauge is unbounded at infinity;
see `radia_mcp.radia_ngsolve` cln_3d).  Kameari accumulation reaches the
Stoll Cauer ladder to 0.000% ONLY for the SPHERE (the special
chi-susceptibility case) and for AXISYM bodies -- NOT for the general 3D
tensor.  So the physical-tensor LTI is NOT obtained by eigen-accumulation.

## The verified route: AAA + NNLS sample fit of the per-frequency tensor

The verified PHYSICAL tensor is the per-frequency 3D HCurl solve
`examples/levitation/ellipsoid/ellipsoid_alpha_tensor_3d.py` (gauged
complex HCurl + a FINE AIR SHELL that resolves the air reaction dipole =
the lift / Re[alpha] part; CompactAMS + COCR; ~2-3% vs the analytic
sphere).  rom_fit.py turns those frequency samples into a passive, stable
LTI

    alpha(s) ~ alpha_inf + sum_k g_k/(1 + s tau_k),  g_k >= 0, tau_k > 0,

(one RC relaxation state per pole + a feedthrough D = alpha_inf), via:

  1. `scipy.interpolate.AAA(s, alpha)` DISCOVERS the dominant real LHP
     poles -- they land on the physical Stoll decay times
     tau_n = mu0 sigma a^2/(n pi)^2 to ~0.00 %.
  2. pole set = those dominant poles UNION a log-spaced filler over the
     sample band (the filler captures the high-order tail AAA buries in
     Froissart pairs; it is ALWAYS present, so the fit degrades gracefully
     -- one coherent method, not a fallback chain).
  3. `scipy.optimize.nnls` real residues g_k >= 0 -> PASSIVE by
     construction; real negative poles -> STABLE by construction.

CONVENTION (a real gotcha): the fit form g_k/(1+s tau_k), g_k>=0 has
Im[alpha(j omega)] < 0 (causal / passive e^{+j omega t}).  A 3D HCurl eddy
solve returns the PHYSICS convention (Im > 0), so CONJUGATE the FEM samples
(np.conj(alpha)) before fitting -- else NNLS cannot match positive Im, the
residues collapse, and the band fit is ~58 % (a real failure).

Why not the obvious alternatives (probe history, do not re-walk):
  - raw AAA: excellent function fit (6e-10) but places Froissart +
    occasional RHP poles -> UNSTABLE LTI; you cannot re-LS residues on the
    full near-cancelling pole set (ill-conditioned -> garbage residues).
  - dense log-grid NNLS only: passive but grid-smeared (~3% floor; a fixed
    grid cannot land a pole exactly at the physical tau).
  - hand-rolled vector fitting: finicky (column scaling / relaxed
    constraint); not worth debugging vs AAA-discover + filler + NNLS.
  - METRIC TRAP: the sphere has alpha(0)=0, so a POINTWISE relative error
    blows up near DC (a phantom "980%").  Use max|fit-data|/max|data|.

## Verified (analytic sphere Stoll spectrum, pure numpy, FEM-free)

Cu sphere a=5mm, alpha(s) = 4 pi a^3 [ -1/2 + sum_n (3/(n pi)^2)/(1+s tau_n) ],
1 Hz..1 GHz, n_filler=20:
  - 16 states, band fit 1.95e-4, passive (all g_k >= 0), stable.
  - dominant poles 184.62 / 46.155 / 20.513 us == analytic
    mu0 sigma a^2/(n pi)^2 to 0.000 %.
  - alpha_inf = -784.9 mm^3 == the perfect-conductor flux-exclusion limit
    -2 pi a^3 = -785.4 mm^3.
Golden: tests/test_levitation_mixed_galerkin_golden.py
`::test_rom_fit_sphere_stoll_spectrum` + `::test_rom_fit_diagonal_tensor_mimo`.

Verified on the ANISOTROPIC body too (FEM triaxial): Cu ellipsoid 5x3x1.5 mm,
the verified 3D HCurl tensor sampled at 9 frequencies per axis (conjugated) ->
3 ROMs -> a diagonal MIMO LTI (13 states), band fit ~3.6-3.9 % (= the
per-frequency FEM data accuracy plus the sparse 9-point sampling).  Dominant
decay times shape-split: tau_z = 60.2 us > tau_x = tau_y = 34.6 us; D_diag =
[-102.5, -122.3, -218.8] mm^3 = the -V/(1-N_i) ordering |z|>|y|>|x| (short axis
strongest) -- matching the static `ellipsoid_alpha_tensor.py` HF anchors.

## API (radia.levitation.mixed_galerkin)

```python
from radia.levitation.mixed_galerkin import (
    passive_foster_fit, FosterROM, diagonal_tensor_state_space)

# sample the verified per-frequency tensor (or any alpha(s) data) on j omega
rom = passive_foster_fit(s, alpha, n_filler=20)     # -> FosterROM
rom.dominant_tau      # AAA-discovered physical Stoll decay times
rom.tau_n, rom.g_n, rom.alpha_inf, rom.band_fit_relerr
A, B, C, D = rom.state_space()                      # passive scalar LTI

# diagonal (principal-axis) 3D tensor -> one MIMO LTI
A, B, C, D, n = diagonal_tensor_state_space([rom_x, rom_y, rom_z])
```
Example `examples/levitation/physical_tensor_rom.py` (default = analytic
sphere, fast; `--fem` = triaxial ellipsoid per-frequency tensor -> diagonal
MIMO LTI).

## Scope / caveats (honest)

- The DOMINANT AAA poles ARE the physical Stoll decay times; the FILLER
  poles are an approximation basis, NOT individually physical.
- The fit is only as accurate as the per-frequency FEM data it consumes
  (~2-3% for a general body; the analytic sphere is the exact anchor).
- These are FOSTER poles tau_n, not Cauer rungs tau_pair[k] (systematic
  ~-6-8% offset; compare like-with-like -- `radia_mcp.radia_ngsolve`
  cln_3d POLICY).
- This is the EXTERIOR-MATCHED physical tensor.  The mixed-Galerkin bulk
  Foster (bulk_foster_via_eigen / bulk_foster_vector_via_eigen) uses the
  INTERIOR-PEC eigenmodes -- a different, interior model; a SIBC tail
  completes it, but it is NOT the same object as this exterior-matched fit.
- Deferred first-principles alternative: a singular-M free-decay GEP with
  an air shell -> the Stoll eigenmodes directly (route B); higher risk
  (gauge + singular-M + open boundary).  The sample fit is the robust route.

## Cross-references
- topic `cln_mor_control` -- the port-CLN ladder this mirrors for the
  open-boundary polarizability tensor
- topic `force_computation` -- F ~ Re[alpha] grad(B^2) the LTI feeds
- `radia_mcp.radia_ngsolve` cln_3d -- 3D Kameari + the Kelvin-accumulation
  open problem (why the eigen route is avoided here)
- `radia_mcp.mor` mor_cln / mor_cln_multiport -- CLN MOR + matrix-CLN theory
- Refs: Landau-Lifshitz ECM sec. 59 (sphere alpha = 4 pi a^3 G(x));
  Stoll 1974 (Bessel eddy spectrum); Nakata-Berthier / Gustavsen 1999
  (vector fitting context); AAA = Nakatsukasa-Sete-Trefethen 2018.
"""


PM_MAGLEV_ZERO_POWER = r"""
# Passive permanent-magnet levitation (zero on-board power)

## Maxwell-Earnshaw constraint

Earnshaw's theorem (1842): a configuration of charges (or PM
poles, or paramagnetic bodies) in **static** equilibrium cannot
be stable in all directions simultaneously under inverse-square
(or 1/r) potentials.

Concretely for PMs:
- **Axial** repulsion (face-to-face same-pole PMs) -> AXIAL stable,
  RADIAL unstable. Mechanical or active radial constraint required.
- **Radial** repulsion (concentric PM rings, opposing polarity) ->
  RADIAL stable, AXIAL unstable.
- **Full 3D PM-only equilibrium is impossible** without auxiliary
  diamagnet, superconductor, or feedback.

## Routes around Earnshaw

| Route | Mechanism | Notes |
|-------|-----------|-------|
| Diamagnetic stabilisation | Bismuth / pyrolytic graphite supplies negative susceptibility | Small force; lab demos (levitating frog @ 16 T) |
| Superconducting flux pinning | Type-II SC pins flux lines (vortex pinning) | YBCO / GdBCO disks; truly passive |
| Active EM control | Feedback PI on gap sensor + EM coil | Modern industrial bearings |
| Rotating PMs (gyroscopic) | Angular momentum stabilises radial drift | Levitron toy |
| Eddy-current EDS | Motion-induced repulsion | See `eddy_current_maglev` topic |

## Axial PM bearing -- Sumitomo Heavy lineage (lab PDF)

Reference: JP Pub. Patent 7-327337 "Permanent-Magnet Magnetic-Levitation
Rotating Device" (Sumitomo Heavy Industries, 1995).

```
   stator coil-pack
   +-------+
   | back  |
   | core  |
   +---+---+   <-- annular NdFeB, magnetised axially
       |
   * * * gap (~mm)  <-- same-pole face-to-face
       |
   +---+---+   <-- annular NdFeB on rotor
   |       |
   | brake |   <-- Al / permalloy annulus on rotor outer
   | ring  |       (with glass-wool thermal insulation)
   +-------+
        |
   motor shaft
```

Patent claim: axial PM repulsion gives semi-permanent low-friction
support; **Al ring on rotor periphery** acts as overspeed brake via
induced eddy currents. Reaction torque is proportional to slip, so
rotor self-limits to a stable constant speed in the face of drive-side
ripple. Glass-wool layer between magnet and Al ring prevents eddy
dissipation from heating the NdFeB above its irreversible-loss
temperature.

Reported result: speed variation within +-0.5 percent at the no-load
operating point. NdFeB R-Fe-B grade preferred over Alnico / ferrite
for compactness and high BHmax.

## Passive PM linear bearing

For linear guideways the same principle applies in 2D: two NdFeB
strips with N-pole-to-N-pole spacing give a vertical levitation
force, but lateral guidance must come from a mechanical V-rail
or a second pair of magnets perpendicular to the levitation pair.

## Halbach axial bearings

A Halbach-magnetised ring (4-segment or 8-segment) concentrates
flux on one face. Used in:
- Flywheel energy storage (vacuum chamber, low loss)
- Vacuum-pump rotor support
- Cross-link: see `halbach_arrays` topic for k-vector picture.

## Magnetic field expression (axial PM stack)

For two coaxial cylindrical PMs of radius `a`, thickness `t`,
magnetisation `M_s`, separated by axial gap `g`:

```
F_axial(g) ~ (mu_0 / 2) * M_s^2 * pi * a^2
             * [ f(g, t, a) - f(g + 2t, t, a) ]
```

where `f(.)` involves elliptic integrals (closed form in
`radia_mcp.radia_ngsolve.analytical_formulas` topic
`ellipsoid_demag_torque`). Radia (MMM/MSC) computes this exactly
for arbitrary shape.

## Cross-references

- topic `earnshaw_stability` -- the FULL theorem, all 5 loopholes, and
  the stability stiffness matrix; the "Routes around Earnshaw" section
  here is the PM-bearing-specific view of that general result
- `radia_mcp.radia_ngsolve.analytical_formulas` -- PM force closed forms
- `radia_mcp.electromagnet` -- compare to EMS (active attractive)
- `eddy_current_maglev` -- active route around Earnshaw via motion
"""


EDDY_CURRENT_MAGLEV = r"""
# Eddy-current EDS (Electrodynamic Suspension)

## Principle (Arago disk in reverse)

A permanent magnet (or PM array) moves relative to a conductive
plate (Cu, Al). By Lenz's law, eddy currents are induced in the
plate that **oppose** the change in flux. The interaction force
has two components:

- **Drag** (parallel to motion, slowing the magnet/plate)
- **Lift** (perpendicular to plate, repelling the magnet)

The lift-to-drag ratio is essentially the "Q" of the moving
electromagnetic system.

```
F_lift / F_drag ~ (v / v_crit)    for v << v_crit
                ~ const                    for v >> v_crit (saturated)
v_crit = 2 * rho / (mu_0 * d)             (d = plate thickness, rho = resistivity)
```

Below v_crit the response is resistive (mostly drag). Above v_crit
the response is reactive (lift dominates). For Al plates of a few
mm thickness, v_crit is on the order of 10-100 m/s.

## Lab demonstration (lab PDF)

Reference: "Magnetic Levitation by Eddy Currents" (Murata Kenji,
Yamagata HS).

Cheap kit:
1. 100-yen-shop USB fan (2x AA, ~3 V)
2. Replace fan blade with acrylic disk
3. Glue 8 NdFeB ring magnets, alternating N/S, around the disk
4. Hold a 0.5 mm Al sheet under the spinning disk in a petri dish
   -> Al sheet spins (Arago's disk)
5. Hang a small NdFeB on a copper-wire pendulum near the spinning
   Al disk -> NdFeB lifts due to repulsive eddy-current image

The same principle scales up to: induction maglev, magnetic
braking on roller coasters, retarder brakes on trucks, and (with
much larger magnetic wheels) experimental EDS vehicles.

## Magnetic-wheel EDS -- Kansai University lineage (lab PDF)

Reference: Saiki, Ino, Inada (Kansai Univ., 2021), "Optimization
on Pole Pitch of Magnetic Wheels and Thickness of Metal Plate
for Floating and Propulsion System using Permanent Magnets".
Builds on Fujii et al. 1995 / 1997 IEEJ Trans (in Japanese).

System: NdFeB ring of 8 alternating poles on a non-magnetic Al
hub, spun by a brushless DC motor (550 to 2000 KV) up to about
20 000 rpm, held a few mm above a Cu or Al plate.

### 2D analytical model

Travelling magnetic wave at the plate surface:

```
B_g(x, t) = B_g * sin(pi * x / tau + omega_s * t)
omega_s = v_s / r = 2 * pi * f          (electrical angular freq.)
f       = (p / 2) * (N / 60)            (p = pole count, N = rpm)
```

Skin depth into the conductor:

```
delta_s = sqrt( 2 * rho / (omega_s * mu_0) )
```

For Cu (rho ~ 1.7e-8 Ohm.m) at 7000 rpm with 8 poles: delta_s ~ 3 mm.
For Al (rho ~ 2.7e-8 Ohm.m) at the same speed: delta_s ~ 4 mm.

Secondary (induced) magnetic field inside the plate:

```
B_2x(z) ~ B_g * (sqrt(2) * tau / (delta_s * pi))
             * exp(-z / delta_s)
             * cos(pi * x / tau + omega_s * t - z / delta_s + pi/4)
B_2z(z) ~ B_g *  exp(-z / delta_s)
             * sin(pi * x / tau + omega_s * t - z / delta_s)
```

Secondary current density (single component, y direction):

```
J_2y = -(B_g / mu_0) * (2 * tau / (delta_s^2 * pi))
       * exp(-z / delta_s)
       * sin(pi * x / tau + omega_s * t - z / delta_s)
```

Levitation force density (Lorentz, J x B):

```
F_L = -B_g^2 / mu_0
      * exp(-2 z / delta_s)
      * (2 * tau / (delta_s^2 * pi)) * (sqrt(2) * tau / (delta_s * pi))
      * sin(...) * cos(... + pi/4)
```

### Scaling laws (key result)

Lift-per-watt:

```
F_L / P = k_0 * sigma * P_in
sigma   = tau / (v_s * pi * delta_s)          (lift/drag ratio in low-speed regime)
```

Peak lift-per-watt scales as `(F_L / P)_peak ~ tau^0.75`. The
empirical fit (Saiki 2021 Fig. 20) gives 1/4-power-decade
improvement when pole pitch goes from 1 cm to 10 cm.

### Engineering implications (Saiki 2021 summary)

1. **Bigger is better.** Magnetic wheels should be as large as
   structural and mass budgets allow. Diameter 4.5 cm -> 9 cm
   roughly doubled (F_L / P).
2. **Fewer poles is better.** Minimum 2 poles theoretically;
   4 poles practical (vibration trade-off).
3. **Plate thickness must exceed delta_s.** Below 3 mm of Cu the
   lift drops linearly because flux leaks through the plate.
4. **Motor characteristic dominates.** A high-torque low-speed
   motor outperforms a high-speed motor at the same input power,
   because lift saturates well below 20 000 rpm.
5. **F_L/P decays at very high RPM.** Saturation of the secondary
   horizontal field B_2x sets the ceiling, not the skin effect.

## Inductrack (LLNL Halbach EDS)

Reference: R.F. Post, D.D. Ryutov, IEEE Trans. Appl. Supercond.
10 (2000) 901. See `halbach_arrays` topic.

## Cross-references

- `radia_mcp.team_benchmark.force_motion.problem_28` -- TEAM 28
  Electrodynamic Levitation benchmark (force, drag vs. velocity)
- `halbach_arrays` -- Inductrack and Halbach magnetic wheels
- `sumitomo_heavy_industrial` -- planar mover uses same physics
- `radia_mcp.radia_ngsolve.analytical_formulas` topic
  `thin_plate_eddy_current` -- closed-form scaling for plate eddy losses
"""


SUMITOMO_HEAVY_INDUSTRIAL = r"""
# Sumitomo Heavy Industries patents -- industrial maglev

Two complementary patents (both readable in the lab PDF set).
They show the design-engineering pragmatics of PM levitation in
production-grade equipment.

## Patent 1: PM axial bearing with eddy-current brake ring
(JP Pub. Patent 7-327337, 1995)

**Application**: magnetic bearings, flywheels.

**Stack** (from top down):
```
Stator      : back-core / cylindrical core / annular NdFeB (axial mag)
Air gap     : ~ several mm
Rotor       : annular NdFeB (axial mag, SAME pole facing stator)
              cylindrical core / back-core
Brake ring  : Al / permalloy annulus on rotor outer periphery
              with glass-wool thermal isolation from the NdFeB
Drive       : shaft connected to external motor (motor 11)
```

**Operating principle**:
- Same-pole face-to-face NdFeB pair gives axial repulsion -> bearing.
- Drive motor turns rotor; reluctance-free axial gap means very
  low friction.
- If drive ripple or inertia tries to over-speed the rotor, the
  rotor's rotational frequency relative to the stator's PM field
  rises, eddy currents grow in the Al ring, the resulting drag
  torque self-limits the rotor speed.
- Glass-wool insulator prevents eddy-current dissipation from
  heating the NdFeB above its irreversible-loss temperature
  (NdFeB irreversible loss starts around 80-150 C depending on
  grade).

**Reported quantitative result**: rotational speed variation
within +- 0.5 % at the no-load operating point. Dimensions in
the embodiment: 80 mm inner / 99 mm outer / 10 mm thick NdFeB
on the stator; 70 mm inner / 99 mm outer / 7 mm thick NdFeB on
the rotor; 7 mm glass-wool isolator; 7 mm thick Al brake ring.

**Material choices**:
- Stator + rotor magnets: R-Fe-B (NdFeB) preferred over Alnico /
  ferrite for compactness.
- Brake ring: Al preferred (light, non-magnetic, high conductivity).
  Permalloy and pure iron also claimed -- magnetic ring sits on
  PM by attraction (no adhesive needed). Al ring requires glue.

**Variants** (claimed):
- Brake ring on stator outer periphery instead of rotor.
- Brake ring inside the cylindrical core inner bore.

## Patent 2: Eddy-current planar mover
(JP Pub. Patent 2007-215264, 2007)

**Application**: semiconductor handling, vacuum-chamber transport,
clean-room positioning. Air-bearing floats the mover; eddy currents
in the conductive base plate provide horizontal thrust.

**Mover assembly** (one unit):
```
Frame  (flat, light alloy)
+- Air bearing (compressed-air pads on lower face)
+- "Drive head" (laminated magnetic / non-magnetic stack)
   |   internal cavity contains ...
   +-- Rotating PM (2-pole NdFeB on horizontal shaft)
   +-- DC motor that turns the rotating PM
Reluctance paths : 3 to 9 ferromagnetic ribs inside the head,
                   connecting the rotating-PM cavity wall to the
                   underside of the head (facing the base plate).
Spacing of underside ribs ~ pole pitch of the eddy-current image.
```

**Base plate**: flat conductive metal (iron in the main embodiment).

**Operating principle**:
- DC motor spins the small NdFeB inside the head.
- As the PM rotates, each ferromagnetic rib magnetises N or S in
  sequence, producing a travelling magnetic field on the head's
  lower face.
- This travelling field induces eddy currents in the base plate;
  the JxB Lorentz force pushes the head (and therefore the mover)
  in the opposite direction to the travelling-field motion.
- Air bearing decouples vertical support from horizontal drive.

**Omnidirectional version** (4 heads, 90 deg apart):
- 4 DC motors individually driven.
- Cross-pair counter-rotation -> linear translation in one axis.
- All-same-direction rotation -> yaw about the central axis.
- Wireless command + on-board battery + laser positioning =>
  no umbilical, supports clean-room / vacuum operation.

**Variants** (Sumitomo claims):
| Variant | Base plate | Driving force | Operating mode |
|---------|------------|---------------|----------------|
| 1       | Conductive plate (Fe, Al) | Eddy currents in base | Induction-motor analog |
| 2       | Conductive plate + array of fixed PMs | Attractive synchronisation | LSM analog |
| 3       | Magnetic plate with periodic teeth | Reluctance attraction | Stepper analog |
| 4       | Variant 3 + coils on teeth | Active electromagnet | Hybrid PM + EM |

**Why this matters**: it is a complete, deployable PM planar
motor with 4-DOF control and no cables. The same physics
(travelling field over conductive plate -> eddy-current thrust)
underlies the magnetic-wheel EDS in `eddy_current_maglev`, but
here the head is **stationary in rotation** -- only the small
2-pole PM inside spins. Big practical advantage: small DC motor
runs at modest speed, large flux density at the base because the
ferromagnetic ribs amplify the path.

## Cross-references

- `eddy_current_maglev` -- shared physics (skin depth, JxB)
- `radia_mcp.motor` -- rotating PM core inside the drive head; linear-motor / LSM analysis
- `radia_mcp.electromagnet` -- variant 4 (coiled teeth)
"""


KANSAI_RESEARCH = r"""
# Kansai University -- Saiki / Fujii magnetic-wheel lineage

## Lineage

| Year | Authors | Where | Topic |
|------|---------|-------|-------|
| 1995 | N. Fujii, K. Ogawa, T. Matsumoto | IEEJ Trans. 115-D, No. 3, 319 | Revolving Permanent Magnet Type Magnet Wheels (founding) |
| 1995 | N. Fujii, G. Hayashi, Y. Sakamoto | IEEJ Trans. 115-D, No. 3, 319 | Self-driven Magnet Wheel + Maglev Test Car |
| 1997 | N. Fujii, K. Ogawa, M. Chida | IEEJ Trans. 117-D, No. 6, 768 | Magnetic Poles vs. Characteristics of Revolving PM Magnet Wheel |
| 1997 | K. Ogawa, Y. Horiuchi, N. Fujii | IEEE TMAG 33, No. 2, 2069 | EM force calculation for magnet wheels (English) |
| 2008 | J. Bird, T.A. Lipo | IEEE TMAG 44, No. 3, 365 | 2D steady-state FEM for electrodynamic wheel |
| 2018 | M. Flankl, T. Wellerdieck, A. Tuysuz, J.W. Kolar (ETH) | IET EPA 12, No. 3, 357 | Scaling laws for EDS in high-speed transport |
| 2021 | T. Saiki, K. Ino, M. Inada (Kansai Univ.) | (this paper, lab PDF) | Pole-pitch + plate-thickness optimisation |

## Headline results of Saiki et al. 2021

1. **Peak F_L / P proportional to tau^0.75** where tau is pole pitch.
   Validated experimentally between tau = 1.2 cm (small wheel,
   p = 8) and tau = 5 cm (large wheel, p = 4).
2. **Diameter doubling** (4.5 -> 9.0 cm) roughly doubled lift-per-power.
3. **Plate thickness > skin depth required**. Cu plate at 7 kRPM has
   delta_s ~ 3 mm; plates thinner than 3 mm bleed flux through.
4. **Motor characteristic matters more than nominal power**.
   At equal P_in, a low-KV motor with high torque at 4 kRPM
   outperformed a high-KV motor that hits 18 kRPM at the same
   input power. (Lift saturates well before the high-KV motor
   reaches its design rpm.)
5. **Reactance-resistance crossover** at v_crit determines
   lift-vs-drag regime. Below v_crit, the eddy phase is near 0
   (pure drag, Arago disk). Above v_crit, the eddy phase
   approaches +pi/2 (pure lift).
6. **B_2x crosses B_2z** at the rpm where torque starts to fall.
   This is the practical "knee" for designing a maglev pod: pick
   operating rpm just above this knee.

## Where the open-literature math lives

- Bird & Lipo 2008 (FEM)                          -- electrodynamic wheel
- Knowles 1982, Ko & Ham 2007 (transient circuit) -- EDS transient dynamics
- Carlstedt et al. 2016                           -- dimensional / similarity scaling
- Flankl et al. 2018 (ETH Kolar group)            -- high-speed scaling laws
- Post & Ryutov 2000 (Inductrack)                 -- passive Halbach + ladder

## Cross-references

- `eddy_current_maglev` -- the 2D analytical derivation
- `sumitomo_heavy_industrial` patent 2 -- engineering deployment
- `halbach_arrays` -- Inductrack route to higher F_L / P
- `radia_mcp.radia_ngsolve.analytical_formulas` topic
  `thin_plate_eddy_current` -- closed forms for plate eddy losses
"""


SCMAGLEV_EDS = r"""
# SCMaglev (Chuo Shinkansen) -- Superconducting EDS
*(Open literature, not lab PDF.)*

## System overview

**Vehicle**: each car carries 4 bogies; each bogie holds 2
superconducting magnets (SCM) on its sides. Each SCM is a NbTi
race-track coil cooled to 4.2 K by liquid helium (with auxiliary
GM coolers).

**Track**: U-shaped guideway. Each guideway wall has two coil
sets:
- **Levitation / guidance coils** (figure-of-8 / null-flux),
  short-circuited or differentially connected.
- **Propulsion coils** (3-phase LSM stator).

## Levitation -- null-flux figure-8

The figure-8 coil is wound such that at the vertical centre
position no net flux links the coil. Vertical displacement of
the SCM breaks the symmetry, induces a current that opposes the
displacement -> restoring force. This is **self-stable** above
the critical speed (typically 100-150 km/h).

Below the critical speed, auxiliary rubber wheels carry the load
(the Yamanashi test track and the Maglev demonstration line both
use this scheme).

## Propulsion -- ground-coil LSM

3-phase AC fed to the ground coils generates a travelling
magnetic wave. The SCM on the vehicle locks onto this wave like
a rotor of a synchronous motor. Inverter frequency is computed
from vehicle position (laser / inductive sensors).

```
v_max ~ 500 km/h (139 m/s) for Chuo Shinkansen design
```

## Guidance

The two figure-8 coils on opposite walls of the U-guideway are
cross-connected. Lateral displacement of the SCM induces a
restoring current that pushes the vehicle back to centre.

## Why superconductors

NbTi SCM:
- Operate near 250-700 kA-turn (current * turns)
- Field at centre: 4-5 T (compare to NdFeB peak ~ 1.4 T)
- Zero resistance => persistent-mode operation possible
- Cryogen consumption: small (lab on-board cooler + bottom-up LHe)

The huge MMF is what enables the 10 cm levitation gap, vs.
~ 10 mm for EMS (Transrapid). Bigger gap -> easier track tolerance,
easier curve negotiation, safer derailment margin.

## Comparison vs. Transrapid (EMS)

| Property                | SCMaglev (EDS)              | Transrapid (EMS)            |
|-------------------------|-----------------------------|-----------------------------|
| Mover magnet            | Superconducting coil        | Electromagnet (controlled)  |
| Levitation gap          | ~ 10 cm                     | ~ 10 mm                     |
| Active control needed   | No (above v_crit)           | Yes (always)                |
| Auxiliary wheels        | Yes (< v_crit)              | No                          |
| Track complexity        | High (guideway coils)       | Medium                      |
| Cryogenics              | Required                    | None                        |
| Max design speed        | ~ 600 km/h                  | ~ 500 km/h                  |
| Energy at standstill    | Very low (SC persistent)    | High (lift current)         |

## Cross-references

- topic `superconducting` -- the Meissner-vs-flux-pinning FORCE physics
  behind SC levitation (SCMaglev uses the null-flux figure-8 EDS variant;
  HTS-bulk pinning is the stationary-bearing variant)
- `eddy_current_maglev` -- physics common to SCMaglev levitation
- `halbach_arrays` -- PM alternative to SC magnets
- `radia_mcp.fusion` -- adjacent SC magnet design expertise
"""


HALBACH_ARRAYS = r"""
# Halbach arrays and Inductrack
*(Mostly open literature, references the lab PDFs where applicable.)*

## Halbach array

A Halbach array is a periodic arrangement of PMs whose
magnetisation rotates by 90 degrees from one segment to the next.
The result is a magnetic field that is strongly enhanced on one
side of the array and almost cancelled on the other side.

```
Pattern (rotating 90 deg between segments):

   ^   ->   v   <-      ^   ->   v   <-
   |        |        |        |
  N->     N|S      ->S       S|N
   strong side : peak |B| ~ 1.4 * B_r
   weak side   : residual |B| ~ 0.1 * B_r
```

(Direction-arrow diagrams approximate. See Mallinson 1973 IEEE TMAG
9 for the magnetic-tape derivation, Halbach 1980 NIM for the
synchrotron derivation.)

**Fundamental harmonic** of a 4-segment-per-period array on the
strong side:

```
B(y) = B_r * (1 - exp(-k * t)) * exp(-k * y) * (1 - sin(k * pi/(2N)) / (k * pi/(2N)))
```

where:
- `k = 2*pi/lambda` (wavenumber, lambda = spatial period)
- `t` = magnet thickness
- `y` = distance from strong face
- `N` = segments per period (4 is the standard "k=1" array)

For thick magnets (`k * t >> 1`) the factor `(1 - exp(-k t))` -> 1.

## Inductrack (LLNL Halbach EDS)

Reference: R.F. Post, D.D. Ryutov, IEEE Trans. Appl. Supercond.
10 (2000) 901, "The Inductrack: a simpler approach to magnetic
levitation".

**Vehicle**: Halbach array on the underside (strong side down).
**Track**: closely-spaced shorted coils (a "ladder" of loops).

Physics:
- The Halbach array's first spatial harmonic creates a travelling
  wave below the vehicle when in motion.
- The shorted ladder responds with induced currents whose phase
  shifts from in-phase (drag dominant) to 90 deg lag (lift
  dominant) as speed increases.
- Above the transition velocity, lift-to-drag ratio can exceed
  100.

**Advantages over SCMaglev**:
- No cryogenics (room-temperature PMs).
- No on-board power for levitation.
- Track is passive (just coils + crossbars).

**Disadvantages**:
- Halbach PMs are expensive in bulk (NdFeB cost).
- Magnetic field decay length on the strong side limits gap
  (lambda / (2*pi) characteristic decay).
- Drag at low speed is high (auxiliary wheels still needed below
  threshold).

**Test demonstrations**:
- General Atomics test track (San Diego), ~ 100 m, 1990s-2000s.
- 25 g passenger comfort limit -> practical for cargo + transit.

## Halbach planar / rotary applications

- **Magnetic gear** (Atallah / Howe 2001): two Halbach rings with
  different pole counts coupled by ferromagnetic flux modulator
  -- equivalent to gear ratio with no mechanical contact.
- **Halbach axial bearing**: see `pm_maglev_zero_power` topic.
- **NMR + MRI permanent magnet**: cylindrical Halbach with
  high uniformity inside the bore.

## Open lab questions

- The lab has not (in this PDF set) implemented a Halbach
  Inductrack. Saiki et al. 2021 use 8-pole alternating-PM
  wheels, NOT Halbach -- but the analytical machinery in
  `eddy_current_maglev` applies directly with B_g replaced by
  the Halbach first harmonic.

## Cross-references

- `eddy_current_maglev` -- shared physics (skin depth, thrust)
- `pm_maglev_zero_power` -- Halbach axial bearings
- `scmaglev_eds` -- SC alternative to Halbach
- `radia_mcp.radia_ngsolve.analytical_formulas` topic
  `fabri_solenoid` -- closed-form for related coil geometries
"""


INDUCTION_LEVITATION = r"""
# Induction (eddy-current) levitation

A time-varying magnetic source above (or moving relative to) a conductor
induces eddy currents; by Lenz's law those currents oppose the change and
the conductor is REPELLED.  This is the physics of:

- the classroom **jumping ring / Thomson ring** (Elihu Thomson, ~1887),
- **EDS** repulsion (the static-frequency analogue of maglev's
  moving-magnet EDS),
- eddy-current **metal sorting** and **non-contact braking**,
- the lift half of **EML melting** (see topic `eml_melting`).

## Mutual-inductance circuit model (the canonical, correct one)

Model the primary (driven, current I1 at angular frequency w) and the
conductor as a secondary loop with self-inductance L2 and resistance R2,
coupled by mutual inductance M(z) that depends on the separation z.

Induced secondary current (phasor):

    I2 = -j*w*M*I1 / (R2 + j*w*L2)

The instantaneous vertical force from the mutual coupling is
I1*I2*dM/dz; its TIME AVERAGE is

    <F_z> = (dM/dz) * (1/2) * Re[I1 * conj(I2)]
          = -(1/2) * |I1|^2 * (dM/dz) * w^2 * M * L2 / (R2^2 + w^2*L2^2)

Because dM/dz < 0 (coupling weakens as the bodies separate), <F_z> is
POSITIVE -> lift.  Two limits:

- **Low frequency** (w*L2 << R2):  |<F_z>| ~ w^2  (tiny -- resistive,
  the induced current is small and nearly in phase).
- **High frequency** (w*L2 >> R2):  the force SATURATES at

    |<F_z>|_max = (1/2) * |I1|^2 * |dM/dz| * M / L2

  i.e. the **perfect-conductor image limit** -- the conductor behaves
  like a mirror and the lift no longer grows with frequency.

The crossover is the **magnetic Reynolds number** Rm = w*L2/R2 ~ 1
(equivalently the conductor thickness vs skin depth delta = sqrt(2*rho/
(w*mu0*mu_r))).  Lift requires Rm >~ 1; below that the field soaks
straight through.

## Jumping-ring quick numbers

For a single aluminium ring (rho = 2.8e-8 Ohm*m) on a 50/60 Hz coil,
the ring is in the rising part of the curve -- this is why the demo
"jumps" only when the ring is COLD (low R2) and barely moves when warm
(R2 up -> force down, the same R2-in-the-denominator dependence).
Driving the coil at kHz pushes any ordinary metal ring into the
saturated image-limit regime.

## Stability

Induction lift is intrinsically NON-conservative (the field does net
work on the eddy currents -> dissipation), so it is NOT bound by
Earnshaw's theorem -- it can be passively stable, which is exactly why
EDS works (see `earnshaw_stability`).  Vertical lift is generically
stable (force falls with height); LATERAL stability depends on field
shaping -- a flat plate under a single coil is laterally unstable, a
cusped / multi-coil field is stable.

## Computing it in Radia / PEEC

1. **PEEC + SIBC** (`radia_mcp.peec`): mesh the conductor surface, drive
   the coil filaments at frequency f, solve for the surface current,
   then integrate the time-average Lorentz force density <J x B>/... over
   the conductor (see `force_computation`).  This captures M(z), L2, R2
   self-consistently -- no lumped fitting needed.
2. **Frequency sweep** confirms the w^2 -> saturation shape above; use it
   to pick a drive frequency in the saturated regime for maximum lift.
3. **TEAM Problem 28** is the validation case -- a coil set lifting an
   aluminium plate; cross-check your PEEC lift against the published
   force-vs-height curve (`radia_mcp.team_benchmark`).

## Cross-references
- `radia_mcp.peec` -- the eddy-current / SIBC solver
- topic `eddy_current_maglev` -- the MOVING-magnet (velocity-driven)
  version of the same physics
- `radia_mcp.team_benchmark` (force_motion.problem_28) -- benchmark
- topic `eml_melting` -- when the levitated conductor is a molten sphere
"""


EML_MELTING = r"""
# Electromagnetic levitation melting (EML)

EML levitates and melts a metal CONTAINERLESS: a high-frequency coil
induces eddy currents in the sample, the time-average Lorentz force
supports it against gravity, and the same eddy currents Joule-heat it
through the melting point.  No crucible -> no contamination, deep
undercooling, and clean thermophysical-property measurement of reactive
melts (Ti, Zr, Ni superalloys, semiconductors).

Foundational reference: **Okress, Wroughton, Comenetz, Brace & Kelly
(1952), J. Appl. Phys. 23, 545** -- "Electromagnetic levitation of solid
and molten metals."

## The governing tension: lift and heat are COUPLED

Both lift and heating come from the SAME induced current, so on Earth you
cannot set them independently:

    P_heat  ~ proportional to  |H|^2 * f^{1/2}      (surface Joule loss)
    F_lift  ~ proportional to  grad(|H|^2)           (field-gradient force)

To carry 1 g of dense metal you need a strong field gradient, which
dumps so much power that the sample OVERHEATS / boils.  This is the
central terrestrial-EML problem.  Mitigations:

- **Field shaping**: a CONICAL or segmented coil with a counter-wound
  upper turn makes a field MINIMUM (a "magnetic bottle") so the sample
  sits at low |H| (low heating) but a steep surrounding gradient (good
  restoring force).  Classic Fromm/Jehn coil geometry.
- **Microgravity**: with almost no weight to support, only weak
  positioning fields are needed -> lift and heat decouple.  This is why
  the **TEMPUS** facility (parabolic flights, TEXUS sounding rockets)
  and the **ISS-EML / MSL-EML** (Egry, Lohofer et al.) exist; they
  measure surface tension, viscosity, density, and specific heat of
  deeply undercooled melts.
- **Cold crucible / skull melting**: a water-cooled segmented copper
  crucible confines the melt with a thin solid "skull" of its own
  material -- a half-levitation industrial variant for reactive alloys.

## Stability

Vertical: stable at the field minimum.  Lateral and rotational:
controlled by the cusped field; molten samples also have internal MHD
STIRRING (the same Lorentz body force drives recirculation), which is
both useful (mixing, homogenisation) and a nuisance (surface
oscillation, droplet ejection at high power).

## Relation to induction heating

EML is "induction heating where the workpiece is the levitated sample".
Everything in `radia_mcp.ih` about skin depth, SIBC, surface power
density, and Karl/ESIM nonlinear surface impedance applies to the
heating side; this server's `induction_levitation` + `force_computation`
topics supply the lift side.  A full EML model couples:

    PEEC / FEM-SIBC eddy currents  ->  P_heat (radia_mcp.ih)
                                   \->  <J x B> lift (force_computation)

## Computing it in Radia

1. Coil + sphere: `radia_mcp.peec` (PEEC + SIBC) or NGSolve FEM-SIBC at
   the drive frequency (typ. 100 kHz - 1 MHz for cm samples).
2. Lift = volume integral of time-average Lorentz force density (see
   `force_computation`); compare against sample weight m*g to find the
   levitation gap and the minimum coil current.
3. Heat = SIBC surface Joule loss; cross-check sample temperature vs the
   melting point using `radia_mcp.ih` workpiece-power recipes.
4. The conical-coil field minimum is a geometry-optimisation target --
   pair with `radia_mcp.optuna` / `radia_mcp.bayesian_opt` to trade lift
   margin against overheating.

## Cross-references
- `radia_mcp.ih` -- the heating side (skin depth, SIBC, P_wp)
- topic `induction_levitation` -- the lift physics
- topic `force_computation` -- the <J x B> integral
- `radia_mcp.peec` -- coil + sample eddy-current solver
"""


MAGNETIC_BEARINGS = r"""
# Active magnetic bearings (AMB) and magnetic suspensions

An AMB holds a rotor with electromagnets whose currents are servo-
controlled from a position sensor -- the engineering workhorse of
attractive (EMS-type) levitation.  Bible: **Schweitzer & Maslen (eds.),
"Magnetic Bearings", Springer 2009**.

## The reluctance-actuator force law

A single horseshoe electromagnet of N turns, pole-face area A, air gap g,
pulling on a ferromagnetic target (neglecting iron reluctance and
fringing):

    F(i, g) = mu0 * N^2 * A * i^2 / (4 * g^2)

The force is ALWAYS attractive and is NONLINEAR in both i and g -- so a
single magnet cannot push, and the bare law is awkward for control.

## Differential (bias-current) linearization

Use two opposed magnets with a bias current i0 and a control delta ix
(top i0+ix, bottom i0-ix), nominal gap g0, rotor displacement x toward
the top:

    F_net ~ k_i * ix + k_x * x

    k_i = mu0 * N^2 * A * i0 / g0^2     (current stiffness, > 0)
    k_x = mu0 * N^2 * A * i0^2 / g0^3   (position stiffness, > 0 here
                                          meaning DESTABILIZING)

The crucial fact: **k_x is a NEGATIVE (open-loop unstable) stiffness** --
displacing the rotor toward one magnet increases that magnet's pull, so
the open-loop bearing is unstable (Earnshaw again).  A feedback law

    ix = -(kp + kd*s + ...) * x

must supply more than k_x of effective stiffness to stabilize.  This
i0-bias trick turns the quadratic actuator into a clean linear
force-current/force-displacement actuator around the operating point --
the single most important AMB design identity.

## Bearing types

| Type | Geometry | Use |
|------|----------|-----|
| **Radial heteropolar** | alternating N/S poles around rotor | general |
| **Radial homopolar** | one axial bias flux, control flux radial | low rotor loss |
| **Axial (thrust)** | two opposed ring magnets | axial location |
| **Conical** | combined radial + axial in one actuator | compact |
| **Hybrid PM-biased** | PM supplies i0 flux, coil supplies ix | low power, zero-bias-current possible |
| **Passive PM** | repelling/attracting PM rings | ONE axis only (Earnshaw) |
| **Superconducting** | HTS pinning (see `superconducting`) | passive, all axes |

A passive PM bearing can be stable in at most 2 of the 3 axes (and never
all 3) -- the remaining axis must be actively controlled or
superconducting.  This is a direct Earnshaw corollary; see
`earnshaw_stability`.

## Flywheel energy storage

AMBs enable contactless, low-loss flywheels: an evacuated rotor on
magnetic bearings has no mechanical friction, so the only spin-down loss
is residual eddy / windage.  PM-biased or HTS-pinned bearings minimize
the standby power.  Cross-link `radia_mcp.motor` (bearingless
motor) -- a bearingless machine merges the AMB and the motor windings.

## Computing it in Radia

1. **Force map**: sweep the rotor over the gap range with Radia MMM/MSC
   (iron yoke + coil), extract F(i, g) via the Maxwell stress tensor
   (`force_computation`), then fit k_i and k_x for the control model.
2. **PM-biased designs**: model the PM bias flux + coil control flux
   superposition directly (Radia handles PM + soft iron together).
3. **Negative-stiffness check**: confirm dF/dx > 0 at the operating
   point (the value you must beat with feedback gain).

## Cross-references
- `radia_mcp.electromagnet` -- DC magnet / pole-face B design
- `radia_mcp.motor` -- bearingless motor (AMB + motor in one)
- `radia_mcp.magnetic_materials` -- soft-iron B-H, PM bias data
- topic `earnshaw_stability` -- why k_x is destabilizing
- topic `force_computation` -- extracting F(i, g)
"""


SUPERCONDUCTING = r"""
# Superconducting levitation

Superconductors give the only PASSIVE, STATIC, fully-3D-stable magnetic
levitation -- the basis of HTS maglev bearings, demonstrators, and
(in transport form) SCMaglev.  Reference: **Moon, "Superconducting
Levitation", Wiley 2004**.

## Two distinct effects -- do not conflate them

| Effect | Type | Field behaviour | Stability ALONE |
|--------|------|-----------------|-----------------|
| **Meissner** | Type-I, and Type-II below Hc1 | B fully EXPELLED (perfect diamagnet, chi = -1) | UNSTABLE laterally (like balancing on a dome) |
| **Flux pinning** | Type-II, Hc1 < H < Hc2 (mixed state) | flux penetrates as quantized vortices, PINNED at defects | STABLE in all directions, including SUSPENSION (hanging below) |

Pure Meissner repulsion is a perfect-diamagnet mirror -- it lifts but a
magnet over a flat SC slides off sideways.  **Flux pinning** is what
makes the famous "frozen" levitation: vortices are pinned to material
defects, so any displacement (up, down, AND sideways) is opposed -- the
magnet is locked in 3D and can even hang UNDER the superconductor.
Engineering HTS levitation = Type-II bulk (YBCO/GdBCO) with strong
pinning, never pure Meissner.

## Field-cooled (FC) vs zero-field-cooled (ZFC)

- **ZFC**: cool the SC with no magnet present, then approach.  Vortices
  are excluded then forced in on approach -> strong repulsion, weaker
  lateral pinning.
- **FC**: cool the SC with the magnet already in place.  The vortex
  lattice freezes the EXISTING field configuration -> the magnet is
  pinned AT that position with restoring force in every direction (and
  can be suspended below).  FC gives the stiff, stable "locked" levitation
  used in bearings.

## The frozen-image model (for force computation)

Force and its HYSTERESIS are captured by a two-image construction
(Kordyuk 1998, "frozen mirror-image" model, open literature):

- a **mirror image** (opposite sign, mirror position) -> the Meissner
  shielding / repulsion, and
- a **frozen image** (same sign, fixed at the field-cooling position) ->
  the pinning / restoring + suspension force,

each scaled by the SC response.  Levitation force vs gap then shows the
characteristic OPEN HYSTERESIS LOOP (different on approach vs retract) --
a fingerprint of irreversible vortex motion, not present in Meissner-only
or PM-PM systems.

## Applications

- **HTS magnetic bearings / flywheels** -- passive, self-stable, near-zero
  standby power (Hull, Supercond. Sci. Technol. 2000, open literature).
- **HTS maglev demonstrators** -- a YBCO sled over a PM Halbach guideway
  (cross-link topic `halbach_arrays`).
- **Levitated rotors, vibration isolation, lab demos.**

## Computing it in Radia

Radia computes the PM / coil field; the SC response is added via the
frozen-image model above (mirror + frozen image sources mirrored about
the SC surface), then force via the stress tensor (`force_computation`).
For full critical-state (Bean model) vortex dynamics, couple to an
external SC solver -- Radia supplies the applied field map.

## Cross-references
- topics `scmaglev_eds`, `halbach_arrays` -- SC maglev transport
- topic `diamagnetic` -- Meissner is the chi = -1 diamagnetic limit
- topic `force_computation` -- stress-tensor force from the image sources
- topic `earnshaw_stability` -- pinning is loophole #5
"""


DIAMAGNETIC = r"""
# Diamagnetic levitation

A diamagnet (susceptibility chi < 0: graphite, bismuth, water, most
organic matter) is pushed toward WEAKER field.  Because a static field in
free space can have a |B| MINIMUM (but never a maximum -- Earnshaw), a
diamagnet can sit in a passively stable, fully static trap.  This is the
ONE loophole that needs neither motion, feedback, spin, nor
superconductivity.  Reference: **Berry & Geim (1997), Eur. J. Phys. 18,
307**; **Brandt, Science 243, 349 (1989)**.

## Force density

For a linear diamagnet (|chi| << 1) of mass density rho, the body-force
density is

    f = (chi / (2*mu0)) * grad(|B|^2)      [N/m^3]

(using grad(B^2/2) = (B.grad)B for a curl-free field).  Since chi < 0 the
force points DOWN the |B|^2 gradient -- toward the field minimum.  Lift
against gravity requires

    (chi / (2*mu0*rho*g)) * d(|B|^2)/dz  >=  1

i.e. the field must supply

    B * dB/dz  >=  mu0 * rho * g / |chi|      [T^2 / m]

## The famous numbers

| Material | chi (volume, SI) | B*dB/dz to levitate |
|----------|------------------|----------------------|
| Water / living tissue | ~ -9e-6 | ~ 1400 T^2/m (the Geim "flying frog", ~16 T bore) |
| Bismuth | ~ -1.6e-4 | ~ 80 T^2/m |
| Pyrolytic graphite | ~ -4e-4 (anisotropic, strongest perp) | ~ 30 T^2/m (levitates over ordinary NdFeB!) |

Pyrolytic graphite over a checkerboard NdFeB array is the classic
benchtop demo: the PM array makes a periodic |B| with many local minima,
each a stable pocket, and the strong perpendicular diamagnetism of
graphite needs no superconductor and no high field.

## Stability

Diamagnetic equilibrium is genuinely stable in all three axes at a true
|B|^2 minimum -- the diamagnet sees an effective potential
U = -(chi/(2*mu0)) |B|^2 * V which, with chi < 0, has a MINIMUM where
|B|^2 does.  This is why Braunbek (1939) proved that a material with a
region of chi < 0 is NECESSARY and SUFFICIENT for passive static
levitation.

## Computing it in Radia

1. Compute |B|^2 from the PM / coil array (Radia MMM/MSC, `rad.Fld`).
2. The trap is where grad(|B|^2) = 0 with a positive-definite Hessian of
   |B|^2 (a minimum); the levitation height is where the upward
   f_z balances rho*g.
3. For the diamagnetic object's own (tiny) back-reaction, register it as
   a weak linear material -- usually negligible for the field, dominant
   only for its own force.

## Cross-references
- topic `earnshaw_stability` -- diamagnets are loophole #1 (and the
  Braunbek necessity result)
- topic `superconducting` -- the chi = -1 (perfect-diamagnet) extreme
- `radia_mcp.magnetic_materials` -- susceptibility data
"""


EARNSHAW_STABILITY = r"""
# Earnshaw's theorem and its loopholes

Earnshaw's theorem (1842) is the master key to all of magnetic
levitation: it tells you WHY naive "just use a strong magnet" never works
and forces every real scheme into one of five escape routes.

## The theorem

A collection of point charges (or fixed magnetic poles, or paramagnetic
bodies) interacting by inverse-square (1/r^2) forces cannot be held in a
STABLE static equilibrium by those forces alone.

Proof sketch: in free space the force field is divergence-free
(div F = 0, because the potential satisfies Laplace's equation,
del^2 U = 0).  A stable equilibrium needs a potential MINIMUM, i.e.
del^2 U > 0 there -- impossible when del^2 U = 0.  Every equilibrium is a
saddle: stable along some axes, unstable along at least one.

Consequence: **no arrangement of permanent magnets or fixed
electromagnets in free space can levitate another magnet stably.**  A PM
over a PM always slips out sideways or flips.

## The five loopholes (each violates a proof assumption)

| # | Loophole | Which assumption it breaks | Example |
|---|----------|----------------------------|---------|
| 1 | **Diamagnetism (chi < 0)** | material is NOT paramagnetic; U gains a real minimum | graphite/water levitation; topic `diamagnetic` |
| 2 | **Induced / eddy currents** | force is NOT conservative / static (the field does work) | EDS, jumping ring; topic `induction_levitation` |
| 3 | **Active feedback** | sources are NOT fixed (current servo-controlled) | EMS, AMB; topic `magnetic_bearings` |
| 4 | **Spin / gyroscopic** | not truly static (precession stabilizes) | Levitron spinning top (Berry & Geim 1997) |
| 5 | **Flux pinning (Type-II SC)** | not a simple 1/r^2 material; vortices pin | HTS bulk; topic `superconducting` |

**Braunbek's refinement (1939)**: among purely STATIC, PASSIVE schemes
(no motion, no feedback, no spin), a material with relative permeability
mu_r < 1 somewhere -- i.e. a DIAMAGNET (including a superconductor as the
mu_r -> 0 extreme) -- is both NECESSARY and SUFFICIENT.  So static passive
levitation == diamagnetic levitation; everything else needs dynamics or
control.

## The stability stiffness matrix (practical test)

Near an equilibrium x0, expand the levitation force F(x) = -K (x - x0).
The 3x3 stiffness K = -dF/dx must be POSITIVE DEFINITE (all eigenvalues
> 0) for stability.  Earnshaw says: for fixed-source magnetostatics,

    trace(K) = -div F = 0      ->   eigenvalues sum to zero
                               ->   at least one is <= 0  ->  unstable.

So the engineering job is to add a mechanism (a loophole) that makes
trace(K) > 0:
- diamagnet / SC: adds a positive-definite block,
- eddy current: adds velocity-dependent (damping + non-conservative)
  terms that fall outside the static trace argument,
- feedback: actively injects positive stiffness (and damping),
- spin: converts an unstable axis into precessional (gyroscopic) stability.

For a passive PM BEARING this is why you get stability in at most 2 axes
(two positive eigenvalues) and instability in the third
(trace = 0 forces the third negative) -- the third axis MUST be actively
controlled or superconducting (see `magnetic_bearings`).

## Cross-references
- topic `diamagnetic` -- loophole #1 + Braunbek necessity
- topic `induction_levitation` -- loophole #2
- topic `magnetic_bearings` -- loophole #3 + the 2-of-3-axes corollary
- topic `superconducting` -- loophole #5
- topic `pm_maglev_zero_power` -- the passive-PM-bearing application
  (axial-stable / radially-unstable; the 2-of-3-axes corollary in practice)
- topics `eddy_current_maglev` / `scmaglev_eds` -- EDS/EMS maglev
  transport apply the same loopholes
"""


FORCE_COMPUTATION = r"""
# Computing the levitation force in Radia / PEEC / NGSolve

Every topic in this server ultimately needs ONE number: the vertical
force on the levitated body.  There are three rigorous routes; pick by
what your solver gives you.

## 1. Maxwell stress tensor (field-based, the Radia default)

Surround the body with a closed surface S in the air/vacuum and integrate

    F = oint_S  T . dA ,   T_ij = (1/mu0) [ B_i B_j - (1/2) delta_ij |B|^2 ]

- Works for PM, soft-iron, and (with the time-average B) eddy-current
  bodies.
- In Radia: evaluate B on a box enclosing the body (`rad.Fld`) and sum
  the stress over the faces; keep the surface in air, >~ 1 element off
  the body (per the NGSolve evaluation-distance policy).
- Vector difference of fields, not |B| difference, throughout.

## 2. Virtual work / energy gradient (robust for gaps)

    F_z = - dW/dz   (constant flux)     or    F_z = + dW_co/dz (constant current)

- Sweep the body over the gap, get W(z) (co-energy for current-driven),
  differentiate.  Excellent for AMB force maps and PM bearings where you
  want F(i, g) anyway.
- Numerically gentler than the stress tensor near saturated iron.

## 3. Time-average Lorentz force (eddy-current / EML, PEEC)

For induced currents J(t) in a conductor (induction levitation, EDS,
EML), the cycle-average body force is

    <f> = (1/2) Re[ J x conj(B) ]        [N/m^3]

    F_z = integral_V  <f>_z  dV

- This is the natural output of a PEEC + SIBC or FEM-SIBC eddy-current
  solve (`radia_mcp.peec`, `radia_mcp.ih`): you already have the complex
  surface/volume current J and field B at the drive frequency.
- The same <J x conj(B)> field that gives lift also gives the MHD stir in
  an EML melt -- integrate the vertical component for lift, the full
  vector for the flow drive.

## Consistency check

Routes 1 and 3 must agree: the Maxwell stress on an enclosing surface
equals the volume integral of the Lorentz density inside (divergence
theorem on T).  Computing BOTH and comparing is the standard
verify-first check for an eddy-current lift model -- a mismatch flags a
surface-too-close or skin-depth-underresolved error before you trust the
number.

## Recipe pointers

- Stress-tensor + virtual-work rigor: `radia_mcp.differential_forms`
  (em_force_recipe, em_force_extras) -- the differential-forms derivation
  and the sign/gauge pitfalls.
- Eddy-current J, B at frequency: `radia_mcp.peec` (PEEC+SIBC),
  `radia_mcp.ih` (FEM-SIBC workpiece).
- Benchmark the result: topic `benchmarks` (TEAM 28).

## Cross-references
- `radia_mcp.differential_forms` -- rigorous EM force (stress vs Lorentz)
- `radia_mcp.peec` / `radia_mcp.ih` -- where <J x B> comes from
- topic `benchmarks` -- numbers to check against
"""


BENCHMARKS = r"""
# Levitation-force validation benchmarks

Trust a levitation force only after it matches a known case.  Use the
FIRST-QUESTION habit: "given THIS analysis, which closed form / benchmark
is the trusted reference?" (cf. `radia_mcp.radia_ngsolve` analytical
formulas).

## TEAM Workshop Problem 28 -- electrodynamic levitation device

The canonical eddy-current LIFT benchmark: a 3-axis axisymmetric coil
set driven at 50 Hz (10 A rms) lifts an aluminium plate (50 x 50 x 5 mm)
to an equilibrium height of ~5 mm; the published quantity is the
levitation force vs plate height (the plate floats where lift = weight).
(Specs per `radia_mcp.team_benchmark.force_motion` problem_28.)

- Use it to validate a PEEC / FEM-SIBC induction-levitation model
  (topic `induction_levitation`).
- Cross-link `radia_mcp.team_benchmark` (force_motion.problem_28) for the
  geometry, material, and reference force-height curve.
- Pass criterion: levitation height and force within the published band;
  check BOTH the Maxwell-stress and the <J x B> routes agree
  (`force_computation`).

## Jumping-ring analytic

The mutual-inductance circuit model in `induction_levitation` IS an
analytic benchmark: for a single ring you can write L2, R2, M(z) in
closed form (or measure them) and compare the predicted
force-vs-frequency curve (w^2 rise -> saturation) against a PEEC solve of
the same ring.  A clean way to validate the frequency dependence before
trusting a complex EML geometry.

## Sphere / disk over a coil

A conducting sphere (or disk) above an AC coil has semi-analytic eddy-
current lift (image-dipole / series solutions, open literature) -- a good
intermediate check between the lumped ring and a full EML sample.

## EML lift coefficient

For EML, the standard non-dimensional check is the lift efficiency vs the
sample's normalized radius x = a/delta (a = sample radius, delta = skin
depth): published lift- and power-absorption functions G(x), H(x)
(Fromm/Jehn, Lohofer, open literature) give the lift force and heating
for a sphere in a known gradient.  Match your PEEC sphere against G(x),
H(x) before committing to a coil-geometry optimisation.

## Diamagnetic / SC checks

- Diamagnetic: the B*dB/dz threshold table in topic `diamagnetic` is a
  direct check (does your computed grad(|B|^2) cross rho*g/|chi| at the
  observed height?).
- Superconducting: the force-vs-gap HYSTERESIS LOOP shape (topic
  `superconducting`) is the qualitative signature; the frozen-image model
  gives the quantitative curve.

## Cross-references
- `radia_mcp.team_benchmark` -- TEAM 28 and the force-motion suite
- topic `force_computation` -- the two routes that must agree
- topic `induction_levitation` -- the analytic ring model
- `radia_mcp.radia_ngsolve` -- analytical-formula reference layer habit
"""


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch magnetic-levitation topics.

    Topics:
        overview                  - Magnetic levitation landscape + lab research (DEFAULT)
        radia_iem_fem             - Radia IEM <-> reduced-potential FEM weak coupling (Yano)
        cln_mor_control           - Cauer Ladder Network MOR for control-coupled maglev (Yano)
        physical_tensor_rom       - Physical polarizability tensor alpha(s) as a passive LTI (AAA+NNLS)
        pm_maglev_zero_power      - Passive PM levitation, Maxwell-Earnshaw
        eddy_current_maglev       - Eddy-current EDS, Kansai 2D model, Arago
        sumitomo_heavy_industrial - JP 7-327337 PM bearing + JP 2007-215264 planar mover
        kansai_research           - Saiki/Fujii magnetic-wheel lineage
        scmaglev_eds              - SCMaglev (Chuo Shinkansen) -- SC-EDS levitation
        halbach_arrays            - Halbach + Inductrack
        induction_levitation      - Eddy-current (AC) lift: jumping/Thomson ring
        eml_melting               - Electromagnetic levitation melting (ties to IH)
        magnetic_bearings         - AMB force-current-displacement, neg. stiffness
        superconducting           - Meissner vs flux pinning, HTS bulk, frozen image
        diamagnetic               - grad(B^2) levitation, graphite/water-frog
        earnshaw_stability        - Earnshaw's theorem + its 5 loopholes
        force_computation         - Maxwell stress / virtual work / time-avg J x B
        benchmarks                - TEAM 28, jumping-ring analytic, EML lift coeff
        all                       - Everything
    """
    topic = topic.lower().strip()

    if topic in ("overview", "intro", ""):
        return OVERVIEW
    if topic in ("radia_iem_fem", "iem_fem", "iem", "weak_coupling",
                 "reduced_potential", "moving_magnet", "iem-fem"):
        return RADIA_IEM_FEM
    if topic in ("cln_mor_control", "cln", "cauer", "cauer_ladder", "mor",
                 "model_order_reduction", "multiport_cln", "control_coupled"):
        return CLN_MOR_CONTROL
    if topic in ("physical_tensor_rom", "tensor_rom", "polarizability_rom",
                 "alpha_rom", "rom_fit", "aaa_nnls", "foster_rom",
                 "stoll_rom", "physical_tensor"):
        return PHYSICAL_TENSOR_ROM
    if topic in ("pm_maglev_zero_power", "pm", "passive_pm",
                 "zero_power", "pm_bearing"):
        return PM_MAGLEV_ZERO_POWER
    if topic in ("eddy_current_maglev", "eddy", "eds", "arago",
                 "magnetic_wheel", "magnetic_wheels"):
        return EDDY_CURRENT_MAGLEV
    if topic in ("sumitomo_heavy_industrial", "sumitomo", "shi",
                 "planar_mover", "patent"):
        return SUMITOMO_HEAVY_INDUSTRIAL
    if topic in ("kansai_research", "kansai", "saiki", "fujii"):
        return KANSAI_RESEARCH
    if topic in ("scmaglev_eds", "scmaglev", "sc_eds", "chuo_shinkansen",
                 "chuo", "jr_maglev"):
        return SCMAGLEV_EDS
    if topic in ("halbach_arrays", "halbach", "inductrack"):
        return HALBACH_ARRAYS
    # -- levitation FORCE physics --
    if topic in ("induction_levitation", "induction", "eddy_lift",
                 "jumping_ring", "thomson_ring", "thomson"):
        return INDUCTION_LEVITATION
    if topic in ("eml_melting", "eml", "melting", "levitation_melting",
                 "containerless", "tempus", "cold_crucible", "skull_melting"):
        return EML_MELTING
    if topic in ("magnetic_bearings", "amb", "bearing", "bearings",
                 "active_magnetic_bearing", "flywheel"):
        return MAGNETIC_BEARINGS
    if topic in ("superconducting", "sc", "meissner", "flux_pinning",
                 "pinning", "hts", "ybco"):
        return SUPERCONDUCTING
    if topic in ("diamagnetic", "diamagnetism", "graphite", "pyrolytic",
                 "frog", "geim"):
        return DIAMAGNETIC
    if topic in ("earnshaw_stability", "earnshaw", "stability", "loophole",
                 "loopholes", "braunbek"):
        return EARNSHAW_STABILITY
    if topic in ("force_computation", "force", "stress_tensor",
                 "maxwell_stress", "virtual_work", "lorentz", "compute"):
        return FORCE_COMPUTATION
    if topic in ("benchmarks", "benchmark", "validation", "team28",
                 "team_28", "team", "verify"):
        return BENCHMARKS
    if topic == "all":
        return "\n\n".join([
            OVERVIEW,
            RADIA_IEM_FEM,
            CLN_MOR_CONTROL,
            PHYSICAL_TENSOR_ROM,
            PM_MAGLEV_ZERO_POWER,
            EDDY_CURRENT_MAGLEV,
            SUMITOMO_HEAVY_INDUSTRIAL,
            KANSAI_RESEARCH,
            SCMAGLEV_EDS,
            HALBACH_ARRAYS,
            INDUCTION_LEVITATION,
            EML_MELTING,
            MAGNETIC_BEARINGS,
            SUPERCONDUCTING,
            DIAMAGNETIC,
            EARNSHAW_STABILITY,
            FORCE_COMPUTATION,
            BENCHMARKS,
        ])

    return (
        f"Unknown topic '{topic}'. Available: " + ", ".join(TOPICS.keys()) + "."
    )
