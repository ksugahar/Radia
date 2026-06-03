"""Magnetic levitation knowledge (incl. the lab's Radia-IEM + CLN maglev research).

Scope: MAGNETIC LEVITATION for transport / suspension / control.  The
former linear-drive (LIM/LSM, end-effect) material was removed -- this
server is maglev-only; for the levitation-FORCE physics (induction lift,
EML melting, magnetic bearings, superconducting, diamagnetic, Earnshaw)
see the sibling `radia_mcp.levitation`.

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

Distilled from W:/.../99_アプリケーション/07_磁気浮上/
plus open literature for standard maglev topics that the
lab does not have direct PDFs for (those sections are clearly marked
"open literature, not lab PDF" so users know to cross-check).

Primary lab PDFs read for this revision:

  - 07_磁気浮上/research/渦電流で磁気浮上.pdf
      (Murata, Yamagata HS) - 1-page Arago-disk demo, NdFeB array
      over rotating Al disk, eddy currents lift small NdFeB pendulum.

  - 07_磁気浮上/04_関西大学/永久磁石を用いた磁気浮上 原本.pdf
      Saiki, Ino, Inada (Kansai Univ., 2021): "Optimization on Pole
      Pitch of Magnetic Wheels and Thickness of Metal Plate for
      Floating and Propulsion System using Permanent Magnets". 2D
      analytical EDS model with skin depth delta_s, secondary B/J
      fields, F_L = J x B, F_L/P scaling proportional to sqrt(pole
      pitch tau). Builds on Fujii (1995, 1997) IEEJ "Revolving
      Permanent Magnet Type Magnet Wheel" lineage.

  - 07_磁気浮上/01_特許/磁気浮上_住重特許1.pdf
      JP 公開特許 7-327337 (Sumitomo Heavy Industries): PM-repulsion
      axial magnetic bearing with eddy-current brake ring (Al
      annulus) for overspeed prevention. Glass-wool thermal
      insulator between magnet and ring.

  - 07_磁気浮上/01_特許/磁気浮上_住重特許2.pdf
      JP 公開特許 2007-215264 (Sumitomo Heavy Industries): planar
      mover with rotating-PM "drive head" (9 reluctance paths),
      air-bearing levitation, eddy currents in conductive base
      plate provide thrust. 4-head omnidirectional version with
      wireless DC-motor control. Stepper / LSM variants.

NOT read on this pass (size budget):
  - 09_リニアドライブ/*.pdf  (each 10-44 MB)
  - 03_教科書/磁気浮上と磁気軸受.pdf  (155 MB textbook)
  - research/RocketScience_Dec2014.pdf  and other > 2 MB items

Sections marked "(open literature)" are textbook-level summaries
included so the knowledge file is self-contained; cross-check
against the unread lab PDFs when those become accessible.
"""


# Authoritative topic enum for the dispatcher tool (wired into
# `maglev_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "Magnetic levitation landscape (EMS/EDS/PM/SC/Halbach) + the lab's Radia-IEM / CLN maglev research line",
    "radia_iem_fem": "Radia IEM (MMM/MSC) <-> reduced-potential FEM weak coupling for moving-magnet eddy-current levitation force; A-phi/T-Omega; no re-mesh on magnet motion (Yano bachelor, lab research)",
    "cln_mor_control": "Cauer Ladder Network (CLN) model-order reduction for real-time control-coupled maglev: ~1/500 speedup, multiport matrix-CLN, 3D gauge A-phi/T-Omega/A-T, TEAM 28 (Yano master, lab research)",
    "pm_maglev_zero_power": "Zero-power passive PM levitation: Maxwell-Earnshaw constraint, axial PM bearings, halbach diamagnetism",
    "eddy_current_maglev": "Eddy-current EDS: Arago-disk physics, magnetic wheels (Fujii/Kansai 2D model), Inductrack",
    "sumitomo_heavy_industrial": "Sumitomo Heavy patents: PM axial bearing with brake ring (JP 7-327337); planar eddy-current mover (JP 2007-215264)",
    "kansai_research": "Kansai University Saiki/Fujii magnetic-wheel lineage: skin depth, F_L/P proportional to sqrt(pole pitch)",
    "scmaglev_eds": "Superconducting EDS (SCMaglev / Chuo Shinkansen): null-flux figure-8 levitation + guidance",
    "halbach_arrays": "Halbach array maglev / Inductrack: passive track coils, threshold-velocity self-stability",
    "all": "Concatenation of every topic above",
}


OVERVIEW = r"""
# Magnetic levitation landscape

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

These are validated on the standard eddy-current levitation benchmarks
(TEAM Problem 7, TEAM Problem 28 = the 1-axis electrodynamic levitation
device).

## Lab focus

- **Radia-IEM + FEM weak coupling** for moving-magnet eddy-current force
  (Yano & Sugahara) -- topic `radia_iem_fem`
- **CLN model-order reduction** for control-coupled maglev (Yano) --
  topic `cln_mor_control`
- **Magnetic-wheel EDS** -- Kansai-Univ collaboration (Saiki et al. 2021)
- **PM axial bearings** -- Sumitomo Heavy lineage (heavy-industrial)
- **Bearingless motor + WPT** (cross-link `radia_mcp.motor`, `radia_mcp.wpt`)

## Sibling server: levitation-force physics

For the levitation FORCE physics itself -- induction (eddy-current) lift,
EML melting, magnetic bearings (AMB), superconducting (Meissner /
pinning), diamagnetic, Earnshaw's theorem, and how to compute the force
(Maxwell stress / virtual work / time-average Lorentz) -- see
`radia_mcp.levitation`.  Split by intent: a maglev SYSTEM (this server)
vs. the levitation FORCE (levitation server).

## Cross-references

- `radia_mcp.levitation` -- levitation-force physics + bearings + EML
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
Also validated against the eddy-current levitation benchmarks TEAM
Problem 7 and TEAM Problem 28.

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
2017) that produces orthogonal current/field modes directly from the FEM
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
- Refs: Kameari, Ebrahimi, Sugahara, Shindo & Matsuo 2017 (CLN
  representation of eddy-current fields); Sugahara, Tanimoto, Takahashi &
  Matsuo 2023 (CLN with constant basis for conductor movement); Matsuo,
  Fujiwara, Kuriyama & Shindo 2019 (multiport matrix CLN); Tanimoto,
  Yano, Sugahara & Nagamine 2025 (3D CLN gauge comparison).
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

## Axial PM bearing -- Sumitomo Heavy lineage (★ lab PDF)

Reference: JP 公開特許 7-327337 "永久磁石式磁気浮上回転装置"
(Sumitomo Heavy Industries, 1995).

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

## Lab demonstration (★ lab PDF)

Reference: 渦電流で磁気浮上 (Murata Kenji, Yamagata HS).

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

## Magnetic-wheel EDS -- Kansai University lineage (★ lab PDF)

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
(JP 公開特許 7-327337, 1995)

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
(JP 公開特許 2007-215264, 2007)

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


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch magnetic-levitation topics.

    Topics:
        overview                  - Magnetic levitation landscape + lab research (DEFAULT)
        radia_iem_fem             - Radia IEM <-> reduced-potential FEM weak coupling (Yano)
        cln_mor_control           - Cauer Ladder Network MOR for control-coupled maglev (Yano)
        pm_maglev_zero_power      - Passive PM levitation, Maxwell-Earnshaw
        eddy_current_maglev       - Eddy-current EDS, Kansai 2D model, Arago
        sumitomo_heavy_industrial - JP 7-327337 PM bearing + JP 2007-215264 planar mover
        kansai_research           - Saiki/Fujii magnetic-wheel lineage
        scmaglev_eds              - SCMaglev (Chuo Shinkansen) -- SC-EDS levitation
        halbach_arrays            - Halbach + Inductrack
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
    if topic in ("pm_maglev_zero_power", "pm", "passive_pm", "earnshaw",
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
    if topic == "all":
        return "\n\n".join([
            OVERVIEW,
            RADIA_IEM_FEM,
            CLN_MOR_CONTROL,
            PM_MAGLEV_ZERO_POWER,
            EDDY_CURRENT_MAGLEV,
            SUMITOMO_HEAVY_INDUSTRIAL,
            KANSAI_RESEARCH,
            SCMAGLEV_EDS,
            HALBACH_ARRAYS,
        ])

    return (
        f"Unknown topic '{topic}'. Available: " + ", ".join(TOPICS.keys()) + "."
    )
