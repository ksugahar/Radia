"""Magnetic levitation FORCE-physics + stationary/industrial levitation.

Scope split vs `radia_mcp.maglev_linear`:

  - `radia_mcp.maglev_linear`  -> TRANSPORT / PROPULSION (EMS/EDS trains,
    SCMaglev, Halbach/Inductrack tracks, LIM/LSM linear motors, magnetic
    wheels).  Vehicle-scale, motion along a track.

  - `radia_mcp.levitation` (THIS server) -> the LEVITATION-FORCE PHYSICS
    and the STATIONARY / industrial / laboratory uses that Radia, PEEC and
    NGSolve are the right tools to COMPUTE:
        * induction (eddy-current) levitation  -- AC source over a
          conductor; the jumping-ring / Thomson-ring lift
        * electromagnetic levitation MELTING (EML) -- containerless
          molten-metal processing (direct tie to radia_mcp.ih)
        * active magnetic bearings (AMB) -- force-current-displacement,
          negative stiffness, flywheels
        * superconducting levitation -- Meissner vs flux-pinning, HTS bulk
        * diamagnetic levitation -- pyrolytic graphite, water/frog
        * Earnshaw's theorem and its loopholes (the WHY of all the above)
        * how to actually COMPUTE the levitation force in Radia / PEEC /
          NGSolve (Maxwell stress, virtual work, time-average Lorentz)

The two servers cross-link freely; pick by intent (move a vehicle ->
maglev_linear; suspend / melt / bear / compute a force -> levitation).

References are standard open literature unless tagged "(lab PDF)".
Equations are ASCII (cp932 console policy).  Key references:

  - Earnshaw, S. (1842) Trans. Camb. Phil. Soc. 7, 97 -- the theorem.
  - Braunbek, W. (1939) Z. Phys. 112, 753 -- diamagnets are necessary
    AND sufficient for passive static levitation.
  - Brandt, E.H. (1989) "Levitation in Physics", Science 243, 349 --
    the canonical review of every levitation mechanism.
  - Berry, M.V. & Geim, A.K. (1997) "Of flying frogs and levitrons",
    Eur. J. Phys. 18, 307 -- diamagnetic / Levitron analysis.
  - Okress, E.C. et al. (1952) "Electromagnetic levitation of solid and
    molten metals", J. Appl. Phys. 23, 545 -- the EML foundation.
  - Moon, F.C. (2004) "Superconducting Levitation", Wiley.
  - Schweitzer, G. & Maslen, E.H. eds. (2009) "Magnetic Bearings",
    Springer -- AMB design bible.
  - Jayawant, B.V. (1981) "Electromagnetic Levitation and Suspension
    Techniques", Edward Arnold.
  - TEAM Workshop Problem 28 -- "An Electrodynamic Levitation Device"
    (eddy-current lift validation benchmark).
"""


# Authoritative topic enum for the dispatcher tool (wired into
# `levitation_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "Levitation-force taxonomy; where Radia/PEEC/NGSolve fit; scope split vs maglev_linear",
    "induction_levitation": "Eddy-current (AC) lift: jumping/Thomson ring, mutual-inductance circuit model, force vs frequency, perfect-conductor image limit",
    "eml_melting": "Electromagnetic levitation MELTING (EML): containerless molten metal, lift-vs-heat coupling, conical coils, microgravity TEMPUS/ISS-EML (ties to radia_mcp.ih)",
    "magnetic_bearings": "Active magnetic bearings (AMB): F = k_i*i + k_x*x linearization, negative position stiffness, bias current, radial/axial, flywheels; passive PM + SC bearings",
    "superconducting": "Superconducting levitation: Meissner (unstable) vs flux pinning (stable 3D), HTS bulk YBCO, field-cooled vs ZFC, frozen-image model, force hysteresis",
    "diamagnetic": "Diamagnetic levitation: f = (chi/2mu0)*grad(B^2), pyrolytic graphite over PM array, Geim water/frog high-field, the one Earnshaw loophole that is passive + static",
    "earnshaw_stability": "Earnshaw's theorem + its 5 loopholes (diamagnet/eddy/feedback/spin/SC-pinning); Braunbek limit; stability stiffness matrix",
    "force_computation": "Compute levitation force in Radia/PEEC/NGSolve: Maxwell stress tensor, virtual work dW/dz, time-average <J x B>; recipes + cross-links",
    "benchmarks": "Validation references: TEAM 28 electrodynamic levitation, jumping-ring analytic, sphere-over-coil, EML lift coefficient",
    "all": "Concatenation of every topic above",
}


OVERVIEW = r"""
# Magnetic levitation: the force-physics landscape

This server is about the FORCE that lifts, and the stationary /
industrial / laboratory uses of it.  For moving a vehicle along a track
(EMS/EDS trains, SCMaglev, Inductrack, LIM/LSM) see
`radia_mcp.maglev_linear` instead.

## The five levitation mechanisms (Brandt 1989 taxonomy)

| Mechanism | Force | Static? | Passively stable? | Radia/PEEC computes it? |
|-----------|-------|---------|-------------------|--------------------------|
| **Ferromagnetic attraction (EMS)** | attractive, controlled current | yes | NO (needs feedback) | MMM/MSC for the magnetic circuit; AMB force map |
| **Induced eddy current (EDS)** | repulsive, AC or relative motion | no (time-varying) | yes (above threshold) | **PEEC + eddy-current solver** (this server's core) |
| **Diamagnetic** | repulsive, toward |B| minimum | yes | yes | Radia field + grad(B^2) |
| **Superconducting (Meissner + pinning)** | repulsive + pinned | yes | yes (Type-II pinning) | field map + frozen-image model |
| **Spin / gyroscopic (Levitron)** | PM repulsion + precession | quasi-static | yes (narrow window) | Radia PM field + dynamics |

Hybrid / actively-controlled magnetic BEARINGS (AMB) sit on top of the
EMS row: a permanent-magnet bias plus feedback-controlled current.

## Why a separate server from maglev_linear

`maglev_linear` answers "how do I move a vehicle".  This server answers
four different questions the lab actually computes:

1. **"How much lift does my AC coil exert on this conductor?"**
   -> induction_levitation (the jumping ring, EDS lift, sorting magnets)
2. **"Can I melt this metal without a crucible?"**
   -> eml_melting (containerless processing; same physics as IH)
3. **"What is the force-current-displacement law of my bearing?"**
   -> magnetic_bearings (AMB, flywheels)
4. **"Will it be stable, and if not, which loophole do I use?"**
   -> earnshaw_stability

## Where Radia / PEEC / NGSolve fit

| Need | Tool | Notes |
|------|------|-------|
| PM / soft-iron static field + force | **Radia MMM/MSC** | `rad.Fld`, force via stress tensor |
| AC eddy-current lift on a conductor | **Radia PEEC + SIBC** | time-average Lorentz force; `radia_mcp.peec` |
| EML melt + lift coupling | **Radia PEEC / NGSolve FEM-SIBC** | `radia_mcp.ih` for the heating side |
| Diamagnetic / SC force | **Radia field + post** | grad(B^2); frozen-image for SC |
| Bearing force map sweep | **Radia MMM/MSC over gap** | parametric; AMB linearization |

## Cross-references

- `radia_mcp.maglev_linear` -- TRANSPORT (EMS/EDS trains, LIM/LSM, wheels)
- `radia_mcp.ih` -- induction heating (EML is levitation + IH together)
- `radia_mcp.peec` -- the eddy-current / SIBC engine for induction lift
- `radia_mcp.electromagnet` -- DC magnet / pole-face B for EMS & AMB
- `radia_mcp.differential_forms` (em_force_recipe) -- rigorous EM force
- `radia_mcp.team_benchmark` (force_motion.problem_28) -- EDS lift benchmark
- `radia_mcp.magnetic_materials` -- PM / soft-iron constitutive data
"""


INDUCTION_LEVITATION = r"""
# Induction (eddy-current) levitation

A time-varying magnetic source above (or moving relative to) a conductor
induces eddy currents; by Lenz's law those currents oppose the change and
the conductor is REPELLED.  This is the physics of:

- the classroom **jumping ring / Thomson ring** (Elihu Thomson, ~1887),
- **EDS** repulsion (the static-frequency analogue of maglev_linear's
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
- `radia_mcp.maglev_linear` (eddy_current_maglev) -- the MOVING-magnet
  (velocity-driven) version of the same physics
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
the standby power.  Cross-link `radia_mcp.maglev_linear` (bearingless
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
- `radia_mcp.maglev_linear` -- bearingless motor (AMB + motor in one)
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
  (cross-link `radia_mcp.maglev_linear` halbach_arrays).
- **Levitated rotors, vibration isolation, lab demos.**

## Computing it in Radia

Radia computes the PM / coil field; the SC response is added via the
frozen-image model above (mirror + frozen image sources mirrored about
the SC surface), then force via the stress tensor (`force_computation`).
For full critical-state (Bean model) vortex dynamics, couple to an
external SC solver -- Radia supplies the applied field map.

## Cross-references
- `radia_mcp.maglev_linear` (scmaglev_eds, halbach_arrays) -- SC transport
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
- `radia_mcp.maglev_linear` -- EDS/EMS transport apply the same loopholes
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

The canonical eddy-current LIFT benchmark: two concentric coils driven at
50 Hz lift an aluminium plate; the published quantity is the levitation
force vs plate height (the plate floats where lift = weight).

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
    """Dispatch levitation-force topics.

    Topics:
        overview              - Levitation-force taxonomy + Radia fit (DEFAULT)
        induction_levitation  - Eddy-current lift: jumping ring, circuit model
        eml_melting           - Electromagnetic levitation melting (ties to IH)
        magnetic_bearings     - AMB force-current-displacement, negative stiffness
        superconducting       - Meissner vs flux pinning, HTS bulk, frozen image
        diamagnetic           - grad(B^2) levitation, graphite/water
        earnshaw_stability    - Earnshaw + the 5 loopholes, stiffness matrix
        force_computation     - Maxwell stress / virtual work / time-average J x B
        benchmarks            - TEAM 28, jumping-ring analytic, EML lift coeff
        all                   - Everything
    """
    topic = topic.lower().strip()

    if topic in ("overview", "intro", "", "landscape", "taxonomy"):
        return OVERVIEW
    if topic in ("induction_levitation", "induction", "eddy", "eddy_current",
                 "jumping_ring", "thomson_ring", "thomson", "eds_static"):
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
    if topic in ("all", "everything", "full"):
        return "\n\n".join([
            OVERVIEW,
            INDUCTION_LEVITATION,
            EML_MELTING,
            MAGNETIC_BEARINGS,
            SUPERCONDUCTING,
            DIAMAGNETIC,
            EARNSHAW_STABILITY,
            FORCE_COMPUTATION,
            BENCHMARKS,
        ])

    valid = ", ".join(sorted(TOPICS.keys()))
    return (
        f"Unknown topic: '{topic}'.\n\n"
        f"Valid topics: {valid}.\n\n"
        f"Call levitation(topic='overview') for the landscape, or "
        f"levitation(topic='all') for everything."
    )
