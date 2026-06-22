"""
Accelerator magnet design — end-pole chamfer analytical theory,
multipole analysis, rotating coil measurement, Radia case studies.
"""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "end_pole": "Analytical chamfer design (Delferriere)",
    "kolkata": "Radia + TOSCA validation case study",
    "rotating_coil": "Multipole measurement + field reconstruction",
    "isochronous_endpack": (
        "Radial field index design: scaling (k=const) vs isochronous "
        "(rising k_iso), and the saturated nonlinear end-pack reshape"
    ),
    "foliate_perturb": (
        "When can a magnet be designed as a 2-D body + end perturbation? "
        "Leaf coupling scales as ~ gap/L (measured slope -0.95)"
    ),
    "two_plane_design": (
        "The two-plane -> 3-D method: design in transverse (r,z) + azimuthal "
        "(s,z) planes, reflect into a 3-D sector pole; FFAG scaling sector, "
        "reflection validity = the L/g leaf coupling"
    ),
    "beam_referenced_twist": (
        "The beam-referenced equipotential SURFACE as the design primitive + "
        "the TWIST: rotate the surface by phi <=> multipole phase n*phi "
        "(the n-fold law, verified on a twisting quadrupole)"
    ),
    "all": "Everything (all topics concatenated)",
}

END_POLE_DESIGN = """
# Analytical end-pole chamfer for accelerator magnets
(Delferriere-de Menezes-Duperrier, CEA Saclay, SOLEIL design context)

## The problem

Accelerator magnets (dipoles, quadrupoles, sextupoles) need to deliver
a precise field profile along the beam axis with tight multipole
tolerances (~10⁻⁴ relative).  In the 2D cross-section, this is done by
careful pole-tip shaping.  But in 3D, the END of the magnet introduces:

- A fringe field region with mixed-multipole content
- An effective magnetic length that depends on end geometry
- Local saturation at the pole tip corners

A simple 45° chamfer (the historical default) is INADEQUATE for
modern light sources (SOLEIL, ESRF, NSLS-II) where tolerances reach 10⁻⁵.

## The analytical model

Start from the 3D magnetic scalar potential expansion:

    V(r, θ, z) = Σ_n (A z + B_n(z)) · r^n · (a_n cos nθ + b_n sin nθ)

For the desired multipole order n (n=1 dipole, n=2 quadrupole, n=3
sextupole), set L_m = L_f (magnetic length = iron length) at the
design field.  Then the longitudinal pole profile r(z) for z in the
end region II is:

    r(z) = ∆ · (1/2 - z/L_f)^(1/n)

where ∆ is a depth parameter, L_f is the iron length, and n is the
multipole order.

For a dipole (n=1):  r(z) = ∆ · (1/2 - z/L_f)   — linear taper
For a quadrupole (n=2): r(z) = ∆ · √(1/2 - z/L_f)   — square-root
For a sextupole (n=3): r(z) = ∆ · (1/2 - z/L_f)^(1/3)   — cube-root

## The depth ∆ determines magnetic length

For the dipole with half air-gap g:
    ∆ = L_f · (g·B_dyn/(μ₀·NI) - 1)

where B_dyn is the desired uniform-field value and NI is the
ampere-turns.

For the quadrupole with bore radius r_g and gradient G_dyn:
    ∆ = L_f · (r_g²·G_dyn/(2·μ₀·NI) - 1)

## Numerical verification (TOSCA)

Tested in TOSCA 3D for both dipole and quadrupole geometries:
- Sharp end (no chamfer): integrated multipoles A_3, A_5, A_7 ~ 10⁻³
- 45° chamfer: ~10⁻⁴
- Analytical profile (10-slope approximation): ~10⁻⁵ — 3 orders of
  magnitude improvement

The first slope is NOT 45° — it depends on ∆ and the multipole order.
Smooth profile via 10 slopes is a manufacturable approximation.

## Connection to Radia

Radia's strength is rapid evaluation of magnet geometries WITHOUT
iron meshing — perfect for parametric studies of end-pole ∆.

A Radia workflow:
1. Use the analytical formula r(z) = ∆ · (1/2 - z/L_f)^(1/n) to
   generate the pole tip surface
2. Use Radia ObjHexahedron / ObjPolyhedron with magnetization
   M_s from the saturation curve
3. Sweep over ∆ values, evaluate the integrated multipoles via
   rad.Fld() along beam trajectory
4. Choose ∆ that minimizes A_3 + A_5 + A_7 (or the relevant
   multipole budget)
5. Validate with TOSCA / NGSolve (using the radia-mcp.radia_ngsolve
   FEM pipeline)

This is exactly the workflow used at SOLEIL, ESRF, and Kindai's
accelerator partners.
"""

KOLKATA_CYCLOTRON = """
# Radia validation: Kolkata Superconducting Cyclotron case study
(Pradhan et al, VECC Kolkata, Cyclotrons 2007)

## The problem

Kolkata SCC has 2 superconducting coils (α-coil + β-coil) excited
to various current combinations.  Internal field can only be measured
up to r = 673 mm (extraction radius).  Beyond that, simulation must
extrapolate.

## Why Radia + TOSCA

The paper compares two simulation codes:

| Aspect | TOSCA (FEM) | Radia |
|---|---|---|
| Solver | Finite element | Boundary integral |
| Mesh | 450k elements needed | None |
| Geometry input | Complex (FEM mesh) | Simple (Mathematica) |
| CPU time | Long (1-2 hours) | Fast (few minutes) |
| Mathematica integration | None | NATIVE (RADIA exports to MATHEMATICA) |
| Saturation handling | Native | Needs trick (see below) |

The paper explicitly says: "RADIA requires geometry creation using
simple MATHEMATICA code, requires less CPU time of the solver."

This validates the radia-mcp ecosystem's Mathematica+Radia coupling —
it's not a new invention, it's the natural way RADIA has always been
used at major accelerator labs.

## The "saturated pole" trick (when FEM mesh is too expensive)

For complex SC cyclotron pole tips (shims, trim-coil inserts,
splittings):

1. Assume pole-tip surface is FULLY SATURATED at design field
2. Compute the field from saturated surface B_sat(r, θ) from ONE
   measured data set (calibration)
3. Compute the coil field B_coils(r, θ, I_α, I_β) analytically
4. The yoke contribution B_iron(r, θ, I_α, I_β) is computed via
   Radia (just the iron, ignoring detailed pole-tip features)
5. Total:  B_total = B_sat + B_coils + B_iron

This DECOMPOSITION trick avoids meshing the detailed pole-tip
geometry while preserving the dominant saturation physics.

## Results

Average iron field: |measured - simulated| < 0.5% across all coil
current combinations (extrapolation to higher currents than
measured).

3rd, 6th, 9th harmonics: agreement < 0.2%.

This validates Radia as a TRUSTWORTHY tool for SC cyclotron design,
including extrapolation beyond measured operating points.

## Connection to radia-mcp

The Kolkata workflow is exactly what radia-mcp.electromagnet supports:
- Build the coil with CoilBuilder
- Build the iron yoke with Radia primitives
- Apply BH curve with MatSatIsoTab
- Compute field via rad.Fld()
- Verify via Mathematica (radia-mcp.mathematica) symbolic identities
- Use TOSCA/NGSolve for cross-validation when needed

The Pradhan 2007 paper is THE blueprint for using Radia at an
accelerator facility.
"""

ROTATING_COIL_MEASUREMENT = """
# Rotating coil measurement and field reconstruction

The standard technique for measuring multipole content of accelerator
magnets.  Closely tied to Radia simulation for validation.

## The principle

A coil rotates around the magnet axis.  As it rotates, the magnetic
flux through it varies sinusoidally with frequency components:
    n = 1 → dipole
    n = 2 → quadrupole
    n = 3 → sextupole
    ...

Fourier-decompose the induced voltage to extract the multipole
amplitudes.

## What you get

For each multipole order n, you measure:
- The complex amplitude C_n (magnitude + phase)
- The "main" component (the design multipole)
- The "allowed" multipoles (with same symmetry as main)
- The "forbidden" multipoles (manufacturing errors)

## Local field reconstruction

A single rotating-coil measurement gives INTEGRATED multipoles
along the entire magnet length.  But for storage-ring optics, the
LOCAL field profile along z matters.

The reconstruction problem:
- Inputs: rotating-coil C_n(z) at multiple z positions OR
  a single integrated measurement
- Outputs: local B_x(x, y, z), B_y(x, y, z) field map

Approach 1: many rotating coils at different z (expensive)
Approach 2: 1 rotating coil + 3D Radia model fit to measurement
Approach 3: integrated + analytical end-pole model (Delferriere)

## Connection to radia-mcp

For Sugahara lab accelerator partner work:
1. Measure: rotating coil → integrated C_n
2. Simulate: Radia 3D model with proposed end-pole ∆ parameter
3. Compare: integrated C_n_sim should match C_n_measured
4. Iterate: adjust ∆ until match
5. Local field: use the validated 3D model to predict B(x, y, z)
   along beam orbit
6. Cross-validate: use radia-mcp.mathematica for symbolic
   verification of any analytical end-pole formulae

This is the canonical magnet validation cycle, supported entirely
by the radia-mcp tool stack.
"""


ISOCHRONOUS_ENDPACK_DESIGN = """
# Radial field-index design: scaling vs isochronous, and the
# saturated nonlinear END PACK (super-ferric pole reshape)

Beyond the longitudinal end-pole chamfer (the `end_pole` topic), a
CIRCULAR machine (FFAG / cyclotron) has an orthogonal design axis: the
RADIAL field index

    k(r) = d log B_y / d log r

that controls how the gap field grows with radius.  This is a
hodograph-native, single-valued SHAPE design (no topology change),
shipped + golden-tested in
`examples/clebsch_hodograph/scaling_ffag_pole_2d.py`.

## Two achromaticities -- and they are DIFFERENT (mutually exclusive)

"Achromatic" (momentum-independent) means one of two distinct things:

| | scaling FFAG | isochronous (cyclotron / non-scaling FFAG) |
|---|---|---|
| invariant | betatron TUNE | revolution TIME |
| field law | B_y(r) = B0 (r/r0)^k | <B>(r) = B0 gamma(r) |
| field index | k = const (rigid) | k_iso(r) = (beta gamma)^2 = beta^2/(1-beta^2) (RISING) |

with beta(r) = beta0 (r/r0) (beta = omega r / c, linear in r at a
constant revolution frequency).  In u = log r a momentum scaling
r -> lambda r is a TRANSLATION; the scaling field is translation-
covariant (k = const, a straight log-log line), while isochronism
DELIBERATELY breaks that symmetry (rising k_iso).  Relativistically you
cannot have both: scaling fixes the tune, isochronism fixes the time.
The pole gap is g(r) ~ 1/B(r) in either case (thin-gap, B ~ 1/g).

## The super-ferric wall = the NONLINEAR END PACK

With a Froehlich mu(B) iron pole, the high-r edge carries the highest B,
saturates FIRST, and the achieved field index DROOPS there -- degrading
achromaticity at the high-energy edge of the momentum acceptance.  For
the ISOCHRONOUS magnet this is most acute: the high-r end must deliver
the STEEPEST rise exactly where the iron gives out.  That high-r,
highest-B region IS the nonlinear end pack.

## The fix: the SAME hodograph machinery for both targets

A 2-parameter pole reshape in the log / von Mises chart

    g = g0 exp(-k u - gamma/2 u^2 - gamma2/6 u^3),
    local index k_geom(u) = k + gamma u + gamma2/2 u^2

(single-valued: the full 2-variable hodograph FOLDS once mu = mu(q), so
the von Mises single-variable chart is used) drives the SATURATED index
back onto the target via a 2-D Newton on (tilt, curvature):

  * scaling  (`run_step3`): target k = const -> flatten the saturated
    index ~7.2x in one Newton step.
  * isochronous (`run_isochronous`): target the RISING k_iso(r) -> drive
    the saturated <B>(r) back onto B0 gamma(r), restoring isochronism
    ~3.1x.  Measured: field-shape residual |<B>/(B0 gamma) - 1| goes
    2.3% -> 0.73% at B_gap ~ 1.33 T > knee Bk = 1.2 T (k_iso rises
    0.28 -> 1.44 across the aperture).

## Certified INTO saturation (the A/phi complementary bracket)

The same operating point is solved BOTH ways -- phi (Dirichlet on the
poles) and A (Dirichlet on the flux walls, driven to the phi-solve's
median flux so both sit at the SAME saturation state).  Monotone BH =>
convex energy => the energy bracket survives into saturation (Synge
hypercircle / Rikabi-Bryant-Freeman).  k_phi(r) and k_A(r) converge from
discretisation-complementary sides, so a tight gap (~5e-4) certifies the
saturated index is PHYSICS, not mesh -- for the scaling AND the
isochronous reshaped pole.

## Hodograph AS the solver (no remesh)

For the linear pole, `run_pullback` solves on a FIXED computational mesh
with the pole shape entering as a pullback DEFORMATION
(mesh.SetDeformation, weight W = |det J|(J^T J)^{-1}), so a reshape is a
new WEIGHT on the same mesh -- Netgen runs ONCE for the whole shape
sweep (reproduces the physical-remesh k(r) to ~5e-4).

## Honest scope (repo-first)

  * The reshape residual (0.73% isochronous) is the higher-order
    mismatch a 2-PARAMETER quadratic reshape leaves against a ~5x rising
    k_iso -- more shape DOF closes it; the Newton converges in one step.
  * This is the RADIAL <B>(r) isochronism ONLY.  The AVF flutter
    (vertical focusing), the betatron tunes, and the orbit<->field
    self-consistency (the closed-orbit r(p)) are SEPARATE problems, not
    modeled here.
  * The no-remesh pullback is shown for the LINEAR pole; threading it
    through the nonlinear saturated Newton is the next rung.

## Connection to radia / where this sits

This is the radial-index companion to the `end_pole` (longitudinal
chamfer) topic and the `kolkata` SC-cyclotron case study.  It is part of
the electromagnet / accelerator-magnet domain (Clebsch-hodograph
pole-face inverse design, `examples/clebsch_hodograph/`,
`docs/clebsch_hodograph/DESIGN_METHODOLOGY.md` section 3.3).  Goldens:
`tests/feec/test_clebsch_hodograph_research.py`
(`test_scaling_ffag_pole_2d_step1/_step2_saturation/_step3_reshape/
_saturated_bracket/_pullback_solver/_isochronous`).
"""


FOLIATE_PERTURB = """
# Foliate-and-perturb: when can a 3-D magnet be designed as a 2-D body
# plus an end perturbation?  (Leaf coupling scales as ~ gap/L)

A beamline / circular magnet is QUASI-2-D: the long beam direction varies
slowly and only the magnet ENDS are genuinely 3-D.  The natural design
scheme is therefore to FOLIATE the magnet into 2-D cross-section "leaves"
along the beam:

    0th order : STACK the 2-D cross-section (leaf) solution along the beam;
    1st order : CONNECT adjacent leaves by a beam-direction perturbation
                (the ENDS / fringe).

This pays off ONLY if the inter-leaf coupling is small AND localised to the
ends.  Whether it does is governed by the ASPECT RATIO L_iron / gap, and is
MEASURED (no model assumption) in
`examples/clebsch_hodograph/leaf_coupling_perturbation_3d.py` on a real
reduced-Omega + CoilBuilder finite-length C-frame dipole (beam = y, gap = z;
geometry parametrised by L so L/gap can be swept).

## What is measured (a straight, constant-gap magnet)

For a straight magnet the BODY slice y=0 IS the 2-D infinite-long leaf, so:

  * delta(y) = || B_perp(.,y) - B_perp(.,body) || / || B_perp(.,body) ||
        the 0th-ORDER leaf-stacking error (~0 in the body, grows at the ends);
  * eps(y) = (g/2) |dBz/dy| / |Bz_body|
        the local PERTURBATION PARAMETER = (transverse scale)/(beam-variation
        scale).  Use THIS, not an operator-norm ratio ||d^2/dy^2|| / ||grad_perp^2||,
        which is trapped at 1 by the current-free Laplace identity
        grad_perp^2 Omega = -d^2 Omega/dy^2 in air;
  * fringe_excess = (L_eff - L_iron)/L_iron
        the integrated 1st-order (inter-leaf) correction.

## THE RESULT: leaf coupling decays as ~ gap/L

A log-log fit of the fringe excess vs L/gap gives slope -0.95 (i.e. ~ gap/L):

    L/gap        2       3       5       8
    fringe     +180%   +111%   +70%    +48%

Extrapolating, the fringe drops to ~10% near L/gap ~ 40.

## Consequences for design

  * A COMPACT magnet (e.g. an end-study dipole, L/gap = 3) is firmly
    NON-PERTURBATIVE: +111% fringe, the 0th-order leaf stack misses ~40% of
    the integrated dipole, and the 3-D-ness is NOT end-localised (it is the
    whole magnet).  You CANNOT foliate it -- solve it fully 3-D.
  * Foliate-and-perturb LANDS only for LONG magnets (L/gap >> 1, ~10% fringe
    near L/gap ~ 40 -- typical beamline dipoles).  There the BODY is a 2-D
    cross-section design and ONLY the ends need 3-D treatment.
  * That 3-D end treatment is exactly the `end_pole` / equipotential-following
    design: an equipotential-following end (Delferriere r(z)=Delta(1/2-z/L)^(1/n),
    the beam-referenced equipotential surface) removes the fringe's HARMONIC
    contamination, but NOT the fringe itself (L_eff > L_iron is a free-space
    effect).  So: the END fixes the integrated STRENGTH (L_eff), the body 2-D
    design fixes the FIELD QUALITY (harmonics).

## Honest scope

This measures the SCALING of the BARE-end 3-D-ness for a STRAIGHT magnet.
Curved orbits (combined-function, the beam-referenced equipotential surface
twisting along a bent orbit) and the orbit<->field closed-orbit
self-consistency are separate problems.  The point established here is the
DECISION RULE: compute L/gap, and only foliate when the implied fringe (~gap/L)
is below your field-integral tolerance.
"""


TWO_PLANE_DESIGN = """
# The two-plane -> 3-D design method (the FFAG scaling sector cell)

An accelerator (bending / focusing) magnet is designed in TWO orthogonal 2-D
planes, and the two designs are REFLECTED into ONE 3-D pole.  The radial-index
design (isochronous_endpack) and the foliate-and-perturb leaf coupling
(foliate_perturb) are fragments of this single method.

## The two planes (circular / FFAG sector geometry: orbit = azimuthal arc)

  * TRANSVERSE plane (r, z), perpendicular to the orbit:
        sets the field the beam SEES -- the scaling field index
        B_z(r) ~ r^k (k = const => momentum-independent tune), pole gap
        g(r) = g0 (r/r0)^(-k).  Engine: scaling_ffag_pole_2d.py (Plane A).
  * AZIMUTHAL plane (s, z), along the orbit (s = r0 * theta):
        sets how the field TURNS ON/OFF -- the sector ENDS, the fringe, the
        effective magnetic length L_eff = INT B_z ds / B_z(body).  A
        finite-length iron pole of gap g(r0) over the sector arc
        L_sector = r0 * dtheta.  Engine:
        ffag_sector_two_plane.py::solve_azimuthal_end (Plane B).

The 3-D pole is the (r,z) gap profile SWEPT around the sector arc and truncated
at the azimuthal ends shaped by Plane B.  The two 2-D designs are EXACT in the
body and couple only at the ends.

## The reflection's validity = the leaf-coupling aspect L_sector / g

Measured (ffag_sector_two_plane.py, ngsolve only, golden-tested):

  * Plane A: scaling index k ~ 4.88 (vs design 5; the naive pole droop,
    reshaped per isochronous_endpack), A/phi bracket ~2e-6 (physics, not mesh).
  * Plane B: EACH sector end adds ~0.75 g of effective length, so
        L_eff = L_sector + ~1.5 g,
    and the fringe excess (L_eff - L_sector)/L_sector falls as ~ gap/L:

        L_sector/g      2      3.5     6      10
        fringe        +73%    +43%    +25%   +15%

    log-log slope -0.98 -- the SAME leaf-coupling ~ gap/L law as the straight
    magnet (foliate_perturb), now on the AZIMUTHAL plane.

So the aspect L_sector / g(r0) is the design lever: a COMPACT cell (L/g = 3 ->
+50% fringe) is NON-perturbative -- the sector ENDS are a genuine 3-D problem,
not a body-stack correction.  The two-plane reflection is exact only as
L/g -> infinity.

## Rung 2 (the 3-D reflection -- ESTABLISHED, ngsolve, golden-tested)

Sweep g(r) around the sector arc -- REVOLVE the (r,z) gap cross-section about
the bend axis -- into a 3-D iron pole, driven as an IRON-POLE EQUIPOTENTIAL
(upper-half model, median z=0 the up-down antisymmetry plane, the high-mu pole
back at Psi=mmf; the 3-D form of Plane A/B's scalar potential).  Verified
(ffag_sector_two_plane.py --rung2): the orbit RECOVERS B_z(r) ~ r^k -- field
index mean ~4.88 (range [4.6, 5.2], design 5; the naive-pole droop + mesh
scatter) -- so the swept g(r) ~ r^(-k) pole reproduces the designed radial field
in full 3-D; the azimuthal sector ends add L_eff/L_sector - 1 ~ +39% (the
Plane-B fringe, cross-checked in 3-D).  The field index is set by the pole
GEOMETRY (a high-mu equipotential forces B ~ 1/g(r)), so it is DRIVE-INDEPENDENT
-- a CoilBuilder + reduced-Omega coil drive sets the field AMPLITUDE (a further
step), not the index.

## Honest scope

Plane B is the LINEAR-iron azimuthal-END geometry (the magnetic-length excess;
saturation is Plane A's lever, composed orthogonally).  The radial profile is
<B>(r) only -- AVF flutter (vertical focusing), the betatron tunes, and the
orbit<->field self-consistency are separate.  This is the METHOD SCAFFOLD: two
2-D hodograph-native designs + a measured L/g reflection criterion; the
curved-orbit twist (combined-function, the beam-referenced equipotential
surface rotating along a bent orbit) is the next axis.
"""


BEAM_REFERENCED_TWIST = """
# The beam-referenced equipotential SURFACE as the design primitive + the TWIST

Instead of solving a magnet and then multipole-EXPANDING the field, make the
beam-referenced equipotential SURFACE the design SPEC.  In the Frenet frame of
the orbit s the iron pole face is

    Omega(r, theta; s) = sum_n r^n b_n(s) sin(n theta + phi_n(s))

(high-mu => H_tangential = 0 => Omega = const), so the multipole (b_n, phi_n)(s)
IS the surface's angular Fourier mode.  Design = prescribe (b_n, phi_n)(s), sweep
the equipotential surface along the orbit, place iron there -- NOT
solve-then-expand.

## The twist (the curved-orbit / combined-function axis)

The genuinely 3-D content of a CURVED-orbit / COMBINED-function magnet is that
the transverse multipole ROTATES along s (the Frenet frame turns with the bend;
a rotating-gradient magnet turns the pole on purpose).  The key fact is the
N-FOLD LAW:

    rotate the equipotential SURFACE by phi  <=>  multipole PHASE turns n*phi.

For the quadrupole (n=2): rotate the pole by phi <=> (b_2, a_2) ->
|b_2| (cos 2phi, sin 2phi), so a quad twisted by 45 deg becomes a pure SKEW quad.

## What is verified (twisting_quadrupole_pole.py, ngsolve only, golden-tested)

The quad pole face is the hyperbola xy = +-r0^2/2 (the Omega=const equipotential,
accel_pole_design.quad_pole_hyperbola).  A 2-D Laplace solve in the aperture with
the 4 hyperbola poles at alternating +-Omega0 recovers a CLEAN quad:
  - skew a_2/|c_2| ~ 5e-6 (pure normal at phi=0)
  - forbidden n=1,3,5 at the ~5e-5 floor
  - leading allowed spurious = the finite-pole 12-pole b_6 ~ 5.6e-3
Rotating the 4 poles by phi rotates the recovered pole orientation
alpha = -(1/2) atan2(a_2, b_2) by EXACTLY phi (slope 1.000, max error 0.00 deg),
and b_6 is rotation-INVARIANT.  The twist is the surface angular mode, measured.

## Honest scope

This is the per-station (Frenet cross-section) 2-D design -- the SLOW-TWIST
(adiabatic) limit d phi / ds -> 0, where the magnet is a stack of 2-D leaves (the
foliate_perturb picture, now twisting).

The CONFLUENCE with the FFAG sector (rung 1-2) is DONE
(combined_function_frenet_sweep.py): a COMBINED-FUNCTION magnet (dipole b1 + quad
gradient b2 in one tilted-gap cross-section) swept along the CURVED orbit it bends.
In the Frenet frame the cross-section is fixed; the Frenet rotation theta(s) =
s/rho (rho = Brho/b1) twists the lab pole, and the n-fold law gives dipole phase
theta, quad phase 2 theta -- BOTH orientations track theta (slope 1.000, err
0.00 deg), the multipole phase-change ratio psi2/psi1 = 2.000.  Verified ngsolve,
golden test_combined_function_frenet_sweep.

REMAINING: a fast twist / tight bend couples adjacent leaves (a longitudinal-field
correction); the twist rate d phi / ds is a leaf-coupling perturbation parameter
(the next rung -- when does the per-station 2-D break).  A spiral sector (pole
twist phi != orbit bend theta) and an s-ramped (b1, b2)(s) are the other
extensions the n-fold law governs.
"""


def get_accelerator_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "end_pole"             - Analytical chamfer design (Delferriere)
      "kolkata"              - Radia + TOSCA validation case study
      "rotating_coil"        - Multipole measurement + field reconstruction
      "isochronous_endpack"  - Radial field index: scaling vs isochronous,
                               saturated nonlinear end-pack reshape
      "foliate_perturb"      - When can a magnet be a 2-D body + end
                               perturbation?  Leaf coupling ~ gap/L
      "two_plane_design"     - The two-plane -> 3-D method: transverse (r,z)
                               + azimuthal (s,z) -> 3-D sector pole
      "beam_referenced_twist" - The equipotential SURFACE as design primitive
                               + the twist (n-fold law, twisting quadrupole)
    """
    topic = topic.lower().strip()
    if topic in ("end_pole", "chamfer", "delferriere"):
        return END_POLE_DESIGN
    if topic in ("kolkata", "cyclotron", "validation"):
        return KOLKATA_CYCLOTRON
    if topic in ("rotating_coil", "multipole", "measurement"):
        return ROTATING_COIL_MEASUREMENT
    if topic in ("isochronous_endpack", "isochronous", "scaling_ffag",
                 "field_index", "end_pack", "endpack"):
        return ISOCHRONOUS_ENDPACK_DESIGN
    if topic in ("foliate_perturb", "foliate", "leaf_coupling", "quasi_2d",
                 "quasi2d", "perturbation"):
        return FOLIATE_PERTURB
    if topic in ("two_plane_design", "two_plane", "twoplane", "sector",
                 "ffag_sector", "reflection"):
        return TWO_PLANE_DESIGN
    if topic in ("beam_referenced_twist", "twist", "twisting", "design_primitive",
                 "equipotential_surface", "n_fold", "rotating_gradient"):
        return BEAM_REFERENCED_TWIST
    if topic == "all":
        return "\n\n".join([
            END_POLE_DESIGN, KOLKATA_CYCLOTRON, ROTATING_COIL_MEASUREMENT,
            ISOCHRONOUS_ENDPACK_DESIGN, FOLIATE_PERTURB, TWO_PLANE_DESIGN,
            BEAM_REFERENCED_TWIST,
        ])
    return (
        f"Unknown topic '{topic}'. Available: all, end_pole, kolkata, "
        "rotating_coil, isochronous_endpack, foliate_perturb, two_plane_design, "
        "beam_referenced_twist."
    )
