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
    "sector_saturation": (
        "The scaling-FFAG SECTOR body driven into SATURATION: the azimuthal (s,z) "
        "end made NONLINEAR (Froehlich mu_r(|B|)) -> the two sector planes respond "
        "OPPOSITELY -- the azimuthal effective length L_eff is ROBUST (gap-reluctance-"
        "dominated; barely moves even as the high-r iron <mu_r> collapses x0.3), while "
        "the radial field index k(r) is FRAGILE (droops, the achromaticity wall) "
        "(scaling_ffag_sector_saturation.py)"
    ),
    "bending_endpack_saturation": (
        "The 2D BENDING-magnet end pack optimized vs iron saturation: the median-plane "
        "field the BEAM sees is NATURALLY saturation-invariant (normalized profile the "
        "same linear/saturated, even flat-cut); the real problem is the pole-TIP CORNER "
        "hot spot (kappa=peak iron|B|/body ~3.7). Minimize kappa over the end chamfer "
        "(depth, exponent) -> the chamfer RELIEVES the corner while the beam field stays "
        "put (bending_endpack_saturation_opt.py)"
    ),
    "excitation_invariant_field": (
        "'Same flux lines when you turn up the current' (excitation-invariant, NOT a "
        "cyclotron): below the iron knee the magnet is LINEAR, so scaling the current "
        "scales B everywhere -> the flux-LINE pattern (streamlines B/|B|) is IDENTICAL. "
        "Metric = air-region flux-line DIRECTION drift D_dir(I)=rms||bhat(I)-bhat(I_lin)|| "
        "across an excitation sweep; the LINEAR control is exactly 0 at every drive "
        "(scaling current cannot move flux lines), saturation is the SOLE breaker. Even "
        "flat-cut is sub-degree invariant (~1.5 mrad); minimizing saturated D_dir over the "
        "end chamfer keeps flux lines ~6-7x MORE invariant (corner relief, same lever as "
        "bending_endpack_saturation) (docs/clebsch_hodograph/excitation_invariant_field.ipynb)"
    ),
    "hodograph_feasibility": (
        "2D LINEAR hodograph outlook for a BENDING magnet, TWO facts of different kinds. "
        "(1) FEASIBILITY = harmonic analysis, NOT the hodograph: demand a mid-plane field "
        "B_y(x,0) flat with edge width d; the gap field is the unique HARMONIC continuation, "
        "and the tanh continuation's singularity at pi*d/2 must clear the gap h -> d* = (2/pi)h "
        "~ 0.64h (edge no sharper than ~0.64 x gap). Cauchy-Kovalevskaya/analyticity fact, NOT a "
        "'limit line' (passive saturation stays ELLIPTIC: mu falls but |B|=mu|H| still rises, "
        "d(mu q)/dq>0, no type change). ngsolve linear-FEM verified: read-off pole reproduces "
        "g(x) to 0.01%, FEM equipotential matches pole to ~1e-6. "
        "(2) TRANSPARENCY = the hodograph proper (partial von Mises, keep x transform ONE "
        "potential): replace y by s=-phi -> the unknown iron-pole free boundary becomes the "
        "FIXED top edge s=s0 of a rectangle; map y(x,s) solves the von-Mises PDE "
        "y_s^2 y_xx - 2 y_x y_s y_xs + (1+y_x^2) y_ss = 0, Jacobian y_s>0 (no fold); residual "
        "->0 under refinement (2.6e-3->7.2e-4), top edge = read-off pole to ~1e-8. Saturation "
        "enters as a coefficient mu(q) on the SAME fixed rectangle (no new free boundary). "
        "(docs/clebsch_hodograph/hodograph_feasibility_2d.ipynb)"
    ),
    "hodograph_bending_sy": (
        "2D bending-magnet END in the s-y (longitudinal) plane: fringe feasibility + "
        "end-shaping DESIGN. Forward 2D s-y is OPEN-boundary (dipole flux return is out-of-plane) "
        "so the effective-length integral is LOG-DIVERGENT (g~1/s); a phi=0 box gives a "
        "box-DEPENDENT fringe. (1) Fringe feasibility (hodograph inverse, box-free): demand "
        "B_y(s,0)=1/2(1-tanh((s-s0)/d)); nearest singularity pi*d/2 -> same bound d*=(2/pi)h~0.64h, "
        "now capping the LONGITUDINAL fringe (Enge edge / EFB) sharpness; manufactured-solution FEM "
        "verifies to 0.01%. (2) End-shaping DESIGN (forward FEM optimization): the pole END is a "
        "FREE termination, NOT an equipotential (reading {phi=-h} to the end curls into the demand "
        "singularity, |B| diverges). Square end = reentrant 270deg corner -> |B|~r^-1/3 (saturation "
        "hot spot). Optimize a forward quarter-ellipse chamfer (length a, rise b) to MINIMIZE peak "
        "pole-face |B| (a LOCAL, box/mesh-convergent quantity): peak|B|/B0 = 2.14 (near-square) -> "
        "1.42 (round) -> 1.03 (a=3.2h,b=1.0h) -> 1.01 (a=4.4h), numerically RECOVERING the classic "
        "Rogowski electrode (no enhancement) without conformal algebra. Honest: 2D L_eff/EFB is "
        "log-divergent (needs 3D); the END design targets no-enhancement (local). "
        "(docs/clebsch_hodograph/hodograph_bending_sy.ipynb; golden validation_test/feec/test_hodograph_bending_sy.py)"
    ),
    "edge_focusing_tracking": (
        "Vertical EDGE FOCUSING of a tilted dipole end, measured by PARTICLE TRACKING (the s-y "
        "end-shaping companion; here the tilt is HORIZONTAL: pole-face rotation beta about the "
        "vertical axis, x-s bend plane). Hard-edge thin lens |1/f_z|=tan(beta)/rho. CORRECT "
        "measure = the linearized vertical HILL INTEGRAL along the reference orbit "
        "1/f_z=(q/p)INT(u_y dBx/dz - u_x dBy/dz)|z=0 ds -- NOT a field-EFB slope (which is "
        "wrong-sign/blows up on a compact dipole; edge_focusing_efb_slope_negative). VERIFIED on a "
        "genuinely Maxwellian tilted fringe (curl-free AND div-free: the only vacuum linear-in-z "
        "continuation is B_s=+B0 z g'(s); a div-free-but-not-curl-free choice has a spurious edge "
        "current sheet and FLIPS the sign): slope vs the law -> 1 as fringe w->0 (0.84->0.99 for "
        "w 0.08->0.005), beta=0 baseline = -0.5 w/rho (finite-fringe residual ->0), rho*(1/f_z) "
        "collapses onto tan(beta); matches the FULL SCOFF/Enge law tan(beta-psi)/rho, "
        "psi=(K1g/rho)(1+sin^2)/cos, to <=0.7 pct. 3D-FEM chain VALIDATED (parallelogram dipole, "
        "both edges tilted, coil FOLLOWING the pole outline, exactly C2-symmetric -> C2 map "
        "symmetrization): closed-orbit symmetric traversal + window-decomposed Hill integrals give "
        "dK_in FEM/model=0.93, FEM/closed-form=0.94 at beta=20 rho=5 (exit DEfocusing matched 0.96; "
        "dK*rho const to 1 pct over rho 5..40; beta=0 floor +0.0013). CROSS-CHECKED 2026-07-13 with "
        "a SECOND independent engine, the FEEC HDiv-VIM (hdiv_scoff_study: iron-only tet mesh, no "
        "air discretization, batch rad.Fld map): dK_in agrees 0.8 pct at matched edge "
        "mesh (+-3 pct across meshes), dK/model 0.92-0.95 in EVERY configuration -> the -5..-7 pct "
        "model deficit is REAL 3D physics, NOT mesh error (earlier attribution RETRACTED): the "
        "local iso-field tilt near x=0 is only ~0.95-0.96 of the geometric tan(beta), so "
        "hard-edge/SCOFF bookkeeping with the geometric angle overpredicts by ~5 pct -- the "
        "EFFECTIVE edge angle is magnet-specific and this chain measures it. DESIGN LESSONS the "
        "measurement exposed: straight coil front across a tilted edge -> iso-field tilt far below "
        "the geometric angle (dK 0.55x; edge-angle bookkeeping ASSUMES the coil follows the pole "
        "contour); rigid whole-coil rotation breaks MMF topology (B0 3x collapse). RADIA NOTE: "
        "CoilBuilder.mirror('xy') now performs a true geometric mirror and is locked by a "
        "pointwise geometry/field regression; rad.RadiaField on a TrfOrnt-wrapped container "
        "still causes 0xC0000374 heap corruption, so bake transforms into primitives. "
        "(docs/clebsch_hodograph/edge_focusing_tracking.ipynb + edge_focusing_fem_results.json; "
        "golden validation_test/feec/test_edge_focusing_tracking.py incl. coil-cleanliness lock)"
    ),
    "beam_referenced_twist": (
        "The beam-referenced equipotential SURFACE as the design primitive + "
        "the TWIST: rotate the surface by phi <=> multipole phase n*phi "
        "(the n-fold law, verified on a twisting quadrupole)"
    ),
    "endpack_two_plane": (
        "The magnet END PACK in two planes -> 3-D: design the END in the x-y "
        "cross-section (shim zeroes b3) AND the s-y longitudinal plane (a "
        "standalone 2-D Laplace fringe -> Rogowski end chamfer + L_eff), reflect "
        "into one 3-D equipotential pole; the chamfer-depth sweep zeroes the "
        "pole-tip corner over-field (endpack_two_plane.py)"
    ),
    "spectrometer_endpack": (
        "The SPECTROMETER end pack, NONLINEAR: the pole-tip corner is a Chaplygin "
        "saturable THROAT (concentration kappa) that reaches the iron knee FIRST "
        "(corner saturates at B_gap=B_K/kappa, below the bulk knee), drifting the "
        "EFB / edge focusing; the Rogowski chamfer is the throat-width knob. "
        "Nonlinear design map at LINEAR cost + B-input A-formulation truth check "
        "(endpack_spectrometer_saturation.py)"
    ),
    "endpack_cobake": (
        "The two-plane end pack CO-BAKED into ONE pole face "
        "z(x,s)=g/2-delta(x/w)^2+lift(s): the x-y shim (delta) AND the s-y Rogowski "
        "chamfer (ghat) in one gap face -> clean transverse b_3,5 AND a rounded "
        "pole-tip corner at once (endpack_cobake.py); the PRECISION construction is a "
        "smooth OCC ThruSections LOFT (endpack_cobake_loft.py) that meshes baseline & "
        "shim at the SAME density (resolves the x-prism staircase artifact)"
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
`docs/clebsch_hodograph/demos/scaling_ffag_pole_2d.py`.

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
pole-face inverse design, `docs/clebsch_hodograph/demos/`,
`docs/clebsch_hodograph/DESIGN_METHODOLOGY.md` section 3.3).  Goldens:
`validation_test/feec/test_clebsch_hodograph_research.py`
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
`docs/clebsch_hodograph/demos/leaf_coupling_perturbation_3d.py` on a real
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


SECTOR_SATURATION = """
# The saturating sector body -- the two sector planes respond OPPOSITELY to saturation

The scaling-FFAG SECTOR is designed in two planes: the radial (r,z) field index
k(r) ~ r^k (the achromaticity) and the azimuthal (s,z) sector end (the effective
length L_eff = INT B_z ds / B_z(body)).  scaling_ffag_pole_2d.py showed the RADIAL
plane droops k(r) at the high-r edge under iron saturation (Step 2, the achromaticity
wall); the azimuthal end was solved only LINEARLY.  scaling_ffag_sector_saturation.py
drives the SECTOR BODY into saturation by making the azimuthal end solve NONLINEAR
(Froehlich mu_r(|B|), the same knee as Step 2) and solving it at the high-r aperture
edge (smallest gap, highest B) and the low-r body.

## The honest result: a CONTRAST between the two planes (golden-tested)

  azimuthal (s,z) L_eff : ROBUST  -- drift < 0.1% even where the high-r iron <mu_r>
                                     collapses 1574->481 (x0.3); the sector end is
                                     GAP-RELUCTANCE-DOMINATED (the same honest scope
                                     as clebsch_dipole_saturation_3d: a large-gap
                                     magnet's gap field softens only mildly with iron
                                     saturation).  The fringe stays ~1.5 gaps.
  radial   (r,z) k(r)   : FRAGILE -- droops dk_hi ~ -0.26 at the high-r edge (the
                                     achromaticity wall, scaling_ffag Step 2).

So SATURATION degrades the radial field SHAPE (the achromaticity) but NOT the
azimuthal end LENGTH.  Design insight: the high-r achromaticity needs the radial
reshape (scaling_ffag Step 3), while the sector ENDS are saturation-robust and need
no nonlinear end correction.  (A sector where L_eff WOULD drift needs R_iron ~ R_gap
-- a small gap or a necked iron path; the scaling pole's large gap makes the end
robust.)
"""


BENDING_ENDPACK_SATURATION = """
# The 2D BENDING-magnet end pack, optimized vs iron saturation

The natural design question for a bending dipole end: does the field stay put as the
iron saturates?  bending_endpack_saturation_opt.py answers it on a small-gap (24 mm),
near-knee 2D (s, z) end pack (Froehlich mu_r(|B|)), and the answer is a TWO-PART surprise.

## (1) The MEDIAN-plane field the beam sees is NATURALLY saturation-invariant

The normalized longitudinal profile p(s) = B_z(s)/B_body is the SAME curve, linear and
saturated -- ||p_sat - p_lin|| ~ 1e-3, L_eff shift < 0.2% -- EVEN for a hard flat cut.
The gap-reluctance-dominated pole face stays equipotential (its <mu_r> at the saturated
drive is still ~500), so the EFB / homogeneity / integrated-field SHAPE are all
excitation-invariant AUTOMATICALLY (only the overall amplitude softens = bulk saturation,
handled by calibration).  This is the (s,z) twin of scaling_ffag_sector_saturation's
azimuthal-L_eff robustness.

## (2) The real saturation problem is the pole-TIP CORNER

The corner concentrates flux -- a raw grid max hits ~6 T iron |B| (kappa = peak/body
~ 3.7), deep past the 1.2 T knee -- a hot spot that limits the achievable field and
wastes iron, WHILE the beam-plane field is fine.

## The optimization: relieve the corner, not the (already-invariant) beam profile

Minimize kappa = peak iron |B| / B_body at the saturated drive over the end chamfer

    z_face(s) = g/2 + depth * ((|s| - s_body)/(s_pole - s_body))^exponent

The chamfer rounds the tip so the flux is not forced through the corner singularity
(kappa -> ~1.0; a smooth L^10 end-iron proxy ~1.9 -> ~1.0, peak-proxy |B| ~3.0 -> ~1.6 T),
and there is a GENUINE 2D optimum in (depth, exponent) -- too much chamfer re-concentrates
the flux elsewhere (kappa is non-monotone in depth).  Optuna TPE if present, else a
built-in grid + local refine.  Golden-tested (test_bending_endpack_saturation_opt).

## The honest engineering answer

"Flux lines invariant under saturation, is that good?" -> for the BEAM it is largely
AUTOMATIC (the median-plane field shape is gap-reluctance-robust); what you actually
optimize is the CORNER the iron saturates.  The regime where the beam field WOULD drift
needs the whole pole (not just a tip) deep in saturation -- a very-small-gap / necked
design; a normal bending end is beam-invariant by construction.
"""


EXCITATION_INVARIANT_FIELD = """
# Excitation-INVARIANT flux lines -- same field-line shape as the current rises

"Same flux lines when you turn up the drive current" (excitation-invariant) is NOT a
cyclotron (where the field is MEANT to change with radius) -- it is the opposite: hold
the field-line PATTERN fixed as the excitation grows.  The result-bearing docs notebook
docs/clebsch_hodograph/excitation_invariant_field.ipynb (+ its .py helper) makes this a
direct, optimizable metric on the small-gap bending end pack, and ties it to the hodograph
design method (see below).

## Linearity => invariance is AUTOMATIC (the key physics)

Below the iron knee the magnet is a LINEAR magnetostatic system: scale the excitation by
alpha and B scales by alpha EVERYWHERE, so the field-LINE pattern (streamlines
bhat = B/|B|) is IDENTICAL -- only the amplitude grows.  So "same flux lines as the
current rises" is automatic in the linear regime; the ONLY thing that can move the flux
lines is the NONLINEARITY, i.e. iron SATURATION (mu(|B|) dropping non-uniformly, first at
the pole-tip corner), which redistributes the flux.

## The metric: air-region flux-line DIRECTION drift vs excitation

D_dir(I) = rms_x || bhat(x; I) - bhat(x; I_lin) ||   over air sampling points near the
pole end (gap + fringe).  D_dir = 0 while linear, grows once the iron saturates.
NOTE: the MEDIAN-plane profile p(s)=B_z/B_body does NOT show this (gap-reluctance-robust,
~1e-3 for any end shape -- see bending_endpack_saturation); the AIR field-line DIRECTION
near the tip/fringe is where saturation actually moves the pattern, so measure THAT.

Two computed controls (the proof):
  * LINEAR control (mu forced constant): D_dir = 0 EXACTLY at every drive -- scaling the
    current cannot move the flux lines (linearity).
  * SATURATED sweep: D_dir grows monotonically with the drive as <mu_r> falls
    (~1960 -> ~580 for a 0.15 -> 1.70 T drive), so saturation IS the sole breaker.

## Result + the design lever

Even a hard FLAT cut is already nearly invariant: D_dir < 1e-2 rad (sub-degree, ~1.5 mrad
at the saturated drive), because the high-mu pole face stays equipotential.  The residual
drift is pole-tip-corner-dominated, so minimizing the saturated D_dir over the end chamfer
(depth, exponent) keeps the flux lines invariant ~6-7x DEEPER into saturation
(~1.5 -> ~0.2 mrad) -- and it drops the corner kappa (~1.8 -> ~1.0) in LOCKSTEP, the SAME
corner-relief lever as bending_endpack_saturation, now judged directly on the air
flux-line geometry across an excitation sweep.  Optuna TPE if present, else grid+refine.
Golden-tested (test_excitation_invariant_field + test_excitation_invariant_linear_control_is_zero).

## What the HODOGRAPH buys (is it "just linear"?) -- NO, that is only part

The hodograph design method's payoffs, and where this study sits:
 1. DIRECT INVERSE DESIGN (the main payoff; holds even for purely LINEAR problems):
    instead of searching a pole shape whose forward-solve gives the target field, the
    hodograph PRESCRIBES the target field and READS the pole off as an equipotential
    (level set) -- one shot, no forward loop; exact multipole content, the exact geometric
    twist / n-fold law.  THIS is why a hodograph-designed pole is excitation-invariant in
    the linear regime: it IS the level set of a linear potential, and linearity => invariant
    flux lines (this notebook's linear control = 0 proves it).
 2. CHAPLYGIN LINEARISATION of SATURATION (the deeper payoff; nonlinear, 2D only): the
    nonlinear PDE div(nu(|grad A_z|) grad A_z)=0 becomes a LINEAR PDE in hodograph (field)
    coordinates (the gas-dynamics Chaplygin transform) -> the whole saturation design curve
    from ONE linear solve, no Picard loop.  3D does NOT auto-linearise.  So "it becomes
    linear" is real and deep, but it is the SATURATION-specific bonus, not the whole story.
 3. UNIFICATION with the Kelvin open boundary (differential forms): hodograph and Kelvin are
    both pullback / coordinate-transform tools; one machinery serves inverse design AND
    exact open boundary.

## The honest engineering answer

This notebook does NOT use the hodograph -- it is a FORWARD Froehlich-FEM + shape-optimization
study.  Its role is to CHARACTERIZE payoff (1)'s linear-regime byproduct (excitation-invariant
flux lines), show where SATURATION breaks it, and verify the fix is the corner-relief lever that
payoff (2) linearises.  "Same flux lines when the current rises" is largely AUTOMATIC -- it is
just linearity; what you design for is keeping it true DEEP into saturation, which -- again --
means relieving the pole-tip corner.  This is the excitation-sweep complement of
bending_endpack_saturation's iron-kappa view: that one minimizes the iron hot spot, this one
minimizes the air flux-line drift, and the two levers COINCIDE.
"""


HODOGRAPH_FEASIBILITY = """
# 2D LINEAR hodograph -- a feasibility boundary for a bending magnet

The linear (mu=const, Laplace) design outlook the hodograph gives on a bending-magnet
cross-section (docs/clebsch_hodograph/hodograph_feasibility_2d.ipynb, runnable + golden).

## The demand and the catch
Demand the mid-plane field B_y(x,0) = g(x): flat over the good-field half-width x0, rolling
off over an edge width d.  The scalar potential is HARMONIC in the gap, so g must extend
upward as its UNIQUE SMOOTH (analytic) continuation -- the mid-plane field is NOT free.

## (1) The feasibility boundary -- HARMONIC ANALYSIS, not the hodograph
The continuation of a tanh edge has its nearest singularity at height y = pi*d/2.  The field
is realizable by an iron pole at gap h ONLY IF that singularity clears the gap:

    d > d* = (2/pi) h  ~  0.64 h        -- the edge is no sharper than ~0.64 x gap.

Below d* the continuation is singular INSIDE the gap and no single smooth iron pole can make
the field.  This is a Cauchy-Kovalevskaya / analyticity fact about the DEMAND -- the hodograph
plays NO role, and it is NOT a 'limit line': ordinary passive saturation stays ELLIPTIC (mu
falls but |B|=mu|H| still RISES, d(mu q)/dq>0), so there is no type change / shock here (unlike
transonic gas dynamics -- do NOT frame this as a limit-line precursor).  The 2/pi is
tanh-edge-specific; the UNIVERSAL statement is the gap scaling ("fringe scale ~ gap", here
PROVEN as a feasibility boundary from analyticity).

## Inverse design + FEM verification
For a feasible demand the iron pole is READ OFF as the equipotential {phi=phi0} (no
forward-solve loop); phi has a closed form (sinh/cosh/sin/cos/atan2).  Verified by an
INDEPENDENT ngsolve linear-FEM manufactured-solution solve: (a) interior B_y(x,0) reproduces
g(x) to 0.01% of B0, (b) FEM phi == analytic phi (L2 ~1e-8), (c) FEM equipotential == read-off
pole to ~1e-6.  Golden: validation_test/feec/test_hodograph_feasibility_2d.py.

## (2) What the hodograph ACTUALLY buys -- the partial von-Mises chart
The Sugahara-lab hodograph is a PARTIAL transform: keep one coordinate, transform ONE
potential.  Keep x, replace y by the scalar-potential coordinate s = -phi (monotone in y since
B_y=-phi_y>0; s=0 on the mid-plane, s=s0 at the pole).  The UNKNOWN iron-pole SHAPE -- a FREE
BOUNDARY in physical space -- becomes the FIXED top edge {s=s0} of the rectangle
[-xr,xr]x[0,s0].  The physical map y(x,s) solves the (quasi-linear) von-Mises PDE

    y_s^2 y_xx - 2 y_x y_s y_xs + (1 + y_x^2) y_ss = 0        (mu = const)

with Jacobian y_s>0 everywhere (single-valued, no fold).  Verified on the closed-form map:
PDE residual -> 0 under refinement (2.6e-3 -> 7.2e-4, ~O(h^2)), min y_s ~ 0.78, and the top
edge reproduces the read-off equipotential to ~1e-8.  THIS is 'design in field space': the
pole is a coordinate line on a FIXED domain, and saturation enters as a COEFFICIENT mu(q) on
the SAME rectangle (no new free boundary).

## Where it sits
(1) is harmonic analysis; (2) is the genuine hodograph transparency.  The nonlinear (x,s) /
Chaplygin design (saturation as a coefficient on the fixed chart) is the natural next build
from this linear baseline.
"""


HODOGRAPH_BENDING_SY = """
# 2D bending-magnet END in the s-y (longitudinal) plane: fringe + end-shaping design

s = beam direction, y = gap (half-model, mid-plane symmetry, half-gap h).  The LONGITUDINAL
plane is where the magnet ENDS and the FRINGE lives -- the home of pole-end / Rogowski design
(docs/clebsch_hodograph/hodograph_bending_sy.ipynb, runnable + golden).

## Forward 2D s-y is OPEN-boundary
A dipole's flux return is OUT of this plane (through the x-y yoke), so a pure s-y slice leaks
flux and the effective-length integral is LOG-DIVERGENT (g(s)~1/s, INT B_y ds ~ INT du).  A
phi=0 air box gives a box-DEPENDENT fringe (verified: EFB set-back keeps growing with box).  So
split the work:

## (1) Fringe feasibility (hodograph inverse, box-free)
Demand B_y(s,0)=g(s)=1/2(1-tanh((s-s0)/d)); the continuation's nearest singularity at y=pi*d/2
gives the SAME bound d>d*=(2/pi)h~0.64h, now capping the LONGITUDINAL fringe (Enge edge / EFB)
sharpness.  A manufactured-solution FEM verifies the demanded fringe + pole-FACE equipotential
to 0.01% / ~1e-6 (box-free).

## (2) End-shaping DESIGN (forward FEM optimization) -- recovers Rogowski
The pole END is a FREE termination, NOT an equipotential: reading {phi=-h} all the way to the
end curls it INTO the demand's singularity (|B| diverges) -- the pole FACE is an equipotential,
the pole END is not.  A SQUARE end has a reentrant 270deg air corner -> |B|~r^-1/3 (saturation
hot spot).  Parametrize the end as a forward quarter-ellipse chamfer (length a, rise b) and
MINIMIZE the peak pole-face |B| by a forward FEM sweep -- the peak is a LOCAL quantity, hence
box/mesh-convergent even though L_eff is log-divergent.  Result:
    peak |B|/B0 = 2.14 (near-square) -> 1.42 (round) -> 1.03 (a=3.2h, b=1.0h) -> 1.01 (a=4.4h),
numerically RECOVERING the classic Rogowski electrode (|B|=B0, no enhancement) WITHOUT conformal
algebra.  The FEM optimization generalizes to arbitrary gaps and tilted ends where the closed
form is awkward.  Golden: validation_test/feec/test_hodograph_bending_sy.py.

## Honest scope
The 2D effective length / EFB is log-divergent (a 3D / finite-magnet quantity); the s-y END
design targets NO field enhancement (local).  The numerical chamfer optimum IS the (low novelty,
correct) Rogowski profile.  Pairs with the x-y cross-section note (HODOGRAPH_FEASIBILITY).
"""


EDGE_FOCUSING_TRACKING = """
# Vertical EDGE FOCUSING of a tilted dipole end -- measured by PARTICLE TRACKING

Companion of the s-y end-shaping note (HODOGRAPH_BENDING_SY shapes the LONGITUDINAL end so the
pole face stays at B0).  Here the tilt is HORIZONTAL: rotate the pole face about the VERTICAL
axis by beta (in the x-s bend plane).  That turns the edge into a thin VERTICAL lens
    |1/f_z| = tan(beta)/rho ,   rho = p/(qB0)  (bend radius).
(docs/clebsch_hodograph/edge_focusing_tracking.ipynb, runnable + golden.)

## Measure it by tracking, NOT by a field-EFB slope
Vertical edge focusing is a SECOND-ORDER, off-mid-plane property of the fringe.  The mid-plane
|B| effective-field-boundary (EFB) slope CANNOT recover it -- it is wrong-sign / blows up on a
compact dipole (lab finding edge_focusing_efb_slope_negative).  The correct measure is the
linearized vertical HILL INTEGRAL along the reference orbit:
    1/f_z = (q/p) INT ( u_y dB_x/dz - u_x dB_y/dz )|_{z=0} ds ,
with u the mid-plane orbit tangent (RK4).  By Maxwell (curl B=0 in air) dB_x/dz=dB_z/dx and
dB_y/dz=dB_z/dy, so the kernel can also be read off the clean mid-plane B_z gradients.

## Verified on a closed-form field
Build a GENUINELY MAXWELLIAN tilted hard-edge fringe (curl-free AND div-free).  Sign trap: with
mid-plane B_z(s,0)=B0 g(s), edge-normal s=(y-y_edge)cos b + x sin b, the ONLY vacuum linear-in-z
continuation is B_s=+B0 z g'(s) -> dB_x/dz=+B0 g'(s) sin b.  A div-free-BUT-NOT-curl-free choice
(B_s=-B0 z g'(s)) is a different field with a spurious current SHEET at the edge and FLIPS the
focusing sign.  Using the curl-free field, the tracker reproduces the hard-edge law:
  - slope of tracked 1/f_z vs tan(beta)/rho -> 1 as the fringe narrows: 0.84 -> 0.99 for
    w = 0.08 -> 0.005 (a finite fringe just spreads the delta-edge);
  - beta=0 baseline = -0.5 w/rho EXACTLY -- an O(w/rho) finite-fringe residual that vanishes;
  - rho*(1/f_z) collapses onto tan(beta) (thin lens, 1/rho scaling).
The magnitude tan(beta)/rho is the invariant; the sign is orientation-dependent (a rectangular
magnet's edges DEFOCUS; this curl-free entrance-edge orientation FOCUSES).  With the classical
first-order fringe correction psi = (K1g/rho)(1+sin^2 beta)/cos beta (K1g = INT g(1-g) ds along
the edge normal; = w/2 for the tanh fringe) the tracked values match the FULL law
tan(beta-psi)/rho to <=0.7 pct at w=0.02, and the beta=0 baseline IS -K1g/rho^2 exactly.
Golden (pure-numpy optics plus a small Radia coil check):
validation_test/feec/test_edge_focusing_tracking.py.

## 3D-FEM chain VALIDATED -- parallelogram dipole, SCOFF/Enge + closed orbit
A NAIVE Hill integral on a plain tilted-edge dipole is NOT enough (the thick 3D fringe + any
coil/edge mismatch pollute it -- do NOT re-attempt the wide-pole/difference variant).  The
protocol that works (PART B of docs/clebsch_hodograph/edge_focusing_tracking.py,
fem_scoff_study; committed results edge_focusing_fem_results.json):
  - TESTBED: PARALLELOGRAM dipole (both edges tilted beta -- the spectrometer configuration)
    with the COIL FOLLOWING the pole outline (rounded-parallelogram loop pair at fixed 20 mm
    normal offset, sides threading the iron circuit).  Parallelogram poles + centered legs +
    the C2-symmetric coil make the WHOLE magnet exactly C2-symmetric -> the C2-odd part of the
    sampled mid-plane map is pure solver error, removed EXACTLY by map symmetrization.
  - MEASUREMENT: closed-orbit SYMMETRIC traversal (2-pass shooting, |x|max ~ sagitta), Hill
    integrals window-decomposed (entrance/flat/exit x T1 transverse-gradient + T2
    edge-crossing), x-uniform model B0*Gin(s_in)*Gout(s_out) from the FEM's OWN mid-line
    profile, SCOFF/Enge closed form with measured K1g and orbit-corrected beta_eff.
  - RESULT (beta=20, rho=5): dK_in FEM +0.0690 vs model +0.0739 (0.93) vs closed form +0.0732
    (0.94) vs hard edge +0.0728; exit DEfocusing matched to 0.96; dK*rho constant to 1 pct over
    rho 5..40; beta=0 spurious floor +0.0013 (1.8 pct of signal), vanishing faster than 1/rho.

## Engine cross-check: FEEC HDiv-VIM -- and the REVISED error budget (2026-07-13)
The SAME chain re-run with the field engine swapped to the FEEC HDiv-VIM (hdiv_scoff_study:
radia.vim.MeshSoftIron on an IRON-ONLY tet mesh -- no air discretization, exact open boundary --
rad.Solve auto dispatch RT1, mid-plane map by ONE batch rad.Fld):
  - dK_in agrees with reduced-Omega to 0.8 pct at MATCHED edge-mesh density (+0.06839 vs
    +0.06896); absolute-dK scatter across engines/meshes is ~+-3 pct (finer HDiv edge mesh
    gives +0.07085); dK_in/model stays 0.92-0.95 in EVERY configuration.
  - ERROR BUDGET (REVISED -- the earlier "iron-mesh error" attribution is RETRACTED): the
    -5..-7 pct deficit vs the x-uniform model is ENGINE-INDEPENDENT real 3D physics.  Mechanism
    measured directly on both engines' maps: the LOCAL ISO-FIELD TILT of the entrance fringe
    near x=0 is only ~0.95-0.96 of the geometric tan(beta) (corner arcs + coil side bars +
    finite pole width).  Hard-edge/SCOFF bookkeeping with the GEOMETRIC edge angle therefore
    overpredicts the vertical edge focusing by ~5 pct for this geometry -- the EFFECTIVE edge
    angle is a magnet-specific quantity that this measurement chain extracts.
  - FRINGE-SHAPE MESH SENSITIVITY: DIAGNOSED AND CURED (2026-07-14 separation runs, JSON
    `mesh_dependence_diagnosis`).  The culprit is the piecewise-constant-M ripple of BULK-size
    GAP-FACE elements at 20 mm standoff (rad.Fld write-back collapses the RT1 solution to
    per-element constant M) -- NOT the edge lines (edge 1.5 mm does not remove the flat-top
    overshoot).  PRESCRIPTION: refine the gap-facing pole faces to face_maxh ~ standoff/3
    (0.006 here) -> overshoot g_max 1.0245->1.0001, K1g 7.84->8.80 mm = reduced-Omega 8.85 to
    0.6 pct (the 12 pct cross-engine K1g gap was entirely this ripple).  LIMITS: keep the
    face/bulk size contrast <= ~2.5x (4 mm faces at 14 mm bulk push the mass-Riesz CG past its
    4000-iter cap, fail-loud); B0 stays RELATIVE-ONLY (+-3.5 pct near-field-evaluation scatter
    across all configurations while the global demag probe is constant to 5 digits; it
    first-order-cancels in the B0-normalized dK).
  - Committed: edge_focusing_fem_results.json `hdiv_vim_cross_check` (runs, agreement, tilt
    probe, mesh_dependence_diagnosis); engines swap via
    fem_scoff_study(solve_midplane=hdiv_solve_midplane); hdiv_* take face_maxh.

## Magnet-design lessons the measurement exposed (all measured failures)
  - A STRAIGHT coil front across a tilted iron edge leaves the iso-field tilt FAR below the
    geometric edge angle (dK deficit 0.55x at beta=20: coil-front misalignment up to
    x_side*sin(beta) ~ 29 mm vs a ~60 mm fringe).  Edge-angle bookkeeping ASSUMES the coil
    follows the pole contour -- and this measurement is sensitive enough to expose it.
  - Rigidly rotating the WHOLE coil to match the entrance edge breaks the MMF linkage topology
    (conductor escapes past the return leg; B0 collapses 3x).  Shear the loop, do not rotate it.
  - A fully-buried winding shunts its MMF in local iron loops (B0 collapse); a long conductor
    overhang leaves a direct-field plateau ("the magnet has no outside").

## Radia implementation status
  - CoilBuilder.mirror() is a true geometric mirror: path points map pointwise, current is
    unchanged, and the magnetic moment follows the axial-vector rule.  Its geometry and field
    symmetry are locked by tests/test_coil_builder_mirror.py.  This study still builds both
    loops explicitly so its reference geometry does not depend on that helper.
  - rad.RadiaField on a TrfOrnt-wrapped CONTAINER crashes the process with 0xC0000374 heap
    corruption (rad.Fld on the same container is fine; TrfOrnt on PRIMITIVES inside a container
    is fine -- CoilBuilder.to_radia relies on that).  Bake transformations into segment
    geometry instead.
Regression locks: tests/test_coil_builder_mirror.py and
validation_test/feec/test_edge_focusing_tracking.py::test_fem_coil_pair_is_clean.
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

THE FAST-TWIST LEAF COUPLING is DONE (twist_rate_leaf_coupling.py): the exact
helical multipole Phi_n = I_n(n k r) sin(n(theta - k s)) (Laplace + helical
symmetry => modified Bessel) deviates from the 2-D stack as eps ~ (ka)^2
(transverse focusing error, slope 2.04) with a longitudinal B_s ~ ka (slope 0.97),
ka = 2 pi a / P the twist per aperture.  As k -> 0 it IS the per-station 2-D stack.
Threshold eps = 1% at pitch/aperture ~ 46 -- the SAME "longitudinal >> transverse
by ~40x" rule as the straight magnet's foliate_perturb (L/gap ~ 40), with the
twist replacing gap/L by a/P.  So the twist axis closes: the twist (n-fold law) ->
the combined-function confluence on a curved orbit -> its validity threshold.

REMAINING: a spiral sector (pole twist phi != orbit bend theta) and an s-ramped
(b1, b2)(s) are the extensions the n-fold law governs.
"""


ENDPACK_TWO_PLANE = """
# The magnet END PACK in two planes -> 3-D (x-y cross-section + s-y end)

This is the two-plane -> 3-D method (see two_plane_design) applied to the magnet
END PACK -- the longitudinal TERMINATION of a STRAIGHT magnet -- as opposed to a
curved FFAG sector cell.  Convention: beam = y (the arc length s of a straight
magnet), gap = z (pole faces z = +-g/2), width = x.  So the design planes are the
literal (x-y) transverse cross-section and the (s-y) longitudinal plane.

NOTE: this is a DIFFERENT "end pack" from isochronous_endpack (which is the RADIAL
field-index end pack of an FFAG, scaling vs isochronous k(r)).  This one is the
longitudinal magnet-END termination.

## The method: two cheap 2-D DESIGN solves, then a 3-D REFLECTION

Unlike the 3-D end study (accel_pole_ends_fem, which reads the end equipotential
OUT of the finished 3-D solve), here BOTH planes are standalone 2-D DESIGN solves
done FIRST; the 3-D solve only reflects + verifies.

  Plane 1 (x-y cross-section) -- the transverse multipole:
    a finite flat pole droops (b3/b1 ~ -3.6e-5); the concave shim
    z = g/2 - delta (x/w)^2 with delta ~ 0.41 mm zeroes it
    (accel_pole_dipole_body_2d.solve).

  Plane 2 (s-y longitudinal) -- the end chamfer:
    a STANDALONE 2-D Laplace fringe (the parallel-plate end, Dirichlet pole
    terminating at the iron end, midplane = 0) -> the interior-gap equipotential
    Phi = Phi_ref bows up past the pole edge = the ROGOWSKI end chamfer shape
    ghat(s); and L_eff ~ 151 mm, a +26% excess over the 120 mm iron (each end
    ~0.75 g) -- endpack_two_plane.solve_sy_endpack.

  Reflection (3-D) -- equipotential-pole drive (PURE LAPLACE, no coil):
    an upper-half pole at Psi = mmf, midplane Psi = 0 (the same high-mu
    equipotential drive as the FFAG rung-2; B = -grad Psi is a plain CF, so the
    whole solve + the integrated multipole analyzer are CHEAP -- no RadiaField).
    Sweeping the chamfer DEPTH with the ghat(s) shape drives the hard-cut
    pole-tip corner over-field (+11%) THROUGH ZERO at ~2.1 mm, while the
    integrated transverse b3,5 stays ~0.3% -- body/Plane-1 dominated (the END
    shape is the wrong lever for it, the honest two-lever split).

## Verified (endpack_two_plane.py, ngsolve only, golden-tested ~8 s --fast)

  Plane 1 : b3/b1 ~ -3.6e-5 (flat droop), shim delta ~ 0.41 mm (zeroes b3).
  Plane 2 : L_eff ~ 151 mm (+26% excess), a real Rogowski bow-out ~21 mm/30 mm.
  3-D     : flat pole-tip corner over-field ~ +11%; depth sweep [0, 2.4, 5, 8 mm]
            -> tip over-field [+11, -2, -15, -26]% (monotone), designed depth
            ~2.1 mm (zero over-field); integrated b3,5 ~0.3% throughout.
  Cross-check: the cheap 2-D s-y chamfer SHAPE predicts the 3-D end
            equipotential bow-out to ~7% rms.

## Why equipotential-pole drive (not a coil)

A CoilBuilder + RadiaField (Biot-Savart) source makes the reduced-Omega assembly
and the per-point field readout SERIAL and ~1000x slower (the analyzer + contour
do thousands of RadiaField evals).  The high-mu iron pole IS an equipotential, so
driving the pole face at Psi = mmf (pure Laplace) gives the SAME gap-field SHAPE
(the field index / multipole content is set by the pole GEOMETRY, drive-
independent -- same argument as two_plane_design rung 2) at a fraction of the
cost.  The coil only sets the field AMPLITUDE (Tesla calibration), not the end
geometry the chamfer designs.

## Honest scope

The 3-D reflection carries the s-y chamfer at a FIXED body pole width; the x-y
shim (delta) is verified as the transverse lever but not yet co-baked into the
same 3-D loft.  A single 3-D pole surface carrying BOTH z = g/2 - delta (x/w)^2
(x-y shim) AND + lift(s) (s-y chamfer) -- a tensor-product loft -- is the next
refinement.  The transverse b3,5 being body-dominated (and barely moved by the
end shape) is the same two-lever split established in accel_pole_ends_fem /
accel_pole_dipole_body_2d.
"""


SPECTROMETER_ENDPACK_SATURATION = """
# The SPECTROMETER end pack, NONLINEAR -- the pole-tip corner is a saturable throat

A large BENDING SPECTROMETER dipole runs near the iron knee, so its end pack must
be designed WITH saturation.  This composes the linear end pack (endpack_two_plane)
with the iron-saturation lever (clebsch_dipole_saturation_2d/3d) on the magnet END.

## The physics (why a spectrometer end needs the nonlinear design)

The linear (high-mu equipotential) end pack found the hard-cut pole END concentrates
flux at the tip CORNER by

    B_corner = kappa * B_gap ,   kappa = tip_enhancement ~ 1.11   (geometry-only).

In a SATURATING iron, kappa > 1 means the CORNER reaches the iron knee FIRST: it
saturates at a gap field

    B_gap_knee = B_K / kappa ~ 1.5 / 1.13 ~ 1.33 T,

~12% BELOW the bulk iron knee B_K = 1.5 T.  Above that the corner mu_r collapses, the
flux can no longer follow the pole edge, and the EFB (effective field boundary
~ L_eff) DRIFTS with excitation.  That is fatal for a spectrometer: the pole-edge
angle's EDGE FOCUSING (the vertical focusing tan(beta)/rho) depends on the EFB, so a
drifting EFB means the optics CHANGE with the field setting.

## The design (the Rogowski chamfer is the corner-throat width knob)

The corner is exactly a Chaplygin saturable THROAT (clebsch_dipole_saturation_2d):
kappa is its inverse cross-section, and it saturates first.  The Rogowski end chamfer
(endpack_two_plane's s-y design) WIDENS that throat -- lowers kappa -- so the corner
knee B_K/kappa RISES toward the bulk knee.  The SAME chamfer that zeroes the LINEAR
corner over-field (the cosmetic field-quality lever) REMOVES the PREMATURE corner
saturation (the hard engineering lever) and keeps the EFB -- and the edge focusing --
STABLE up to the bulk limit.  Linear and nonlinear levers POINT THE SAME WAY;
saturation gives the chamfer its hard justification.

## Verified (endpack_spectrometer_saturation.py, golden-tested, ngsolve only)

1. The MAP at LINEAR cost: reuse the linear equipotential corner concentration
   kappa(chamfer) (endpack_two_plane's depth sweep) and overlay the Froehlich iron BH
   (B_K=1.5 T, mu_r0=2000, from clebsch_dipole_saturation_3d): the corner knee
   B_K/kappa(chamfer) -- flat 1.33 T -> 2.4 mm chamfer 1.55 T -> 5 mm 1.75 T.  To
   clear a B_op = 1.45 T operating field without corner saturation needs kappa <= 1.034
   (~1.4 mm chamfer; the linear cosmetic optimum kappa=1 is ~1.9 mm -- the same lever).
   The whole nonlinear end-pack map = 4 linear equipotential solves + a BH overlay --
   the Chaplygin "nonlinear analysis done linearly" applied to the END corner.

2. DESIGN-GRADE, components validated: the lumped kappa-throat overlay is the same
   lumped-magnetic-circuit class as clebsch_dipole_saturation_2d (~10% vs FEM), and
   its two ingredients are independently verified elsewhere:
     - the corner concentration kappa = the LINEAR equipotential tip_enhancement
       (endpack_two_plane.py, golden-tested) -- geometry-only, drive-agnostic;
     - the Froehlich BH + the well-conditioned B-input A-formulation that backs the
       iron saturation = clebsch_dipole_saturation_3d.py (the documented cure: the
       reduced-Omega mu(|H|) Picard STALLS at high mu, the A-formulation does not).
   So the composition B_corner = kappa*B_gap until B_K is well-founded.

## Honest scope

A fully coil-driven 3-D corner-saturation FEM is the documented expensive extension,
NOT run in the example: the equipotential/MMF drive forces flux across the gap (a
UNIFORM applied field does NOT reproduce the corner concentration kappa -- measured),
and the coil's Biot-Savart B_s projection is the serial bottleneck.  The corner kappa
softens before the hard knee (a real FEM is the truth).  The curved/rotated-EFB edge
focusing (the horizontal x-s edge contour) and the fully-saturating sector body are
the remaining spectrometer extensions (see endpack_two_plane + two_plane_design).
"""


ENDPACK_COBAKE = """
# The two-plane end pack CO-BAKED into ONE pole face

The completion of endpack_two_plane: that example's 3-D reflection carried the s-y
chamfer at a FIXED body width (the x-y shim was verified as the transverse lever, not
baked into the same 3-D pole).  endpack_cobake.py bakes BOTH into one gap face

    z_face(x, s) = g/2 - delta (x/w)^2  +  lift(s)
                    \\___ x-y shim ___/    \\_ s-y _/   (a tensor-product face),

and shows BOTH levers act AT ONCE.  The 4 cases (the same equipotential-pole drive +
integrated analyzer as endpack_two_plane._solve_3d_endpack):

  baseline (flat cut):  corner tip 1.16 (over-fields),  transverse b_3,5 ~ 0.7%
  shim only   (delta):  corner tip 1.23,                transverse b_3,5 ~ 0.07%  (x-y lever)
  chamfer only (ghat):  corner tip 0.96 (rounded),      transverse b_3,5 ~ 0.3%   (s-y lever)
  BOTH (co-baked):      corner tip 1.02 (rounded)  AND  transverse b_3,5 ~ 0.07%

i.e. ONE pole face delivers a clean integrated transverse harmonic AND a rounded
pole-tip corner -- the two cleanly-separated two-plane levers (delta ~ 0.41 mm from
Plane 1; ghat the Rogowski shape from Plane 2) composed in 3-D.

## Honest scope of the staircase build (endpack_cobake.py)

The exact delta(x/w)^2 shim needs an x-VARYING face, built there as an x-prism
STAIRCASE (per-slab shim offset), so the no-shim cases mesh coarser than the shim
cases -- the per-case ABSOLUTE numbers are research-grade, not precision.  The locked
claim is the CO-EXISTENCE of both levers in the (well-resolved) BOTH pole; the per-lever
CAUSATION is golden-locked separately (the x-y shim zeroes b_3: accel_pole_dipole_body_2d
/ endpack_two_plane Plane 1; the s-y chamfer drives the corner over-field through 1:
endpack_two_plane's depth sweep).

## Precision construction (endpack_cobake_loft.py, OCC ThruSections, golden-tested)

The clean construction the staircase pointed to.  The gap face z(x,s)=g/2-delta(x/w)^2
+lift(s) is a SMOOTH OCC ThruSections LOFT through per-x-station cross-section wires
(each carrying its shim offset delta(x_i/w)^2 + the chamfer lift(s)), so the surface is
smooth in x (no facets).  The headline is MESH CONSISTENCY:

    smooth LOFT       ne(shim)/ne(baseline) ~ 0.97  (same density)
    x-prism STAIRCASE ne(shim)/ne(baseline) ~ 36    (merges delta=0 slabs -> coarse;
                                                     steps delta>0 -> fine)

so the loft RESOLVES the documented staircase artifact, and the co-baked pole's b_3,5 +
corner become a PRECISION claim.  On the consistent mesh both levers still act: the
chamfer rounds the corner (tip 0.99), the shim REMOVES the transverse content the
chamfer introduces (both 0.47% < chamfer-only 0.84%) and returns it to the baseline
mesh-noise floor (~0.5%).  (Absolute b_3,5 differs from the staircase table above
because the meshes differ; the loft's value is the precision one.)

## Rotated-EFB EDGE FOCUSING -- characterized NEGATIVE, not shipped

The rotated-edge EDGE FOCUSING extension was retried this session and is recorded as a
CHARACTERIZED NEGATIVE (not a shipped example).  A first attempt read the INT B
effective-field-boundary angle out of the EQUIPOTENTIAL drive and it attenuated to
~0.47*beta_cut (unattributed).  The retry used the genuine forward engine -- a rigidly
rotated WHOLE magnet (iron AND the CoilBuilder coil rotated together by beta) driven by
reduced-Omega + Biot-Savart, with the EFB read both as INT B_z dy / B_z(body) and as the
half-field crossing y_half(x).  Across all variants (parallelogram vs rigid rotation,
integral vs half-field EFB, narrow vs wide pole) the result is robust: at beta=0 the EFB
slope is cleanly ~0 (unbiased), but for beta>0 the field-EFB slope does NOT recover the
geometric edge angle (wrong sign, many times tan(beta); at larger beta the per-line
B_z(body) normalization passes through zero).

Attribution: the per-beam-line field integral INT B_z dy through a tilted FINITE magnet is
NOT a local edge tracker -- the compact fringe is fully 3-D and the field-EFB slope simply
is not the edge angle (the prior ~0.47 was the same surrogate failing, not the drive).  The
genuine edge-focusing strength is a TRAJECTORY quantity (the vertical kick INT (v x B) along
particle orbits through the fringe) and needs PARTICLE TRACKING, not a field-EFB slope.
Per the repository honest-results policy it is kept as a characterized open problem, not
shipped as a working example.
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
      "sector_saturation"    - The saturating SECTOR body: azimuthal L_eff ROBUST
                               (gap-dominated) vs radial k(r) FRAGILE -- planes differ
      "beam_referenced_twist" - The equipotential SURFACE as design primitive
                               + the twist (n-fold law, twisting quadrupole)
      "endpack_two_plane"    - The magnet END PACK in two planes (x-y cross-
                               section + s-y end) -> 3-D equipotential reflect
      "spectrometer_endpack" - The SPECTROMETER end pack NONLINEAR: the pole-tip
                               corner = saturable throat, B_K/kappa knee, EFB drift
      "endpack_cobake"       - The two-plane end pack CO-BAKED into one pole face
                               (x-y shim + s-y chamfer): clean b_3,5 AND rounded corner
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
    if topic in ("sector_saturation", "saturating_sector", "sector_sat",
                 "ffag_saturation", "gap_reluctance", "leff_robust"):
        return SECTOR_SATURATION
    if topic in ("bending_endpack_saturation", "bending_endpack", "bending_end",
                 "endpack_opt", "corner_relief", "saturation_invariant",
                 "endpack_saturation"):
        return BENDING_ENDPACK_SATURATION
    if topic in ("excitation_invariant_field", "excitation_invariant", "excitation_sweep",
                 "flux_line_invariant", "flux_line_invariance", "same_flux_lines",
                 "invariant_flux_lines"):
        return EXCITATION_INVARIANT_FIELD
    if topic in ("hodograph_feasibility", "feasibility", "limit_line", "edge_sharpness",
                 "field_space_design", "hodograph_2d_linear", "inverse_design"):
        return HODOGRAPH_FEASIBILITY
    if topic in ("hodograph_bending_sy", "bending_sy", "sy_plane", "end_shaping",
                 "pole_end_design", "chamfer_optimization", "rogowski_recover",
                 "longitudinal_fringe"):
        return HODOGRAPH_BENDING_SY
    if topic in ("edge_focusing_tracking", "edge_focus", "edge_tracking",
                 "tilted_edge", "vertical_focusing", "hill_integral", "beta_edge",
                 "edge_focusing_tan"):
        return EDGE_FOCUSING_TRACKING
    if topic in ("beam_referenced_twist", "twist", "twisting", "design_primitive",
                 "equipotential_surface", "n_fold", "rotating_gradient"):
        return BEAM_REFERENCED_TWIST
    if topic in ("endpack_two_plane", "end_pack_two_plane", "magnet_end",
                 "end_chamfer", "rogowski", "endpack_2plane"):
        return ENDPACK_TWO_PLANE
    if topic in ("spectrometer_endpack", "spectrometer", "nonlinear_endpack",
                 "corner_saturation", "saturable_endpack", "efb", "edge_focusing"):
        return SPECTROMETER_ENDPACK_SATURATION
    if topic in ("endpack_cobake", "cobake", "co_bake", "tensor_pole",
                 "shim_chamfer", "both_planes", "endpack_cobake_loft", "cobake_loft",
                 "tensor_loft", "loft", "thrusections"):
        return ENDPACK_COBAKE
    if topic == "all":
        return "\n\n".join([
            END_POLE_DESIGN, KOLKATA_CYCLOTRON, ROTATING_COIL_MEASUREMENT,
            ISOCHRONOUS_ENDPACK_DESIGN, FOLIATE_PERTURB, TWO_PLANE_DESIGN,
            SECTOR_SATURATION, BENDING_ENDPACK_SATURATION, EXCITATION_INVARIANT_FIELD,
            HODOGRAPH_FEASIBILITY, HODOGRAPH_BENDING_SY, EDGE_FOCUSING_TRACKING,
            BEAM_REFERENCED_TWIST, ENDPACK_TWO_PLANE, SPECTROMETER_ENDPACK_SATURATION,
            ENDPACK_COBAKE,
        ])
    return (
        f"Unknown topic '{topic}'. Available: all, end_pole, kolkata, "
        "rotating_coil, isochronous_endpack, foliate_perturb, two_plane_design, "
        "sector_saturation, bending_endpack_saturation, excitation_invariant_field, "
        "hodograph_feasibility, hodograph_bending_sy, edge_focusing_tracking, "
        "beam_referenced_twist, endpack_two_plane, spectrometer_endpack, endpack_cobake."
    )
