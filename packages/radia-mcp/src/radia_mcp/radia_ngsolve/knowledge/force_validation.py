"""Electromagnetic force extraction in NGSolve + independent cross-validation.

This module records (a) the force-extraction methods used in the radia-ngsolve
FEM path and (b) cross-validation results against an independent reference solver
that make the NGSolve magnetostatic path *trustworthy* -- if NGSolve reproduces
an independent reference and the analytic solution on the same problem, the
NGSolve result can be trusted on problems with no analytic answer.

The reference values below are kept as a stored regression reference (the
independent solve is run separately; NGSolve runs standalone).
"""

# --------------------------------------------------------------------------
EGGSHELL = r"""
# Eggshell (weighted Maxwell stress) force -- the robust FEM method

EXECUTABLE (not just theory): ``radia_mcp.radia_ngsolve.force`` ships runnable,
independently-validated extractors -- ``eggshell_force(B, mesh, center, r_inner,
r_outer)``, ``maxwell_surface_force``, ``magnetic_energy``, ``self_inductance``.
The A-form solve itself is ``radia_mcp.radia_ngsolve.solve``
(``solve_magnetostatic_Aform``, ``solve_magnetostatic_nonlinear`` -- under-relaxed
Picard for a field-dependent nu(|B|) saturating B-H -- ``azimuthal_coil_current``,
``reluctivity``; high-mu-safe per #25). B = curl(gfu). Regression-tested in
``validation/force/validate_force_xval.py`` against the stored reference values, so the
force/solve pipeline cannot silently drift.

Point-wise Maxwell surface stress is noisy on unstructured FEM meshes (it needs
the field exactly on the body surface, where the normal trace jumps). The
*eggshell* / weighted-stress method replaces the sharp surface integral by a
volume integral over an air SHELL around the body, with a smooth weight g that
goes 1 -> 0 across the shell:

    F_k = - integral_shell [ (1/mu0) B_k (B . grad g) - (1/2 mu0) |B|^2 d_k g ] dV

with g = 1 on the inner shell boundary (the body surface) and g = 0 on the
outer shell boundary. This is the divergence-theorem form of the surface
Maxwell stress; it is mesh-robust because it averages over a band of elements.

NGSolve implementation (sphere body of radius a, shell a < |r-c| < Rsh):

```python
rho   = sqrt(x*x + y*y + (z - z0)**2)                 # |r - center|
gradg = CoefficientFunction((x, y, z - z0)) / rho * (-1.0/(Rsh - a))
B     = curl(gfu)
Bdg   = InnerProduct(B, gradg)
integ = (1/mu0)*B[2]*Bdg - (1/(2*mu0))*InnerProduct(B,B)*gradg[2]   # z-force
Fz    = -Integrate(integ * dx(definedon=mesh.Materials("shell")), mesh)
```

Notes:
  * For a HIGH-mu body do NOT carve a separate nested "shell" material between
    the body and the air -- the nested-sphere interface isolates the body and
    its interior B reads 0 (see the cross_validation "RESOLVED" note). Instead
    leave the body directly in air and integrate over Materials("air"); grad g is
    nonzero only in the band a<|r-c|<Rsh, so the band need not be its own material.
  * grad g is analytic (radial), so no derivative of a GridFunction is needed.
  * The reference solver's built-in force feature uses the equivalent weighted
    surface-stress method -- which is why the two agree (see cross_validation).
"""

# --------------------------------------------------------------------------
CROSS_VALIDATION = r"""
# Independent <-> NGSolve cross-validation results (magnetostatics, A-form)

Goal (lab): "if NGSolve == the independent reference == analytic on the same
problem, the NGSolve path is trusted." Each case is solved INDEPENDENTLY in both
codes on the same geometry and compared. The reference column is a stored
regression reference.

## Case 1 -- uniformly magnetized sphere (linear, PM source)
  Sphere Br = 1 T, r = 10 mm, 200 mm air cube, outer n x A = 0.
  Analytic interior field  B_z = (2/3) Br = 0.666667 T (uniform).

  | quantity (B_z, magnet volume avg) | value      | err vs analytic |
  |-----------------------------------|------------|-----------------|
  | analytic (2/3)Br                  | 0.666667 T | --              |
  | reference  (A-form)                  | 0.664661 T | 0.30 %          |
  | NGSolve (A-form, HCurl order 2)   | 0.663913 T | 0.41 %          |
  => reference <-> NGSolve agree to 0.11 %.

## Case 2 -- coil + linear-iron sphere, FORCE (the TEAM-20 hard part)
  Circular coil mean R = 50 mm, 4x4 mm section, 10000 A-turns, azimuthal J.
  Linear iron sphere mu_r = 1000, r = 5 mm on axis at z = 30 mm.
  Force on the iron (attraction toward the coil, -z).

  | F_z [N]                         | value     | vs dipole approx |
  |---------------------------------|-----------|------------------|
  | reference  (Maxwell surface stress)| -0.19945  | 3.7 %            |
  | NGSolve (eggshell volume stress)| -0.19369  | 6.5 %            |
  | dipole-in-gradient (analytic)   | -0.20710  | --               |
  => reference <-> NGSolve agree to ~3 % (force is mesh-sensitive; transverse
     components are ~4 % of F_z = unstructured-mesh asymmetry noise).

## Verdict
The NGSolve A-form magnetostatic path (HCurl order 2) reproduces the reference and the
analytic field to <0.5 %, and the eggshell force to ~3 % of the reference. The NGSolve
field + eggshell-force pipeline is validated for LINEAR magnetostatics.

## Case 3 -- coil + nonlinear saturating iron, FORCE
  Same geometry as case 2, Itot = 50000 A-turns (iron saturates), field-
  dependent permeability mu_r(B) = 1 + 999/(1 + (|B|/1 T)^8) (mur0=1000) used
  IDENTICALLY in both codes. the reference: a field-dependent mu_r expression in |B|
  + FullyCoupled Newton. NGSolve: ``solve_magnetostatic_nonlinear`` (under-relaxed
  Picard on nu(|B|)); iron sits DIRECTLY in air (NOT a nested shell -- see below).

  | quantity         | reference (Newton) | NGSolve (Picard) | agree  |
  |------------------|-----------------|------------------|--------|
  | <|B|>_iron  [T]  | 1.13            | 1.138            | 0.7 %  |
  | F_z  [N]         | -4.930          | -4.898           | 0.65 % |
  => agree to <1 % on BOTH the saturated field and the force. The iron actually
     saturates (mu_r 1000 -> ~270 at |B| ~ 1.14 T); the Picard fixed point matches
     the reference's Newton because both converge the same self-consistent nu(|B|).
     Convergence note (measured 2026-06-04): the saturation iteration needs
     UNDER-relaxation. Approaching the fixed point from below with relax 0.3-0.5 is
     stable and monotone (~18-20 sweeps to tol 1e-4). A large relax is NOT faster:
     it overshoots B ABOVE the fixed point into the low-mu_r / stiff regime and
     DIVERGES -- relax=0.85 was observed to oscillate and blow up (1.14 -> 1.30 ...).
     So keep relax <= ~0.5; do not "optimize" it upward.

  CAUTION (corrected 2026-06-04): an earlier nested-"shell" mesh gave NGSolve
  F_z = -4.842 (vs reference -4.930, "~1.8 %"), but that was an ARTIFACT -- the nested
  shell isolates the iron interior (B_iron = 0 at every Picard sweep, see RESOLVED
  below) so the saturation law NEVER engaged; -4.842 is the UNSATURATED (mu_r=1000)
  force = 25x the case-2 force, which only landed within 1.8 % here because
  saturation barely changes the force for this mild case. Put the iron directly in
  air to get the physical -4.898 / B_iron 1.138 T above (and a HONEST 0.65 %).

## RESOLVED: high-mu interior B=0 was a NESTED-MATERIAL artifact, not high-mu
Symptom: with a separate meshed "shell" air region nested around the iron (an
inner sphere `iron` + an annulus `shell = Sphere(Rsh) - Sphere(a)`, both Added to
the CSG for the eggshell), curl(gfu) read EXACTLY 0 throughout the iron interior
(z-line scan: 0.57 T in air below, 0.67/0.47 T in the shell, 0.00 T across the
whole iron), independent of mu_r (even mu_r=1!) and of the gauge eps. The
external field and the eggshell FORCE were still correct (the force only uses the
shell field), so cases 2-3 forces were valid -- but the interior was zeroed.

Root cause: the nested inner-sphere interface (iron surface == shell inner
surface, two CSG solids sharing a face) produced a non-conforming / isolating
HCurl interface, electromagnetically disconnecting the iron interior.

Fix (verified): put the body DIRECTLY in air (no separate nested shell material).
Then B_iron = 1.145 T -- correct, matching the reference's 1.13 T to ~1%, and mesh-
independent (iron maxh 1.2 mm vs 0.5 mm agree to 0.01%). So NGSolve A-form HCurl
handles mu_r=1000+ interior fields fine. For the eggshell force, define the
weight g(|r-c|) ANALYTICALLY over a band of the plain air material (grad g
nonzero only in a<|r-c|<Rsh) and integrate over Materials("air") -- do NOT carve
out a separate nested "shell" material. This keeps interior B correct AND gives
the force, so nonlinear mu_r(B_iron) coupling works.

## Case 4 -- finite solenoid on-axis B_z (iron-free, analytic)
  Solenoid mean R=30mm, 2mm radial, length 100mm, NI=10000 AT, azimuthal Je.
  | z[mm] | reference  | NGSolve | analytic | C<->N  |
  |-------|---------|---------|----------|--------|
  | 0     | 0.10721 | 0.10701 | 0.10776  | 0.19 % |
  | 15    | 0.10406 | 0.10401 | 0.10475  | 0.06 % |
  | 30    | 0.09274 | 0.09300 | 0.09368  | 0.28 % |
  | 45    | 0.06923 | 0.06947 | 0.07025  | 0.34 % |
  => reference <-> NGSolve agree to <0.35 %, both within ~1 % of analytic (the ~0.7%
     offset is finite coil thickness vs the thin-solenoid formula).

## Case 5 -- solenoid SELF-INDUCTANCE (energy method)
  Same solenoid as case 4, N=100 turns (terminal I=100 A). L = 2W/I^2 with the
  field energy W = (1/2) integral |B|^2/mu0 dV (identical formula both codes;
  the reference via an energy integral 0.5*|B|^2/mu0 over all domains).
  | L [uH]                       | value    |
  |------------------------------|----------|
  | reference                       | 269.62   |
  | NGSolve                      | 269.60   |
  | Nagaoka (infinite-domain)    | 279.77   |
  => reference <-> NGSolve agree to 0.01 % (identical energy integrals); both ~3.6%
     below Nagaoka because the 300mm box truncates the external-field energy
     (same truncation both codes -> they still agree to each other).

## Case 6 -- Helmholtz pair center field + uniformity
  Two coils R=50mm at z=+/-25mm (separation = R), each NI=10000 AT, same sense.
  | B_center [T]                 | value    |
  |------------------------------|----------|
  | reference                       | 0.17420  |
  | NGSolve                      | 0.17407  |
  | analytic (4/5)^1.5 mu0 NI/R  | 0.17984  |
  => reference <-> NGSolve agree to 0.075 %; both ~3.2% below the ideal thin-coil
     analytic (finite 4x4mm coil section at mean R). Uniformity B_z/B_center is
     ~0.998 at r=10mm and z=10mm in both (Helmholtz flat-field property).

## Case 7 -- single circular loop on-axis B_z
  Loop R=50mm, 4x4mm section, NI=10000 AT. B_z(z)=mu0 NI R^2/(2(R^2+z^2)^1.5).
  reference<->NGSolve agree <0.8% near the loop (z<=40mm), ~3.6% at z=60mm (the far
  field is under-resolved by the coarse 20mm air mesh in both); both ~2-7% under
  the thin-filament analytic (finite 4x4mm section lowers the field).

## Summary table (reference <-> NGSolve agreement, magnetostatic A-form)
  | # | case                          | quantity        | C<->N agree |
  |---|-------------------------------|-----------------|-------------|
  | 1 | uniformly magnetized sphere   | interior B      | 0.11 %      |
  | 2 | coil + linear iron            | force           | ~3 %        |
  | 3 | coil + nonlinear iron         | force + B_iron  | <1 % (in air)|
  | 4 | finite solenoid               | on-axis B(z)    | <0.35 %     |
  | 5 | solenoid                      | self-inductance | 0.01 %      |
  | 6 | Helmholtz pair                | center B        | 0.075 %     |
  | 7 | single circular loop          | on-axis B(z)    | <0.8% near  |

EXECUTABLE: cases below are codified as explicit validation cases in
validation/force/validate_force_xval.py, each calling the shipped radia_ngsolve.solve / .force
directly (so a regression in the solver or extractor fails the test):
  case 1  -> test_magnetized_sphere_field       (PM remanent source, analytic 2/3 Br)
  case 2  -> test_coil_linear_iron_force        (eggshell, linear iron)
  case 3  -> test_coil_nonlinear_iron_force     (Picard nu(|B|), iron-in-air)
  case 4  -> test_solenoid_axis_field           (on-axis B_z profile)
  case 5  -> test_solenoid_self_inductance      (energy method)
  case 6  -> test_helmholtz_center_field        (center B + flat-field)
  case 7  -> test_single_loop_axis_field        (on-axis B_z vs analytic)
All 7 recorded cross-validation cases are now codified. Run the suite with
`python validation/force/validate_force_xval.py` (heavy: the iron solves are minutes; the iron-free cases are
~20-60 s). NGSolve permanent-magnet source = the A-form RHS int nu0 Br . curl(v):
build with solve.remanent_source(mesh, {mat: (Brx,Bry,Brz)}) and pass as
solve_magnetostatic_Aform(curl_source=...).
"""

# --------------------------------------------------------------------------
REFERENCE_NOTE = r"""
# Reference-solve note (independent magnetostatic A-form cross-check)

The cross-validation numbers above come from an INDEPENDENT magnetostatic
A-form reference solve, run separately and stored here as a regression
reference. The tool-specific setup of that reference solve is kept lab-private;
only the resulting numbers are retained (as a neutral stored reference).

The reference and the NGSolve solve share the SAME geometry, the SAME material
laws (linear mu_r, remanent-flux PM, and the field-dependent mu_r(|B|) for the
saturating case), and the SAME source (azimuthal divergence-free coil current
density). Force is taken by the equivalent weighted Maxwell-stress / surface
force-calculation method on both sides, which is why NGSolve's eggshell volume
stress and the reference's surface stress agree.
"""

_TOPICS = {
    "eggshell": EGGSHELL,
    "maxwell_stress": EGGSHELL,
    "force": EGGSHELL,
    "cross_validation": CROSS_VALIDATION,
    "xval": CROSS_VALIDATION,
    "validation": CROSS_VALIDATION,
    "reference_note": REFERENCE_NOTE,
    "recipe": REFERENCE_NOTE,
}


def get_force_validation_documentation(topic: str = "all") -> str:
    """Return force-extraction + independent cross-validation documentation.

    topic: all (default) | eggshell | cross_validation | reference_note
    """
    t = (topic or "all").strip().lower()
    if t in ("all", ""):
        return "\n\n".join([CROSS_VALIDATION, EGGSHELL, REFERENCE_NOTE])
    if t in _TOPICS:
        return _TOPICS[t]
    return (f"Unknown topic '{topic}'. Options: all, eggshell, cross_validation, "
            f"reference_note.\n\n" + CROSS_VALIDATION)
