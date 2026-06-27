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
FORCE_METHOD_MAP = r"""
# Electromagnetic force method map for radia-ngsolve

This is the index page for force/torque work in ``radia_mcp.radia_ngsolve``.
Use it before adding a new force extractor: most practical cases reduce to one
of the methods below, with an analytic sanity check available.

## Which method to use

| Situation | Primary method | radia-ngsolve entry point |
|---|---|---|
| Body in air, global 2D/3D force | Weighted Maxwell stress | ``force.eggshell_force*`` |
| Rotating machine torque from field solution | Weighted Maxwell stress torque | ``force.eggshell_torque*`` |
| 2D sector torque to whole-machine torque | Stack/symmetry scaling | ``machine_scaling.MachineScaling``, ``machine_scaling.torque_scaling_summary`` |
| Simple closed integration surface in air | Maxwell surface stress | ``force.maxwell_surface_force*`` |
| One surface patch / sign convention teaching | Local Maxwell traction | ``force.maxwell_stress_tensor_air`` / ``force.maxwell_traction_summary`` |
| First-order `.vol` boundary triangle force trace | P1 triangle traction load / Maxwell traction / boundary pressure/vector-traction rows | ``force.surface_triangle_constant_traction_load_summary``, ``force.surface_triangle_maxwell_traction_summary``, ``NetgenTriTetVolMesh.boundary_pressure_force_moment_rows``, ``NetgenTriTetVolMesh.boundary_traction_force_moment_rows`` |
| CAD-side analytic box face loads | build123d face pressure / vector traction rows | ``build123d.modeling.box_face_pressure_moment_rows``, ``build123d.modeling.box_face_traction_moment_rows`` |
| Current-carrying conductor force | Lorentz volume integral | ``force.lorentz_force_2d``, ``force.planar_lorentz_force_summary``, ``force.parallel_wire_lorentz_force_summary`` |
| Discrete force rows to net force/torque | Resultant and pivot moment sum | ``force.force_moment_resultant_summary`` |
| Axisymmetric actuator / coil force | Axisymmetric weighted stress | ``force.eggshell_force_axi`` |
| Uniform air-gap holding force | Magnetic pressure | ``force.air_gap_*`` and ``solve.magnetic_circuit_gap_force`` |
| Uniform / sampled cylindrical air-gap torque | Maxwell shear stress | ``force.air_gap_shear_*``, ``force.air_gap_shear_torque_from_angle_samples`` |
| RF beam / waveguide radiation pressure | Time-averaged Poynting momentum flux | ``force.radiation_pressure_*``, ``force.poynting_patch_force_summary``, ``force.time_average_maxwell_*`` |
| Energy/inductance-derived checks | Field energy / virtual work | ``force.magnetic_energy*``, ``force.inductance_*``, ``force.virtual_work_force_*`` |
| dq motor operating torque | Lumped dq model | ``solve.dq_torque`` and companions |
| Synchronous/induction machine torque curves | Circuit-level machine model | ``solve.synchronous_power_angle_torque`` / ``solve.induction_machine_torque`` |
| MEMS/electrostatic force | Electric Maxwell stress / capacitance gradient / weighted stress | ``electrostatics.parallel_plate_capacitor_energy_force``, ``electrostatics.capacitance_gradient_force_summary``, ``force.electrostatic_traction_summary``, ``force.electrostatic_eggshell_force*`` |

## Core identities

Magnetostatic Maxwell stress in air:

```text
T = (1/mu0) (B tensor B - 0.5 |B|^2 I)
F = integral_S T n dS
```

Weighted-stress / eggshell form replaces a sharp surface integral by a
volume-band integral in air:

```text
F_k = - integral_band [(1/mu0) B_k (B.grad g)
                       - (1/(2 mu0)) |B|^2 d_k g] dV
```

Use this when possible: it averages over elements, avoids surface-trace jumps,
and is the most robust extractor for unstructured FEM meshes.

Lorentz force on an imposed-current conductor:

```text
F = integral J x B dV
```

For 2D out-of-plane current ``Jz`` and in-plane ``B=(Bx, By)``:
``Fx = -int Jz By dA`` and ``Fy = int Jz Bx dA``.  Pass the total field; the
self-field integrates to zero by symmetry for the net conductor force.

For two long parallel wires, put wire 1 at the origin and wire 2 at separation
vector ``r`` in the xy plane.  With positive current along ``+z``:

```text
B_1(r) = mu0 I1 / (2 pi |r|) (zhat x rhat)
F_2/L = I2 zhat x B_1 = -mu0 I1 I2 / (2 pi |r|) rhat
```

Like currents attract, so the right-hand wire feels a force back toward wire 1.
This is the signed/directional gate behind
``force.parallel_wire_lorentz_force_summary`` and the scalar
``solve.two_wire_force_per_length``.

Discrete force rows to force/torque:

```text
F = sum_i F_i
M_p = sum_i (r_i - p) x F_i
```

This is the final common reduction for force rows from Maxwell-stress patches,
Lorentz elements, pressure faces, and nodal loads.  A force couple has zero net
force and a pivot-independent moment.

P1 constant surface traction load:

```text
F_e = A t,      f_i = F_e / 3
```

The equal nodal load preserves force and moment because the centroid of a
linear triangle is the mean of its three vertices.  Use this as the readable
boundary-load gate before replacing ``t`` with Maxwell traction on `.vol`
surface triangles.

Planar CAD face vector traction:

```text
F_face = A_face t,      M_p = (c_face - p) x F_face
```

Use this when checking a mesh or CAD handoff before the traction is generated
from fields.  Scalar pressure is tied to the oriented normal
(``F = p n A``); vector traction is already a global vector, so it does not
change direction when it is applied to another box face.

Uniform air-gap pressure:

```text
p = B_gap^2 / (2 mu0),     F = p A
```

This is the fast sanity check for solenoids, relays, magnetic circuits, and
pole-face holding force.  It is the local Maxwell stress evaluated in the gap.

Sampled cylindrical air-gap torque:

```text
tau(theta) = Br(theta) Bt(theta) / mu0
T = r^2 L integral tau(theta) dtheta
```

This is the machine-torque post-processing path for exported air-gap samples.
The uniform helper is a closed-form special case; the sampled helper integrates
segment contributions so sector scaling and harmonic signs stay visible.

2D sector-to-machine torque scaling:

```text
T_whole = T_2d * length_unit_m^2 * stack_length_m * symmetry_factor
```

Use ``machine_scaling.torque_scaling_summary`` when a 2D torque is per unit
axial depth and was solved on a rotational sector.  The helper makes the mesh
unit, active stack length, and modelled-sector multiplier visible in the JSON
record, which prevents comparing a sector N/m value against a whole-machine
N m target by accident.

Electrostatic normal pressure:

```text
E = V / d
p = 0.5 eps |E|^2
F = p A = 0.5 eps A V^2 / d^2 = W / d
fixed voltage coordinate force: F_x = 0.5 V^2 dC/dx
fixed charge coordinate force:  F_x = 0.5 Q^2 C^-2 dC/dx
```

This is the electric Maxwell-stress twin of magnetic air-gap pressure.  Use the
parallel-plate identity as the first MEMS/electrostatic force gate, then compare
the same pressure against ``force.electrostatic_traction_summary``.  For shape
or displacement sweeps, use ``electrostatics.capacitance_gradient_force_summary``:
if the coordinate is a gap/height, ``dC/dx`` is usually negative and the signed
force points toward smaller gap.  Then move to a full electric weighted-stress
extraction.

RF/time-harmonic radiation pressure at normal incidence:

```text
p = (A + 2 R) I / c,     F = (A + 2 R) P / c
one-sided scattering:    F = (1 + R - T) P_inc / c
```

where ``I`` is time-average Poynting intensity, ``P`` is integrated incident
power, ``A`` is absorptance, ``R`` is reflectance, and ``T`` is transmittance.
The two factors are identical because ``A = 1 - R - T``.  Use this for absorber,
short, mirror, transmitted through-line, and waveguide-port force sanity checks
before attempting a full time-averaged Maxwell-stress surface integral.

For a vector Poynting patch calculation with propagation direction ``k`` and
surface normal ``n`` pointing out of the illuminated side:

```text
P_inc = |S| A_patch max(0, -k.n)
F_abs = absorptance * P_inc k / c
F_ref = -2 reflectance * P_inc max(0, -k.n) n / c
```

This reproduces the oblique-incidence ``cos^2(theta)`` normal force and the
absorbed tangential momentum term.

Complex phasor Maxwell stress:

```text
<T_ij> = q Re(eps E_i E_j* + mu H_i H_j*)
         - q/2 (eps |E|^2 + mu |H|^2) delta_ij
q = 1/2 for peak phasors, q = 1 for RMS phasors
```

For a normally incident plane wave, ``<T> n = -(I/c) n`` on a control surface
whose normal points along propagation.  The receiving surface force is the
opposite sign, with the reflected field doubling the pressure for an ideal
mirror.

Energy/coenergy:

```text
constant current:      F = dW_co/dx
fixed flux/field:      F = -dW/dx
linear energy check:   L = 2 W / I^2
```

For nonlinear B-H force, prefer a weighted-stress or coenergy virtual-work
method.  ``force.virtual_work_force_from_displacement_samples`` is the readable
post-processor for displacement sweeps: use ``energy_kind="coenergy"`` for
fixed-current samples and ``energy_kind="stored_energy"`` for the negative
stored-energy derivative.  Finite-difference energy is useful as an independent
check, but it requires a stable geometry perturbation and matched meshes or
careful remeshing.

## Validation anchors to remember

- Uniform magnetized sphere: average interior ``B`` should approach ``2 Br / 3``.
- Two long parallel wires: ``F/L = mu0 I1 I2 / (2 pi d)``.
- Force couple: two opposite tangential forces at radius ``r`` give
  ``T = 2 r F`` and zero net force.
- Wire above high-permeability plane: image-current limit gives
  ``F/L = mu0 I^2 / (4 pi h)``.
- Uniform air gap: ``p(1 T) = 1 / (2 mu0) = 397887.35772973835 Pa``.
- RF normal incidence: perfect absorber ``F = P/c``; perfect reflector
  ``F = 2 P/c``.
- Round-wire internal inductance at DC: ``L_int = mu0 / (8 pi)`` per length.
- dq PM torque: ``T = (3/2) p (lambda_m iq + (Ld-Lq) id iq)``.

## Practical pitfalls

- Keep the Maxwell stress integration surface or eggshell band entirely in air.
- Do not integrate on material interfaces where the magnetic field trace jumps.
- Curve circular/curved geometry before force extraction; polygonal geometry can
  create a systematic force error even when the algebra is correct.
- For high-permeability or nonlinear bodies, place the body directly in the air
  region and use an analytic eggshell band.  Do not insert a nested shell material
  that electrically isolates the body from the surrounding HCurl region.
- For harmonic complex fields, use the time-averaged Maxwell stress with the
  ``1/2`` and ``1/4`` factors, not the static quadratic expression.
- For 2D planar results, forces are per unit out-of-plane length and torque is
  per unit length.  Axisymmetric helpers include the ``2 pi r`` weight and return
  full 3D force/inductance.

## Related radia-ngsolve topics

- ``ngsolve_usage("lorentz_force")``: conductor force, image force, busbars.
- ``ngsolve_usage("air_gap_force")``: magnetic-circuit holding force.
- ``validation_test/fem_readable/validation_surface_maxwell_force_trace.py``:
  first-order `.vol` boundary triangles to Maxwell traction and P1 nodal force loads.
- ``validation_test/fem_readable/validation_surface_triangle_constant_traction_load.py``:
  constant P1 surface-triangle traction to equivalent nodal loads, with
  force/moment preservation.
- SHOWCASE NOTEBOOK ``docs/fem_readable/fem_readable_validation_archive.ipynb``:
  the 12 readable FEM validation scripts, their summary JSON payloads, full
  source, and SHA-256 hashes, executed and synchronized with
  ``docs/fem_readable/fem_readable_validation_archive_results.json``.
- ``examples/cubit_mesh_export/validation_vol_boundary_pressure_moment.py``:
  named `.vol` sidesets to pressure force/moment rows and generic resultant reduction.
- ``examples/cubit_mesh_export/validation_vol_boundary_traction_moment.py``:
  named `.vol` sidesets to constant vector-traction force/moment rows and
  generic resultant reduction.
- ``validation_test/build123d_netgen_gmsh_flow/validation_build123d_cubit_pressure_moment.py``:
  build123d analytic box pressure moments checked against named `.vol` rows.
- ``validation_test/build123d_netgen_gmsh_flow/validation_build123d_cubit_traction_moment.py``:
  build123d analytic box vector-traction moments checked against named `.vol`
  rows.
- ``validation_test/electric_machine/validation_planar_lorentz_block_force.py``:
  planar ``Jz x B`` block force checked against the directional two-wire
  Lorentz summary.
- ``validation_test/electric_machine/validation_machine_torque_scaling.py``:
  2D sector torque to whole-machine N m scaling with stack length, mesh unit,
  and symmetry factor visible.
- ``validation_test/electric_machine/validation_virtual_work_force_displacement_sweep.py``:
  displacement energy/coenergy samples to force, including the fixed-current
  versus fixed-flux sign gate.
- ``validation_test/electric_machine/validation_sampled_air_gap_shear_torque.py``:
  sampled cylindrical air-gap ``Br``/``Bt`` to torque, including uniform and
  harmonic analytic gates.
- ``validation_test/electrostatics/validation_parallel_plate_electrostatic_force.py``:
  capacitance-energy, Maxwell-pressure, and traction-equivalence gate.
- ``validation_test/electrostatics/validation_capacitance_gradient_force.py``:
  fixed-voltage / fixed-charge capacitance-gradient force with signed
  gap/closing coordinates.
- SHOWCASE NOTEBOOK ``docs/electrostatics/electrostatics.ipynb``: the five
  electrostatic identities above, executed and rendered (parallel-plate,
  coaxial, capacitance-gradient, sphere-over-plane, Maxwell traction).
- ``examples/rf_waveguide/validation_scattering_radiation_force.py``:
  one-sided normal-incidence scattering force from reflectance/transmittance,
  equivalent to the absorber/reflector radiation-pressure factor.
- ``ngsolve_usage("dq_torque")`` / ``ngsolve_usage("mtpa")``: machine torque maps.
- ``ngsolve_usage("electrostatic_force")``: MEMS electric force.
- ``force_validation("cross_validation")``: stored regression-reference cases.
- SHOWCASE NOTEBOOK ``docs/force_validation/force_validation.ipynb``: the
  coenergy/torque, Maxwell-contour, two-wire Lorentz, torque-waveform, and
  virtual-work force identities, executed and rendered.
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
    "method_map": FORCE_METHOD_MAP,
    "map": FORCE_METHOD_MAP,
    "overview": FORCE_METHOD_MAP,
    "electromagnetic_force": FORCE_METHOD_MAP,
    "em_force": FORCE_METHOD_MAP,
    "torque": FORCE_METHOD_MAP,
    "method_selection": FORCE_METHOD_MAP,
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

    topic: all (default) | method_map | eggshell | cross_validation | reference_note
    """
    t = (topic or "all").strip().lower()
    if t in ("all", ""):
        return "\n\n".join([FORCE_METHOD_MAP, CROSS_VALIDATION, EGGSHELL, REFERENCE_NOTE])
    if t in _TOPICS:
        return _TOPICS[t]
    return (f"Unknown topic '{topic}'. Options: all, eggshell, cross_validation, "
            f"reference_note, method_map.\n\n" + FORCE_METHOD_MAP)
