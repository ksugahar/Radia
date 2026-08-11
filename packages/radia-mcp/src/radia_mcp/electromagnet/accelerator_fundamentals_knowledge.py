"""Curated accelerator fundamentals for ``radia_mcp.electromagnet``.

The material in this module is a public-safe synthesis of a 12-volume
accelerator textbook corpus.  It records engineering rules, equations, scope
boundaries, and source locators; it does not reproduce the source PDFs.
"""

from __future__ import annotations

TOPICS: dict[str, str] = {
    "accelerator_fundamentals": (
        "Accelerator families, magnetic rigidity, reference coordinates, and "
        "the lattice-to-magnet requirements contract"
    ),
    "beam_optics_contract": (
        "Transfer matrices, Twiss parameters, emittance, dispersion, tune, "
        "chromaticity, and the field quantities a magnet solver must deliver"
    ),
    "enge_fringe_integrals": (
        "Enge reference-orbit form factors I1/I2, effective pole-face maps, "
        "and momentum-indexed soft-edge diagnostics for FFAG magnets"
    ),
    "accelerator_magnet_types": (
        "Dipole, quadrupole, sextupole, combined-function, corrector, solenoid, "
        "kicker, septum, undulator, and wiggler roles"
    ),
    "accelerator_magnet_design": (
        "Iron-dominated magnet first estimates, nonlinear 2-D/3-D design, "
        "coils, cooling, tolerances, field quality, and design gates"
    ),
    "rapid_cycling_magnets": (
        "Ramped-magnet eddy currents, laminations, stranded conductors, "
        "hysteresis, power-supply ripple, grounding, and circuit-field coupling"
    ),
    "superconducting_accelerator_magnets": (
        "Critical-surface margin, multifilament conductors, coupling loss, "
        "persistent-current harmonics, stability, quench, and protection"
    ),
    "accelerator_magnet_measurement": (
        "Rotating/search coils, Hall and NMR probes, magnetic center, "
        "integrated multipoles, current history, and commissioning closure"
    ),
    "accelerator_model_boundaries": (
        "What a magnet model does not settle: RF acceleration, space charge, "
        "collective effects, impedance, beam loss, radiation, and tracking"
    ),
    "accelerator_sources": (
        "Searchable guide to the 12 accelerator textbooks used for this "
        "knowledge layer"
    ),
}


ACCELERATOR_FUNDAMENTALS = r"""
# Accelerator Fundamentals for Electromagnet Design

## What accelerates and what guides

- Electric fields and RF cavities change particle energy.  A static magnetic
  field changes direction but does no work: `q v x B` is perpendicular to the
  velocity.
- Magnets nevertheless determine whether the accelerated beam remains inside
  the available aperture.  They define the reference orbit, linear focusing,
  chromatic correction, orbit correction, injection/extraction, and insertion
  devices.
- Linacs use a single-pass sequence of RF and focusing elements.  Cyclotrons,
  synchrocyclotrons, synchrotrons, FFAGs, and storage rings reuse magnetic
  elements over many turns, so field reproducibility, tune, resonances, ramp
  history, and integrated errors become central.

## Main accelerator families and their magnet implications

| Family | Synchronism or trajectory principle | Main magnet implication |
|---|---|---|
| Electrostatic accelerator | A static voltage accelerates the beam in one or more passes | Transport magnets see a mostly fixed rigidity; voltage and insulation usually set the energy ceiling. |
| Linear accelerator | RF structures accelerate a single-pass beam along an open trajectory | Quadrupole or solenoid focusing is distributed between accelerating sections; alignment and wakefield sensitivity accumulate along the line. |
| Classical cyclotron | Constant dipole field and fixed RF frequency rely on the nonrelativistic cyclotron frequency | Pole shaping must provide vertical focusing while preserving the average field; relativistic phase slip limits energy. |
| Synchrocyclotron | RF frequency is swept as the particle becomes relativistic | The main magnet can remain nearly static, but acceleration is pulsed rather than continuous. |
| Isochronous cyclotron | The average field rises with radius to preserve revolution frequency | Sector flutter and spiral angle supply focusing while a deliberately nonuniform field maintains isochronism. |
| Synchrotron | Dipole field and RF frequency/phase follow the increasing rigidity on a nearly fixed orbit | Dipoles and quadrupoles ramp reproducibly; eddy currents, hysteresis, tracking between families, and power-supply ripple become first-order design inputs. |
| Storage ring | The beam circulates at nearly fixed energy for many damping or collision times | Integrated multipole errors, magnetic-center stability, orbit correction, tune, chromaticity, dynamic aperture, and insertion-device feed-down dominate. |
| FFAG | Fixed fields guide beams over a momentum range using alternating-gradient optics | Large good-field aperture and controlled nonlinear field law replace magnet ramping; orbit excursion must be included in the field-quality domain. |
| Beam line or spectrometer | A single-pass lattice transports, focuses, selects, or analyzes momentum | Magnet acceptance, fringe fields, dispersion, collimation, and an end-to-end field-map convention are essential. |

Betatrons and microtrons are useful historical and conceptual families:
betatrons accelerate through transformer induction and require the betatron
flux condition, while microtrons reuse an RF cavity with orbit timing chosen
so successive passes remain synchronous.  They should not be treated as
ordinary fixed-rigidity transport lines.

## Reference coordinates

Use the accelerator Frenet frame when translating optics into a field model:

- `s`: path length along the reference orbit;
- `x`: horizontal displacement, normally in the bend plane;
- `y`: vertical displacement;
- `h(s) = 1/rho(s)`: reference-orbit curvature.

For a straight Radia model, the repository convention may represent `s` by
the global `z` axis.  Record the transform explicitly.  Never compare a field
map, multipole phase, or skew/normal sign until charge sign, beam direction,
coordinate handedness, and the complex-field convention are fixed.

## Magnetic rigidity: the first magnet input

For particle momentum `p` and charge `q = Z e`,

    B rho = p / |q|
    B rho [T m] = 3.33564095198152 * p [GeV/c] / |Z|

The rigidity, not beam energy by itself, sets magnetic deflection.  For a
dipole on the reference trajectory,

    theta = integral(B_perp ds) / (B rho)

with the sign supplied by the charge and coordinate convention.  Therefore a
dipole requirement is fundamentally an integrated field requirement.  A
center-field target without magnetic length and end-field convention is
incomplete.

For a quadrupole with body gradient `G = dB_y/dx`, define the normalized
strength

    k1(s) = G(s) / (B rho)

and, in the thin-lens approximation, `1/f = integral(k1 ds)`.  One transverse
plane focuses while the other defocuses; alternating-gradient cells obtain
net confinement in both planes.

For a sextupole, a common convention is

    B_y + i B_x = (B''/2) * (x + i y)^2
    k2(s) = B''(s) / (B rho)

but codes differ over factorials and signs.  Exchange integrated normal/skew
multipoles with an explicit convention rather than a bare value called
"sextupole strength."

## Accelerator requirement handoff to Radia

Before CAD or meshing, the lattice/beam owner supplies:

1. particle species, charge state, momentum range, and beam direction;
2. reference trajectory and coordinate transform;
3. magnet family and required integrated main field or gradient;
4. physical and good-field apertures, reference radius, and allowed orbit
   excursion including dispersion;
5. allowed integrated normal and skew multipoles over the operating range;
6. current waveform, repetition rate, flat-bottom/flat-top timing, and field
   history/pre-cycle;
7. alignment, magnetic-center, reproducibility, ripple, and stray-field limits;
8. required field-map domain and the downstream tracking convention.

Radia then owns the electromagnetic realization: coil/yoke geometry,
materials, saturation and hysteresis, open boundary, end fields, losses where
applicable, and checked field artifacts.  Beam optics/tracking closes the loop
against the original lattice requirements.
"""


BEAM_OPTICS_CONTRACT = r"""
# Beam Optics Contract for Magnet Analysis

## Linear optics

Near the reference orbit, uncoupled transverse motion is represented by Hill
equations and 2 x 2 transfer matrices.  A valid linear element map is
symplectic; in one uncoupled plane its determinant is one.  For each plane the
Twiss parameters satisfy

    beta * gamma - alpha^2 = 1

and the Courant-Snyder invariant is

    gamma*x^2 + 2*alpha*x*x' + beta*x'^2 = epsilon.

The phase advance accumulated around a ring gives the tune.  Integer,
half-integer, coupling, and higher-order resonance conditions make small
systematic or random field errors accumulate coherently.  Magnet acceptance
therefore cannot be reduced to peak `|B|` error at one point.

## Dispersion and chromaticity

Particles with relative momentum error `delta = Delta p/p` have a displaced
closed orbit.  To first order,

    x(s) = x_beta(s) + D_x(s)*delta
    sigma_x^2 = beta_x*epsilon_x + D_x^2*sigma_delta^2

so the physical aperture and the field-quality region must include both the
betatron envelope and dispersive orbit.  The tune shift

    Delta Q_x,y = xi_x,y * delta

defines chromaticity.  Sextupoles correct chromaticity only when placed where
dispersion is nonzero, and they also introduce amplitude-dependent tune and
nonlinear resonances.  A better local sextupole field is not automatically a
better ring unless the lattice-level dynamic aperture is checked.

## Field quantities to exchange

For each current or ramp state, export at least:

- the reference-orbit field and `integral(B_perp ds)`;
- integrated normal and skew multipoles at the declared reference radius;
- body and integrated transfer strengths (`integral(k1 ds)`,
  `integral(k2 ds)`, and higher terms as required);
- effective magnetic length and the longitudinal distribution of the main
  and error harmonics;
- fringe-field maps for edge focusing and injection/extraction studies;
- current, field-history, material-state, temperature, geometry, and
  coordinate provenance.

Use a Maxwell-consistent 3-D map for particle tracking through ends.  A local
effective-field-boundary slope is not a universal substitute for the Lorentz
kick through a finite fringe.  Evaluate edge focusing from the trajectory or
the appropriate field integral in the declared frame.

For the Enge reference-orbit displacement and effective edge map, including
the distinction between Enge profile coefficients and the form factors
`I1`/`I2`, use `electromagnet_usage("enge_fringe_integrals")`.

## Field-quality interpretation

The useful distinction is:

- **body multipoles** diagnose the nominal 2-D cross-section;
- **integrated multipoles** determine the leading beam kick;
- **longitudinal harmonic profiles** diagnose where end packs, saturation,
  coil ends, or assembly errors create the integral.

For an ideal normal multipole with the usual mid-plane symmetries, leading
allowed families are `n = 1,3,5,...` for a dipole, `n = 2,6,10,...` for a
quadrupole, and `n = 3,9,15,...` for a sextupole.  Manufacturing and alignment
errors generate forbidden and skew terms.  Always state the complex-field,
normal/skew, normalization, reference-radius, and integration conventions.
"""


ENGE_FRINGE_INTEGRALS = r"""
# Enge Fringe Integrals for Reference Orbits and Edge Maps

## Do not confuse two different Enge parameter families

1. **Enge profile coefficients** `a_i` describe a fitted longitudinal falloff,
   for example

       S(z) = 1 / (1 + exp(2*p(z/epsilon)))
       p(u) = a1 + a2*u + ... + ak*u^(k-1).

2. **Enge form factors** `I1` and `I2` are integrals of the realized fringe
   field.  `I1` changes the reference trajectory; `I2` changes the effective
   vertical edge focusing.  They can be evaluated from an Enge fit, but they
   are not the polynomial coefficients and should normally be evaluated from
   the computed or measured field profile.

In this topic `g` is the **full vertical pole gap**.  A code using `HGAP` may
instead define `HGAP = g/2`.  Some references use a dimensional coordinate
normal to the pole face rather than the dimensionless `sigma` below; factors
of `g` then move between the coordinate and the integrals.  Never exchange an
`I1` or `I2` value without its coordinate, gap, field-normalization, entry/exit,
charge-sign, and pole-angle conventions.

## Effective field boundary and dimensionless coordinate

For a rectangular dipole exit, let `beta` be the pole-face rotation, `rho` the
body reference-orbit radius, `(xi,zeta)` a rectilinear frame at the face, and

    sigma = (zeta*cos(beta) + xi*sin(beta)) / g.

Let `By(sigma)` be the soft-edge mid-plane bending field, `B0` its body value,
and `By0(sigma)` the equivalent sharp-cutoff field.  Place `sigma=0` so the
soft and sharp profiles have the same field integral:

    integral(By d sigma) = integral(By0 d sigma).

With bounds in the constant-field and zero-field regions, define

    I1 = integral integral ((By0 - By)/B0) d sigma d sigma_star
    I2 = integral (By*(B0 - By)/B0^2) d sigma.

These are dimensionless under the `sigma` convention above.  The first-order
pure-dipole exit result is

    delta_x_Enge = g^2 * I1 / (rho*cos(beta)^2)
    tan(beta_v) = tan(beta)
                  - (g/rho)*(1 + sin(beta)^2)/cos(beta)^3 * I2

and the reduced map is

    x  = x0 + delta_x_Enge
    x' = x0' + x0*tan(beta)/rho
    y  = y0
    y' = y0' - y0*tan(beta_v)/rho.

Thus `I1` is the Enge reference-orbit correction: an extended pure-dipole
fringe gives the same net bend as the equal-integral hard edge but shifts the
reference trajectory parallel to it.  `I2` modifies the effective vertical
pole-face angle and therefore the `R43`-type focusing term.  Apply the entry
map in its own local frame and sign convention; do not copy the exit signs by
inspection.

## Relation to pole-face geometry

For a rectangular-bend geometry represented as an `SBEND`, keep the orbit
geometry separate from the fringe correction.  With bend angle `phi`, pole
angles `E1,E2`, and straight magnetic length `Lmag`, the convention used in
CERN-ACC-NOTE-2018-0059 gives

    E1 + E2 = phi
    rho = Lmag / (sin(E1) + sin(E2))
    Larc = rho * phi.

`E1/E2` or `beta` describe the pole-face geometry.  `I1/I2` describe the
soft-edge field about that geometry.  `Larc` is the design-orbit path length,
not the physical yoke length and not automatically the straight magnetic
length.

## FFAG use: make the form factors orbit dependent

Bell and Abell's non-scaling FFAG example uses a mid-plane combined-function
field

    By(x,0,z) = T(x)*S(z),       T(x) = B0 + G*x,

over 31--250 MeV proton closed orbits with about 0.2 m radial excursion.
Entrance/exit angles change with momentum, adjacent magnet fringes can overlap,
and the soft-edge model materially changes tunes and resonance-crossing loss.
Consequently a single magnet-wide `I1,I2,beta` tuple is not an adequate FFAG
contract.  For momentum samples `p_j`, record at every entrance and exit:

    orbit_j, beta_j, rho_j, B0_j, g, I1_j, I2_j,
    delta_x_Enge_j, beta_v_j, and the local frame.

The recommended design loop is:

1. recover the closed orbit for every required momentum;
2. evaluate the field along the local pole-face normal and locate the
   equal-integral effective field boundary;
3. integrate `I1_j` and `I2_j` as orbit-indexed diagnostics/reduced maps;
4. compute the full Lorentz trajectory and transfer map through the combined
   3-D field, including overlapping fringes;
5. optimize orbit closure, tune/transfer-map targets, aperture, and resonance
   margins jointly across momentum; then repeat the closed-orbit solve.

For a combined-function field the Enge displacement moves the particle into a
different `T(x)`, so the final deflection need not remain parallel to the
hard-edge orbit.  In that case `I1/I2` are valuable initialization and
diagnostic quantities, but detailed ray tracking through the realized field is
the acceptance calculation.

## HDiv-MMM connection and the air-mesh-free claim

HDiv-MMM is well matched to this FFAG loop when the design variables are
ferromagnetic yoke/pole elements:

- mesh only the magnetic material with an HDiv space (BDM1 for production
  accelerator work); do not create a vacuum/air volume mesh or a Kelvin shell;
- represent coils by the audited Biot-Savart source and evaluate the resulting
  field at closed-orbit and fringe quadrature points in otherwise unmeshed air;
- reuse the compressed Laplace interaction operator while changing active iron
  elements or deforming the material mesh;
- batch the many momentum/orbit observation rows and adjoint contractions.

This is a real advantage for FFAGs because the good-field region spans the full
orbit excursion and because repeated pole/yoke topology changes would otherwise
force repeated remeshing of a large 3-D air volume.  The accurate claim is
**air-volume-mesh-free magnet optimization**, not "no vacuum computation":
field observations, closed-orbit recovery, Lorentz tracking, and tune/map
evaluation in the vacuum remain mandatory.

Do not overclaim the reduced map.  If particles are tracked through the full
soft-edge field, its Enge shift and focusing are already present; multiplying
an additional `I1/I2` edge map double counts them.  Also retain a full-field
check for saturation, combined-function gradients, overlapping fringes, large
incidence angles, and off-mid-plane motion.

## Executable Python PoC

The reduced target family is available without an air volume mesh:

    from radia.ffag_topopt import (
        FFAGSoftEdgeCellSpec,
        build_ffag_cell_target_family,
    )

    spec = FFAGSoftEdgeCellSpec.bell_abell(full_gap_m=0.10)
    family = build_ffag_cell_target_family(
        [31.0, 140.0, 250.0], spec=spec, n_segments=256,
        bend_field_band=2e-3, transfer_matrix_band=2e-3)

`family.references` records the periodic reduced orbit, rigidity, transverse
offset, sampled soft-edge field, and cell map at every energy.
`family.fringe_integrals` records the declared-gap `I1/I2` diagnostic.
`family.objective` is a `MultiMomentumTransferMatrixObjective` and supplies the
analytic `transform` and `transform_jacobian` consumed by the binary master.

For a BDM1 iron candidate mesh, build all physical vacuum rows with
`build_multi_orbit_field_response_matrix(...)`, then call
`optimize_hdiv_mmm_magnet_from_transfer_matrices(...)`.  The latter sends one
fused multi-momentum contract to the existing whole-element Schur,
ACA--QR--TSVD, and 0--1 LP driver.  The optimizer does not finite-difference a
design variable.  A completed design must still recover the full-field closed
orbits and recompute the maps around those recovered orbits.

## Sources

- G. I. Bell and D. T. Abell, "Fringe field simulations of a non-scaling FFAG
  accelerator," arXiv:1202.0805.
- R. Molloy and S. Blitz, "Fringe Field Effects on Bending Magnets, Derived for
  TRANSPORT/TURTLE," arXiv:1310.8630.
- O. E. Berrig, "Review of path-length calculations for rectangular bending
  magnets with arbitrary pole-face angles," CERN-ACC-NOTE-2018-0059.
"""


ACCELERATOR_MAGNET_TYPES = r"""
# Accelerator Magnet Types and Selection

| Element | Leading field | Beam purpose | Dominant engineering concern |
|---|---|---|---|
| Dipole | approximately uniform transverse `B` | Bend/reference orbit | integrated field, ends, saturation, ripple |
| Quadrupole | `B` linear in transverse position | Focus one plane, defocus the other | gradient integral, magnetic center, allowed `n=6,10,...` |
| Sextupole | `B` quadratic in position | Chromatic correction | nonlinear errors, alignment feed-down, small coil space |
| Higher multipole | higher transverse order | resonance correction or nonlinear optics | very small error budget and positioning |
| Combined-function | dipole plus gradient (and sometimes higher terms) | Bend and focus in one yoke | coupled excitation and reduced tuning freedom |
| Corrector/steerer | weak dipole, often two-axis | Closed-orbit correction | reproducibility, cross-talk, low-current hysteresis |
| Solenoid | longitudinal field | Coupled focusing, low-energy transport | transverse coupling and fringe rotation |
| Kicker | fast pulsed transverse field | Turn-selective injection/extraction | rise/fall time, flatness, impedance, timing jitter |
| Septum | strong field beside a field-free channel | Separate circulating/extracted trajectories | leakage field, thin septum heating/stress |
| Undulator/wiggler | alternating periodic field | Radiation source or damping | phase error, trajectory integral, end compensation |

## Technology selection

- **Iron-dominated normal-conducting magnet**: preferred when moderate field,
  precise pole shaping, straightforward service, and controllable DC/ramped
  operation are more important than compact peak field.
- **Air-core or transmission-line pulsed magnet**: preferred for very fast
  fields where iron loss, saturation, and pulse propagation would dominate.
- **Permanent-magnet or hybrid magnet**: avoids continuous coil power but needs
  temperature, radiation, aging, sorting, shimming, and tuning provisions.
- **Superconducting magnet**: preferred when aperture-field product or power
  makes normal conduction impractical, provided cryogenics, stored energy,
  quench protection, conductor magnetization, and training are accepted.

Selection is driven by the full operating cycle, not only the maximum field.
A storage-ring DC quadrupole and a rapid-cycling synchrotron quadrupole with the
same gradient require different yokes, conductors, power systems, and QA.
"""


ACCELERATOR_MAGNET_DESIGN = r"""
# Accelerator Magnet Design Workflow

## First magnetic-circuit estimates

For an unsaturated iron-dominated dipole with one dominant air gap `g`,

    N I approximately B_gap * g / mu0.

For an ideal quadrupole with pole-tip radius `r_p`, body gradient `G`, and pole
face `x*y = +/- r_p^2/2`,

    B_pole approximately G*r_p
    N I per pole approximately G*r_p^2/(2*mu0).

With `B_y = (B''/2)(x^2-y^2)` for a sextupole,

    N I per pole approximately B''*r_p^3/(6*mu0).

These are sizing estimates, not acceptance results.  Finite permeability,
leakage, finite pole width, coil windows, holes, joints, three-dimensional
ends, and saturation move the required current and harmonics.

## Design sequence

1. Convert lattice requirements into rigidity, aperture, integrated main
   field/gradient, reference radius, and multipole limits.
2. Choose topology and symmetry without deleting physically required current
   return paths, end fields, holes, or asymmetric services.
3. Size the gap, pole, yoke, coil window, and ampere-turns with a magnetic
   circuit.  Use the actual `B-H` curve; a rule such as "keep soft iron near
   1.5-2 T" is only an early warning, not a material model.
4. Optimize the 2-D body with nonlinear material data.  Pole width and shims
   trade good-field radius against coil/yoke space.  Check low-field remanence
   and high-field saturation, not only the nominal current.
5. Solve the full 3-D magnet with real coil ends, end packs/chamfers, joints,
   holes, and open boundary.  Match integrated field, effective length, fringe
   shape, peak conductor field, and integrated multipoles.
6. Design the electrical/thermal system: turns, current, conductor area,
   resistance versus temperature, inductance, voltage, stored energy,
   insulation, cooling, pressure loss, and fault cases.
7. Propagate manufacturing, assembly, alignment, material, and current errors
   into magnetic center and integrated multipoles.
8. Define pre-cycle, measurement, fiducial, acceptance, and field-map
   procedures before releasing drawings.

## Coil and cooling checks

At minimum evaluate `P_cu = I^2 R(T)`, coolant heat balance
`m_dot = P/(c_p Delta T)`, hydraulic pressure drop, maximum conductor and
insulation temperature, voltage to ground, radiation environment, and loss of
flow.  Hollow conductors are common at high DC loss; rapid-cycling conductors
may require stranded/transposed construction to control circulating currents.
Cooling velocity cannot be increased without limit: erosion, corrosion,
vibration, water chemistry, and local flow changes become lifetime constraints.

## Required gates

- **Geometry/labels**: checked solver mesh and material/boundary ownership.
- **Magnetostatic**: low/nominal/high current, saturation, leakage, force.
- **3-D field quality**: integrated main field, multipoles, ends, stray field.
- **Electrical/thermal**: resistance, voltage, power, cooling, stored energy.
- **Transient**: mandatory for ramps, pulses, conductive structures, or a
  field-history-dependent material.
- **Beam closure**: reference trajectory, optics strength, fringe tracking,
  aperture, and dynamic consequences of accepted errors.
"""


RAPID_CYCLING_MAGNETS = r"""
# Rapid-Cycling Magnets, Eddy Currents, and Power Supplies

## Why the DC model is insufficient

A changing field induces currents in the yoke, pole ends, clamps, vacuum
chamber, cooling structures, and coil conductor.  These currents can:

- delay and phase-shift the field relative to current;
- distort body and end harmonics during the ramp;
- change effective length and leakage field;
- generate local heat and mechanical force;
- make a current-ripple spectrum differ from the field-ripple spectrum.

The first screening length for a linear material is

    skin_depth = sqrt(2/(omega*mu*sigma)).

It is useful for choosing a lamination or strand scale, but nonlinear
permeability, hysteresis, pole-end geometry, and electrical contacts require a
transient field solve for final acceptance.

## Construction rules

- Laminate the yoke with electrically insulated sheets oriented so desired
  flux remains easy while large eddy-current loops are interrupted.
- Select electrical steel from the actual frequency, peak induction, loss,
  permeability, coercivity, isotropy, stacking factor, punching stress, and
  reproducibility requirements.  "Thinner is always better" ignores cost,
  stacking, stiffness, and material tradeoffs.
- Split or insulate end plates, clamps, vacuum chambers, and structural loops
  that would otherwise act as shorted turns.
- Use stranded or transposed conductors when conductor eddy loss matters.
  Strand insulation and transposition must survive impregnation and service.
- Model the true waveform.  A resonant sinusoidal supply, arbitrary pattern
  supply, flat-bottom hold, and pulsed kicker are different problems even at
  the same repetition frequency.

For classical lamination loss, the eddy component scales qualitatively as
`B_peak^2 * f^2 * thickness^2 / resistivity`; hysteresis loss scales roughly
with cycle rate times loop area.  Treat those as design trends and use measured
material loss data for prediction.

## Circuit and field are one system

An accelerator electromagnet is not just `R + L`.  Include mutual coupling,
stray capacitance to ground, cables, bus topology, shunts/damping, supply
output impedance, and grounding.  Common-mode currents can create unwanted
multipoles depending on pole-to-pole wiring.  A symmetric connection should be
checked electrically and magnetically, not inferred from the schematic name.

Field reproducibility also depends on the `B-I` path.  Define injection/flat
bottom, ramp, flat top, reset, and pre-cycle.  Compare simulation and
measurement at the same current history and time, and report transfer
functions from current to integrated field/harmonics over the frequency band
of interest.
"""


SUPERCONDUCTING_ACCELERATOR_MAGNETS = r"""
# Superconducting Accelerator Magnets

## Operating margin

The admissible current is set by the conductor critical surface, commonly
represented as `J_c(B,T,strain)`, not by a single critical-current number.
Evaluate the peak field and temperature on the conductor, including coil-end
enhancement and self field.  Define current-sharing, temperature, field, and
mechanical margins for the complete operating cycle.

## Multifilament and cable physics

- Fine filaments reduce flux-jump energy and persistent-current magnetization.
- Twisting interrupts long coupling loops; matrix transverse resistivity and
  twist pitch set coupling-current time constants.
- Rutherford, cable-in-conduit, or other cables add interstrand coupling,
  current redistribution, contact resistance, and mechanical constraints.
- Ramped operation has superconductor hysteresis loss, interfilament and
  interstrand coupling loss, and eddy loss in stabilizer/support structures.

These mechanisms are also field-quality mechanisms.  Persistent-current and
coupling-current multipoles can dominate at low field even when the geometric
coil is perfect; iron saturation and coil deformation dominate elsewhere in
the cycle.  A static main-field calculation alone cannot certify a ramped
superconducting accelerator magnet.

## Stability, quench, and protection

A defensible design couples electromagnetic, thermal, circuit, and mechanical
models and checks:

1. credible local disturbances and enthalpy/temperature margin;
2. minimum propagation zone and normal-zone propagation;
3. detection threshold and validation delay;
4. current decay through dump, heaters, CLIQ-like systems, or quench-back as
   applicable;
5. hot-spot temperature, turn-to-turn and terminal voltage, stored-energy
   extraction, and force/stress during the fault;
6. failed sensor, failed heater, and loss-of-cooling cases.

Radia can own the magnetic geometry and field-quality pieces.  Conductor-level
AC loss, quench propagation, cryogenic hydraulics, and protection remain named
coupled analyses; do not silently replace them with a magnetostatic surrogate.
"""


ACCELERATOR_MAGNET_MEASUREMENT = r"""
# Accelerator Magnet Measurement and Commissioning

## Match the instrument to the observable

- **Rotating/harmonic coil**: normal and skew multipoles, magnetic center, and
  integrated field when the coil spans the magnet.  Calibrate effective area,
  radius, angular encoder, integrator, and bucking ratios.
- **Search/flip coil**: flux and integrated field changes; useful for transfer
  functions and effective length.
- **Hall probe**: local vector-field mapping.  Control probe temperature,
  angular alignment, planar/transverse sensitivity, active-area averaging,
  and positioning.
- **NMR probe**: accurate absolute `|B|` in a sufficiently uniform field;
  unsuitable as a universal local vector or gradient probe.

## Centering and feed-down

An offset measurement axis converts a higher multipole into apparent lower
orders.  For example, a quadrupole produces an apparent dipole when the coil is
off center; a sextupole feeds down to quadrupole and dipole.  Solve for the
magnetic center and roll before interpreting forbidden harmonics, and keep the
mechanical-to-magnetic fiducial transform in the result.

## Acceptance sequence

1. Verify geometry, coil resistance/inductance, insulation, cooling, and
   sensor polarity before excitation.
2. Apply the declared pre-cycle and measure the full `B-I` transfer, including
   low-field remanence and high-field saturation.
3. Measure body and integrated harmonics at the specified reference radius and
   current states; report uncertainty and convention.
4. Map longitudinal fields/ends where effective length or edge focusing is a
   requirement.
5. For ramped/pulsed magnets, measure time alignment, current-to-field transfer,
   flatness, ripple, repeatability, and thermal steady state.
6. Transfer the magnetic center/axis to external fiducials and verify polarity
   and family consistency before installation.
7. Compare model and measurement only after matching current history,
   temperature, coordinates, integration length, and normalization.

The final acceptance object is the beam-relevant field integral and its
uncertainty/provenance, not a screenshot of a visually plausible field map.
"""


ACCELERATOR_MODEL_BOUNDARIES = r"""
# Boundaries of an Accelerator Magnet Model

## Questions the electromagnet model answers

- What coil/yoke/material geometry produces the requested field integral?
- How do saturation, hysteresis, eddy currents, temperature, ends, and
  tolerances change the main field and multipoles?
- What are conductor/yoke fields, losses, forces, voltages, and checked field
  artifacts for downstream analysis?

## Questions that require another model

- RF phase stability, bucket area, longitudinal capture, and acceleration;
- space-charge tune shift/spread, beam-beam effects, wake fields, and
  collective instabilities;
- vacuum-chamber impedance, beam loading, electron cloud, and feedback;
- beam loss, activation, radiation damage, and shielding;
- full-ring closed orbit, tune, chromaticity, resonance driving terms,
  dynamic aperture, and momentum acceptance;
- injection/extraction efficiency and septum/kicker timing in the lattice.

The textbook corpus is explicit that high-intensity beams are not a collection
of independent particles: space charge changes incoherent tune and can create
tune spread, while coherent modes and wall impedance add separate limits.
Likewise, a magnet can meet local multipole specifications and still degrade a
ring through placement, systematic accumulation, fringe kicks, or nonlinear
resonances.

## Required coupled loop

    lattice/beam requirements
        -> Radia electromagnetic design
        -> checked 3-D field map and integrated multipoles
        -> particle/optics tracking
        -> revised aperture, ends, correctors, or tolerances

Keep this loop explicit.  Radia is the electromagnetic engineering platform in
the loop; it should not claim that a field solve alone is an accelerator
performance proof.
"""


SOURCE_CATALOG: tuple[dict[str, object], ...] = (
    {
        "id": "eddy_currents_lammeraner_stafl",
        "filename": "EddyCurrents.pdf",
        "title": "Eddy Currents",
        "authors": "Jiri Lammeraner and Milos Stafl",
        "pages": 231,
        "focus": (
            "Skin effect, laminations, solid cores, coil-winding current "
            "redistribution, special conductors, approximate methods, and loss."
        ),
        "locators": "Contents on PDF pp.4-5; especially Chapters 2, 8, and 13.",
        "keywords": "eddy current skin depth lamination winding conductor loss",
    },
    {
        "id": "oho_2010_electromagnet_power",
        "filename": "OHO_OHO10_nakamura_20100822.pdf",
        "title": "Magnets and Vacuum: Electromagnets and Power Supplies",
        "authors": "OHO 2010 lecture note",
        "pages": 16,
        "focus": (
            "Dipole/quadrupole/sextupole basics, B-I history, R-L-stray-C "
            "equivalent circuits, wiring symmetry, ripple, resonance, and supplies."
        ),
        "locators": "PDF pp.1-6 for magnets; pp.7-15 for wiring, ripple, and supplies.",
        "keywords": "power supply wiring ripple resonance B-I common mode quadrupole",
    },
    {
        "id": "oho_1984_circular_accelerators",
        "filename": "OHO_txt-1984-Ⅱ.pdf",
        "title": "Principles of Circular Accelerators and Storage Rings",
        "authors": "OHO 1984 lecture note",
        "pages": 112,
        "focus": (
            "Reference orbit, magnetic rigidity, transfer matrices, strong "
            "focusing, tune, resonances, chromaticity, synchrotrons, and rings."
        ),
        "locators": "PDF p.13 for B rho; p.42 for chromaticity; p.56 for resonances.",
        "keywords": "beam optics rigidity matrix tune resonance chromaticity storage ring",
    },
    {
        "id": "oho_1984_accelerator_magnets",
        "filename": "OHO_txt-1984-Ⅲ.pdf",
        "title": "The Role and Design of Magnets in Accelerators",
        "authors": "OHO 1984 lecture note",
        "pages": 29,
        "focus": (
            "Magnet roles, field expansion, dipole/quadrupole design, leakage, "
            "coils/cooling, superconducting conductors, and field measurement."
        ),
        "locators": "PDF pp.1-16 for optics and normal magnets; pp.24-29 for superconducting magnets and measurement.",
        "keywords": "dipole quadrupole hyperbola leakage coil cooling superconducting measurement",
    },
    {
        "id": "oho_1986_beam_instability",
        "filename": "OHO_txt-1986-C-1.pdf",
        "title": "Beam Instability and Space-Charge Limits",
        "authors": "OHO 1986 lecture note",
        "pages": 43,
        "focus": (
            "Incoherent space charge, tune shift/spread, image effects, wall "
            "impedance, wake fields, head-tail and transverse instability."
        ),
        "locators": "PDF pp.1-8 for space charge; later chapters for impedance and instability.",
        "keywords": "space charge tune spread impedance wake instability head-tail",
    },
    {
        "id": "oho_2000_space_charge",
        "filename": "OHO_txt-2000-Ⅲ.pdf",
        "title": "Space-Charge Effects in High-Intensity Accelerators",
        "authors": "OHO 2000 lecture note",
        "pages": 25,
        "focus": (
            "Incoherent/coherent motion, Laslett tune shift, distributions, "
            "resonance interpretation, mitigation, aperture, and simulation."
        ),
        "locators": "PDF pp.1-7 for definitions and tune shift; remaining pages for stability and mitigation.",
        "keywords": "high intensity space charge coherent incoherent Laslett aperture loss",
    },
    {
        "id": "oho_2001_proton_synchrotron_magnets",
        "filename": "OHO_txt-2001-9.pdf",
        "title": "Magnets for High-Intensity Proton Synchrotrons",
        "authors": "OHO 2001 lecture note",
        "pages": 19,
        "focus": (
            "Rapid-cycling versus high-field magnets, sizing, electrical steel, "
            "lamination/skin effect, stranded aluminum conductors, heat, and materials."
        ),
        "locators": "PDF pp.1-4 for sizing; pp.5-8 for steel/laminations; pp.8-19 for conductors and thermal design.",
        "keywords": "rapid cycling synchrotron electrical steel lamination stranded conductor thermal",
    },
    {
        "id": "oho_beam_transport",
        "filename": "OHO_ビーム輸送の基礎.pdf",
        "title": "Beam Transport Fundamentals: Beam Optics and Electromagnets",
        "authors": "OHO lecture note",
        "pages": 28,
        "focus": (
            "Lorentz motion, element matrices, fringe/edge focusing, multipole "
            "fields, magnet construction, pulsed magnets, cooling, and measurement."
        ),
        "locators": "PDF pp.1-15 for optics/magnet fields; pp.16-21 for pulsed magnets; pp.22-27 for cooling/measurement.",
        "keywords": "beam transport transfer matrix edge focusing kicker cooling rotating coil Hall NMR",
    },
    {
        "id": "meot_zgoubi",
        "filename": "ZGOUBI_Understanding the Physics of Particle Accelerators.pdf",
        "title": "Understanding the Physics of Particle Accelerators",
        "authors": "Francois Meot",
        "pages": 651,
        "focus": (
            "Numerical ray tracing and worked accelerator physics for electrostatic "
            "systems, cyclotrons, synchrotrons, FFAGs, beam lines, and spectrometers."
        ),
        "locators": "Chapters 8-9 strong/weak focusing, 10-11 FFAG, 12-14 beam lines and optical elements.",
        "keywords": "ZGOUBI tracking cyclotron synchrotron FFAG beam line spectrometer optics",
    },
    {
        "id": "kamei_kihara_accelerator_science",
        "filename": "パリティ物理学_加速器科学.pdf",
        "title": "Accelerator Science",
        "authors": "Tohru Kamei and Motohiro Kihara",
        "pages": 293,
        "focus": (
            "Comprehensive accelerator history and fundamentals: RF acceleration, "
            "betatron motion, strong focusing, rings, linacs, insertion devices, and applications."
        ),
        "locators": "Contents on PDF pp.7-11; Chapters 2, 5, 8, 9, 10, and 16 are central to magnet/beam coupling.",
        "keywords": "accelerator science RF betatron strong focusing storage ring linac undulator",
    },
    {
        "id": "riken_accelerator_science_history",
        "filename": "加速器科学.pdf",
        "title": "Accelerator Science: RIKEN Historical and Application Overview",
        "authors": "RIKEN historical overview",
        "pages": 52,
        "focus": (
            "Development of cyclotrons, heavy-ion linacs, ion sources, the RIBF "
            "program, and the scientific/application context of accelerator facilities."
        ),
        "locators": "Read as historical/application context rather than a magnet-design formula source.",
        "keywords": "RIKEN history cyclotron heavy ion RILAC RIBF application",
    },
    {
        "id": "funaki_sumiyoshi_multifilament",
        "filename": "多芯線と導体.pdf",
        "title": "Fundamentals of Superconducting Engineering: Multifilamentary Wires and Conductors",
        "authors": "Kazuo Funaki and Fumio Sumiyoshi",
        "pages": 345,
        "focus": (
            "Critical-state conductor physics, twist/coupling currents, AC loss, "
            "current distribution, thermal/magnetic stability, quench, protection, and measurement."
        ),
        "locators": "Contents on PDF pp.5-9; Chapters 2-7 and 9 cover coupling, loss, stability, quench, and measurement.",
        "keywords": "superconducting multifilament twist coupling AC loss stability quench protection",
    },
)


def get_accelerator_source_guide(query: str = "") -> str:
    """Return a portable, public-safe source guide filtered by keywords."""
    terms = [term for term in query.casefold().split() if term]
    selected = []
    for entry in SOURCE_CATALOG:
        haystack = " ".join(str(value) for value in entry.values()).casefold()
        if all(term in haystack for term in terms):
            selected.append(entry)

    if not selected:
        return f"No accelerator textbook sources match query: {query!r}"

    total_pages = sum(int(entry["pages"]) for entry in SOURCE_CATALOG)
    lines = [
        "# Accelerator Textbook Source Guide",
        "",
        (
            f"Curated corpus snapshot: {len(SOURCE_CATALOG)} PDFs, "
            f"{total_pages} pages. The PDFs are not distributed by radia-mcp; "
            "this guide stores only bibliographic and topical locators."
        ),
    ]
    if query:
        lines.extend(["", f"Filter: `{query}` ({len(selected)} matches)"])

    for entry in selected:
        lines.extend(
            [
                "",
                f"## {entry['title']}",
                f"- ID: `{entry['id']}`",
                f"- Source file: `{entry['filename']}`",
                f"- Authors/context: {entry['authors']}",
                f"- Length: {entry['pages']} pages",
                f"- Use for: {entry['focus']}",
                f"- Locator: {entry['locators']}",
            ]
        )
    return "\n".join(lines)


DOCUMENTS: dict[str, str] = {
    "accelerator_fundamentals": ACCELERATOR_FUNDAMENTALS,
    "beam_optics_contract": BEAM_OPTICS_CONTRACT,
    "enge_fringe_integrals": ENGE_FRINGE_INTEGRALS,
    "accelerator_magnet_types": ACCELERATOR_MAGNET_TYPES,
    "accelerator_magnet_design": ACCELERATOR_MAGNET_DESIGN,
    "rapid_cycling_magnets": RAPID_CYCLING_MAGNETS,
    "superconducting_accelerator_magnets": SUPERCONDUCTING_ACCELERATOR_MAGNETS,
    "accelerator_magnet_measurement": ACCELERATOR_MAGNET_MEASUREMENT,
    "accelerator_model_boundaries": ACCELERATOR_MODEL_BOUNDARIES,
    "accelerator_sources": get_accelerator_source_guide(),
}
