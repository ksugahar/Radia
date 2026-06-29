"""Public-safe CAE loop learning rules for radia-ngsolve.

This module distills recurring validation-loop lessons into MCP-readable
guidance.  It deliberately avoids private paths, commercial solver provenance,
and benchmark values; keep those in internal cross-validation logs.
"""

TOPICS = {
    "overview": "How to close the loop from validation artifact to MCP learning",
    "dual_lane": "How one validation artifact teaches public and source-tool MCP lanes",
    "mesh_geometry_vol": "Geometry, Cubit/build123d mass properties, and Netgen .vol gates",
    "force_moment": "Force, moment, Maxwell traction, Lorentz, and coenergy gates",
    "motor_airgap_torque": "Motor air-gap Maxwell shear torque from Br/Bt harmonics",
    "electrostatic_layered_dielectric": "Layered dielectric stack capacitance, D-continuity, and energy gates",
    "acoustic_impedance_power": "Acoustic impedance reflection, absorption, and boundary power",
    "rf_acoustic_passivity": "Acoustic/RF passivity and power-balance identities",
    "mcp_closure": "How to decide whether an MCP server has actually learned",
}


OVERVIEW = r"""
# CAE loop learning overview

A validation artifact is evidence, not learning by itself.  A loop is closed
only when the evidence is converted into at least one durable MCP artifact:

1. a knowledge entry that an agent can retrieve,
2. a small test or validation script,
3. a lint/policy rule that prevents a repeated mistake, or
4. a reusable helper with focused verification.

Use this order at every slot boundary, before advancing to the next solver or
tool in the rotation:

1. Read each JSON/Markdown artifact and identify the physical identity or API
   contract that was checked.
2. Classify it as public-safe, private-only, or not stable enough.
3. Split the lesson into two lanes: public/open learning and source-tool
   learning.
4. Encode public-safe lessons in radia-mcp without private paths or solver
   provenance.
5. Encode source-tool lessons in the corresponding private MCP/converter lane
   when the artifact exposed API, parser, session, or workflow behavior.
6. Verify with the narrowest meaningful test for every lane that changed.
7. Say "learned" only after the MCP update and verification both exist.

Good loop artifacts teach students as well as agents: they name the governing
identity, state the tolerance, record the failure mode, and explain the next
gate to run.

Session diagnostics are a valid learning artifact when they unblock a solver
slot, but keep them separate from physics validation.  Record which existing
session was reused, whether direct MCP discovery failed, which fallback path
worked, and whether any solver process was started.  Do not turn a healthy
session-reuse result into a physics claim.

Do not wait until a full loop is over to learn.  A full-loop summary is only a
roll-up of slot-level learning that should already have been attempted.
"""


DUAL_LANE = r"""
# Dual-lane loop learning

The CAE loop is strongest when one artifact teaches twice.

Public/open lane:

* Extract solver-independent math, physics, geometry, meshing, and validation
  rules.
* Put those rules in public-safe radia-mcp knowledge, tests, lint, notebooks,
  or reusable helpers.
* Remove private paths, commercial solver provenance, and benchmark values.

Source-tool lane:

* Capture the tool-specific API or workflow behavior that made the artifact
  possible or caused the failure.
* Examples include session discovery, attach/reuse rules, file export
  preconditions, parser edge cases, unit conventions, table-column
  interpretation, and clearer failure messages.
* Keep private/commercial provenance in the owning private MCP or converter
  lane, not in public radia-mcp.

Both lanes can be useful at once.  A passivity artifact can become a generic
power-balance rule for radia-mcp and also a private session or export rule for
the source tool.  A force artifact can become a solver-independent Lorentz or
coenergy gate and also improve a private parser or automation message.
"""


MESH_GEOMETRY_VOL = r"""
# Mesh / geometry / .vol loop lessons

The reusable mesh lesson is: do not trust a mesh-export file merely because it
exists.  Validate the semantic inventory.

For Netgen `.vol` used as FEM/BEM input:

* Accept only triangle surface elements and tetrahedron volume elements in the
  first-order education path.
* Reject quad/hex/wedge/pyramid instead of silently converting them.
* Check `volumeelements > 0`; a boundary-only `.vol` is not a volume FEM mesh.
* Check boundary triangles against adjacent tetrahedron faces.  Orphan boundary
  triangles indicate an open or incorrectly exported surface.
* Preserve one-based node ids for readable FEM/BEM trace views.
* Validate boundary areas, vector areas, normals, pressure resultants, and
  moment resultants on simple boxes before using a complex model.

For CAD mass properties:

* Compare build123d/OCC volume and surface area against analytic boxes,
  cylinders, and spheres before using generated geometry downstream.
* For Cubit/Coreform exports, register material volume blocks before exporting
  solver-facing `.vol` files; otherwise the downstream volume inventory may be
  empty even if the surface inventory looks plausible.
* Sum per-surface areas when the goal is total boundary area; avoid assuming a
  similarly named volume API returns the boundary-area quantity you need.
* Run Cubit headless and wait for process completion before reading generated
  files.

Role split:

* Netgen/OCC is enough for tet-only meshes, especially the readable H1/HCurl
  and FEM/BEM teaching path.
* Cubit/Coreform slots should spend their budget on hex-led and mixed
  hex+pyramid+tet routes, because that is where Cubit adds unique value.
* For a mixed Cubit `.vol`, first run a semantic inventory gate that recognizes
  hex, pyramid, wedge, tet, quad, and triangle records.  Do not feed it to the
  tri/tet education parser and do not silently split pyramids into tets unless
  a downstream solver contract explicitly asks for that conversion.
"""


FORCE_MOMENT = r"""
# Force / moment loop lessons

Force gates should compare independent descriptions of the same quantity:

* Lorentz force on parallel conductors: force-per-length scales as
  `mu0*I1*I2/(2*pi*d)` and changes sign when either current is reversed.
* Virtual work/coenergy: compare force or torque against a finite-difference
  derivative of coenergy, but use absolute tolerance near zero crossings.
* Maxwell pressure/traction: integrate vector area and pressure/traction on a
  simple closed box first; uniform pressure on all faces must cancel.
* Moment resultants: always state the pivot.  A nonzero moment can vanish when
  the pivot is moved to the line of action.
* Torque waveforms: check periodicity, sign convention, and amplitude scaling
  before trusting a sampled table.

Do not use a single solver output as its own proof.  Each gate should include a
closed form, a conservation identity, a symmetry/antisymmetry check, or an
independent discretization identity.
"""


MOTOR_AIRGAP_TORQUE = r"""
# Motor air-gap torque loop lesson

For rotating-machine checks, a compact public validation gate is the cylindrical
air-gap Maxwell shear identity:

* `tau(theta) = Br(theta)*Bt(theta)/mu0`
* `T = r^2*L*integral tau(theta) dtheta`

For one harmonic pair,
`Br = Br0*cos(n theta)` and `Bt = Bt0*cos(n theta + phi)`, the average shear is
`0.5*Br0*Bt0*cos(phi)/mu0`.  That gives three useful checks:

* `phi = 0`: positive torque.
* `phi = pi/2`: zero torque, so use an absolute tolerance.
* `phi = pi`: negative torque with the same magnitude as the in-phase case.

Use this as a motor slot sanity gate before trusting a heavier FE torque
extraction.  It checks sign convention, phase convention, sector scaling,
radius/stack-length scaling, and the difference between mesh-independent
harmonic torque and mesh-sensitive weighted-stress extraction.

In radia-ngsolve, the executable helper is
`air_gap_shear_torque_from_angle_samples`: feed angle samples and Br/Bt samples,
then compare with the closed form above.
"""


ELECTROSTATIC_LAYERED_DIELECTRIC = r"""
# Electrostatic layered-dielectric loop lesson

For a parallel-plate stack with layers normal to the field, the normal electric
displacement is constant through all layers:

* `C = eps0*A/sum(d_i/eps_ri)`
* `D = eps0*V/sum(d_i/eps_ri)`
* `E_i = D/(eps0*eps_ri)`
* `Delta V_i = E_i*d_i`

This is a compact public gate for dielectric assignment, interface continuity,
terminal charge, and energy-density integration.  It is stronger than checking
capacitance alone: if a solver accidentally leaves every domain as vacuum, the
capacitance, interface potential, layer fields, and energy split all fail in a
diagnostic way.

In radia-ngsolve, use `layered_parallel_plate_stack_summary` to record the
analytic values and residuals from a solver artifact.
"""


ACOUSTIC_IMPEDANCE_POWER = r"""
# Acoustic impedance power loop lesson

For a planar acoustic impedance boundary, use this solver-independent gate:

* `R = (Zs - Z0)/(Zs + Z0)`
* `absorption = 1 - |R|^2`
* `P_boundary = 0.5*Re((1+R)*conj((1-R)/Z0))` for unit peak incident pressure

The gate catches three common mistakes:

* A purely reactive impedance should absorb zero power; use an absolute
  tolerance for this zero target.
* A matched impedance has `R=0`, absorption one, and boundary power equal to
  the incident power.
* Passive lossy impedances must have nonnegative absorption; active/negative
  resistance should be reported as a passivity violation, not silently accepted.

In radia-ngsolve, use `acoustic_impedance_reflection_summary` for single cases
and `acoustic_impedance_reflection_sweep_summary` for sweeps.  Keep the
reflection coefficient, absorption, boundary power, and residual in the
artifact so the next agent can see whether the failure is a sign convention,
phasor convention, or passivity issue.
"""


RF_ACOUSTIC_PASSIVITY = r"""
# RF / acoustic passivity loop lessons

For impedance, scattering, and radiation-pressure workflows, passivity is the
first sanity gate.

Acoustic impedance boundary:

* With normalized impedance `Zs/Z0`, reflection is
  `R = (Zs - Z0)/(Zs + Z0)`.
* Absorption is `1 - |R|^2` for passive boundaries.
* Purely reactive impedance should have zero absorption; use an absolute
  residual for this zero target and relative residuals for nonzero cases.

Two-port S-parameters:

* Reciprocity gate: `S12 == S21` when the modeled network is reciprocal.
* Passivity gate: the largest eigenvalue of `S^H S` must be no larger than 1.
* Power balance gate: for each unit incident port excitation,
  outgoing power plus absorbed power must equal one.
* Keep return loss, insertion loss, absorbed power, and passivity residual in
  the artifact so later agents can diagnose why a sweep failed.
* Treat one-port match quality as its own row: `S11` gives `|Gamma|`, VSWR,
  return loss, mismatch loss, reflected power, and transmitted power.  MATLAB
  teaching notebooks can use the same scalar gate as an optimization objective
  or constraint, but it should not be merged with `S21` insertion loss.
"""


MCP_CLOSURE = r"""
# MCP closure rule

Use these labels precisely:

* `collected`: a validation artifact exists.
* `distilled`: the artifact has been reviewed and turned into a public-safe or
  private-only lesson.
* `encoded`: the lesson has been added to MCP code, knowledge, tests, lint, or
  a reusable helper.
* `verified`: a focused test/lint/selftest has passed after the encoding.
* `learned`: encoded and verified.

If only cross-validation files were written, say "collected", not "learned".
This keeps the server honest and prevents repeated overclaiming.

Apply the labels per lane.  A public lesson can be learned while the
source-tool lesson is still only a candidate, and the report should say that
plainly.

Apply the labels per slot.  Advancing the rotation without at least recording
the MCP learning status makes later review harder and weakens the loop.
"""


_TOPIC_TEXT = {
    "overview": OVERVIEW,
    "dual_lane": DUAL_LANE,
    "mesh_geometry_vol": MESH_GEOMETRY_VOL,
    "force_moment": FORCE_MOMENT,
    "motor_airgap_torque": MOTOR_AIRGAP_TORQUE,
    "electrostatic_layered_dielectric": ELECTROSTATIC_LAYERED_DIELECTRIC,
    "acoustic_impedance_power": ACOUSTIC_IMPEDANCE_POWER,
    "rf_acoustic_passivity": RF_ACOUSTIC_PASSIVITY,
    "mcp_closure": MCP_CLOSURE,
}


def get_loop_learning_documentation(topic: str = "overview") -> str:
    """Return public-safe loop-learning guidance."""
    key = (topic or "overview").strip().lower()
    if key == "all":
        return "\n\n".join(_TOPIC_TEXT[k] for k in TOPICS)
    if key in _TOPIC_TEXT:
        return _TOPIC_TEXT[key]
    available = ", ".join(sorted([*TOPICS, "all"]))
    return f"Unknown topic {topic!r}. Available topics: {available}"
