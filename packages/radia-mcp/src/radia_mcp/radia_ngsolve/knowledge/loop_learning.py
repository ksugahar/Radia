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
    "fem_bem_trace_orientation": "FEM/BEM trace packages with normal-flux orientation evidence",
    "fem_bem_solver_report": "FEM/BEM coupled solves with solver-report identity",
    "bem_demag_source_mesh": "BEM demag source-balance surface mesh identity",
    "electrostatic_layered_dielectric": "Layered dielectric stack capacitance, D-continuity, and energy gates",
    "acoustic_impedance_power": "Acoustic impedance reflection, absorption, and boundary power",
    "rf_acoustic_passivity": "Acoustic/RF passivity and power-balance identities",
    "geometric_time_integration": "Energy-drift checks for geometric time integration teaching gates",
    "source_native_seed_queue": "How to start multi-tool loop slots from source-native examples without leaking provenance",
    "autonomous_basic_learning": "How to process a full source-native queue into basic learning rows and solver-ready follow-ups",
    "artifact_feedback": "How cross-validation JSON and notebook/result artifacts become MCP knowledge",
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

For force rows replayed from an external 2D solver, bind the numeric result to
the input model artifact as well as the solution and postprocess artifacts.
The reusable package should carry a model-input artifact id, digest, and path
next to the loaded-solution id and output-table digest.  That prevents a copied
force table from passing an analytic value check while it is actually tied to
stale geometry, block labels, or source definitions.
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

When a sampled air-gap torque scalar is promoted from a solver export to a
notebook or optimizer, keep the upstream model input package with it.  The
torque-result package should repeat the field-table artifact, sample-grid
artifact, integration method, torque-output artifact, and the project/model
input artifact id, digest, and path.  This prevents a plausible Maxwell-shear
torque value from being joined to stale project geometry, material, or current
definitions.

In radia-ngsolve, the executable helper is
`air_gap_shear_torque_from_angle_samples`: feed angle samples and Br/Bt samples,
then compare with the closed form above.
"""


FEM_BEM_TRACE_ORIENTATION = r"""
# FEM/BEM trace orientation lesson

A first-order tri/tet FEM/BEM handoff is not complete just because the trace
matrix is one-hot and the boundary rows have the expected ids.

Keep `trace_basis_schema_id` with the trace package.  It binds the volume H1
nodal basis, boundary/surface P1 basis, and compact trace-row ordering before
the row is reused by a BEM kernel, notebook result, or optimizer-visible
postprocess table.

Keep the normal-flux evidence as a separate artifact:

* `normal_flux_artifact_id`: the report that checked stored triangle orientation
  against the outward-from-volume convention.
* `normal_flux_digest`: the actual orientation/sign table or flux-balance
  report consumed by the notebook.
* `normal_flux_convention`: for example `outward_from_volume`.

This matters for scalar potential, acoustic, and low-frequency BEM teaching
examples.  A stale normal-flux sign report can leave trace rows, kernel family,
and assembly/quadrature ids looking valid while the Neumann/source sign is
wrong.  Use `netgen_vol_first_order_fem_bem_trace_package_handoff` with
`require_normal_flux_artifact=True` before reusing normal derivatives, flux
integrals, or surface-source rows.
"""


FEM_BEM_SOLVER_REPORT = r"""
# FEM/BEM solver-report identity lesson

A FEM/BEM teaching package is still not reusable as a solved result just
because the trace matrix, BEM kernel manifest, assembly/quadrature ids, normal
orientation, and coupled-system digest are present.

Keep the linear solve report as its own evidence:

* `linear_solver_report_artifact_id`: which solve report was generated.
* `linear_solver_report_digest`: the concrete report content consumed by the
  notebook/result package.
* `linear_solver_name`: for example a minimum-norm rank-deficient teaching
  solve, direct factorization, or Krylov method.
* `linear_solver_tolerance`: the target residual/solve tolerance.
* `linear_solver_residual_norm`: the measured residual norm.
* `linear_solver_iteration_count`: iteration count or `1` for a direct solve.
* `result_artifact_id`: the solved-result package or executed notebook result
  that consumed the solve report.
* `run_started_at` or `created_at_utc`: an ISO-like timestamp for when the
  result package was produced.
* `tool_version`: the solver or teaching environment version that produced the
  result package.
* `notebook_source_artifact_id`, `notebook_source_digest`, and
  `notebook_source_path`: the notebook or source script revision that produced
  the visible result package.
* `parameter_set_artifact_id`, `parameter_set_digest`, and
  `parameter_set_path`: the initial values or design variables used by the
  notebook/result package.
* `objective_observable_id` and `objective_observable_family`: the scalar
  objective or teaching observable that an optimizer, panel, or replayed
  notebook will consume.
* `timing_breakdown_s`: a compact timing ledger with about four dominant stages
  so later notebooks can see whether mesh read, trace assembly, solve, or JSON
  write dominated the teaching run.

This keeps a valid coupled-system artifact from being paired with a stale solve
report, stale notebook/source script, stale parameter defaults, or an
optimization objective from a different observable family.  Use
`netgen_vol_first_order_fem_bem_trace_package_handoff` with
`require_linear_solver_report=True` and, when notebooks or optimization reuse
the row, `require_parameter_set_artifact=True`.  Then keep the result artifact
id, timestamp, version, notebook/source identity, parameter-set identity,
objective-observable identity, and compact timing ledger with archived notebook
or JSON results.
"""


BEM_DEMAG_SOURCE_MESH = r"""
# BEM demag source-balance surface mesh lesson

For PM demagnetization and magnetic-charge BEM workflows, the source-balance
row is not only a scalar residual.

Keep these identities separate:

* `surface_mesh_id`: which closed PM surface mesh was used.
* `surface_mesh_digest`: the actual mesh artifact, node/face order, and
  surface-row set used by the BEM source balance.
* `surface_row_count`: a quick row-count guard against truncated surface
  ledgers.
* `source_balance_artifact_id` and `source_balance_digest`: the computed
  source-balance evidence.
* `source_convention`: for example `sigma_m = M dot n`, with the normal
  convention recorded explicitly.

The mesh identity and the source-balance result identity are twins, but they
are not the same artifact.  A stale surface mesh can produce a plausible
near-zero source-balance residual, especially for symmetric PM teaching cases.

Before promoting a demag-margin notebook or comparing against an open solver,
run `pm_demag_margin_screening_package_gate` with the expected BEM surface mesh
id, mesh digest, row count, source-balance artifact, source-balance digest, and
source convention.  Negative controls should include stale mesh digest, wrong
row count, stale source-balance digest, and wrong source convention.
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
* Keep the frequency-axis identity with S-parameter rows: `frequency_grid_id`,
  `frequency_grid_digest`, row count, selected row index, and selected
  frequency should travel together before passivity, equivalent-circuit, BEM,
  or notebook reuse.
* Keep the project/model input artifact identity with solver-ready
  Touchstone/S-parameter manifests: `model_input_artifact_id`,
  `model_input_digest`, and `model_input_path` should travel with the port,
  grid, row, output, and timing evidence.  A valid S-parameter table is not the
  same evidence package after geometry, materials, ports, or solver setup were
  regenerated.
* Keep the exact Touchstone export recipe with the same package:
  `export_recipe_artifact_id`, `export_recipe_digest`, and
  `export_recipe_path` identify the macro/script/postprocess recipe that
  produced the raw file or selected row.  The recipe identity is part of the
  operator/source ledger, just like reference planes, port basis, and current
  conventions in electromagnetic papers.
* Keep the postprocessed Touchstone output table schema with the same package:
  `touchstone_output_schema_id`, `touchstone_output_columns`, and
  `touchstone_output_units` distinguish a full two-port table from scalar-only
  rows such as a lone `S21` objective export.  A raw/output artifact digest can
  still be correct while the table layout is stale, reordered, or in different
  units.
* Keep the port-mode-basis schema with the port-mode-basis value:
  `touchstone_port_mode_basis_schema_id` records how single-ended versus
  mixed-mode power waves, reference orientation, and port order are interpreted.
  The value `single_ended_power_wave_modes` alone is not enough provenance for
  notebook, BEM, optimization, or equivalent-circuit reuse.
* Pair solver-ready Touchstone/S-parameter manifests with
  `solver_result_artifact_provenance_timing_gate`: record parseable
  `created_at_utc` / `run_date_utc`, solver and MCP versions, and about four
  dominant `timing_breakdown_s` stages before notebook, equivalent-circuit,
  BEM, or optimization reuse.
* Use the same discipline for solver-derived result tables from frequency,
  design, or parameter sweeps: `sweep_axis_id`, `sweep_axis_digest`, and
  `sweep_axis_row_count` belong beside the solution-data and exported-table
  artifacts.  A copied table with the wrong parameter grid is a different
  evidence package even when the dataset, solution tag, columns, and units
  still look plausible.
* When the table is produced by a parameter study, notebook panel, or
  optimization objective, keep `parameter_set_artifact_id`,
  `parameter_set_digest`, `parameter_set_path`, `objective_observable_id`, and
  `objective_observable_family` with the result-table package.  Column names,
  units, solution data, and sweep axis can all be correct while the row still
  belongs to a stale design-variable set or a different objective family.
* Solver-derived result tables also need solver-configuration identity:
  `solver_configuration_artifact_id`, `solver_configuration_digest`,
  `solver_sequence_tag`, `linear_solver`, and `relative_tolerance`.  A copied
  table produced with a different solver setup is different evidence even when
  the solution tag, sweep axis, output artifact, and numeric rows still look
  plausible.
* Replayable result tables should keep schema identities separate:
  `result_table_schema_id` for layout, `physics_convention_schema_id` for
  physical meaning, `result_postprocess_row_convention_schema_id` for row
  reduction/objective semantics, and `result_component_basis_schema_id` for
  component columns, complex representation, and basis/normalization ordering.
  A stale component basis can keep the same values and row convention while
  changing how the columns should be interpreted.
* Keep return loss, insertion loss, absorbed power, and passivity residual in
  the artifact so later agents can diagnose why a sweep failed.
* Treat one-port match quality as its own row: `S11` gives `|Gamma|`, VSWR,
  return loss, mismatch loss, reflected power, and transmitted power.  MATLAB
  teaching notebooks can use the same scalar gate as an optimization objective
  or constraint, but it should not be merged with `S21` insertion loss.
"""


GEOMETRIC_TIME_INTEGRATION = r"""
# Geometric time-integration loop lesson

For teaching dynamics, use a Hamiltonian toy problem before a heavy field solve.
The harmonic oscillator is enough to expose a real numerical-analysis lesson:

* explicit Euler is a useful negative control because its energy drifts upward
  on a fixed time grid;
* symplectic Euler does not conserve the exact energy pointwise, but its energy
  error remains bounded and oscillatory;
* implicit midpoint is symplectic and, for the linear oscillator, preserves the
  quadratic Hamiltonian to roundoff.

Record every method on the same time grid: method name, step size, step count,
omega, initial energy, final relative drift, and maximum relative energy drift.
Then run `geometric_integrator_energy_drift_gate` before claiming that a MATLAB
or notebook dynamics example demonstrates a geometric integrator.  The gate is
also a compact way to teach why preserving structure can matter more than only
reducing the local truncation error.
"""


SOURCE_NATIVE_SEED_QUEUE = r"""
# Source-native seed queue lesson

A multi-tool validation rotation should not be seeded with generated toy cases
when the goal is to improve tool-aware MCP behavior.  Every source-tool slot
should begin from a real source-native example:

* official tutorials, model pages, training decks, or application-library
  models from the owning ecosystem;
* community or manual examples for open/community tools;
* upstream package examples for CAD/mesh libraries;
* lab examples only when they are directly derived from those source-native
  workflows and are kept in the private provenance lane.

Generated scripts are still useful, but only as replay harnesses, reduced
public analogues, schema gates, or negative controls.  They are not the lesson
seed for the source-tool lane.

For each queued slot, record:

* `source_native_example`: a scrubbed id or internal pointer to the real source;
* `source_native_type`: local artifact, public URL, training model, manual
  example, or upstream example;
* `lesson_axis`: the API, file-format, physics, meshing, or solver behavior to
  extract;
* `intended_validation`: the later solver or open-reference gate;
* `learning_lanes`: public/open and source-tool states, using candidate until
  a real MCP code/knowledge/test change is verified.

Public radia-mcp should store only the generic rule and the open validation
candidate.  Private provenance, commercial model names, file-system paths, and
benchmark values stay in the owning internal artifact.  This lets a 20-round
queue be useful immediately as a review ledger while keeping the word
`learned` reserved for encoded and verified MCP changes.
"""


AUTONOMOUS_BASIC_LEARNING = r"""
# Autonomous basic learning pass

Use an autonomous basic-learning pass when a source-native queue is large
enough that manual per-slot review becomes the bottleneck.  The pass should
process every queued slot and leave a concrete machine-readable record for
each one:

* `lesson_family`: mesh/geometry, FEM/BEM, force/torque, RF/acoustic,
  optimization, thermal/eddy, session/API, source-MCP policy, or general
  source-native.
* `source_present` and `required_fields_present`: the seed is usable before
  solver-ready promotion begins.
* one lightweight `computed/reference/tolerance/pass` row per slot, using a
  public analogue when possible and a metadata-contract row for session/API or
  policy slots.
* `learning_lanes.public`: verified only when the public analogue row and
  queue gates pass.
* `learning_lanes.source_tool`: candidate for commercial/private source-tool
  lanes until the owning MCP or converter receives a focused edit and test.
* `next_action`: the solver-ready or private-MCP follow-up produced by the
  autonomous pass.

This makes a 160-slot queue actionable without pretending all 160 commercial
or live solvers have run.  The basic pass answers "is every slot classified,
checked, and ready for the next stage?"  Heavy source-tool execution remains a
separate solver-ready queue.

In radia-mcp, the public helper is
`loop_autolearn.build_autonomous_basic_learning_artifact`; the validation CLI
is `validation/loop_learning/autonomous_basic_learning.py`.  Pair its output
with the computed/reference row gate, `source_native_seed_queue_gate`, and the
MCP feedback artifact gate.
"""


ARTIFACT_FEEDBACK = r"""
# Cross-validation and notebook artifact feedback

Cross-validation JSON, executed notebooks, and result sidecars are not MCP
knowledge until the distilled lesson can be retrieved or checked by the MCP
server.

Use this promotion sequence:

1. Gate the result provenance: schema, run date, versions, result artifact id,
   and the dominant `timing_breakdown_s` stages.
2. Gate the output-table metadata: schema id, columns, units, independent axis,
   row convention, physics convention, component basis, and output artifact.
3. If a notebook produced or consumed the result, bind the notebook source:
   `notebook_source_artifact_id`, `notebook_source_digest`, and
   `notebook_source_path`.
4. Distill one public-safe lesson: the physical identity, API/file-format
   contract, or failure mode that should shape the next MCP answer.
5. Record the MCP target: knowledge topic, helper, lint, validation gate, or
   notebook/panel policy that was updated.
6. Record focused verification: pytest, policy lint, notebook execution, or a
   solver-free gate that passed after the MCP update.
7. Set `learning_lanes.public` to `verified` only after steps 4-6 exist.

Use `cross_validation_artifact_to_mcp_feedback_gate` as the final feedback
check.  It composes with `solver_result_artifact_provenance_timing_gate` and
keeps a hard line between collected evidence and learned MCP knowledge.

Recommended artifact fields:

* `public_lesson` or `mcp_feedback.public_summary`
* `learning_targets` containing a public MCP target such as `radia-mcp`
* `verification.public` with the focused command/result
* `result_artifact_id` or `notebook_result_artifact_id`
* `notebook_source_artifact_id`, `notebook_source_digest`,
  `notebook_source_path` when a notebook is involved
* `versions`, `created_at_utc`, `execution.run_date_utc`,
  `timing_breakdown_s`, and output schema/columns/units

If these fields are absent, say the artifact is collected or distilled, not
learned.  This is the mechanical rule that turns validation archives into
MCP behavior.
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

When a direct MCP connector misses an already-running shared solver session,
close the slot with session evidence rather than starting another process.
Record the direct-discovery status, `matlab.engine.find_matlab()` engine list,
selected shared engine name, successful shared-engine eval status, and
`started_new_process=false` / `killed_process=false`.
Use `shared_solver_session_health_gate` to keep a reusable session-health check
separate from physics residuals or FEM/BEM result values.  A pure MATLAB or
notebook validation can still use the selected shared engine when that is the
session of record; the artifact must say which engine was used and that no
solver process was started or killed.  If the selected engine name is absent
from the recorded `find_matlab()` list, keep the slot at `needs_attention` even
when passive diagnostics found a solver-owned MATLAB process.

For an external solver slot, a visible shared engine name is not enough to call
the slot live.  Close the discovery false-negative only after code executed in
the selected shared engine reports both a successful solver-session attach and
a solver-native preflight verdict.  Record those two payloads as session-health
evidence, not as physics validation.
"""


_TOPIC_TEXT = {
    "overview": OVERVIEW,
    "dual_lane": DUAL_LANE,
    "mesh_geometry_vol": MESH_GEOMETRY_VOL,
    "force_moment": FORCE_MOMENT,
    "motor_airgap_torque": MOTOR_AIRGAP_TORQUE,
    "fem_bem_trace_orientation": FEM_BEM_TRACE_ORIENTATION,
    "fem_bem_solver_report": FEM_BEM_SOLVER_REPORT,
    "bem_demag_source_mesh": BEM_DEMAG_SOURCE_MESH,
    "electrostatic_layered_dielectric": ELECTROSTATIC_LAYERED_DIELECTRIC,
    "acoustic_impedance_power": ACOUSTIC_IMPEDANCE_POWER,
    "rf_acoustic_passivity": RF_ACOUSTIC_PASSIVITY,
    "geometric_time_integration": GEOMETRIC_TIME_INTEGRATION,
    "source_native_seed_queue": SOURCE_NATIVE_SEED_QUEUE,
    "autonomous_basic_learning": AUTONOMOUS_BASIC_LEARNING,
    "artifact_feedback": ARTIFACT_FEEDBACK,
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
