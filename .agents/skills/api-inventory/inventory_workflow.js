// api-inventory: bundled fan-out workflow (read-only).
//
// Run via:  Workflow({scriptPath: ".agents/skills/api-inventory/inventory_workflow.js"})
// Requires the multi-agent opt-in (ultracode, or an explicit user request).
//
// It reads (never edits) the radia API surface, categorizes every API FAMILY per
// the "Reduce Proprietary API Surface" policy into one of five buckets, and
// returns a ready-to-commit markdown audit report.
//
// MAINTENANCE: the AREAS file lists below reflect the src/radia layout as of
// 2026-06-26. After running Step 1 (scout) in SKILL.md, refresh these lists if
// modules were added/removed/renamed. AREAS are HARDCODED here on purpose (a
// prior run showed Workflow `args` not propagating to agents) -- do not move them
// into `args`.

export const meta = {
  name: 'radia-api-inventory',
  description: 'Read-only inventory of the radia proprietary Python+pybind API surface, categorized per the "Reduce Proprietary API Surface" policy (plumbing/method/user-intent/deprecated)',
  phases: [
    { title: 'Inventory', detail: 'parallel readers categorize each API area by family' },
    { title: 'Synthesize', detail: 'merge into one inventory + recommendations + markdown report' },
  ],
}

const ROOT = 'S:/Radia/01_GitHub'

const FRAMEWORK = `
You are taking inventory (stock-take) of the **radia** package's PROPRIETARY API surface,
to apply the repo policy "Reduce Proprietary API Surface -- Plumbing to netgen/ngsolve,
Methods Stay" (CLAUDE.md). Categorize each API FAMILY (group by name prefix / role,
do NOT list 200 individual functions) into exactly one bucket:

- **plumbing-delete**: netgen/ngsolve (or MKL/OCC/GMSH/Cubit) already provides this.
  Delete from radia and delegate. Examples the policy names: mesh generation &
  representation, mesh I/O, geometry/CAD kernels, visualization & mesh export,
  generic linear algebra, and the geometry PRIMITIVES (ObjHexahedron/ObjRecMag/...)
  AS THE USER's hand-built-mesh API (replaced by .vol -> soft_iron_from_mesh + intent objects).
  Give the netgen/ngsolve replacement.
- **method-keep**: a genuine numerical METHOD NGSolve lacks (Radia's reason to exist):
  rad.Fld analytic open-boundary field, HDiv-VIM, axifem (Henrotte),
  DtN/FEM-Kelvin, PEEC, BEM, sparsesolv, HACApK, analytical_formulas, coil_builder
  (mesh-free Biot-Savart source), levitation/ECB, stream-function. KEEP.
- **method-demote**: a method that should be KEPT but DEMOTED from a user-facing pybind
  API to an internal C++/representation detail over time (e.g. ObjHexahedron as the
  internal element representation behind SoftIron). Keep != expose.
- **user-intent**: the intended USER layer (intent-based): SoftIron, Magnet, CoilBuilder,
  rad.Fld, rad.Solve, materials. KEEP/PROMOTE.
- **deprecated-drop**: already-removed shells, back-compat shims, dead/legacy, or
  duplicates to retire (e.g. CndLoop/FldVTS already removed; esim_vtk_export kept
  "not re-exported"; _b3d_shim; veriloga/legacy).

For each family return: group_name, members_sample (a few representative names),
approx_count, bucket, one-line rationale, replacement_or_note.

IMPORTANT: This is READ-ONLY inventory. Do NOT edit any file. HACApK
(src/ext/HACApK/, src/core/rad_hacapk.*) is **method-keep AND under active development
by another agent** -- categorize it but flag "do-not-touch (active dev)".
Read the actual files under ${ROOT}; ground every call in what the code exposes.
`

const ITEM = {
  type: 'object', additionalProperties: false,
  properties: {
    group_name: { type: 'string' },
    members_sample: { type: 'array', items: { type: 'string' } },
    approx_count: { type: 'integer' },
    bucket: { type: 'string', enum: ['plumbing-delete', 'method-keep', 'method-demote', 'user-intent', 'deprecated-drop'] },
    rationale: { type: 'string' },
    replacement_or_note: { type: 'string' },
  },
  required: ['group_name', 'members_sample', 'approx_count', 'bucket', 'rationale', 'replacement_or_note'],
}
const AREA_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    area: { type: 'string' },
    groups: { type: 'array', items: ITEM },
    notes: { type: 'string' },
  },
  required: ['area', 'groups', 'notes'],
}

const AREAS = [
  { key: 'pybind-core', focus:
    'The C++ pybind surface re-exported by `import radia`. Read src/lib/radia_pybind.cpp (~172 .def), src/lib/radentry.cpp, src/lib/rad_cln_api.cpp, src/lib/rad_peec_matrices_api.cpp. Group by family: Obj* (geometry primitives & containers), Fld*/field eval, Mat* (materials), Rlx*/Solve/SolverConfig, Uti*, Drw*/drawing, RadiaField CF, IMA, background, PEEC api, CLN api. Decide bucket per family.' },
  { key: 'py-core-methods', focus:
    'Top-level python METHOD/solver modules under src/radia/. Read: kelvin_solver.py kelvin_source.py kelvin_material.py kelvin_geometry.py kelvin_validate.py kelvin_identify_ngsolve.py scalar_potential_solver.py vector_potential_solver.py dielectric_solver.py equivalence_source.py cohomology.py cohomology_cut.py infinite_element.py biot_savart.py analytical_magnet.py cylindrical_magnet.py round_bodies.py ima_field.py clebsch_potential.py em_material.py energy_play_model.py hysteresis_io.py netgen_mesh_import.py gmsh_post_export.py soft_iron.py step_mesh_builder.py scalar_bie_sibc.py. Categorize each module group.' },
  { key: 'py-peec-esim-coil', focus:
    'PEEC / BEM / ESIM / coil python modules under src/radia/. Read: peec_topology.py peec_coupled.py peec_matrices.py peec_msc_schur.py peec_hacapk_solver.py peec_mesh_import.py peec_proximity.py peec_shield.py peec_shielded.py peec_bundle.py prima_hacapk.py lanczos_reduction.py bem_sibc_solver.py fasthenry_parser.py filament_bundle.py veriloga_generator.py workpiece_surface.py esim_cell_problem.py esim_workpiece.py esim_coupled_solver.py esim_anderson.py esim_hantila.py esim_multiport.py esim_vtk_export.py coil_builder.py coil_from_cad.py coil_from_step.py coil_geometry.py coil_profile.py coil_profile_occ.py coil_spec.py coil_topology.py coils.py round_bodies.py. Categorize.' },
  { key: 'py-apps-panels', focus:
    'Application / panel / notebook / analysis python under src/radia/. Read: analysis.py em_design.py em_notebook.py ih_design.py ih_notebook.py ih_optimize.py ih_pipeline.py ih_claude_proposer.py motor_design.py motor_notebook.py pcb_design.py pcb_notebook.py streamfunction_design.py streamfunction_notebook.py streamfunction_volume.py stream_function.py notebook_workbench.py panel_design_common.py radia_gui_base.py radia_ih.py radia_em.py radia_motor.py radia_pcb.py radia_streamfunction.py radia_ngsolve.py install_panels.py setup_cubit.py _heat_panel.py _b3d_shim.py. Categorize (most are application layer, not core proprietary methods).' },
  { key: 'subpackages-pyd', focus:
    'radia subpackages + compiled submodule .pyds. Read src/radia/__init__.py (the re-export hub + 2-layer wrappers Solve/ObjCnt/SolverConfig/UtiDelAll/SoftIron/set_demag_backend) and survey subpackages: src/radia/vim/ (HDiv-VIM), src/radia/bem/, src/radia/open_boundary/, src/radia/maglev/, src/radia/analytical_formulas/, src/radia/tools/. Plus the .pyd modules: axifem.pyd, cln_core.pyd, peec_matrices.pyd, sparsesolv_ngsolve.pyd. Categorize each subpackage/submodule as a family.' },
]

phase('Inventory')
const areaResults = await parallel(AREAS.map((a) => () =>
  agent(`${FRAMEWORK}\n\nYOUR AREA: ${a.key}\n${a.focus}`,
    { label: `inv:${a.key}`, phase: 'Inventory', effort: 'medium', schema: AREA_SCHEMA })
))
const valid = areaResults.filter(Boolean)
log(`inventory areas done: ${valid.length}/${AREAS.length}; total groups=${valid.reduce((n, r) => n + (r.groups?.length || 0), 0)}`)

phase('Synthesize')
const SYNTH_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    exec_summary: { type: 'string' },
    counts_by_bucket: { type: 'object', additionalProperties: { type: 'integer' } },
    top_deletion_candidates: { type: 'array', items: { type: 'string' } },
    top_demotion_candidates: { type: 'array', items: { type: 'string' } },
    recommended_user_layer: { type: 'array', items: { type: 'string' } },
    safe_now_vs_blocked: { type: 'string' },
    markdown_report: { type: 'string' },
  },
  required: ['exec_summary', 'counts_by_bucket', 'top_deletion_candidates', 'top_demotion_candidates', 'recommended_user_layer', 'safe_now_vs_blocked', 'markdown_report'],
}
const synth = await agent(
  `${FRAMEWORK}\n\nThe parallel inventory of all areas produced:\n"""${JSON.stringify(valid, null, 2)}"""\n\n` +
  `Synthesize ONE coherent API inventory (stock-take). Produce:\n` +
  `- exec_summary: 4-6 sentences on the shape of the surface and the biggest reduction opportunities.\n` +
  `- counts_by_bucket: {plumbing-delete, method-keep, method-demote, user-intent, deprecated-drop} family counts.\n` +
  `- top_deletion_candidates: concrete plumbing/deprecated families to remove first (lowest risk, clear netgen/ngsolve replacement), each "name -- why -- replacement".\n` +
  `- top_demotion_candidates: methods to keep but un-pybind / move behind the intent layer, each "name -- why".\n` +
  `- recommended_user_layer: the lean intent-based API the user SHOULD see (SoftIron/Magnet/CoilBuilder/Fld/Solve/materials + the core method entry points).\n` +
  `- safe_now_vs_blocked: which removals are safe to do NOW (read-only areas) vs must wait -- HACApK is under active dev by another agent (do-not-touch); flag anything else risky.\n` +
  `- markdown_report: a COMPLETE ready-to-commit markdown audit doc with a per-bucket table (group | members | count | rationale | replacement/action), an exec summary, and a "phased reduction plan" section. Use ASCII (cp932-safe). This is the deliverable.\n`,
  { label: 'synthesize:inventory', phase: 'Synthesize', effort: 'high', schema: SYNTH_SCHEMA },
)
return synth
