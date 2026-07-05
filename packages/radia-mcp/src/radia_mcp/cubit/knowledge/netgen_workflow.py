"""
Netgen/NGSolve high-order curving workflow knowledge base.

This is the most important knowledge module. It covers:
- Export path: APREPRO command (C++) `export netgen`
- Workflow: Cubit geometry -> mesh -> export netgen -> NGSolve Mesh()
- CallbackGeometry and ACIS curving (compact_netgen)
- .vol file as sole interface between Cubit and NGSolve
- Troubleshooting guide
"""

WORKFLOW_OVERVIEW = """
# Netgen High-Order Curving: Workflow Overview

## Two Export Paths (produce identical results)

### Path A: APREPRO Command (recommended, fast)
```
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
# -> mesh.vol (with curvedelements section) + mesh.vol.json (CAD reference)
```
Uses NetgenCurver (compact_netgen C++ static link). No Python, no DLL dependency.

### Path B: Python (reference, deprecated)
Path B (`extract_curved_mesh`) has been removed. Use Path A (`export netgen`) for all workflows.

Both paths use **CallbackGeometry** to delegate surface/edge projection
to Cubit's ACIS kernel via `closest_point_trimmed`. No STEP files,
no OCC geometry, no SetGeomInfo needed.

## .vol as Sole Interface

```
Cubit (ACIS geometry) -> export netgen -> .vol (self-contained)
                                           |
                              NGSolve: Mesh("mesh.vol")
                              (no Cubit, no STEP needed)
```

The .vol file contains: mesh points, volume elements, surface elements,
material labels, boundary labels, and curvedelements section (high-order
curving coefficients). NGSolve reads it without any geometry file.

### Continuous-loop O-grid hex sphere gate

For Coreform/Cubit-led hex validation, a compact live gate is:
`create sphere radius 1`, `volume 1 scheme sphere`, `volume 1 size 0.25`,
`mesh volume 1`, `block 1 add volume all`, then export Netgen `.vol` at
orders 1, 2, and 3.  The expected inventory is 56 hexes and no tet/wedge/
pyramid elements.  The CAD volume is exactly `4*pi/3`; NGSolve should load
each `.vol` directly and integrate `CF(1)` with high-order quadrature.  Typical
2026-06-29 live results were order-1 rel err 0.2336, order-2 rel err 0.00211,
and order-3 rel err 0.00131.  Do not call `mesh.Curve()` after reading this
high-order `.vol`; the curving is already baked into the file.

On Windows PowerShell, prefer `coreform_cubit.com -nographics -batch script.py`
or the lab launcher background path when you need to wait for batch completion
and capture logs.  `coreform_cubit.exe` can behave as a GUI stub and return
immediately without running the Python batch script.

### Continuous-loop O-grid hex sphere eigenvalue gate

The same O-grid sphere is also a compact PDE-spectrum gate, not only a volume
gate.  Load the order-3 `.vol` as-is in NGSolve, build an H1 Dirichlet
generalised eigenproblem with `laplace_dirichlet_eigenvalues(mesh, n_modes=2,
order=3, shift=-1.0)`, and compare against the ball identities
`lambda_1=(pi/R)^2` and `lambda_l1=(4.493409457909064/R)^2`.  The verified
2026-06-29 live gate at `R=1` returned `lambda1=9.86268248871533` and
`lambda_l1=20.178703158644527`, with relative errors `7.01e-4` and `5.96e-4`.
This catches export/curving issues that pure volume integration can miss, and
it confirms that Cubit high-order hex `.vol` files can drive a real scalar
Dirichlet eigenvalue solve.

### Continuous-loop mapped hex brick volume/area gate

For a planar all-hex sanity check, use `create brick x 2 y 3 z 4`,
`volume 1 scheme map`, `volume 1 size 0.5`, `mesh volume 1`, and
`block 1 add volume all`, then `export netgen "... .vol" order 1 overwrite`.
The expected inventory is 192 hexes and no tet/wedge/pyramid elements.  Cubit's
CAD volume and summed `get_surface_area(surface_id)` values should be exactly
24 and 52.  NGSolve should load the exported `.vol` and integrate volume and
boundary area to machine precision; typical 2026-06-29 live results were
volume rel err `7.25e-15` and surface-area rel err `2.05e-15`.

Batch-script gotcha: `coreform_cubit.com -nographics -batch script.py` can play
Python line by line.  Avoid multi-line `dict(...)`, loops, or parenthesized
blocks in quick validation scripts; write one assignment per physical line or
use a normal Python process to generate the script.

### Continuous-loop Cubit mass-property sidecar gate

Before Cubit exports a mesh or hands CAD rows to build123d/CST/radia-ngsolve
cross-checks, save a compact mass-property sidecar: volume name, Cubit
`get_volume_volume`, summed `get_surface_area(surface_id)` over the volume's
surfaces, and bounding-box size from `get_total_bounding_box`.  Replay it with
`cubit_mass_property_sidecar_gate`.

This is deliberately earlier than `.vol` inventory.  Volume is the common CAD
currency, but volume alone can hide a wrong scale or clipped face; the sidecar
should also carry total surface area and bbox dimensions whenever they are
available.  For a simple `1.5 x 2.0 x 0.75` mapped brick the expected values
are volume `2.25`, surface area `11.25`, and bbox size `[1.5, 2.0, 0.75]`.
Slot210 adds the unit contract for Cubit/build123d/CST CAD cross-checks:
record `length_unit`, `area_unit`, and `volume_unit` (or a `units` mapping)
with the sidecar, then pass the expected units into
`cubit_mass_property_sidecar_gate`.  A perfect numeric volume is still
`needs_attention` if one CAD lane reports `mm^3` while the replay gate expects
`m^3`.

Slot274 adds material/block label identity for mixed routes.  If a hex core,
pyramid transition, and tet region are exported as one solver-ready package,
the mass-property sidecar should keep rows named like `hex_core`,
`pyramid_transition`, and `tet_region`, and those names should match the `.vol`
inventory material labels.  A transition material may report zero CAD volume in
the sidecar, but only when that exact row name is passed as an allowed zero
measurement; missing or stale material names stay `needs_attention`.

### Continuous-loop Cubit export package identity gate

Before docs notebooks, panel notebooks, or solver-ready validation consume a
Cubit export, keep the files as one named package.  The `.vol`, companion
`.vol.json`, raw Coreform batch result, and optional mass-property sidecar
should carry a stable `export_id`, `geometry_id`, order, and routing hint.
Replay that bookkeeping with `cubit_export_package_identity_gate`.

This gate catches a quiet but expensive mistake: a valid `.vol` can be paired
with an old sidecar or a raw JSON from a different geometry.  The package gate
requires the `.vol.json` path to pair with the `.vol`, requires raw result
presence, checks shared `export_id`/`geometry_id`, and can verify the inventory
source plus `cubit_hex_or_mixed_path` routing hint.

Slot306 adds the observable identity to the same package.  Keep
`export_observable_id` and `export_observable_family` next to the emitted
`export_output_artifact_id`/digest/path.  The output artifact says which file
was emitted; the observable identity says what that file is being used to
measure or prove, such as `netgen_vol_inventory`, `quality_distribution`, a
sidecar material map, or a solver-ready mesh contract.  This prevents a fresh
`.vol` file from being paired with a stale inventory/quality/result-table
interpretation.

### Continuous-loop headless batch quality package gate

For hex-led Cubit work, keep the headless batch result and quality replay rows
as one package too.  The raw batch JSON should record `export_id`,
`geometry_id`, Cubit version, `-nographics -batch` command line, output paths,
and `pass=true`; the quality replay row should carry the same
`export_id`/`geometry_id`, element type, count, and status.  Replay the pair
with `cubit_headless_batch_quality_package_gate` before a docs notebook or
solver-ready run consumes the quality evidence.

This gate catches a second stale-row class: a valid quality distribution can be
reused with the wrong mesh, or a GUI command line can sneak into a headless
automation slot.  Treat missing headless command evidence, mismatched
`geometry_id`, and zero quality count as `needs_attention`.

When the run also exported a `.vol`, pass the parsed inventory into
`cubit_headless_batch_quality_package_gate` as `export_inventory`.  The gate
then binds the inventory source path, positive `volumeelements`, volume-kind
count such as `hex: 64`, and `cubit_hex_or_mixed_path` routing hint to the same
headless package.  This is the slot194 lesson: a batch raw JSON and a quality
row are not enough once the export file is handed to a notebook or solver; the
export inventory must agree with the same mesh count and route.

Slot234 adds the route-separation check for the MATLAB/Gypsilab `.vol` policy:
if a quality package claims a Cubit hex-led route, the parsed inventory must
not be tri/tet-only and must contain the quality element kind being replayed.
This prevents a valid hex quality row from being paired with a Netgen
tri/tet-only export, and keeps Cubit hex/mixed evidence out of the educational
tri/tet `.vol` reader.

Slot242 adds process-evidence checks for headless Cubit automation.  When the
raw batch JSON records process metadata, `cubit_headless_batch_quality_package_gate`
now requires `process_mode=headless_batch`, `-nographics -batch`, a disabled GUI
daemon, a batch script path that appears in the command line, a recorded process
exit code, and either `exit_code=0` or an explicit headless startup/warning note
for the known valid-artifact/nonzero-exit path.  This keeps batch/process health
separate from mesh quality while preventing a persistent GUI daemon or stale
script path from becoming solver-ready evidence.

Slot314 tightens the noisy-exit path.  If a headless Cubit batch ends nonzero
after writing valid `.vol`/quality artifacts, record
`process_exit_policy=artifact_evidence_over_process_exit` and
`solver_ready_claimed=false`.  This policy does not make a nonzero process exit
healthy by itself; it only says that archived artifact evidence can be replayed
while solver-ready promotion remains blocked until inventory, quality, and
process evidence agree.  A nonzero exit with `solver_ready_claimed=true` or a
missing/wrong `process_exit_policy` should fail before downstream notebooks or
LLM-driven solver imports reuse the mesh.

Keep the batch process exit status separate from the archived raw JSON.  A
headless run can produce a valid mesh-quality JSON and still end with a startup
or plugin freshness warning.  Record that warning explicitly, and do not promote
plugin-specific export evidence until the export plugin freshness has been
checked; simple Cubit API quality replay can still be useful when the raw JSON
and package identity gates pass.

Slot369 adds a mesh-quality ledger identity gate.  After
`cubit_quality_distribution_gate` proves that a scaled-Jacobian/Jacobian list is
healthy, replay `cubit_mesh_quality_ledger_identity_gate` before reusing the
row in notebooks or solver-ready imports.  The ledger should carry
`mesh_quality_artifact_id`, `mesh_quality_digest`, `quality_metric_set_id`,
`mesh_artifact_id`, `mesh_digest`, `export_id`, `geometry_id`, routing hint,
element-type counts, `min_scaled_jacobian`, and `negative_jacobian_count`.
Negative controls should include a stale quality digest, a nonzero
negative-Jacobian count, and a tri/tet-only inventory paired with a Cubit
hex-led quality row.  This turns Cubit 2026.6-style quality-metric learning and
high-order mesh literature into a reusable artifact contract instead of a loose
table of nice-looking minimum values.

Slot376 extends that ledger to execution evidence.  When a quality ledger is
promoted beyond a raw Cubit replay, require `created_at_utc`, Cubit/Coreform
`version`, nonnegative `elapsed_s`, and a compact `timing_breakdown_s` with the
dominant stages.  Reject stale versions, non-parseable timestamps, sparse timing
breakdowns, and timing totals that cannot fit inside the recorded elapsed time.

Slot390 extends the same quality ledger to parametric meshing and optimization
evidence.  When a Cubit quality row is produced from a mesh-size, smoothing,
scheme, or design-variable parameter set, carry `parameter_set_artifact_id`,
`parameter_set_digest`, `parameter_set_path`, `objective_observable_id`, and
`objective_observable_family` with the ledger and replay
`cubit_mesh_quality_ledger_identity_gate(..., require_parameter_set_artifact=True)`.
Reject stale parameter-set digests, missing parameter-set paths, and wrong
objective families before a quality row is reused by notebooks, panels, or
solver-ready mixed-mesh imports.

Slot425 adds postprocess-row convention schema identity to the quality ledger.
`quality_metric_set_id` says which metric family was computed, but it does not
say how a distribution was selected, aggregated, reduced to a minimum, or handed
to a notebook objective.  Carry
`mesh_quality_postprocess_row_convention_schema_id` and replay
`cubit_mesh_quality_ledger_identity_gate(...,
require_quality_postprocess_row_convention_schema=True)` before using quality
rows as panel defaults, solver-ready evidence, or optimization objectives.
Reject stale scalar-row conventions and missing postprocess-row convention
schemas even when `quality_metric_set_id`, mesh digest, and parameter/objective
identity all look plausible.

Slot432 adds component-basis schema identity to that same quality ledger.
`quality_metric_set_id` names the metric family and
`mesh_quality_postprocess_row_convention_schema_id` names the row reduction, but
neither says which element family, quality component, coordinate basis, or
normalization the row represents.  Carry
`mesh_quality_component_basis_schema_id` and replay
`cubit_mesh_quality_ledger_identity_gate(...,
require_quality_component_basis_schema=True)` before a Cubit quality row becomes
a notebook default, panel slider target, or optimization objective.  Reject
stale scalar-value component bases and missing component-basis schemas even
when metric-set, postprocess-row convention, mesh digest, and parameter/objective
identity all pass.

Slot330 binds that headless/process evidence to the mixed solver-ready package.
When `cubit_mixed_solver_ready_package_gate` is used for a hex+pyramid+tet
handoff, pass the verified `cubit_headless_batch_quality_package_gate` as
`headless_batch_quality_gate` whenever process evidence is available.  The
package then requires the same export/geometry identity plus
`process_mode=headless_batch`, `-nographics -batch`, disabled GUI daemon, batch
script identity, and a successful or explicitly documented process exit before
the mixed mesh is promoted.  This prevents a good mixed topology package from
being accidentally paired with a GUI run, stale batch script, or undocumented
nonzero exit.

Slot354 tightens installed-version evidence.  When replaying
`cubit_headless_installation_route_gate`, record the actual console binary path
as `coreform_cubit.com`, not the GUI `coreform_cubit.exe` stub, and require the
`.com -version` probe command to use that same recorded binary.  Release-note
knowledge such as Coreform Cubit 2026.6 stays a watchlist until that installed
binary reports the same version.  This keeps live headless evidence separate
from documentation-only release learning and from GUI daemon paths.

Slot346 adds the downstream solver-route manifest.  A mixed Cubit
hex+pyramid+tet package should say how each element family will be consumed:
hex cells are the primary volume-FEM region, pyramid cells are explicit
transition bridge cells, and tet cells are compatibility or subregion cells.
Record `solver_route_package_id`, `route_policy`, `downstream_solver`,
`tet_only_owner=netgen_tri_tet_path`, `no_implicit_tetization=true`,
`volume_routes`, and `surface_routes`, then replay the contract with
`cubit_mixed_solver_route_manifest_gate` and pass it into
`cubit_mixed_solver_ready_package_gate`.  Pyramid cells are not display-only
mesh noise and should not be silently split or tetized before a solver-ready
claim.

Slot383 adds downstream solver-reader contract identity to that same route
manifest.  When a mixed hex+pyramid+tet route is promoted as solver-ready, the
manifest should also carry `solver_contract_artifact_id`,
`solver_contract_digest`, and `solver_contract_path` (or the
`downstream_solver_contract_*` aliases) for the NGSolve/radia-ngsolve reader
contract that actually accepts those element families.  A fresh Cubit route
manifest with a stale solver-reader digest, or no contract path, remains
`needs_attention`; do not infer solver readiness from element counts alone.

Slot418 adds solver-route convention schema identity to the same manifest.
`route_policy` and the element-route rows say what the package claims, while
`solver_route_convention_schema_id` records the versioned meaning of hex
primary, pyramid transition, tet compatibility/subregion, surface trace, and
no-implicit-tetization roles.  Pass
`expected_solver_route_convention_schema_id` plus
`require_solver_route_convention_schema=True` to
`cubit_mixed_solver_route_manifest_gate`.  A value-only route convention or
missing schema id should fail even when element counts, route policy, and the
downstream solver-reader contract digest look plausible.

Slot404 binds the emitted Netgen `.vol` file to its Cubit `.vol.json` sidecar.
The export package already carries export id, geometry id, order, observable id,
and output digest, but a stale sidecar can still have the right filename and
wrong counts.  When sidecar count metadata is available, pass
`require_vol_sidecar_inventory_counts=True` to
`cubit_export_package_identity_gate` and require the sidecar `n_elements`,
`n_points`, and `order` to match the parsed `.vol` inventory.  For the
`01_Tet_Hex_Pyramid_order1.vol` fixture this means 12 volume elements, 13
points, and order 1.  A stale sidecar element count should fail before a mixed
mesh is promoted to a panel, notebook, or solver-ready route.  Also do not infer
order from the filename: in a headless Cubit 2025.12 smoke run, omitting the
`order` argument while writing an `*_order1.vol` filename produced a sidecar
with `order=2`; adding `order 1` made the sidecar record `order=1`.

Slot411 adds schema identity for the same sidecar.  Counts and order protect
against stale data, but an old `.vol.json` layout can still carry plausible
`n_elements`, `n_points`, and `order` fields while meaning "material volume
table" instead of "Netgen `.vol` inventory sidecar".  Store
`vol_sidecar_schema_id` on the `.vol.json` artifact row and pass
`expected_vol_sidecar_schema_id` plus `require_vol_sidecar_schema=True` to
`cubit_export_package_identity_gate` before notebooks or solver-ready routes
consume the package.  A stale legacy sidecar schema or a missing schema id
should fail even when the sidecar filename and counts still match the parsed
`.vol`.

### Continuous-loop Coreform Cubit 2026.6 release routing gate

Coreform Cubit 2026.6 was released on 1 June 2026.  The official release notes
list anisotropic tetrahedral meshing, cohesive element generation, higher-order
Tetra10/Tri6 Jacobian and scaled-Jacobian metrics, improved triangle/tet
robustness, lower-memory Sculpt refinement, expanded solver/file compatibility
including 64-bit Exodus IDs, namespaces for tracking names through operations,
GNN feature extraction, and an included Python 3.12 runtime:
`https://coreform.com/coreform-cubit/release-notes/v2026-6/`.

Use `cubit_release_feature_routing_gate` to turn that public release knowledge
into lab routing before treating it as validation evidence.  The CAE-AI Lab
policy remains: Cubit is the hex-led and mixed hex+pyramid+tet lane; tet-only
education stays on Netgen/OCC unless a slot explicitly needs Coreform's advanced
tet controls.  For 2026.6, map higher-order Tetra10/Tri6 metrics into archived
quality replay, cohesive elements into explicit interface/block identity
examples, anisotropic tet meshing into an advanced reference lane, and solver
I/O improvements into format-compatibility notes rather than `.vol` parser
relaxation.

For loop learning, close a Coreform/Cubit slot only when the headless batch
command, `result_artifact_id`, `result_output_schema_id`, top
`timing_breakdown_s`, mesh-quality evidence, and `.vol` or Exodus export
inventory are recorded together.  A GUI-only screenshot or a mesh file without
sidecar identity is a candidate, not a verified MCP lesson.

### Continuous-loop third-party curvilinear handoff manifest gate

The new literature on high-order curvilinear mesh generation from third-party
meshes reinforces a practical Cubit rule: high-order export is not solver-ready
until the imported mesh, CAD/geometry association, curved export order, routing
hint, and quality metrics travel together.  Use
`cubit_curvilinear_handoff_manifest_gate` for this package.  A healthy package
records a `third_party_mesh` or `imported_mesh` source, hex or mixed volume
kinds, preserved boundary ids, a CAD projection/association policy, Netgen
`.vol` or other curved export order at least 2, no implicit tetization or
element splitting, a bounded CAD projection error such as
`projection_quality.max_distance <= projection_quality.tolerance`, a
scaled-Jacobian/Jacobian quality minimum, and `negative_jacobian_count = 0`.

Do not relax the first-order tri/tet `.vol` parser for this.  Tet-only
education remains on Netgen/OCC; Cubit owns the hex-led or mixed curvilinear
handoff lane.  The manifest is the bridge between a mesh-recovery paper idea
and a readable CAE-AI Lab workflow: it says what geometry the high-order nodes
are curved to, what topology was preserved, how far the projected boundary
nodes are from the intended CAD entity, and which quality evidence allows the
export to move toward NGSolve or another solver.

### Continuous-loop mapped hex quality replay gate

For an all-hex quality baseline, use a mapped brick such as
`create brick x 3 y 2 z 1`, `volume 1 scheme map`, `volume 1 size 0.25`,
`mesh volume 1`, `block 1 add hex all`, then export Netgen `.vol` orders 1
and 3.  The verified 2026-06-29 live slot produced 384 hexes, no
tet/wedge/pyramid elements, scaled-Jacobian min `0.9999999999999999`, CAD
volume 6, external area 22, and NGSolve `.vol` volume/BND relative errors
`1.78e-14` and `1.52e-14`.

Use `cubit_hex_quality_gate` on the archived Cubit quality list before trusting
a hex-led export.  Inventory and NGSolve integration still matter: the stock
Netgen export may write `surfaceelementsuv`, and order-3 adds `curvedelements`
without changing the first-order topology.

For slot-level learning, also run the quality distribution replay gate on the
same archived list.  `cubit_quality_distribution_gate` records p05/p50/p95
style quantiles and a compact histogram, so the MCP server can distinguish a
globally weak mesh from a mostly good mesh with a small damaged tail.  This is
the preferred evidence format before promoting Cubit hex meshes into heavier
validation or solver-ready examples.

### Continuous-loop .vol label metadata gate

Before solver setup, confirm that the Netgen `.vol` export carries material
labels and boundary names, not just element topology.  The replay helper
`cubit_vol_label_metadata_gate` checks the parsed `materials` and `bcnames`
sections, required material names, required boundary names, and the presence of
volume/surface elements.  Use it before assigning sources, boundary conditions,
or notebook panels from a Cubit/Coreform export; a volume/area/quality pass is
not solver-ready if labels were lost.

### Continuous-loop mixed hex+pyramid+tet order-series gate

For hex-led mixed workflows, replay the small Coreform/Cubit fixture exported
as Netgen `.vol` orders 1 through 5.  Use
`cubit_mixed_order_series_inventory_gate` to keep the expected topology
invariant across orders: 1 hex, 1 pyramid, 10 tets, 6 quad boundary faces, and
10 triangle boundary faces.  The files grow with order because
`curvedelements` data is added, but the routing class remains
`cubit_hex_or_mixed_path`.  This is why `cubit_vol_inventory` routes from
`surfaceelements` and `volumeelements`, not from file size or sidecar material
volumes.  A pyramid transition block can report zero material volume in the
`.vol.json` sidecar while still being present as a real topology record in the
`.vol`.

### Continuous-loop mixed transition metadata gate

For hex-to-tet handoff, do not treat the pyramid bridge as a bookkeeping
curiosity.  It is the explicit transition topology between the hex-led Cubit
lane and the tet region.  Replay the `.vol` inventory with
`cubit_mixed_transition_metadata_gate`: require hex, pyramid, and tet volume
kinds, require quad and triangle surface families, require the routing hint
`cubit_hex_or_mixed_path`, and require a pyramid transition block label such
as `pyramid_transition` or `pyram`.

This is intentionally checked from `.vol` arity and labels rather than from the
companion `.vol.json` material-volume table.  The sidecar can report zero
volume for the pyramid transition block while the `.vol` still contains the
pyramid element that makes the mixed mesh conformal.  Surface arity matters
too: a hex-led mixed handoff should expose quad faces as well as triangle
faces before boundary labels or NGSolve BND rows are trusted.

Slot282 adds the interface-adjacency ledger.  After the mixed-transition gate
passes, replay `cubit_mixed_interface_adjacency_gate` with rows such as
`hex_to_transition` on a quad face touching `hex_core` and
`pyramid_transition`, plus `transition_to_tet` on triangle faces touching
`pyramid_transition` and `tet_region`.  This catches a stale interface ledger:
the `.vol` can still contain the right hex/pyramid/tet and quad/triangle
counts while the solver-ready boundary row has swapped or forgotten which
surface is the hex-pyramid interface and which surface is the pyramid-tet
interface.

### Continuous-loop live mixed hex+pyramid+tet NGSolve BND gate

The stock Cubit `export netgen` path can write boundary records as
`surfaceelementsuv` rather than plain `surfaceelements`.  Treat
`surfaceelementsuv` as the same surface inventory for routing; otherwise a
valid live export looks like it has no boundary faces.  A verified 2026-06-29
headless batch used the split `2 x 1 x 1` brick fixture (`1` mapped hex on one
side, `1` pyramid transition, `10` tets on the other side), scaled by `0.001`,
and exported orders 1 and 3.  NGSolve loaded both mixed `.vol` files directly:
volume matched `2e-9`, and `Integrate(1, mesh, BND)` matched `11e-6`, not the
external brick area `10e-6`, because the split material interface of area
`1e-6` is included once.  This is the same material-interface rule as the
multi-block hex gates: compare BND to external area plus material-interface
area, and do not use an empty surface inventory as evidence that BND integration
will fail.

### Continuous-loop mixed solver-ready package gate

After the mixed transition, export identity, BND-area, and quality gates pass
individually, bind them with `cubit_mixed_solver_ready_package_gate` before a
notebook or solver-ready validation consumes the mesh.  A healthy mixed Cubit
package records:

* `.vol` inventory with `hex`, `pyramid`, and `tet` volume kinds;
* `cubit_mixed_transition_metadata_gate` status `ok`;
* `cubit_export_package_identity_gate` status `ok` and routing hint
  `cubit_hex_or_mixed_path`;
* `cubit_ngsolve_bnd_area_includes_material_interfaces_once` status `ok`;
* `cubit_quality_distribution_gate` status `ok` with positive quality count.
* optional `cubit_mixed_interface_adjacency_gate` status `ok` to freeze the
  hex-pyramid and pyramid-tet interface roles.
* optional `cubit_curvilinear_handoff_manifest_gate` status `ok` to bind
  third-party/CAD association, boundary id preservation, projection tolerance,
  curved export order, no implicit element conversion, and zero
  negative-Jacobian evidence to the same mixed package.
* optional export output artifact identity: `export_output_artifact_id`,
  `export_output_digest`, and `export_output_path` bind the actual emitted
  `.vol`/sidecar package consumed by notebooks or solver-ready steps.

This is the Cubit counterpart to build123d CAD handoff gates: Cubit owns the
mixed mesh evidence, while tet-only `.vol` education remains on the Netgen/OCC
route.

Slot338 adds the literature-driven curvilinear handoff row to this package
gate.  The high-order third-party mesh lesson is discovery followed by strict
validation: an imported mesh must record CAD/source association and projection
quality before the mixed solver-ready package can reuse its curved `.vol`.

Slot362 binds Cubit scheme traces to the actual exported mesh artifact.  When
`cubit_meshing_scheme_trace_gate` is promoted to solver-ready evidence, record
`export_output_artifact_id`, `export_output_digest`, and `export_output_path`
beside the command digest, volume schemes, and export order.  A fresh journal
trace with an old `.vol` digest is a stale mesh package and should fail before
NGSolve/radia-ngsolve consumes it.

### Continuous-loop submodel boundary handoff mesh package gate

When a Cubit hex-led or mixed `.vol` is used as a local/zoomed submodel, bind
the mesh inventory to the parent-to-local boundary handoff metadata with
`cubit_submodel_boundary_handoff_mesh_package_gate`.  A healthy package records
the `.vol` source, `volume_kind_counts` with at least the expected `hex` family,
`cubit_hex_or_mixed_path` routing, boundary labels from `bcnames`,
`parent_model_id`, `submodel_region_id`, `zoom_boundary_id`,
`boundary_transfer_quantity`, `boundary_transfer_error_estimate`,
`boundary_transfer_error_unit`, `local_refinement_rule`,
`transition_policy` when a pyramid bridge is present, and
`target_observable_id`.  This is the mesh-side companion to the generic
submodel gate: local hex refinement is not enough unless the inherited boundary
condition, mixed-mesh family inventory, transition policy, and error budget are
attached to the same artifact.

Slot202 adds the routing-policy row to this package gate.  When
`cubit_release_feature_routing_gate` has encoded the lab rule, pass it into
`cubit_mixed_solver_ready_package_gate(..., routing_policy_gate=...)`.  The
package is then rejected if Cubit is mislabeled as the tet-only default route or
if Netgen/OCC is not recorded as the ordinary tet-only owner.

### Continuous-loop two-block hex interface gate

For material-interface bookkeeping, create two adjacent mapped bricks
(`1 x 2 x 3` and `2 x 2 x 3`), `imprint`/`merge`, assign one block per
volume, and export Netgen `.vol` orders 1 and 2.  The expected inventory is
144 hexes, no tet/wedge/pyramid elements, CAD volume 18, external area 42, and
one shared material interface of area 6.  NGSolve `Integrate(1, mesh, BND)` on
this `.vol` returns 48, not 42: it includes the material interface once.
Therefore do not compare `BND` area blindly with the external CAD area when
multi-material blocks share an internal face.  Use Cubit's
`get_relatives("surface", sid, "volume")` or equivalent inventory to separate
external boundary surfaces from material interfaces before FEM/BEM coupling.

The 2026-06-29 slot90 replay keeps the same rule but records the Cubit version
explicitly before making any feature claim: the installed headless executable was
Coreform Cubit `2025.12`, while public Coreform 2026.6 release notes list newer
features such as anisotropic tetrahedral meshing, cohesive element generation,
higher-order Tetra10/Tri6 Jacobian metrics, improved triangle/tet robustness,
and expanded solver/file compatibility.  For the CAE-AI Lab loop this means:
keep Cubit as the hex-led and mixed hex+pyramid+tet lane; route ordinary tet-only
education meshes to Netgen/OCC; use 2026.6 tet/cohesive features only when that
version is actually installed and the artifact records the version.

Slot98 turns the 2026.6 higher-order quality-metric lesson into a replayable MCP
gate: archive the raw Tetra10/Tri6 metric list and pass it through
`cubit_element_quality_gate(element_type="Tetra10"|"Tri6", metric="scaled
Jacobian"|"Jacobian", min_value=...)` before routing the mesh downstream.  This
does not claim that 2026.6 is installed on INTEL11; it records that higher-order
tet/tri Jacobian metrics are a lower-bound quality contract, while live hex-led
Coreform work on this machine targets the installed 2025.12 headless executable.

Slot226 rechecked the installed-version lane directly with
`coreform_cubit.com -version`.  The synchronous console probe reported
`status: ValidStudent` and `Coreform Cubit Version 2025.12 Build 3d8d3af7`.
Record this as installation/session evidence only: it proves the local headless
route can report a valid license and installed version, but it is not a mesh
quality result.  The public 2026.6 release-note features remain a watchlist
until a 2026.6 executable is installed and replayed.

Slot90 used two mapped blocks with dimensions `1.25 x 1.5 x 0.75` and
`1.75 x 1.5 x 0.75`, exported Netgen orders 1 and 3.  The live inventory was
216 hexes, 270 quad surface records, no tet/wedge/pyramid records, CAD volume
3.375, external area 15.75, and one material interface of area 1.125.  NGSolve
loaded both orders and integrated volume `3.374999999999978` and BND area
`16.87499999999997`, confirming that BND equals external area plus the shared
interface once.  The batch-script lesson was also repeated: avoid top-level
multi-line dict literals in `coreform_cubit.com -nographics -batch` scripts;
write result dictionaries with one assignment per physical line.

### Continuous-loop three-block hex quality gate

For a slightly heavier hex-led bookkeeping gate, create three adjacent mapped
bricks with widths `1.0, 1.5, 2.0`, common `2 x 3` cross-section, `imprint`/
`merge`, assign one block per volume, and export Netgen `.vol` orders 1 and 3.
The expected inventory is 216 hexes, no tet/wedge/pyramid elements, CAD volume
27, external area 57, and two material interfaces whose total area is 12.
Record Cubit's scaled-Jacobian quality in the raw artifact; the live 2026-06-29
gate had min scaled Jacobian `0.9999999999999999` over all 216 hexes.  NGSolve
loads both `.vol` files and integrates volume 27 while `Integrate(1, mesh, BND)`
returns 69, again confirming that BND includes material interfaces once.
Replay this convention with
`cubit_bnd_area_interface_gate(external_area=57, material_interface_area=12,
ngsolve_bnd_area=69)` so the expected comparison is executable and not just a
prose warning.

### Continuous-loop curved hex cylinder order-series gate

For curved hex validation, a cylinder should use Cubit's hex-led auto/sweep
route rather than forcing every webcut quarter volume to `scheme map`.  The
failed map route can produce "Trouble finding logical box" and zero exported
mesh; `volume all scheme auto` selected a sweep-compatible all-hex mesh for the
live slot.

The verified live gate used `create cylinder height 1.25 radius 0.5`, webcut
with x/y planes, `volume all scheme auto`, `volume all size 0.18`, and exported
Netgen `.vol` orders 1, 3, and 5.  The inventory stayed 224 hexes with no tet,
wedge, or pyramid elements.  NGSolve direct `.vol` integration showed the
curved-geometry convergence clearly:

| order | volume rel err | BND area rel err |
|---:|---:|---:|
| 1 | 0.0255 | 0.00816 |
| 3 | 4.08e-5 | 1.88e-5 |
| 5 | 1.57e-7 | 6.99e-8 |

For multi-volume cylinder webcuts, compare NGSolve `BND` area with external
CAD area plus material-interface area, not external area alone.  The planar
cut interfaces are part of NGSolve's boundary integration set.

### Continuous-loop annular hex tube capacitance field gate

For Cubit's main hex-led route, validate one actual field solve, not only CAD
volume.  A robust live fixture is an annular tube: outer cylinder radius `1.0`,
inner cylinder radius `0.6`, length `1.5`, subtract inner from outer, webcut by
the x/y planes, `volume all scheme sweep`, and export Netgen `.vol` orders 1
and 3.  The heavier order-5 closure repeats the same fixture at orders 1, 3,
and 5.  The sweep mesh should be all hex; the verified slot had 864 hexes, no
tet/wedge/pyramid elements, and min scaled Jacobian `0.9951847266721953`.

In NGSolve, load the `.vol` directly.  The exact checks are:

* volume `pi*(b^2-a^2)*L`,
* BND area = external CAD area + planar webcut material-interface area,
* capacitance per length `2*pi/log(b/a)` for the annular Laplace problem.

The order-3 live result had volume rel err `2.63e-6`, BND area rel err
`1.43e-6`, and capacitance-per-length rel err `1.00e-6`.  This is the
important upgrade over a geometry-only gate: curved hex export preserves a field quantity through the `.vol` path.

The order-5 closure tightened the geometry gate to volume rel err `4.11e-9`
and BND-area rel err `2.27e-9`, while the capacitance-per-length rel err stayed
at `1.02e-6`.  For webcut annuli, do not select electrodes only by
`normal[2] == 0`; that also catches the planar material-interface cuts.  Use a
radial-normal check such as `abs(n dot rhat) > 0.9` before classifying an outer
or inner cylindrical electrode.  If Cubit warns that the C++ plugin is out of
date, record it in the slot artifact and refresh the installed plugin before
using the result as production evidence.

## Choose Your Workflow

1. **Is your geometry planar (no curved surfaces)?**
   -> Use any export format. Curving is not needed.
   -> Simplest: `export netgen "mesh.vol" order 1` + `Mesh("mesh.vol")`

2. **Do you only need 2nd order (not 3rd+)?**
   -> Use `export netgen`:
   ```
   Cubit -> mesh -> export netgen "mesh.vol" order 2 -> Mesh("mesh.vol")
   ```

3. **Do you need 3rd order or higher?**
   -> Use `export netgen` with higher order:
   ```
   Cubit -> mesh -> export netgen "mesh.vol" order 3 -> Mesh("mesh.vol")
   ```
   Works for ANY geometry shape — cylinder, sphere, torus, cone,
   Boolean operations, freeform surfaces, etc. Supports order 1-5.

## Accuracy: p-Convergence Results (Verified 2026-04-02)

All shapes tested with ACIS CallbackGeometry + edge snapping:

| Shape | Surfaces | Curves | p=2 V err | p=3 V err | p=5 V err | p=5 A err |
|-------|----------|--------|-----------|-----------|-----------|-----------|
| Sphere | 1 | 0 | -0.023% | +0.002% | -0.000003% | -0.000003% |
| Cylinder | 3 | 2 | -0.003% | +0.001% | -0.000001% | -0.000001% |
| Frustum (cone) | 3 | 2 | -0.009% | +0.002% | -0.000006% | -0.000009% |
| Torus | 1 | 0 | -0.026% | +0.004% | -0.000018% | -0.000011% |
| Box with hole | 7 | 14 | +0.004% | -0.001% | +0.000003% | +0.000003% |

Key: p=5 achieves 10^-5 to 10^-6 % error for ALL shapes, matching OCC native accuracy.
Runnable OCC-native companion (no Cubit dependency): `docs/ngsolve_user_meeting/volume_area_convergence.ipynb`,
with durable numeric results in `docs/ngsolve_user_meeting/volume_area_convergence_results.json`.

| Method | p=2 Error | p=5 Error | Max Order | Complexity |
|--------|-----------|-----------|-----------|------------|
| export netgen (ACIS) | ~0.003-0.03% | ~1e-6% | 5 (tet/hex/wedge) | Low |
| OCC mesh.Curve() | ~0.003-0.03% | ~1e-6% | 5+ (tet) | Low |
| 1st order (no curving) | ~1.4% | N/A | 1 | None |

## Key Principle

`export netgen` uses Cubit's ACIS kernel for surface projection via
CallbackGeometry (compact_netgen C++ static link). The mesh curving is done
entirely in the C++ plugin — no STEP files, no OCC geometry, no Python
dependency. The .vol file contains the curved mesh ready for NGSolve.
"""

WORKFLOW_EXPORT_CURVED = """
# export netgen APREPRO Command Reference

## Signature

```python
import tempfile
from ngsolve import Mesh
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `order` | int (1-5) | Polynomial order for mesh curving (1=linear, 2=quadratic, ..., 5=quintic) |
| `overwrite` | flag | Overwrite existing file |

## Returns

`.vol` file with curvedelements section + companion `.vol.json` with CAD reference values.

## How It Works Internally

1. Reads mesh topology (nodes, elements, surface tris/quads) from Cubit
2. Extracts **1D segment elements** on geometry curves (edges between surfaces)
   - Each segment has surfnr1/surfnr2 (0-based FD indices of adjacent surfaces)
   - Arc-length normalized dist parameter for each endpoint
3. Creates a `netgen.meshing.Mesh` with volume + surface + segment elements
4. Uses **CallbackGeometry** with ACIS callbacks (C++ static link):
   - `project_func`: Project point onto Cubit surface (ACIS closest_point_trimmed)
   - `normal_func`: Surface normal at a point (ACIS normal_at)
   - `edge_project_func`: Project point onto Cubit curve (ACIS curve.closest_point_trimmed)
5. Calls `BuildCurvedElements(order)`:
   - Surface nodes: projected via `PointBetween` (surface callback)
   - Edge nodes (on curves between surfaces): projected via `PointBetweenEdge` (curve callback)
   - Requires segments with correct surfnr1/surfnr2 for `use_edge` flag
6. Writes the curved mesh to `.vol` file with curvedelements section

### Critical Implementation Details (Edge Snapping)

- **1D segments are required**: Without them, `BuildCurvedElements` does not set
  `use_edge[edgenr]=1` and edge curving is skipped entirely.
- **surfnr1/surfnr2 must be 0-based** FD indices (not -1). BuildCurvedElements passes
  these directly to `PointBetweenEdge` as surfi1/surfi2.
- **CallbackGeometry receives 0-based surfnr**: The `PointBetweenEdge` implementation
  must check `surfi >= 0` (not `> 0`) and convert to 1-based before calling
  `edge_project_func` (which uses 1-based Cubit surface/curve indices).
- **Without edge snapping**: Sphere (1 surface, no edges) works perfectly, but
  cylinder (3 surfaces, 2 edges) area does not p-converge (stuck at -0.4%).
- **With edge snapping**: Both sphere and cylinder converge to 1e-6% at p=5.

## Why ACIS (Not OCC)?

- **No seam lines**: ACIS represents a cylinder as 1 surface; OCC splits
  it at the seam into 2 faces. ACIS has no parametric discontinuities.
- **No STEP needed**: Direct access to Cubit's geometry kernel, no file exchange.
- **Universal**: Works for ANY surface type — analytic, BSpline, freeform,
  Boolean results. Not limited to cylinder/sphere/torus/cone.
- **Exact**: ACIS surface projection is exact for the underlying geometry.

## Example: Basic Usage

```python
import tempfile
from ngsolve import Mesh

# Create and mesh geometry in Cubit
cubit.cmd("create cylinder height 2 radius 0.5")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Export with curving
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Verify
from ngsolve import Integrate, CF
import math
R, H = 0.5, 2.0
expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Example: Any Geometry Shape

```python
import tempfile
from ngsolve import Mesh

# Complex Boolean geometry — no special handling needed
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("create cylinder height 4 radius 0.3")
cubit.cmd("subtract volume 2 from volume 1")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Works for any geometry — cylinder holes, fillets, chamfers, BSplines...
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```

## Requirements

- Coreform Cubit 2025.12+ with Radia plugin installed (`cubit-plugin-install`)
- NGSolve 6.2.2603+ (curvedelements Load, hex/prism curving)
"""

WORKFLOW_CYLINDER = """
# Example: Cylinder

## Step-by-Step

```python
import tempfile
from ngsolve import Mesh, Integrate, CF
import math

R = 0.5   # Radius
H = 2.0   # Height

# Step 1: Create geometry in Cubit
cubit.cmd(f"create cylinder height {H} radius {R}")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "domain"')

# Step 3: Export with curving (no STEP, no OCC, no SetGeomInfo!)
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Step 4: Verify
expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Note

No STEP export/reimport needed. No OCCGeometry. No set_cylinder_geominfo().
`export netgen` handles everything via Cubit's ACIS kernel.
"""

WORKFLOW_SPHERE = """
# Example: Sphere

## Step-by-Step

```python
import tempfile
from ngsolve import Mesh, Integrate, CF
import math

R = 0.5

# Step 1: Create geometry in Cubit
cubit.cmd(f"create sphere radius {R}")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Step 3: Export with curving
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Step 4: Verify
expected_vol = 4/3 * math.pi * R**3
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```
"""

WORKFLOW_TORUS = """
# Example: Torus

## Step-by-Step

```python
import tempfile
from ngsolve import Mesh, Integrate, CF
import math

R_MAJOR = 1.0
R_MINOR = 0.3

# Step 1: Create geometry in Cubit
cubit.cmd(f"create torus major {R_MAJOR} minor {R_MINOR}")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.08")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Step 3: Export with curving
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Step 4: Verify
expected_vol = 2 * math.pi**2 * R_MAJOR * R_MINOR**2
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Half-Torus Coil (IH Sample)

For induction heating, a half-torus coil with source/sink terminals:

```python
# Native torus + webcut = clean curving (1 toroidal surface)
cubit.cmd("create torus major radius 0.11 minor radius 0.01")
cubit.cmd("webcut volume 1 with plane xplane noimprint nomerge")
cubit.cmd("delete volume 2")  # Keep half-torus

cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size auto factor 5")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "coil"')
cubit.cmd('sideset 1 add surface 2')
cubit.cmd('sideset 1 name "source"')
cubit.cmd('sideset 2 add surface 3')
cubit.cmd('sideset 2 name "sink"')
```

NOTE: Do NOT use `sweep` to create torus geometry for high-order meshing.
`create torus` (native ACIS) produces a single toroidal surface that
curves correctly. `sweep` splits the surface at z=0, which can cause
cross-projection issues in ACIS closest_point_trimmed at high order.
Use `create torus` + `webcut` instead.
"""

WORKFLOW_COMPLEX = """
# Complex Geometry (Boolean Operations)

With `export netgen`, complex geometries require NO special workflow.
Boolean operations, multiple curved surfaces, freeform surfaces — all
are handled automatically by the ACIS kernel.

## Example: Brick with Cylindrical Hole

```python
import tempfile
from ngsolve import Mesh, Integrate, CF
import math

BRICK_SIZE = 2.0
R_HOLE = 0.3

# Step 1: Create geometry with Boolean operations in Cubit
cubit.cmd(f"create brick x {BRICK_SIZE} y {BRICK_SIZE} z {BRICK_SIZE}")
cubit.cmd(f"create cylinder height {BRICK_SIZE*2} radius {R_HOLE}")
cubit.cmd("subtract volume 2 from volume 1")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Step 3: Export with curving — works for any geometry
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Step 4: Verify
expected_vol = BRICK_SIZE**3 - math.pi * R_HOLE**2 * BRICK_SIZE
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Why No Special Workflow?

The old workflow required:
- Creating geometry in OCC (not Cubit)
- Calling name_occ_faces() for face name mapping
- STEP export from OCC
- STEP reimport in Cubit with 'noheal'
- export_netgen_with_names() for name-based mapping
- set_*_geominfo() for each curved surface type

With `export netgen`, ALL of this is eliminated. The ACIS kernel handles
surface projection for any geometry directly, regardless of complexity.
"""

WORKFLOW_GMSH_2ND_ORDER = """
# Alternative: Netgen .vol 2nd Order Workflow (Simplest)

If you only need 2nd order elements and don't need 3rd order or higher,
the APREPRO export netgen command provides a simple workflow.

## Step-by-Step

```python
from ngsolve import Mesh, Integrate, CF

# Step 1: Create geometry and mesh in Cubit
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.2")
cubit.cmd("mesh volume 1")

# Step 2: Register blocks
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "sphere"')

# Step 3: Export to Netgen .vol with order 2
cubit.cmd('export netgen "mesh.vol" order 2 overwrite')

# Step 4: Read into NGSolve
mesh = Mesh("mesh.vol")
# Done! No geometry reference needed at compute time.
```

## Advantages
- No geometry reference needed at compute time
- Very simple workflow (APREPRO command)
- Good accuracy (~0.003%)

## Limitations
- Supports order 1-5 (arbitrary order via ACIS CallbackGeometry)
"""

WORKFLOW_ACCURACY = """
# Accuracy Guide: Choosing the Right Order

## Volume Error by Method and Order

| Method | Order 1 | Order 2 | Order 3 | Order 4 | Order 5 |
|--------|---------|---------|---------|---------|---------|
| No curving | ~1.4% | - | - | - | - |
| export netgen | ~1.4% | ~0.003% | ~0.0004% | ~0.00005% | ~0.000006% |

## When Higher Order Matters

- **Structural/thermal FEM**: Order 2 usually sufficient
- **Electromagnetics (curl-curl)**: Order 2-3 recommended
- **BEM inductance extraction**: Order 3 recommended (surface accuracy critical)
- **Acoustic/wave propagation**: Order 3-5 for dispersion control
- **Geometry verification only**: Order 2 is fine

## Mesh Size vs Order Trade-off

For a target accuracy, you can either:
- **h-refinement**: More elements, keep order low
- **p-refinement**: Fewer elements, increase order

High-order curving (order 3+) is most beneficial when:
- Geometry has high curvature
- Coarse meshes are needed (computational cost)
- High accuracy is required on curved boundaries

## Verification Pattern

Always verify accuracy after curving:

```python
import math
from ngsolve import Integrate, CF, BND

# Volume check
expected_vol = math.pi * R**2 * H  # Exact volume
computed_vol = Integrate(CF(1), mesh)
vol_error = abs(computed_vol - expected_vol) / expected_vol * 100

# Surface area check (optional)
expected_area = 2 * math.pi * R * H + 2 * math.pi * R**2
computed_area = Integrate(CF(1), mesh, VOL_or_BND=BND)
area_error = abs(computed_area - expected_area) / expected_area * 100

print(f"Volume error: {vol_error:.4f}%")
print(f"Area error: {area_error:.4f}%")
```
"""

WORKFLOW_FIELD_SOLVE = """
# High-Order Hex in an Actual FIELD Solve (not just volume)

Curving improves not only the VOLUME integral but a solved FIELD quantity. Verified on the
coaxial-capacitor Laplace problem: a hex-meshed annular TUBE (inner radius a, outer b,
length L), exported `order N`, carried into an NGSolve Laplace solve, gives the per-length
capacitance against the exact closed form

    C/L = 2 pi eps / ln(b/a)        (eps = 1 here).

Because the exact solution is purely radial V(r)=ln(b/r)/ln(b/a), the flat ends are exact
Neumann (dV/dz = 0) -- NO end effects -- so the only error is the hex geometry's chord error
on the CURVED cylindrical electrodes, which curving removes.

## Result (a=1, b=2, L=2; 707 hexes; FIELD order fixed, GEOMETRY order varied)

| geom order | volume err | C/L err  |
|------------|------------|----------|
| 1 (faceted)| +0.049%    | +0.854%  |
| 2 (curved) | +0.001%    | +0.002%  |
| 3 (curved) | +0.000%    | +0.000%  |

Faceted hex mis-states the capacitance by ~0.85% (chord error on the round electrodes);
curved hex (order 2-3) nails it. Curving matters for FIELDS, not just geometry checks.

## Boundary conditions without exported sideset names

`export netgen` may not carry your Cubit sidesets as named NGSolve boundaries. Robust trick:
select boundaries by GEOMETRY in NGSolve and impose Dirichlet by a penalty. For the tube,
pick the cylindrical electrodes by their RADIAL normal (|n_z| < 0.5) -- this excludes the
flat ends (|n_z| ~ 1), which must stay Neumann. Do NOT select the electrodes by radius alone:
the flat-end corner rings also have r ~ a or b and pinning them to V=1/0 over a finite width
wrongly distorts the field (gave a stuck +7.5% error until switched to the normal selector).

```python
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction, grad, dx, ds,
                     Integrate, IfPos, sqrt, x, y, specialcf)
mesh = Mesh("tube_o3.vol")                 # high-order hex from Cubit
r = sqrt(x*x + y*y); n = specialcf.normal(3)
cyl    = IfPos(0.25 - n[2]*n[2], 1.0, 0.0) # |n_z|<0.5 -> cylindrical face only
on_in  = cyl*IfPos(0.5*(a+b) - r, 1.0, 0.0)
mask   = cyl                               # penalise both electrodes
fes = H1(mesh, order=3); u, v = fes.TnT(); alpha = 1e7
A = BilinearForm(fes, symmetric=True); A += grad(u)*grad(v)*dx + alpha*mask*u*v*ds
f = LinearForm(fes); f += alpha*on_in*v*ds  # target V=1 inner, 0 outer
A.Assemble(); f.Assemble()
gfu = GridFunction(fes)
gfu.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")*f.vec
C_per_L = Integrate(grad(gfu)*grad(gfu), mesh, order=10) / L   # = 2W/V^2 / L, V=1
```

Same pattern (swap eps->mu/sigma/k, or curl-curl for magnetics) carries any field
solve onto curved high-order hex. Loader rule from `troubleshooting` still applies: read the
curving with `Mesh()` + high `order=` quadrature, never `m.Curve()` (no CAD ref in the .vol).

## EIGENVALUE solves on the high-order hex mesh (not just source-driven fields)

The loaded curved-hex `.vol` also supports EIGENVALUE problems, where curving matters most
(eigenvalues are sensitive to the boundary shape). The Laplace-Dirichlet spectrum
`-nabla^2 u = lambda u, u=0 on the surface` of a meshed BALL has the exact lowest eigenvalue
`(pi/R)^2` (radial s-mode `sin(pi r/R)/r`), then `(4.493409/R)^2` (l=1, triply degenerate).

Build a hex SPHERE (`volume 1 scheme sphere` O-grid), `block 1 add volume all`,
`export netgen ... order 3`; in a separate ngsolve process load the `.vol` AS-IS (no
`mesh.Curve`) and call
`radia_mcp.radia_ngsolve.waveguide.laplace_dirichlet_eigenvalues(mesh, n, order=3)`. On a coarse
56-hex order-3 sphere (R=0.5) this gives lambda_1 = 39.450 vs exact `4 pi^2 = 39.478` (rel err
7e-4) and lambda_2 = 80.71 vs the l=1 mode 80.76 -- the curved hex recovers the round-domain
spectrum on very few elements. These eigenvalues are also the modal DECAY rates of the transient
heat equation (`multiphysics.solve_heat_transient`, `T_n ~ exp(-alpha lambda_n t)`).

The same works on a hex CYLINDER (`create cylinder height Lz radius R`): the Dirichlet spectrum
is `lambda = (j_mn/R)^2 + (p pi/Lz)^2` (Bessel-zero radial + axial); on a 544-hex order-3 cylinder
(R=0.5, Lz=1) lambda_1 = 33.0022 vs the exact `(j01/R)^2+(pi/Lz)^2 = 33.0024` (rel err 4e-6), with
lambda_2/lambda_3 matching the (0,1,2) and (1,1,1) modes -- the curved lateral surface is carried
by the order-3 export.
"""

TROUBLESHOOTING = """
# Troubleshooting High-Order Curving

## export netgen() Fails or Produces Wrong Results

### Symptom: RuntimeError during Curve()

**Cause 1**: Missing boundary blocks
```
Fix: Ensure both volume and surface blocks are registered:
     cubit.cmd("block 1 add tet all")
     cubit.cmd("block 2 add tri all")
```

**Cause 2**: Mesh quality too poor for high-order curving
```
Fix: Reduce element size or improve mesh quality:
     cubit.cmd("volume all size 0.05")  # Smaller elements
     cubit.cmd("smooth volume all")     # Improve quality
```

## Empty Mesh (0 Elements)

**Cause**: No blocks registered
```
Fix: cubit.cmd("block 1 add tet all")
     cubit.cmd("block 2 add tri all")
```

## Missing Boundary Elements

**Cause**: Only volume element block, no surface element block
```
Fix: cubit.cmd("block 2 add tri all")  # Add boundary elements
```

## Volume Error > 1%

**Cause**: 1st order mesh without curving
```
Fix: Use export netgen(cubit, order=2) or higher
     Or use Gmsh 2nd order alternative
```

## Netgen Import Error: "No module named 'netgen'"

```
Fix: Use system Python with CUBIT_PATH environment variable.
     Cubit's bundled Python cannot import ngsolve.
     System Python with CUBIT_PATH can access BOTH Cubit API and NGSolve.

     # Step 1: Set CUBIT_PATH
     set CUBIT_PATH="C:/Program Files/Coreform Cubit 2025.12/bin"

     # Step 2: Run with system Python (which has NGSolve installed)
     python my_script.py
```

In the script:
```python
import sys, os

# CRITICAL: Import NGSolve BEFORE Cubit to avoid DLL conflicts
import ngsolve
from ngsolve import Mesh

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])
```

**Key insight**: By using system Python with `CUBIT_PATH`, scripts can access both
the Cubit API and NGSolve/Netgen simultaneously. This is essential for the
export netgen() workflow.

## NGSolve/Cubit DLL Conflict: Import Order Matters

When using system Python with both NGSolve and Cubit, **NGSolve MUST be imported
BEFORE Cubit**. If Cubit is imported first, its bundled DLLs (VTK, etc.) conflict
with NGSolve's Netgen library, causing `ImportError: initialization failed` on
`from netgen import libngpy`.

**Correct import order** (NGSolve first):
```python
import ngsolve                    # MUST be first - loads Netgen DLLs cleanly
from ngsolve import Mesh

import sys, os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
import cubit                      # Safe: Netgen DLLs already loaded
cubit.init(['cubit', '-nojournal', '-batch'])
```

**Wrong import order** (causes DLL conflict):
```python
import sys
sys.path.append("C:/Program Files/Coreform Cubit 2025.12/bin")
import cubit                      # Loads Cubit's bundled VTK DLLs
import ngsolve                    # FAILS - Netgen can't initialize
```

**Root cause**: Cubit bundles its own versions of VTK and other shared libraries.
When `import cubit` executes, these DLLs are loaded into the process. When NGSolve
subsequently tries to load Netgen's `libngpy`, the already-loaded Cubit DLLs
conflict with Netgen's expected library versions, causing initialization failure.

**Rule of thumb**: Always `import ngsolve` (and any `from netgen...` imports) at
the very top of the script, before adding Cubit to `sys.path` or importing `cubit`.

## Standalone cubit.init() Segfaults in the Radia Panel (verified 2026-06)

**Symptom**: From system Python, `import cubit; cubit.init([...])` prints the banner,
auto-plays the generated Radia startup shim
(`%ProgramData%/Radia/Cubit/radia_startup.py` for all-users installs),
then crashes with exit code -1073741819 (0xC0000005 access violation) BEFORE your first
`cubit.cmd()` runs. Happens with or without `-nographics` and regardless of the ngsolve
import order above. Cause: when the Radia Cubit *panel* plugin is installed, its
startup shim loads the PySide6 toolbar under headless embedded Python (no GUI main
window). The single-process recipe above only works when that panel plugin is absent.

**Fix -- robust two-process pattern** (use this whenever the Radia panel is installed):

  1. EXPORT in a child process via the real Cubit executable -- it degrades the panel
     gracefully ("Cubit main window not found -- Radia Export menu not installed", then
     continues) instead of crashing:

         coreform_cubit.exe -nographics -batch -nojournal mesh_export.py

     In a script launched this way, `cubit` is PRE-INJECTED -- do NOT `import cubit` or
     call `cubit.init()` (that re-inits and can crash). Just use `cubit.cmd(...)`.

  2. LOAD the resulting `.vol` in a SEPARATE ngsolve-only process (no `import cubit`).

  This split also sidesteps the NGSolve/Cubit DLL conflict entirely -- the two libraries
  never share a process, so import order stops mattering.

**`-batch` .py playback is LINE-ORIENTED**: multi-line Python compound statements
(for/if/try blocks, multi-line parenthesized tuples) break with
`SyntaxError: '(' was never closed`. Keep every statement on ONE physical line.
`print()` output may be swallowed -- write results to a file and read it back.

**Exit code 2 / "ERROR: Errors found during session." is benign teardown noise** from the
panel's "main window not found"; verify success from the OUTPUT FILE, not the exit code.

Verified 2026-06-14: hex cylinder (R=0.5, H=2), order-3 `export netgen` (572 hexes,
770 -> 25368 curved nodes) -> NGSolve volume 1.57082 vs pi*R^2*H = 1.57080, error
0.0017% (vs ~0.4% for a 1st-order hex cylinder: curving cuts the error ~250x).
"""

DELETED_APIS = """
# Deleted APIs (Replaced by export netgen APREPRO command)

The following APIs have been **completely removed**.
Do NOT use them — they no longer exist.

## Removed Functions

| Function | Replacement |
|----------|-------------|
| `export_NetgenMesh()` | `export netgen` APREPRO command |
| `export_netgen()` (alias) | `export netgen` APREPRO command |
| `export_netgen_with_names()` | `export netgen` APREPRO command |
| `extract_curved_mesh()` | `export netgen` APREPRO command |
| `name_occ_faces()` | Not needed (no OCC geometry) |
| `set_cylinder_geominfo()` | Not needed (ACIS handles curving) |
| `set_sphere_geominfo()` | Not needed |
| `set_torus_geominfo()` | Not needed |
| `set_cone_geominfo()` | Not needed |
| `compute_cylinder_uv()` | Not needed |
| `compute_sphere_uv()` | Not needed |
| `compute_torus_uv()` | Not needed |
| `compute_cone_uv()` | Not needed |

## Why They Were Removed

The old workflow required multiple steps:
1. STEP export from Cubit (ACIS -> STEP)
2. STEP reimport into Cubit (to match OCC face topology)
3. OCC geometry loading (OCCGeometry(step_file))
4. Mesh export with geometry reference (export_netgen(geometry=geo))
5. Per-shape UV computation (set_*_geominfo())
6. Manual mesh.Curve(order)

`export netgen` replaces ALL of these steps with a single APREPRO command.
It uses CallbackGeometry to delegate surface projection to Cubit's ACIS
kernel directly, without any STEP file exchange or OCC geometry.

## Migration Guide

### Old Code (REMOVED)
```python
# This code NO LONGER WORKS — all these functions are deleted
geo = OCCGeometry("cylinder.step")
ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)    # DELETED
cubit_mesh_export.set_cylinder_geominfo(ngmesh, radius=R, height=H)  # DELETED
mesh = Mesh(ngmesh)
mesh.Curve(3)
```

### New Code
```python
import tempfile
from ngsolve import Mesh
# Single APREPRO command replaces everything
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```

### Old Complex Workflow (REMOVED)
```python
# This code NO LONGER WORKS
cubit_mesh_export.name_occ_faces(shape)
shape.WriteStep("geometry.step")
geo = OCCGeometry("geometry.step")
cubit.cmd('import step "geometry.step" noheal')
ngmesh = cubit_mesh_export.export_netgen_with_names(cubit, geo)
cubit_mesh_export.set_cylinder_geominfo(ngmesh, radius=R_HOLE, height=H)
mesh = Mesh(ngmesh)
mesh.Curve(3)
```

### New Code
```python
import tempfile
from ngsolve import Mesh
# Create geometry directly in Cubit, no OCC needed
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("create cylinder height 4 radius 0.3")
cubit.cmd("subtract volume 2 from volume 1")
# ... mesh and blocks ...
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```
"""


def get_netgen_documentation(workflow: str = "overview") -> str:
	"""Return Netgen workflow documentation by topic."""

	KELVIN_AUTO = """
# Kelvin Auto-Add in Cubit Workflow (2026-04-14)

## Overview

Kelvin open-boundary transformation is automatically added when the user
clicks "Radia-NGSolve" -> OK in the Cubit GUI. Kelvin is added automatically
-- there is no separate "Kelvin Transform" menu item.

## How It Works (register_toolbar.py)

1. User creates physical geometry (coil + air + optional workpiece hole)
2. User clicks "Radia-NGSolve" -> selects analysis mode -> OK
3. register_toolbar.py checks: is "kelvin" block already present?
   - YES -> skip, proceed to export
   - NO  -> auto-detect R, symmetry, and call add_kelvin_cubit()
4. export netgen -> .vol (with Kelvin + periodic identification)
5. Launch analysis window (calc_fem.py reads .vol with Kelvin)

## Auto-Detection Logic

### Sphere Radius (R)
- Find the "air" block volumes
- Find the largest-area surface on those volumes -> outer boundary
- R = max vertex distance from origin on that surface
- Works for full sphere, hemisphere (1/2), quarter sphere (1/4), octant (1/8)

### Symmetry Planes
- Check all vertices of air volumes
- If ALL vertices have x >= 0 AND some are at x = 0 -> "x" symmetry
- Same for y, z
- Result: [] (full), ["z"] (1/2), ["x","z"] (1/4), ["x","y","z"] (1/8)

### Offset Direction
- auto_offset_direction(symmetry) picks a free axis
- ["z"] -> offset in x; ["x","z"] -> offset in y

## What add_kelvin_cubit() Does

1. Creates exterior sphere (same R, at offset position)
2. Webcuts for symmetry planes
3. Copies mesh from interior sphere surface -> exterior sphere (1:1 nodes)
4. Tet-meshes the exterior sphere volumes
5. Assigns blocks ("kelvin"), sidesets ("kelvin_int", "kelvin_ext")
6. Creates GND vertex + nodeset at Kelvin center

## CRITICAL: Kelvin Domain Must Be Tet

Spherical geometry is best approximated by triangular high-order elements.
Hex elements on sphere surfaces introduce systematic geometry error.
Always use `scheme tetmesh` for Kelvin volumes.

## User Requirements

The user only needs to provide:
- Coil geometry (block "coil" with source/sink sidesets)
- Air sphere (block "air") containing the coil
- Optional: workpiece hole (sideset "sibc") -- subtracted from air, NOT meshed

The user does NOT need to:
- Create Kelvin geometry manually
- Know the Kelvin sphere radius or offset
- Specify symmetry planes
- Add Kelvin manually (it is auto-added on Solve)

## .jou Example (Minimal)

```python
# User only creates physical geometry:
reset
create sphere radius 0.06     # air sphere
sweep surface 1 axis ...      # coil (inside air sphere)
subtract volume <coil> from volume <air> keep_tool
imprint all; merge all
mesh volume all
block 1 add volume <coil>; block 1 name "coil"
block 2 add volume <air>;  block 2 name "air"
sideset 1 add surface <gap1>; sideset 1 name "source"
sideset 2 add surface <gap2>; sideset 2 name "sink"
# That's it. Kelvin is added automatically on "Solve".
```
"""

	topics = {
		"overview": WORKFLOW_OVERVIEW,
		"export netgen": WORKFLOW_EXPORT_CURVED,
		"simple_cylinder": WORKFLOW_CYLINDER,
		"simple_sphere": WORKFLOW_SPHERE,
		"simple_torus": WORKFLOW_TORUS,
		"complex": WORKFLOW_COMPLEX,
		"complex_named": WORKFLOW_COMPLEX,  # Alias for backward compat
		"accuracy": WORKFLOW_ACCURACY,
		"field_solve": WORKFLOW_FIELD_SOLVE,
		"capacitance": WORKFLOW_FIELD_SOLVE,
		"kelvin_auto": KELVIN_AUTO,
		"gmsh_2nd_order": WORKFLOW_GMSH_2ND_ORDER,
		"troubleshooting": TROUBLESHOOTING,
		"deleted_apis": DELETED_APIS,
		# Legacy aliases that redirect to new content
		"setgeominfo_api": DELETED_APIS,
		"seam_problem": DELETED_APIS,
		"tolerance_tuning": DELETED_APIS,
		"uv_math": DELETED_APIS,
		"multi_surface": WORKFLOW_COMPLEX,
		"freeform": WORKFLOW_COMPLEX,
		"simple_cone": WORKFLOW_EXPORT_CURVED,
	}

	workflow = workflow.lower().strip()
	if workflow == "all":
		main_topics = [
			WORKFLOW_OVERVIEW, WORKFLOW_EXPORT_CURVED,
			WORKFLOW_CYLINDER, WORKFLOW_SPHERE, WORKFLOW_TORUS,
			WORKFLOW_COMPLEX, WORKFLOW_GMSH_2ND_ORDER,
			WORKFLOW_ACCURACY, WORKFLOW_FIELD_SOLVE, TROUBLESHOOTING, DELETED_APIS,
		]
		return "\n\n".join(main_topics)
	elif workflow in topics:
		return topics[workflow]
	else:
		return (
			f"Unknown workflow: '{workflow}'. "
			f"Available: all, {', '.join(k for k in topics if k not in ('complex_named', 'setgeominfo_api', 'seam_problem', 'tolerance_tuning', 'uv_math', 'multi_surface', 'freeform', 'simple_cone'))}"
		)
