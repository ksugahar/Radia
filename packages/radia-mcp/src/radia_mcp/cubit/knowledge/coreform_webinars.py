# -*- coding: utf-8 -*-
"""Coreform Cubit webinar / tutorial knowledge, synthesized from the official
Coreform YouTube channel auto-captions (https://www.youtube.com/@Coreform,
"Coreform Cubit Tutorials" playlist).

PUBLIC + attributed + PARAPHRASED (no verbatim transcript reproduction; ASR
captions cleaned and corrected against Cubit domain knowledge). This is an
incremental ingestion; topics are added per batch.

Batch 1 (20 videos): Python scripting & automation, hex-meshing fundamentals
& strategy, and third-party-solver front-end workflows. Complements the
sibling knowledge modules (cpp_sdk, custom_toolbar, scripting, mesh_diagnostics,
netgen_workflow) -- those stay authoritative for their narrow topics; this
module is the broader technique knowledge mined from the tutorial corpus.
"""

PYTHON_AUTOMATION = r"""
## Cubit Python Scripting & Automation

Coreform Cubit exposes everything through two layers: the **Cubit Command
Language (CCL)** -- the line-oriented commands the GUI emits and that you write
to `.jou` journal files -- and a **Python API** wrapping that command engine plus
the ACIS geometry kernel. CCL is what most users journal; Python is for power
users who need data structures, control flow, and external libraries. Python is
a true superset: CCL is reachable from Python as a single `cubit.cmd(...)` call,
so you lose nothing by moving to Python and gain lists/dicts/loops,
numpy/scipy/sympy/sklearn, file I/O, and subprocess control.

**Key rule (stated across videos):** issue geometry/mesh *operations* through
`cubit.cmd(...)`, NOT the direct C++ object-mutation methods. The Python API is
SWIG-generated; direct mutators touch the kernel's raw C++ memory and Python
manages those poorly. Use `cubit.cmd` to CREATE/MODIFY; reserve the object
methods for QUERIES.

### Core working functions (the ~90% you use)
```python
import cubit
cubit.cmd("reset")
cubit.cmd("brick x 10 y 5 z 1")
x, y, z = 10, 5, 1
cubit.cmd(f"brick x {x} y {y} z {z}")          # f-string splices computed values
```
Query / introspection (build loops, avoid hard-coded IDs):
- `cubit.parse_cubit_list("volume", "with x_coord > 0")` -> tuple of matching IDs
  (mirrors any CCL `list ... with ...`; prototype the `list` query first, then wrap).
  Works on mesh entities too: `parse_cubit_list("node", "in hex 1")`.
- `cubit.get_last_id("volume")` -> ID of the most recently created entity
  ("I just made it with cubit.cmd, give me its ID") -- the relative-indexing
  pattern that immunizes scripts against ID renumbering.
- `cubit.get_entities("volume")` -> all current IDs of a type (reflects deletions).
- `cubit.get_id_string(list_of_ids)` -> a CCL-parsable ID string to splice many
  IDs into one command.
- `cubit.get_geometry_owner("node 595")` -> the geometry entity owning a mesh entity (per video).
- `cubit.is_merged("surface 6")` -> whether a surface is merged (per video).
Entity objects -- `cubit.volume(id)`, `.surface(id)`, `.curve(id)`, `.vertex(id)` --
expose safe QUERY methods: surface type (planar/cylindrical), principal curvatures,
normals at a location, nearest-point-on-surface; curve `position_from_u(0.5)` and
`tangent(...)` (method spellings per video; verify in the Python Interface docs).

### Record / replay: GUI -> CCL -> Python
1. Keep the **GUI and Journal Editor side by side**; every GUI action is journaled
   as the equivalent CCL command (you always get correct syntax this way).
2. Copy commands from the **History** tab into a `.jou`.
3. **Convert journal -> Python in one click** (Journal Editor "Python"/wrap button):
   imports `cubit` and wraps each line as `cubit.cmd("...")`.
4. **Parameterize**: replace numeric literals with variables (aprepro in CCL, or
   Python variables via f-strings). Rewriting a recorded session into a parametric
   script with `get_last_id`/relative indexing collapses size dramatically (one
   example ~500 CCL lines -> ~160 Python lines) and survives ID renumbering.
Robustness habit: begin scripts with `reset`/`undo` (and `delete mesh` before the
meshing section) so a mid-session mistake is fixed by rerunning from a clean state.

### aprepro vs Python parameterization
- **aprepro** (`#{name = value}` / `#{name}` substitution) is the native CCL way to
  declare/reuse variables and to template solver input files (`{var}` slots).
- **Python f-strings/`.format`** are preferred once in Python: same substitution plus
  full expression evaluation, computed values, and ecosystem integration. Use aprepro
  for pure-CCL journals / text templating; Python variables when any computation,
  looping, or data structure is involved.

### Extended (criteria-based) selection
`... with <criterion>` acts on entities by geometry instead of ID -- invaluable for
imported geometry with unknown IDs:
```text
delete surface all with x_coord < 0
mesh curve with length < cyl_radius and y_coord < height
```
From Python these criteria strings are the 2nd argument to `parse_cubit_list`.

### Headless / batch / external execution
1. **Cubit driving Python (GUI session):** the scripting tab is a live Python
   interpreter inside Cubit (Tools -> Options -> Layout shows the script tab; set
   Python 3 under General and reboot Cubit).
2. **Python driving Cubit (external/headless):**
   ```python
   import sys
   sys.path.append(r"<cubit_install>/bin")   # so `import cubit` resolves
   import cubit
   cubit.init(["cubit", "-nojournal", "-noecho"])   # MUST precede any cubit.cmd
   cubit.cmd("brick x 1")
   ```
   - `cubit.init([...])` is mandatory; forgetting it is the most common "crash on
     first command" error. Useful flags: `-nojournal`, `-noecho`, banner/warning/info
     suppression. (Machine note: a Radia panel present at startup can segfault the
     standalone `cubit.init()`; a 2-process split avoids it -- environment-specific,
     not from these videos.)
   - External Python's **version must match Cubit's bundled Python** (any distro --
     e.g. a conda env -- as long as the version matches).
   - Add modules with the *packaged* pip:
     `<install>/bin/python3/python -m pip install --user numpy scipy sympy`.

### Docs path for the API
Online/packaged help -> **Appendix -> Python -> Python Interface**; entity classes and
the **Cubit Interface** function listing. Many methods are documented in C/C++
signature form (SWIG) -- another reason to prefer `cubit.cmd`.

### Parametric sweeps & convergence/refinement studies (workflow structure)
The recurring loop is **mesh -> solve -> query result -> decide -> remesh**, all in Python:
- **Convergence study:** start coarse; run solver; read a result at a FIXED spatial
  probe (use the API to find the node nearest a fixed location like `0 0 0.1` so the
  same physical point is sampled after every remesh); `delete mesh`; remesh finer;
  compare against the previous value; loop until the change is within tolerance. Append
  values to a list for plotting.
- **Solution-adaptive refinement:** after a solve, read element/integration-point
  results; flag elements with too-large within-element variation; refine just those; repeat.
- **Geometry-parameter optimization:** vary a parameter (e.g. hole radius via
  `tweak surface ... offset ...`, previewable), remesh, solve, read stress at a fixed
  probe node, `while`-loop until under an allowable or max iterations.
- **Refinement study via external solver (openCFS):** wrap the geometry/mesh script as a
  Python class using `get_last_id` + `for`-loops; sweep a resolution variable
  (e.g. `skin_interval` 10->25); per value remesh, return element count, export Exodus,
  invoke the solver binary as a black box, read back a monitored quantity (magnetic
  energy), stop when it stabilizes.
- **Geodynamics remesh loop (batch mode, no GUI):** generate an initial sphere/shell,
  deform node coords externally by editing the Exodus `.exo` in place; adaptive
  refinement writes a per-node **sizing variable** into the Exodus, re-imports, `delete
  mesh`, remeshes; also **untangle** negative-Jacobian cells and build fault interfaces
  by `subtract`-ing projected geometries. Batch efficiency is what makes 36-remesh
  inverse problems feasible.

### DAKOTA + MOOSE optimization (how Cubit regenerates the mesh per design point)
- **DAKOTA** (Sandia, open-source) = outer optimizer/workflow manager (global/local,
  surrogate/EGO, GA, Newton; DOE/list studies; per-eval working-dir create/cleanup;
  parallel/async; failure handling -> substitute a large objective; restart files).
  **MOOSE** = the FEM solver. **Python** = glue. **Cubit** = per-evaluation mesher.
- **Black-box file interface:** a `.in` study defines variables/responses; DAKOTA writes
  a params file matching a `params.template` (`{var}` slots); a `dakota_interface.py` +
  `qoi.py` regex parser (boilerplate from the Dakota GUI wizard) read params/parse
  results; a master `evaluate_iteration.py` does **build_model (runs Cubit) ->
  submit_moose (subprocess/slurmpy) -> compute objectives/constraints -> write results back**.
- **Per design point the mesh is regenerated:** `build_model` rebuilds geometry from the
  current design variables and meshes it (building geometry from parameters and naming
  blocks/sets at creation sidesteps ID instability). Swapping solvers changes only
  `build_model`/`submit_*`.
- **Sizing functions via Exodus (the unique Cubit/Exodus capability):** write an arbitrary
  scalar field (numpy, or a prior solution's von-Mises/error field) into the Exodus mesh
  with the **SEACAS `exodus.py`** module (`import exodus`; add a nodal variable array;
  close to save), then in Cubit `import sizing function` and map min/max field -> min/max
  element size for a graded mesh (docs: Mesh Operations -> Adaptivity and Sizing
  Functions). Exodus is Cubit's native format AND MOOSE's preferred format -> no conversion.

### Custom mesh import (unsupported formats, e.g. Gmsh .msh)
Parse the file in Python and feed elements into the **`MeshImport` class**; then
`cubit.create_mesh_geometry_tet(...)` to convert to a mesh-based (facet) geometry, and
**delete the mesh-import object afterward** (the C++ wrapper holds the previous import and
crashes the NEXT import otherwise). Result: a facet geometry you can still add sets to.

### Packaging scripts for reuse
- **Custom toolbar button** (Tools -> Custom Toolbar Editor): add a script button (icon
  optional) pointing at a `.py`; clicking runs it. (See the `custom_toolbar` knowledge.)
- **Custom command panel:** build the UI in Qt Designer, wire via XML, adapt the script
  into an "apply" routine (docs: Customizing the User Experience -> Adding Command Panels).

### Reproducibility gotchas
- Pipe solver/Cubit `stdout`/`stderr` to a log so you can track progress and errors.
- CCL is one-command-per-line (why copy-from-history and journal->Python wrapping are clean).
- `-nojournal`/`-noecho` + banner/warning/info suppression keep batch logs readable.
- Prefer relative IDs (`get_last_id`/`get_entities`) over hard-coded IDs; name
  blocks/nodesets/sidesets at creation; build geometry from parameters.

### NOT covered (don't infer)
- No PyPI distribution of the Cubit Python API (it wraps proprietary ACIS; the free
  **Cubit Learn** edition, 50,000-element export cap, is the accessible route).
- No built-in sketcher (one presenter solved constraints with sympy, fed via cubit.cmd).
- Exact spellings of some entity query methods are as spoken -- verify in the docs.

### Sources
Introduction to Python scripting in Coreform Cubit; Coreform Cubit: Automate Workflows
with Python; Automated optimization workflow using Cubit's Python API, MOOSE & DAKOTA;
How to perform a refinement study with Python scripting (openCFS); 2025 update: Creating
a mesh with Python; Using Cubit scripting to study the geodynamics of a planetary satellite.
"""

MESHING_STRATEGY = r"""
## Hex Meshing Fundamentals & Strategy in Coreform Cubit

### Why hex meshing is hard (and when to bother)
For tets a valid surface *triangulation* suffices (robust auto algorithms). For hexes
there is **no general algorithm** to fill an arbitrary bounding quad surface mesh -- an
open problem -- so your job is to **decompose geometry into blocks Cubit can fill**.
When you need hex "depends on your physics": nonlinear/plasticity strongly favors hex (or
quadratic tets); pure linear-elastic with a safety factor may not. Hex payoff: higher
accuracy per DOF + better conditioning (easier nonlinear convergence). Hex cost: more
structure -> you **cannot vary element size as freely**; a fine region forced next to a
coarse one yields high-aspect/skewed elements through the swept direction.

### The schemes and when each applies
- **map** -- fully structured; topologically box-like (logical cube: 6 sides, 8 corners).
- **submap** -- structured on a union of mappable sub-blocks (L-shape; many-to-one on a
  face). Can be rejected by a tiny spurious curve ("does not admit submap") -> fix via
  collapse/composite/detect-small-features.
- **sweep** -- mesh source surface(s), sweep to a target. Supports one-to-one and
  many-to-one, **NOT one-to-many**. The workhorse.
- **multisweep** -- many-source-to-one-target (per video).
- **pave** -- all-quad surface scheme (surface workhorse); often "pave-then-sweep".
- **polyhedron** -- semi-structured; decomposes a volume AS IF a convex polyhedron and
  fills with HEXES (not polyhedral elements). Good for a cube with a spherical corner
  removed. Auto-detection rarely picks it -> assign manually (`volume <id> scheme polyhedron`).
- **sphere** -- primitive scheme for spherical bodies (O-grid).
Recognize a sweepable volume: clear source/target with a sweep direction; mappable/pavable
side surfaces. Set manually when auto fails: `volume <id> scheme sweep source <surf...>
target <surf>` (Ctrl-multiselect sources, then the single target).

### Decomposition strategy (the heart of semi-automatic hex meshing)
Per part: **clean -> decompose (web cuts) -> composite -> imprint & merge -> assign
schemes/sizes -> mesh most-constrained volumes first.**
- **Web cut = scissors:** completely disconnects resulting volumes (non-conforming until
  imprint+merge). Types: `plane`; plane-from-curve (pick fraction along a curve -> tangent
  plane); plane from 3 vertices; sheet **extended from surface**; **sweep surface**
  (perpendicular, inward/outward -- preferred when available); sweep curve along a vector;
  `cylinder radius`; loop; or cut with a custom **tool body** (skin curves -> cover ->
  unite -> `webcut volume <id> tool body <id>`, tool auto-deleted).
- **Where to cut (heuristics):** break the body into recognizable regions first; **convert
  tangencies into perpendicularities** (a fillet cut by a plane leaves a knife-edge that
  gives terrible hexes -- cut so the boundary meets the surface perpendicularly);
  avoid cutting through spline surfaces (push cuts onto analytic geometry); add
  mesh-quality cuts to make near-rectangular map regions; cut at/near maximum-curvature so
  the worst angle improves (not degrades) under refinement.
- **Imprint & merge -> conforming meshes:** *imprint* copies topology of one volume onto a
  neighbor's face (volumes stay distinct); *merge* fuses the co-located duplicate surfaces
  so they SHARE the surface and its mesh (this is what makes it conforming). `imprint merge`
  / `merge volume <a> <b>`; `draw surface with is_merged` shows merged interfaces;
  **tolerant imprint** closes small gaps (use cautiously). Verify each volume meshes
  independently BEFORE imprint+merge (merging forces shared sizes -> harder).
- **regularize** (`regularize volume all`): removes redundant topology (interior imprinted
  curves, extra vertices on a straight edge) -- frequently the fix when a volume won't
  imprint/mesh.
- **Mesh order rule:** mesh the **most-constrained volumes first** (polyhedron / fillet /
  shoulder regions), then loose cylinders/bulk; restrict updates with `... except with
  is_meshed` to avoid re-meshing finished volumes.

### Virtual / composite geometry -- repairing dirty CAD without changing the B-rep
- **Real (ACIS) ops:** `tweak` (blunt a tangency; tweak curve remove), `remove surface ...
  extend` (recover a sharp edge / rebuild a fillet), boolean `remove overlap`, `heal`,
  `regularize`.
- **Composite (key virtual tool):** merges adjacent surfaces into ONE macro surface the
  mesher treats as a single face (**dashed curves = the mesher ignores this boundary**).
  Uses: make a non-sweepable topology sweepable (many targets -> one); remove
  sliver/over-constraining faces from web cuts; enable map/sweep on a many-faced source.
  Selection helpers (right-click): select similar / select cavity / select continuous /
  blend chain; Tab cycles stacked selections.
- Other repairs: `collapse` (curve/surface) a tiny entity into a neighbor; `remove surface
  ... extend`; **split periodic surface** (a wrap-around cylindrical surface often won't
  map until split adds the seam); create-surface-from-bounding-curves / skin / stitch /
  unite (cover gaps, make IGES sheet bodies watertight); **detect small features**
  (ITEM -> remove small features) finds zero-area surfaces / tiny curves that silently break meshing.

### Building a mesh from DIRTY CAD
Dirty CAD comes from lossy neutral-format translation and the math limit of trimmed/spline
surfaces. First steps: `vertex vis on` (reveals double-curves); Properties page shows
analytic vs **spline** (a Python helper can color splines red -- avoid web-cutting through
them); **heal** with the analyze option (diagnose bad vertices/co-edges; auto-heal can
convert exactly-analytic splines back to analytic); manage gaps/overlaps (boolean
`remove overlap`, often scripted over the volume list).
**Absolute tolerance & scaling (critical):** Cubit is unitless; `set geometry accuracy`
default **1e-6**. Strategy: find the smallest feature you care about and **scale the model
so that feature ~ 1** (nearest power of 10): e.g. a 0.002 layer -> `volume all scale 1000`.
Sub-features near 1e-3 become ignorable noise; a well-scaled, mostly-analytic model often
heals to ZERO bad geometry. Recover units on export with `transform mesh output scale
<1/factor>` (`transform mesh output reset` before changing; cumulative).

### Sizing control
`volume <id> size <val>` / approximate-size slider (preview, then read the chosen value);
per-curve interval/scheme (e.g. 20 elements along a curve, or a curvature scheme);
**`set maximum arc span <deg>`** forces minimum segments per arc (default ~400 = ignore;
set 10 / 2.5 for finer curve capture); tet **gradation** (1.3 fast/fewer; 1.05 slow/higher
quality) and **deviation angle** (e.g. 5deg); **minimum tet layers** / `tet mesh proximity
layers on ... value 3` forces N tets through thin regions (NOT a true boundary layer);
`volume ... refine` and mesh-scaling for large meshes; sizing can be driven from an
external mesh's sizing function (stress/density field). `set default element tri/tet`
makes auto-meshing produce tets for the session (handy for big assemblies).

### Quality, diagnostics, smoothing
**Scaled Jacobian** is the recommended metric (accounts for size; inversion -> 0). Targets:
0..1, **> ~0.3 is "pretty darn good"**; well-handled fillet/sweep ~0.7-0.9; left tangencies
~0.01. Caveat: the *real* metric is whether your SOLUTION results make sense (always do
V&V). Inspect by coloring by scaled Jacobian, click an element to read it. **Smoothing:**
Laplacian (`smooth ...`, include boundary nodes) recovers inverted/poor elements from
spline projections; tighten `set smooth tolerance` for finer meshes (reset to default 0.05).

### ITEM wizard & meshability
The ITEM tab gives a guided pipeline: build geometry -> preview size for an element budget
-> prepare/build meshable topology (web cuts) -> mesh. **Check mesh ability** lists each
volume meshable/not; **Decompose volume** proposes candidate web cuts (cycle with arrows,
preview, execute -- the executed command appears in the panel as a learning aid).

### Sculpt (overlay/Cartesian) -- the dirty/organic fallback
For very dirty or **organic** geometry (bones, CT isosurfaces, cast parts) that can't be
decomposed: **sculpt** overlays a background Cartesian grid, keeps interior cells, and
sculpts to the boundary; runs in **parallel** (choose cores) to an element budget. Interior
hexes are near-perfect; **surface** hexes have Cartesian-origin skew/high-aspect (the known
weakness; thin regions may be missed). For organic models sculpt (or immersed/overset) is
often the only hex route; a quadratic **tet** is the pragmatic alternative.

### Decision heuristics & common failures
- Clean feature-based CAD: decompose (web cut + imprint/merge + composite) -> swept/mapped/
  polyhedron hexes. Very dirty/organic: sculpt or tets. Mixed/time-constrained: hex where
  cheap, tet the hard sub-volumes, single **pyramid** transition layer (mesh hex first).
- Won't imprint/mesh -> `regularize`. "does not admit submap" / tet collapses to sliver ->
  `detect small features` then collapse/composite/remove. Periodic surface won't map ->
  split it. Non-sweepable (many targets/pinching links) -> composite. Tangency knife-edges
  -> recut perpendicular. Gaps/overlaps -> fix scale+tolerance, `remove overlap`, tolerant imprint.

### Productivity tips
Journal as you go (History -> .jou with `#` comments; `play selected` to reorder/parametrize);
group parts into `part_<id>` groups (survive web cuts) for unit sizing/imprint; exploit
symmetry (mesh an octant -> `copy reflect` -> `merge` collocated nodes); `?` / `set ?` print
command syntax (`set developer commands on` for some); speed scripts with `set echo/info/
warning/journal off` and coarse `graphics tolerance`; STL ("mesh-based") is a weaker kernel
than ACIS -- prefer native/STEP/SAT, or import STL with feature angle then re-mesh on facet
surfaces / `tet mesh tri all` / sculpt.

### Sources
Coreform Cubit Basics: Hex Meshing Fundamentals; Hex-meshing in Coreform Cubit; Hex meshing
deep dive (automotive bushing); Semi-automatic hex meshing of complex parts; Strategies for
solving tricky meshing problems; Building a mesh from dirty CAD in Coreform Cubit.
"""

SOLVER_WORKFLOWS = r"""
## Cubit as a Front-End Mesher for Third-Party CAE/CFD Solvers

Cubit is a meshing / pre-processing tool, NOT a solver and NOT a parametric CAD modeler.

### Core mental model: Exodus shapes everything
Cubit's native format is **Exodus** (Sandia/DOE, open-source); its conventions
(blocks/sidesets/nodesets/element types) ARE Exodus conventions. Exodus is solver-agnostic
and stores BOTH mesh and results -> Cubit has no built-in post-processor (use **ParaView**/
VisIt on Exodus). Native Exodus consumers include **MOOSE** (INL) and lab codes; the broader
**SEACAS** project provides an Exodus Python module, binary<->ASCII and Exodus<->MATLAB
converters, a `decomp` parallel-partitioner, parallel I/O, and an Exodus `diff` for
regression/CI. Geometry kernel is **ACIS** -> importing native `.sat` avoids a translation.

### The block / sideset / nodeset model (mesh -> solver mapping)
| set | holds | solver meaning |
|---|---|---|
| **Block** | elements (a material region/part) | material, element type, section |
| **Sideset** | element faces/sides (surfaces) | surface BCs: pressures, fluxes, contact, CFD patches |
| **Nodeset** | nodes | nodal BCs: fixed displacements, loads, constraints |
- A **block** is where you choose the **element type** (hex8, hex20/27, tet4/tet10, wedge,
  pyramid, beam, shell, ...). Exodus also stores per-type **attributes** (beam area/moments,
  shell thickness) that map to solver section cards (see the Exodus manual shipped with Cubit).
- Blocks survive geometry edits: web-cutting a volume after assigning it keeps the sub-volumes
  in the same block. Assign blocks early.
- **Imprint+merge** makes neighbors share nodes (conformal) instead of tied contact.
- Cubit deliberately has **limited generic BC/contact/solver-card** support -- every solver
  names BCs differently; use a solver plugin, edit the deck externally, or build a custom
  pre-processor on the API.

### Export formats (only those the videos name)
| target | export | notes |
|---|---|---|
| MOOSE / DOE (Sierra, Nalu, Truchas) | **Exodus** (.e/Genesis) | native, lossless; IGA -> Exodus with Bezier extraction |
| Abaqus | `.inp` | often edit element formulation in the deck |
| Calculix | Abaqus-style `.inp` (plugin) | element names follow Abaqus (C3D8, C3D20) |
| LS-DYNA | `.k` | supported |
| Nastran / ANSYS / I-deas / Patran | respective | listed as supported (verify extensions in docs) |
| OpenFOAM | **polyMesh** dir | Cubit writes points/faces/owner/neighbour/boundary |
| Fluent | `.msh` | native |
- **meshio** (open-source Python) fills export gaps (consumes Exodus/.inp). **Polyhedral
  export is NOT supported** (even though Exodus allows it).
- **Critical gotcha -- Exodus stores topology, not formulation.** Element type in Exodus =
  node count/ordering only; quadrature (reduced/full), hybrid, incompatible modes are the
  solver's job. Abaqus/Calculix encode formulation in the NAME (`C3D8`, `C3D8R` reduced
  (verify), `C3D8I` incompatible modes, `C3D8IH` hybrid, `C3D20` quadratic) -> after export
  you frequently **edit the element-type string** in the deck (the Calculix plugin can
  override Cubit's element type with the Calculix type directly).

### Canonical workflow (order matters)
1. **Import/clean CAD** -- prefer the most direct translator; prefer ACIS `.sat` over STEP
   (Cubit's kernel IS ACIS). 2. **heal**. 3. **decompose** (web cuts along geometric
   features, not hard coords). 4. **defeature** (Cubit proposes blend-removal solutions to
   preview). 5. **composite** (virtual topology) -- merge surfaces so many-to-many sweeps
   become one-to-one, or ignore over-constraining curves; composite BEFORE imprint/merge.
   6. **imprint + merge** (after decomposition/defeature/composite). 7. **mesh** from the
   most-constrained volume outward (hint sweep direction/source-target where auto fails).
   8. **assign blocks/sidesets/nodesets**. 9. **export** in the solver's format. 10. **solve**
   externally; **post in ParaView**.
- **Quality before export:** scaled-Jacobian analysis (negative = inverted, shown red);
  smooth (Laplacian) to recover an initially-poor spline/fillet mesh -- don't give up on it.
  `?` lists valid commands for the current keyword. `sculpt` is the parallel hex mesher
  (voxel -> snap, akin to snappyHexMesh) but doesn't guarantee a valid all-hex mesh; classic
  Cubit meshing is largely serial (tet mesher multithreaded).

### MOOSE specifics
Exodus always; block/sideset/nodeset NAMES map into MOOSE's Mesh/BCs/Materials blocks.
Use enough through-thickness layers on thin members (~5 node layers, per video; verify).
MOOSE has runtime **h-adaptivity** (accurate but expensive -- run 1-2 iters, view in
ParaView, then hand-refine). MOOSE supports **CLI variable overrides** (block-path syntax)
and parametric-study blocks -> highly scriptable. **IGA in MOOSE:** build a Bezier
extraction in Cubit, write to Exodus, load in the deck (libMesh/ParaView extended for Bezier;
exact CAD even at coarse resolution; fewer DOFs). Coreform is building a MOOSE input editor
inside Cubit (roadmap/alpha).

### Calculix (via the Coreform plugin)
Open-source solver; Abaqus-like `.inp`. The plugin (GitHub; Tools -> Plugins -> restart)
adds a second model tree. Flow: Cubit prep (web-cut axles into clean cylinders, composite
sheet faces for a sweep direction, imprint+merge) -> Material dialog -> Section (solid) ->
Displacements/BCs -> Surface interaction (contact stiffness ~5-50x E) -> Contact pairs on
sidesets (**master = the COARSER mesh**, slave = finer) -> Field/history outputs (U/E/S) ->
Step (static; initial/min/max increment) -> Job (live monitor; converter threads parallel on
Linux). Post: load `.frd` back into Cubit, convert to ParaView; add **`PARTIAL`** on convert
for contact (values don't exist at every node); multiblock inspector + threshold show
in-contact surfaces. Gotchas: **unitless** (keep a consistent system; mm-N -> MPa); **no true
shells** (2D expanded to 3D "thick shell"); element naming follows Abaqus; can override Cubit
element type with the Calculix type.

### CFD generally + OpenFOAM
- **No single best CFD mesh** -- it depends on the (unknown) flow field; the dominant error is
  **false (numerical) diffusion** from upwind schemes when **flow direction is misaligned with
  mesh lines** -- it does NOT vanish with refinement (affects hex AND tet equally), only when
  mesh lines parallel the streamlines. So **align mesh with the flow**; build stream-aligned
  meshes rather than just refining. **Boundary layers:** Cubit's boundary-layer tool (pick wall
  surfaces, set layers/growth) grows a graded near-wall mesh; bias via `curve <id> scheme bias
  factor 1.25` + `propagate curve bias volume <id>`; combine with an O-grid. Resolution is
  memory-bound (~16 GB ~ 200 cells/direction uniform) -> design a mesh that reproduces the
  physics at realistic size. Hybrid (tet core + hex/boundary-layer near walls) is fine.
- **OpenFOAM (polyMesh):** for industrial geometry, Cubit generates `constant/polyMesh`
  (points/faces/owner/neighbour/boundary). The BC names you assign in Cubit become the **patch
  names** in the `boundary` file. You still hand-author `0/` (fields p/T/U with dimensions and
  per-patch types fixedValue/zeroGradient), `constant/` (transport/thermo), `system/`
  (controlDict/fvSchemes/fvSolution) -- **start from the closest bundled tutorial**, don't write
  from scratch. 2D via `empty` patches; axisymmetric via **wedge** BCs on a single sector.

### Python automation of the export step
Anything the GUI can do is in Python/CCL; convert journal->Python, then extend it to
rename/move the exported polyMesh into the OpenFOAM `run/` case dir, or to export Exodus named
by a test number and launch MOOSE via CLI overrides. A master CSV of parameters drives a
sweep: read table -> push part numbers into Cubit -> mesh (the same script often works across
similar parts) -> export -> launch solver. Mix Cubit/Python in one journal via `#!python` /
`#!cubit` shebangs (e.g. scipy.optimize to find max-curvature -> cubit.cmd web-cut there).
Performance: pre-fetch an entity object once, then call its query methods (don't re-traverse
`cubit.curve(id)` repeatedly).

### Building a full pre-processor on Cubit (OEM/API)
Two APIs: the **Python API** (scripting/automation; web back-ends can drive Cubit headless)
and a **C++ + Qt SDK** (OEMs add GUI widgets/logic; each custom action also emits a journal
command). Example: CAE Fidesis added custom material/block/section/BC dialogs and analysis
types on stock Cubit; their SaaS front-end is a browser talking to a Cubit-Python back-end;
journal replay enables simulation-driven design (tweak CAD -> re-import -> replay -> re-mesh/
re-solve).

### Solver-agnostic mesh-quality guidance (fatigue/structural)
Don't blindly defeature stress-concentrating features (cracks start there) -- but defeaturing
is a legitimate tradeoff ("a simulation you can't run is useless"). **Linear tets: avoid** for
accuracy (checkerboarding is a red flag); **quadratic tets:** OK but expensive; **linear hexes:
best compromise**. Resolve small radii (coarse facets = artificial stress concentrators). Use
the result gradient as a convergence proxy. True singularities (sharp bonded edges) never
converge -- model an explicit crack and study tearing-energy-vs-crack-size, or compare life at
a fixed small offset. Membrane skinning captures free-surface plane-stress where cracks
initiate. Tactics: exploit symmetry (copy-with-mesh + reflect), mesh critical+challenging
regions first, always try smoothing before abandoning a mesh, journal+script everything,
composite to drop over-constraining curves, imprint+merge for conformal meshes.

### Limitations to design around
No robust generic BC/contact/solver-card generation (by design); no built-in post-processor;
no sketcher / not a parametric modeler; no local coordinate systems (no native per-element
material orientation); no polyhedral export; meshing largely serial (sculpt is the parallel
hex option, no validity guarantee); free **Cubit Learn** capped at 50,000 exported elements.

### Sources
How A Systems Engineer Uses Coreform Cubit (MOOSE); Why engineers choose Coreform Cubit for
third-party solvers; Coreform Cubit: a powerful front end for third-party CAE solvers;
Improving MOOSE workflows through Coreform Cubit; Boosting Calculix with Better Meshes; Using
Coreform Cubit in a CFD workflow; Using Cubit and Python to develop OpenFOAM CFD analysis;
Devil in the details: accurate fatigue calculations with Cubit and Endurica CL.
"""

_TOPICS = {
    "python_automation": PYTHON_AUTOMATION,
    "meshing_strategy": MESHING_STRATEGY,
    "solver_workflows": SOLVER_WORKFLOWS,
}

_ALIASES = {
    "python": "python_automation", "automation": "python_automation",
    "scripting_api": "python_automation", "dakota": "python_automation",
    "api": "python_automation", "optimization": "python_automation",
    "meshing": "meshing_strategy", "hex": "meshing_strategy",
    "decomposition": "meshing_strategy", "webcut": "meshing_strategy",
    "dirty_cad": "meshing_strategy", "sculpt": "meshing_strategy",
    "sweep": "meshing_strategy", "composite": "meshing_strategy",
    "solvers": "solver_workflows", "exodus": "solver_workflows",
    "moose": "solver_workflows", "calculix": "solver_workflows",
    "openfoam": "solver_workflows", "cfd": "solver_workflows",
    "abaqus": "solver_workflows", "export": "solver_workflows",
}

_INDEX = """# Coreform Cubit webinar/tutorial knowledge (synthesized from @Coreform YouTube)

Batch 1 topics (get via cubit_docs "coreform_<topic>"):
  python_automation  - Cubit Python API + automation: cubit.cmd, query fns,
                       journal->Python, aprepro vs f-strings, headless init,
                       convergence/refinement loops, DAKOTA+MOOSE, sizing
                       functions via Exodus, custom mesh import.
  meshing_strategy   - Hex meshing: schemes (map/submap/sweep/pave/polyhedron),
                       web-cut decomposition, imprint&merge, composite/virtual
                       geometry, dirty-CAD repair + tolerance scaling, sizing,
                       scaled-Jacobian quality, ITEM wizard, sculpt fallback.
  solver_workflows   - Cubit as front-end for MOOSE/Abaqus/Calculix/OpenFOAM/
                       Fluent: block/sideset/nodeset model, Exodus, export
                       formats + formulation gotcha, CFD false-diffusion, OEM API.

Aliases: python/automation/dakota -> python_automation; hex/decomposition/
sculpt/dirty_cad -> meshing_strategy; moose/calculix/openfoam/exodus -> solver_workflows.
More batches (getting-started, advanced meshing, ML features, domain demos) pending.
"""


def get_coreform_webinar_documentation(topic: str = "index") -> str:
    """Coreform Cubit tutorial-corpus knowledge (synthesized from @Coreform).

    Args:
        topic: "index" (default) for the topic list, "all" to concatenate every
            section, a topic key (python_automation / meshing_strategy /
            solver_workflows), or an alias (python, hex, moose, ...).
    """
    topic = (topic or "index").lower().strip()
    if topic in ("index", "", "list", "help"):
        return _INDEX
    if topic == "all":
        return "\n\n".join(_TOPICS[k] for k in
                           ("python_automation", "meshing_strategy", "solver_workflows"))
    resolved = _ALIASES.get(topic, topic)
    if resolved in _TOPICS:
        return _TOPICS[resolved]
    return (f"Unknown coreform topic {topic!r}. Available: index, all, "
            f"{sorted(_TOPICS)}; aliases: {sorted(_ALIASES)}.")
