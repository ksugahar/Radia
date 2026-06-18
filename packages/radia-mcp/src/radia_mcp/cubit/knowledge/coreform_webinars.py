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

GETTING_STARTED = r"""
## Getting Started: GUI, the ITEM wizard, and the Command Panel

### GUI layout
Center = graphics window (select entities, watch updates). Left = model tree (volumes, then
surfaces/curves; blocks/sidesets/nodesets appear here). Properties page = info for the
selected entity (some fields editable: size, meshed-vs-true area check, block element type).
Command line = history of issued commands (the GUI emits commands; you can type your own).
Right tabbed panels switch between **Power Tools** (holds the **ITEM** wizard + diagnostics;
best for new users) and the **Command Panel** (graphical command front-end; intermediate+).
Missing panels: re-enable from the **View** menu (ITEM = wizard-hat icon). Command Panel
hierarchy (Tools>Options): classic `mode>entity>action` vs new `mode>action>entity`, plus a
breadcrumb option for small screens.

### Why a workflow (tet vs hex)
Tet meshing is ~solved (valid surface triangulation -> guaranteed fill). Hex is NOT a solved
problem -> a human decomposes geometry into Cubit's recipes: **map** (logically 4-sided ->
structured grid), **pave** (quad surface with holes/non-4-sided), **sweep** (pave source(s) ->
march through mapped linking surfaces -> target); plus sphere/polyhedron specials. Use hex for
nonlinear/dynamic accuracy & efficiency per DOF; tet for speed/ease. The prep loop iterates:
import -> heal -> defeature -> decompose (web cut) -> mesh -> (decompose more for quality).

### Beginner first mesh (Command Panel, clean part)
File>Import (e.g. .sat). Mode>Mesh, Entity>Volume. Action>Intervals -> Automatic sizing ->
select volume -> Apply Size. Action>Mesh -> Automatically calculate (scheme auto) -> Apply
Scheme and Mesh. Quality action -> pick metric (Shape) -> Apply (elements color-coded). Done.

### ITEM wizard (Immersive Topology Environment for Meshing)
Top-to-bottom guided pipeline in Power Tools; issues the same commands as manual work but
**auto-imprints and auto-merges** for you; non-modal (exit/re-enter freely). Steps:
1. Import/create geometry (prefer native CAD; STEP/IGES are dirtier).
2. Set up FEA model: pick hex/tet; drag size slider OR enter an **element budget** (e.g.
   10000) and Apply -> Cubit back-solves a target size.
3. Prepare geometry -> **Run Check Diagnostics** (ordered, recommended sequence): fix invalid
   topology (Auto Heal), remove small features, connect volumes (gaps/overlaps), build
   meshable topology. Green != "all removed" -- right-click a feature to **Mark as OK** (keep
   a physically relevant small feature). Remove small features previews a solution (often
   `remove surface <id> extend` -> a fillet becomes a sharp corner); right-click>Execute.
   Connect volumes: detect overlaps -> draw in red -> reduce/remove from the larger volume.
   Build meshable topology -> **Check Meshability** (meshable vs not) -> select a non-meshable
   volume -> **Decompose Volume** offers candidate web cuts; cycle with **up/down arrows**,
   preview the blue web, pick the EFFICIENT cut (one cut can split several features), Execute.
4. Mesh: type `all` -> Generate Mesh (failed-mesh troubleshooting is ITEM's weakest part).
5. Validate: define quality metric (scaled Jacobian for implicit/nonlinear; element size for
   explicit; Shape = combined heuristic) with min/max thresholds; color-code worst elements;
   refine/smooth/delete.
6. Boundary conditions: make blocks/sidesets/nodesets from CAD entities (Ctrl-click multi-
   select); meshing auto-pulls the corresponding faces/nodes into the set.
7. Export: native Genesis/Exodus (broadest downstream, incl. MOOSE) or Abaqus/Nastran/Patran/
   LS-DYNA/Fluent/OpenFOAM/...; option to force-overwrite.

### Blocks / sidesets / nodesets (Exodus model)
Block = material/element-type/section (element type e.g. HEX8 linear vs HEX20 quadratic set on
Properties page). Sideset = element faces with connectivity/orientation (pressures, integrated
BCs, contact). Nodeset = bag of nodes (nodal BCs, output tracking). Make sets from CAD
entities; many solvers treat sidesets/nodesets equivalently (some users use sidesets for all
BCs). Cubit's own material/BC support is basic -> name the sets, apply materials/BCs in the
solver deck.

### ITEM vs Command Panel vs journals/Python
ITEM = new/intermittent users (guided, auto imprint/merge). Command Panel = full control
(schemes, sizing, bias, virtual geometry). Command line/journals = every GUI action emits a
command (History tab) -> copy into a .jou (Tools>Journal Editor), clean (drop `preview`/undo
lines), replay. Python: `#!python` shebang in the command line, or external `import cubit`;
prefer `cubit.cmd("...")` strings; Journal Editor has a journal->Python button. (Deeper:
cubit_docs "coreform_python_automation".)

### Beginner decision points
Auto-mesh when diagnostics show no issues / Check Meshability passes. Ask ITEM to decompose
when a volume is "not meshable / no scheme set". Defeature small irrelevant fillets/curves
(they force tiny elements / kill the explicit time step); Mark-as-OK when a small feature
matters. Heal non-native imports. Imprint&merge for conforming meshes (skip merge only for
intentional contact). Hex refinement PROPAGATES through merged neighbors (structured), and
Cubit does conforming refinement (no hanging nodes); tets refine locally more easily.

### Sources
Your first 15 minutes; Getting the most out of Coreform Cubit; Introduction: Hex meshing for
beginners with the ITEM wizard; Introduction to Coreform Cubit; Coreform Cubit Basics; 2025
update: Your first 15 minutes; 2025 update shorts (first/second mesh, ITEM wizard, Command Panel).
"""

ADVANCED_MESHING = r"""
## Boundary layers, hybrid meshes, fluid regions, and large assemblies

### Boundary-layer (BL) meshing
A BL = thin, highly-resolved, graded layer hugging walls (CFD velocity, but also thermal/
stress) growing out into a coarse core. Cubit "extrudes" wall entities (surfaces in 3D)
inward; where BLs meet it resolves a **boundary-layer intersection** by internal angle (deg):
**End** (~right-angle convex), **Corner** (inverse), **Side** (share a side; use for near-flat
AND sharp internal corners), **Reversal** (trailing-edge / C-grid / O-grid recirculation).
**internal continuity** is ON by default -> forces all intersections to Side and refuses other
types; turn OFF to set per-junction types. For many-surface CFD, Cubit often can't auto-assign
intersection types well -> assign BLs surface-by-surface (labor-intensive).
Parameters: first-row depth (CFD: your y+-derived first height; Cubit doesn't compute y+),
growth/bias (e.g. 1.2), number of layers (e.g. 10). **Assign BL to wall SURFACES, not curves**
(curves are tedious and drift); BL grows normal to surfaces, lateral size inherited from the
volume/surface mesh size. Selection: "all surfaces in volume 1", or "all except surfaces 1 to
6" (exclude an outer box for external flow). BL previews yellow (selected->orange).
**Single-layer BL = structured collar** around bolt/stress holes (turn internal continuity
off). Limitation: only linear (corner) nodes follow curvature; mid-side interior nodes are NOT
projected onto curved walls -> smooth, or limit edge length on high-curvature walls.
**Sizing order: volume -> surface -> curve -> vertex** (a volume size OVERRIDES a prior
surface size). Inline counts: right-click>measure a distance, then APREPRO distance/size.

### Hybrid (hex-pyramid-tet) meshes
Mesh wall SURFACES first (quads by default), assign the volume a **tetmesh** scheme -> the tet
mesher respects the BL: hex layers at the wall, **a single pyramid layer** bridging hex quad
faces to tet triangles, tets in the core. Use when you want hex on boundaries (contact,
surface stress, BL) but an all-hex interior is impractical. Element conversion: `thex`
(tet->hex), `htet` (hex->tet); pyramid->hex is unsupported/unverified. No automatic all-hex
unstructured fill with BLs (needs manual decomposition); Sculpt is the parallel hex option
(its BL definition lives outside Cubit, per video).

### Fluid region from solid CAD
Internal flow (pipe): clean wetted geometry (`remove surface <id> extend` to close a
non-watertight step), select one inner surface -> **select continuous** -> copy&transform the
inner wall, **create surface from bounding curves** to cap each opening, then **volume create
from bounding surfaces** (optional stitch) = watertight fluid volume. External flow: make a
bounding box (extend 100%), `subtract volume <body> from <box>` -> external void; grow BLs off
the body surfaces only (exclude box faces). Conjugate/coupled: **imprint and merge** fluid+
solid for a contiguous interface; web cut (sweep surface perpendicular) to make regions
sweepable; assign blocks per region (solid vs fluid).

### Large-assembly workflow (~4000-part, ~130M tets)
Pipeline: clean (heal precision artifacts; remove large overlaps; defeature irrelevant detail
-> if a part is just wrong, fix the CAD) -> imprint&merge -> blocks/sidesets -> mesh ->
evaluate/repair quality -> export. **Work in small bites**: group structurally-related parts,
imprint&merge BY GROUP (a bug localizes to "group 2" not "volume 3607"). Commands:
`group "vacuum_vessel" add volume in selection` (push GUI pick via `in selection`);
`imprint volume in vacuum_vessel`; `merge volume vacuum_vessel` (reports surface-pair count =
sanity check). **Tolerant imprint** for sliver gaps/overlaps that break standard imprint:
`imprint tolerant volume in <group>` (default tol ~0.05mm; override `tolerance <v>`; slower;
per-group tuning is the long pole). Undo a bad imprint: `regularize volume in <group>`
(reverts to clean CAD). Diagnose: hide a part (press **I**) to see neighbor imprints;
View>Power Tools volume-overlaps -> draw the sliver; validate with a quick gravity/vibration
run (part "falls through floor" = missed imprint). Selection vocab: groups; **select similar
volumes** (all copies of a repeated part; also builds big sidesets); `in` (containment);
`except` (exclusion); `with is_merged = false` (un-merged = exposed/air faces, for BCs/audit);
id ranges `1 to 6`. Quality at scale: `mesh metric` (per-volume tet quality -> find sliver
parts), then **smoothing** (relocate nodes, no remesh -- remeshing a giant model is expensive).
Default to **tet for the bulk** of huge assemblies (all-hex not worth the decomposition).
Script it (Python/C++) once part counts hit the hundreds; headless on HPC for high-memory
nodes; **HDF5-format Exodus** for huge meshes (plain Exodus can't hold them; reads into MOOSE
faster). 2nd-order tet solves are less stable -> tune sizing/quality to avoid negative Jacobians.

### Sources
How to Build Boundary Layers and Hybrid Meshes; How to create a hybrid mesh; How to create a
fluid region from CAD geometry; Useful commands when meshing a large assembly; Using "imprint
and merge" to mesh a massive 4000-part nuclear assembly; Demo: meshing a massive nuclear
assembly; How to mesh a 4000-part nuclear assembly.
"""

ML_AND_GUI = r"""
## Machine-learning features (2023.8) and the custom-GUI intro

All ML runs OFFLINE (training data ships with Cubit). Entry point: the **Geometry
(Diagnostics) power tool** -> Options -> **Load ML Models** (once per session) -> Done ->
Analyze. Results are expandable lists; **Show Solutions** previews candidate ops (double-click
to execute).

### ML-enabled defeaturing
Rule-based diagnostics flag geometry that LOOKS bad (small curves/surfaces, bad angles, blend
chains) against thresholds (manual or Auto Size). **ML instead PREDICTS the meshing outcome
without meshing** (trained on many tet-meshed CAD models). Under "Tet Mesh Poor Quality
Metrics" it predicts per-feature: **scaled Jacobian** (angle quality, default flag <0.2),
**(scaled) in-radius** (element size, <0.1), **deviation** (curvature following). Show
Solutions proposes ~a dozen CAD ops (Remove Surface, Composite, Collapse Angle, Blend/
Tangency, Rebuild Topology) and **predicts the resulting metric** for each. Batch: select
several flagged entities -> right-click **"Execute ML Recommended Solutions"**. Limits: **tet
only** (hex training data is hard to generate); Remove edits REAL geometry (irreversible
without Undo -> checkpoint to .sat/.cub); B-rep only (not mesh-based); in 2023.8 you can't add
your own tet-mesh-ML training data.

### ML part classification & reduction
Classification = shape recognition: predict a volume's category (bolt/gear/pin/spring/washer/
nut/insert/screw). Geometry power tool -> Load ML Models -> **Part Classification** -> Analyze
`all` -> lists each category with members ("48 bolts"). Right-click: Draw / Zoom / **Select
Similar Volumes**. **Reduction** = Show Solutions gives per-category recipes (Preview/Simplify
bolt = strip threads/fill holes; **Reduce Bolt (spider/wagon-wheel)** = replace solid with FE
edges + central rebar; Core; **Reduce Bolt and Fit** = align/fit to hole, web-cut, hex-mesh,
auto-assign blocks). Select many parts -> apply one recipe to all (e.g. 8 bolts -> 3 blocks
head/shank/plug with materials). **Custom training IS supported & persistent** (unlike tet-ML):
`classify volume <id>` (predicted category), `... confidence` (per-category 0-1),
`classify volume <id> <category>` (add a training example + rebuild),
`reclassify volume <id> <newcat>` (override / new category -> unload+reload models to see it),
`classify list` (U marks user-trained). Command/Python use: `draw volume with category bolt`,
`group "bolts" add volume with category bolt`; `cubit.get_ml_classification(...)` (API; verify
name). Algorithm: **random forest (scikit-learn, bundled)** over a fixed feature vector
(volume, bbox, genus, surface area, area/volume ratio...). Best on B-rep (STEP/SAT); STL works
if imported as mesh-based geometry. No semantic context (a flat annulus -> "washer"); correct
via confidence + user training.

### Custom GUI (intro)
2023.8 bundles **PySide6** (Qt for Python) -> easiest way to extend the GUI. **Custom toolbars**
save actions as buttons (a Cubit-command string, or a Python/PySide6 script); exportable/
shareable (Coreform GitHub examples = DAGMC, tire toolbars; MIT-licensed code, PySide6 is
LGPL). Python is the easy path; for proprietary logic write a **C++ component** (deeper:
cubit_docs "cpp_sdk_*"; Python-toolbar depth: the cubit_toolbar_guide tool). Steps: Tools>
Custom Toolbar Editor (hammer icon) -> Add toolbar (name+file) -> Add button (simple-command /
Python-script / command-panel button). PySide6-in-Cubit conventions (per video): start with
`#!python` shebang; `cubit` namespace is ALREADY imported (don't re-import); inside Cubit
`__name__ == "__coreform_cubit__"` (not "__main__"); Cubit's importer dislikes parenthesized
imports (use backslash continuations); build a QDialog (QLineEdit rows -> read text ->
try/except float-convert -> run cubit.cmd commands). A `cubit_util` helper parents dialogs to
the main window (find the QMainWindow named "Claro"). Import shared toolbars: right-click the
toolbar panel -> Import a packaged `.tar.gz`.

### Sources
Machine learning in Coreform Cubit 2023.8: Part classification and reduction; Machine
Learning-enabled Defeaturing in Coreform Cubit 2023.8; How to create a custom GUI (Intro).
"""

NEUTRONICS_FUSION = r"""
## Cubit for neutronics / fusion multiphysics (DAGMC + Exodus)

(Domain-demo corpus -- captions = narration only, on-screen syntax not captured; exact
commands marked "(per video; verify)".)

### One Cubit model, two consumers
A single Cubit geometry serves BOTH the CAD-based Monte-Carlo neutron transport (**DAGMC**
surface mesh, `.h5m`) and the FE solver (**Exodus** volume mesh -> MOOSE/Cardinal/Open FUSION).
Generate the DAGMC tri surface mesh, then tet-fill the SAME surfaces so the FE boundary
triangles coincide with the DAGMC facets -> conformal, mass-conserving coupling (Cardinal runs
OpenMC tallies + a MOOSE heat solve on the SAME Exodus file, no geometry remapping).

### DAGMC export
- **Integrated DAGMC export** (native, recent Cubit) meshes surfaces with the REAL tri-mesher
  -> facets are **watertight on generation** (no separate sealing). The legacy plugin used
  graphics facets (not topology-aware, not guaranteed watertight, needed a sealing step).
- Steps: import -> **`imprint body all` + `merge body all`** (CRITICAL for robust particle
  tracking across shared surfaces; MCNP import imprints automatically, OpenMC-adapter journals
  do NOT) -> assign materials/BCs -> TriMesh `surface all` -> File>Export>DAGMC (`.h5m`).
- Faceting controls (replace the old faceting-tolerance knob): **deviation angle** (max
  triangle-vs-surface angle; smaller=finer; 1 coarse, 0.5 dense, 10 for unimportant air
  boundaries) and **anisotropic ratio** (~100 default). K-effective fidelity tracks faceting
  tolerance -> curved boundaries need enough edges for VOLUME conservation (triangle aspect
  doesn't matter for transport, but DOES matter if the tris also seed a tet mesh -> aim ~unity
  there + "split over-constrained edges"). A post-hoc script reports max triangle->surface
  distance to recover the true faceting tolerance.
- **Metadata** DAGMC/OpenMC read: modern = mesh **blocks** for materials (+ named material
  objects) and **sidesets** for BCs named `boundary:vacuum|reflecting|transmission`; legacy =
  group names `mat:<id>/rho:<density>`. A "graveyard" outer shell (material `graveyard`)
  terminates particles in MCNP-style decks (OpenMC uses sideset BCs instead).
- Select many surfaces by predicate: `sideset 1 add surface all with y_coord < 1e-6`.

### CSG <-> CAD
- **MCNP import** (native; rename input to `.i`): evaluates CSG half-space booleans -> very high
  volume IDs, slow on lattices (disable graphics / import headless via Python -> save `.cub`).
- **`openmc_to_cad`** (openmc-cad-adapter) emits replayable Cubit journals; `--world-size`
  bounds the infinite CSG regions; pass cell IDs to convert a single cell. OpenMC can't export
  back to CAD. Units cm (Cubit unitless -> values import as-entered).

### Volume mesh for coupling
After the DAGMC tri mesh, **`mesh triangles`** (per video; verify) fills selected surface tris
into a **tet** mesh whose boundary == the DAGMC facets exactly -> export Exodus for MOOSE/
Cardinal. Conformal coupling = tet only (track-length estimators); curved/2nd-order transport
elements not yet supported.

### Open FUSION Toolkit topology metadata
Encode multivalued-potential jumps for thin-wall eddy/MHD with **web cuts** (create matching
mesh surfaces) + **nodesets** (one vertex per port-hole; need holes-1 jumps), passed through
Exodus. Sheet (thin-wall) models: copy an outer face (Copy/Transform, no transform) then delete
the parent volume. Toroidal: APREPRO-parameterized cross-section, revolve/clone around the
torus; export **Hex27** for curvature on coarse grids.

### Sources
Neutronics on exact CAD geometry (Cubit/DAGMC); Nuclear multiphysics workflow (Cubit/OpenMC/
MOOSE/Cardinal); Cubit + the Open FUSION Toolkit; OpenMC neutronics via DAGMC and Cubit; Live
Cubit meshing demonstration for MOOSE.
"""

DOMAIN_APPLICATIONS = r"""
## Cubit across domains: sculpt (organic), lattice, tire, geomechanics, student edition

(Domain-demo corpus -- captions = narration only; exact commands "(per video; verify)".)

### Sculpt / overlay hex meshing -- for organic / image-derived geometry
The headline technique for CT/MRI/isosurface/STL geometry where decompose+sweep is impractical
(bio tissue, mineral grains, Mars samples, evolving melt-pool/level-set fronts). Sculpt immerses
the geometry in a **Cartesian background grid**, sculpts boundary hexes to the surface, smooths
-> all-hex, no manual decomposition. **Caveat (decisive):** sharp edges are ROUNDED OFF (esp.
non-axis-aligned), and the BEST elements are interior, WORST near the surface -- the opposite of
what contact/large-deformation wants. So for manufactured parts (fillets/holes) prefer
decomposed hex or tet; sculpt for organic only.
- CLI-first (Greg's pattern): in the GUI choose "do not run sculpt" -> Cubit writes an input
  deck (STL + `.i` + `.diatom` assembly). `sculpt` exe is in `.../Coreform Cubit/bin/sculpt`;
  `sculpt -h` / `sculpt -h <option>` / `sculpt -i <deck>.i -j<N>` (parallel). **Sculpt settings
  are NOT saved in the .cub** -> the `.i` deck is the source of truth (`$` = comment).
- Parallel output = decomposed Exodus (`name.e.8.0`, ...); recombine with **`epu -auto <file>`**
  (ParaView can also read the pieces).
- BCs by mesh-spatial query (sculpt mesh is dis-associated from CAD):
  `nodeset 11 add node with y_coordinate < 0` (highlight first to confirm). Quality: sculpt
  prints min scaled Jacobian; conforming adaptive refinement via `adapt_type`/`adapt_levels`.
- Advanced: supply a custom "mother mesh" (fine in a hot-spot, coarse outside); the `.diatom`
  points to the STL sculpted INSIDE it -> mesh stays identical run-to-run except where geometry
  changes (ideal for evolving fronts).

### Evolving level-set / remesh loop (coupled fluid<->solid)
Per time step export the fluid morphology as a **mesh-based geometry** (read the tet mesh into
Cubit as geometry; NOTE: mesh-based geometry **cannot be webcut** -> use sculpt), delete the
low-quality interface tets, sculpt a fresh hex mesh, loop over hundreds of iterations via
APREPRO (`{a}` index substitution) / sed / Jinja.

### Coreform Lattice GC (3D-printing lattices, add-on)
"Geometry-Compliant" lattices that follow the geometry (not clipped on a grid). Build a
**U-spline** hex-swept mesh; each element's geometric map deforms a chosen **unit cell**
(triangle tessellation) -> supports strut AND TPMS/gyroid cells; U-spline continuity makes
inter-cell interfaces smooth (cells effectively pre-merged). Pipeline: CAD region -> webcut/
partition into sweepable regions -> **imprint+merge** -> composite -> **mesh by sweep** (source
surface first; extrude + redistribute nodes) -> set U-spline (degree 2, continuity 1) -> `build`
-> `build u_spline lattice` (unit cell jack/octet_truss/PVB/truncated_sphere or a tessellation
path). Custom cell: center at 0,0,0 + scale extents to +-0.5 -> `fold` -> export `.json`. Native
parallel slicer -> binary CLI (Materialise Magics); also VTK/OBJ (-> STL for FEA).

### Tire meshing (2D cross-section -> 3D revolve)
Hex for tires (near-incompressible rubber converges better in contact; fewer DOF; layered hex
captures plies/belts for wear). The Endurica plugin (MIT, github.com/coreform-LLC, PySide6)
automates: import line drawing -> make surfaces from lines (stitch with merge tol ~= 1/2 the
smallest curve length, since AutoCAD lines don't close) -> classify tread/belt/ply regions by
ray-firing -> blunt acute tangencies (split/add material/rejoin) -> skew-control cut lines ->
imprint+merge -> composite curves -> mesh (mark belt/ply/chafer surfaces **mapped** so rebar
sits on quad midlines) -> reflect half-section -> create rebar (plies 2 elements thick to find
the midline) -> clean bad triangles toward scaled-Jacobian ~0.3 -> export **Abaqus axisymmetric**.
Pave-and-sweep blocker = triangular linking surfaces (must be map/sub-map) -> composite the
front-face + the small side-surface curves to make it sweepable.

### Geomechanics / Irazu (FDEM)
Surface-mesh the fault/discontinuity network first (imprint+merge, tri scheme, smooth, shape
~0.6), then tet-mesh the volume with a scheme that **respects the existing fault triangles**
(embeds the discontinuities in the conforming tet mesh). Element budget: split into a fine inner
volume + a coarse outer volume with bias on boundary curves (~10x runtime saving). Groups ->
blocks (materials) + nodesets (fault nodes, roller/pin faces) -> export Abaqus `.inp` (Irazu
reads it). Mesh quality sets the explicit time step -> a few tiny elements cost runtime.

### Coreform Cubit Learn / Associate (free non-commercial edition)
The FULL Cubit tool suite (CAD import/cleanup, hex+tet meshing, U-splines, Python API) free for
non-commercial use; the **ONLY cap = 50,000 elements on EXPORT** (mesh/visualize freely).
Academic licenses (sales@coreform.com / resellers) lift the cap. Fits a free pipeline: Cubit
Learn (mesh) -> MOOSE/FEniCS/CalculiX/OpenFOAM (solve) + Dakota (DOE/optimization via Cubit's
Python black-box API) -> ParaView (view), Exodus as the interchange format.

### Cross-domain reusable bits
GUI->journal->automate is the universal habit (play whole or "play selected"); APREPRO for
in-Cubit loops, Python (`cubit.cmd`) for heavier logic / external-solver-driven meshing; batch:
`cubit -batch -nographics -input <jou> -working_directory $PWD`. Absolute-tolerance SCALING
trick: scale the model so the smallest feature ~ order 1, mesh, scale back on export (fixes many
curved-geometry quality issues). Concentric-shell idiom: `volume N copy scale 0.8` repeated.

### Sources
Geomechanics with Irazu; Meshing Mars (paleomagnetics); Automatic hex meshing for biological/
material sciences; Coreform Lattice GC (3D printing); Hex meshing an evolving level-set; How to
mesh a tire tread; Tire model building (Endurica); Coreform Cubit Associate (student pipeline).
"""

GETTING_STARTED_2025 = r"""
## Getting Started (2025 series): UI workflow, journaling, Python, Associate edition

This topic distills Coreform's "Getting started" YouTube playlist (2025.8-era UI refresh plus several
older Trelis/Cubit overview videos). It assumes the generic getting_started topic already covers the
GUI panes, the ITEM wizard concept, and Exodus block/sideset/nodeset basics, so it focuses on what is
distinctive to this playlist: the three parallel "first mesh" paths, journaling/Python scripting,
extended parsing/selection, pillowing, vertex types, smoothing, and the free Associate/Learn edition.
Coreform Cubit is the commercial Sandia-codeveloped hex/tet mesher; attribution to Coreform and YouTube
is fine and required.

### 2025.8 UI refresh and the Power Tools home of ITEM
The 2025.8 release reskinned the toolbar with scalable high-DPI icons (structure unchanged), and
simplified licensing: instead of a license file, you log in once with your Coreform account email and
password and the credentials are remembered for later launches (verify via Help > Product Activation).
The ITEM wizard (Immersive Topology Environment for Meshing, the little wizard-hat-with-wand icon) now
lives as a tab inside the Power Tools dock, alongside the Command Panel. New users are steered to Power
Tools / ITEM; intermediate users shift to the Command Panel; advanced users drive the command line and
then journal/Python. If a panel is missing (Power Tools, ITEM, Properties Page), re-enable it from the
View menu. ITEM presents a vertical top-to-bottom checklist whose steps mirror the real hex-meshing
loop: import/create geometry -> set up FEA model -> prepare geometry (heal invalid topology, remove
small features, connect volumes) -> mesh -> validate quality -> define boundary conditions (sets) ->
export. Each step exposes hyperlinks and "done" buttons; a red exclamation flags unresolved issues, and
you can mark deliberately-kept small features as OK (green check) so a step turns green without removing
them. "Set up FEA model" lets you preview a hex mesh, type an element budget (e.g. 10k or 100k) and have
Cubit back-solve a target mesh size, and toggle the default scheme between hex and tet.

### Three parallel paths to your first mesh
The playlist teaches the same simple.sat / knuckle.sat workflow three independent ways so users can pick
their comfort level. (1) ITEM-in-Power-Tools: walk the wizard checklist, click Run Check Diagnostics,
remove small features, Generate Mesh, Validate (scale Jacobian or the meta "shape" metric), and export.
(2) Command Panel: drive the same operations through mode->action->entity dropdowns -- e.g. mesh mode,
volume entity, the intervals/scheme/mesh icon trio (automatic sizing -> auto-calculate scheme -> mesh),
then a Quality panel. The Command Panel is "training wheels for the command line": every panel is just
a widget builder that emits an ordinary Cubit command, identical to what you would type by hand, so any
edit field accepts the same syntax (including extended parsing) you would type at the prompt. (3) Command
line: the fastest path for frequent users. Learn syntax two ways -- run a GUI command and read the
echoed command in the history, or type a command followed by "?" to list all argument options (e.g.
"brick ?"). Examples taught: brick x 3 y 2 z 1 (single width arg copies to all dims), cylinder height
10 radius 3, subtract 2 from volume 1, mesh volume 1, and quality volume 1 scale jacobian global draw
mesh. Cubit is unitless; the size number just drives intervals.

### Second mesh: decomposition
The second-mesh tutorials add the core hex idea that complex CAD must be web-cut (partitioned) into
shapes Cubit has a recipe for, then re-stitched with imprint+merge. Cubit's recognized hex-meshable
recipes are mapped meshes, swept meshes (mesh one or more source surfaces and sweep layer-by-layer to a
target), the polyhedron scheme (an arrangement of mapped sub-meshes for cube-with-pocket shapes), and a
seven-block sphere scheme. Web cuts taught include: web cut with a sheet extended from a surface, web
cut normal to a curve, and web cut by coordinate plane with an offset. After cutting, imprint copies
neighbor topology across an interface and merge collapses the duplicate collocated entities into one
shared entity -- without imprint+merge the meshes stay disconnected (no shared nodes), and stresses
would only transfer via tied contact. ITEM's "Build Meshable Topology" step (Decompose Volume / Check
Meshability) previews candidate web cuts as a blue "web" surface; arrow keys cycle alternatives, and
choosing a cut that severs several appendages at once is more efficient than the first suggestion. Note
that because hex schemes are structured, a refinement on one merged volume propagates through its
neighbors -- a key hex-vs-tet gotcha.

### Journaling and the Journal Editor
Cubit is fundamentally command-driven: the GUI is an application layered on top that sends commands and
reads data back. Every GUI/Command-Panel action is journaled. Ensure journaling is enabled under Tools >
Options > History (commands then echo to the command line and append to a default .jou journal file).
The built-in Journal Editor (file extension .jou) lets you record, clean, and replay sessions; many
users instead use an external editor (vi/emacs) on a second monitor. You can File > Import the History
tab (or the Command/Script tabs) into a new journal-editor buffer, then prune undo groups and dead
exploration, comment/uncomment lines, and save reusable journals. A common pattern is one journal that
sets up geometry/defeature/decompose and another that sets sizing/schemes/smoothing/BCs/export, since
journals are shareable, replayable, and parameterizable -- Coreform support often asks for a user's
journal, fixes it, and sends it back. In the history tab you can also right-click any line (or a
highlighted portion) and "play selected" to re-run just those commands. The Journal Editor has a one-
click toggle that converts a Cubit-command journal into a Python script (see next section). Tools >
Options > Custom Tools provides ten persistent custom toolbar buttons, each holding a mini-journal of
commands that runs on click and survives between sessions.

### Python scripting intro
Cubit is fully scriptable through a Python API. The Journal Editor's Python button rewraps a journal as
Python where the workhorse is cubit.cmd("...") -- it issues any Cubit command as a string, so f-strings,
numpy, loops, and external data can drive the script. You can run Python inside Cubit (toggle the script
tab between Cubit-command and Python modes) or externally: import sys; sys.path.append(<cubit bin>);
import cubit; cubit.init(...), then cubit.cmd(...). A script beginning with the "#!python" shebang line
is interpreted as Python rather than Cubit commands when played. Useful API helpers shown: cubit.get_last_id
(grab the id of an entity you just created via cmd), cubit.get_id_string(list) (turn an id list into a
parsable string for a command), cubit.parse_cubit_list("surface","in volume 1") (resolve a selector to
ids, combinable with f-strings), cubit.get_entity_name(type,id), and per-entity objects like
cubit.surface(id) exposing geometry queries (surface type plane/cylindrical, closest_point_trim, center
point, etc.). Practical demos: loop volumes sorted largest-to-smallest to assign size so the smallest
volume controls element size; color surfaces by type for a meshing pre-flight; build a mid-surface shell
by copying/uniting the top surface, meshing it, then moving each node to the midpoint between top and
bottom via closest_point_trim; and convergence/optimization loops that remesh, re-solve (e.g. the free
CalculiX/CCX Cubit plugin), read a probe node's Von Mises stress, and refine until an error tolerance is
met or a geometry parameter (hole radius vs allowable stress) converges. Coreform's guidance: prefer
cubit.cmd over direct C/C++-backed object mutators (the API is SWIG-generated; direct memory-mutating
methods exist but are less robust for creating/modifying entities).

### GUI overview tour and setup
The overview videos tour the layout: central graphics window (blue background; keyboard focus there
enables shortcuts), model tree (volumes, sheet bodies, and set assignments), properties page (live
info on selected entities, with some fields editable in place -- e.g. set a requested size, rename an
entity, or compare meshed vs CAD area), the command line/history, and the right-hand Power Tools and
Command Panel docks. All windows are dock windows (stackable, repositionable). First-run housekeeping:
File > Set Directory to move the working dir off the install folder; Tools > Options to set graphics
tolerance, toggle CAD vertices, recolor entities (a common preference is black curves/vertices to feel
CAD-like), choose breadcrumb-trail command panels (saves vertical space), and remap left/middle/right-
drag to rotate/pan/zoom. Help resources: Help > About lists support@ and forum.coreform.com; built-in
HTML docs (also online) with index/search; ~two dozen "Cubit tips"; and a graphics-window keyboard-
shortcut reference. Handy graphics shortcuts: X/Y/Z to slice and view internal mesh along an axis, J/K
to step the cut plane, Q to exit. A power-user habit emphasized throughout: watch the journaled command
each GUI action emits, and you will learn to type commands faster than hunting panels.

### Extended parsing (aprepro-style command expressions)
Cubit's command parser supports topological traversal with the "in" keyword and criteria-based filtering
with the "with" keyword; together they enable powerful one-line selections (documented under Environment
Control > Entity Selection and Filtering). Topological example: draw hex in face in surface 2, or draw
node in face in surface 2. Criteria example: subtract volume all with volume < 0.5 from volume 1 to punch
holes from many small plugs at once; correct the self-subtraction warning with ... with volume > 1.7
except volume 1 from volume 1. The same extended syntax works inside any Command Panel pick widget (the
panel just builds a normal command). Other functions/keywords shown: with x_coordinate, surface area
filters, num_parents (1/2/>2 parent classification of a curve -- mirrors the toolbar curve-valence tool
coloring orange=1, blue=2, white=>2), and is_merged. You can right-click any populated pick widget to
draw/locate/highlight what the expression resolves to before applying. Cubit's parser also handles
aprepro arithmetic expressions in numeric fields. Recommended homework: read the entity-specification
docs and experiment.

### Selection tools: X-ray, extended selection, nodeset creation
X-ray selection (toolbar icon with an X) lets a rubber-band/box drag pick entities hidden behind others
-- without it you only get visible faces; with it a box drag through a brick selects back faces and
interior hexes too, useful on dense assemblies. Extended Selection / Pick Extended (right-click context
menu after an initial selection) opens a dialog driven by Python filter scripts loaded from a folder:
pick all surfaces of a volume, adjacent volumes, curves of a surface, radial selection within a distance,
etc. You can drag entities between source and target columns to chain filters, multi-select to remove,
and copy selected ids/names to the clipboard to paste onto the command line or into a pick widget. Custom
filters are Python classes subclassing cubit_gui.SelectionFilter (class name must match the .py filename
for auto-instantiation), reimplementing display_name(), run_filter() (using self.source_types/
get_source_ids and cubit.get_entity_name to filter, e.g. substring match on names like "hex"/"nut"/
"socket"), and optionally get_ui_file() to attach a Qt .ui widget (built in Qt Creator/Designer, both
free) so the filter takes user input (e.g. a name line-edit referenced by object name via
get_line_edit_value). Reloading edited filters requires a Cubit restart. Extended entity specification
also builds sophisticated node/side sets: e.g. nodeset on surfaces "with not is_core_merged" captures
all exterior nodes; remove nodes "in surface 141" by topology; remove nodes "with x_coordinate <= -4.9";
combine criteria. Working at the geometry level (assign sets to CAD surfaces/volumes; meshed entities are
inherited automatically) is generally easier than picking mesh entities directly.

### Boolean operations
Three booleans, accessed via geometry > volume > boolean action (or command line): intersect (new body
from the shared overlap region), subtract (subtract-body id removed from body id; e.g. subtract 2 from
1), and unite/union (combine bodies into one; works even for non-touching bodies, yielding one body with
multiple volumes). Each has a "keep originals" option. Unite has an "include mesh" option that produces a
united meshed body from already-meshed inputs, but it requires the bodies to be merged first.

### Smoothing, vertex types, and free surfaces
Smoothing improves quality by moving nodes without changing connectivity, and must be applied
progressively curves -> surfaces -> volumes (smoothing a surface leaves its bounding-curve nodes fixed;
smoothing a volume leaves surface/curve nodes fixed). Algorithms shown: Laplacian (curves), mean-ratio
(surfaces, optimization-based), Winslow, and condition-number (volumes). Orthogonal smoothing is a newer
boundary smoother: "near orthogonal" makes boundary-incident edges nearly perpendicular while smoothing
the interior (works on arbitrary mapped surfaces, multiple at once, order may matter); "fully orthogonal"
requires geometry where a ray from a start curve hits the opposite curve exactly once with a spanning
node line -- i.e. a mapped mesh -- as in an elliptical annulus, where you can fix a center curve to
preserve its mesh. Vertex types steer sweeping/sub-mapping: Cubit classifies each vertex (relative to a
surface) by the number of adjacent quads as End=1, Side=2, Corner=3, Reversal=4, normally set
automatically from the corner angle. Display them with "draw surface <id> vertex type" (shows id + a
letter E/S/C/R). When a sweep fails (watch the output window for "could not assign vertex types"
warnings), set them manually via the submap scheme's Advanced fields per surface (and on the opposite
surface), treating a tricky shape like an L-block; copying the journaled vertex-type commands into a
journal is usually faster than the GUI. Using free surfaces to mesh a volume: imprint loose surfaces
placed atop a volume's face onto that face, delete the free surfaces, and the face is now split into
sub-surfaces you can size independently to get the desired mesh.

### Mesh pillowing (boundary-layer-like)
Pillowing inserts one or more layers of hexes along chosen surfaces when web cutting, sources/targets,
vertex types, and smoothing still leave poor quality -- typically where a mapped mesh is forced into a
circular region creating near-180-degree element angles. Access via Command Panel: mesh mode > surface
entity > refinement action > pillowing. Selecting one surface inserts a layer next to it (each top-layer
hex split in two); selecting two adjacent surfaces wraps a continuous layer around the corner. Mechanism:
a "shrink set" of hexes is shrunk away from its neighbors and a pillow layer is formed in the gap; an
optional "through surface" argument controls where the layer exits the volume (single/two-surface picks
auto-designate through surfaces -- the echoed command shows all hexes designated and all surfaces except
the kept one as through surfaces). For complex cases, assemble the shrink set as a group (e.g. all hex
adjacent to surface 31, plus all hex in volume 5) then run "pillow hex in <group> through surface 32",
and follow with progressive smoothing -- this pushes the bad 180-degree angle into the interior where
smoothing can recover quality (demos went from ~0.1 shape to ~0.48, and a marginal scale-Jacobian up to
~0.39).

### Localization, history, install, and the Associate/Learn free edition
Localization: the GUI (menus, labels, tooltips -- not command-window text) can be translated using Qt
Linguist on the provided english.ts file; save as cubit_<lang>_<country>.qm into the bin folder so Cubit
loads the .qm matching the system locale at launch (the video demos a Japanese build and a Polish
translation). Trelis history: some videos say "Trelis" because that was csimsoft's brand; Coreform
acquired csimsoft in 2019 and rebranded Trelis to Coreform Cubit -- Coreform Cubit 2020.2 simply follows
Trelis 17.1, with identical workflows. Download/install (2025): buy/register at coreform.com, follow the
account email's login link to set a password, download the latest build from your account portal, run
the setup-wizard launcher, then log in with your account email/password in the activation window to tie
the seat to your license. The free edition (Coreform Cubit Learn, also called Associate) gives full
access to the entire tool suite for hobbyists/students/researchers with one functional cap: mesh export
is limited to 50,000 elements (some videos quote ~50k; an older slide said 50k for non-commercial use),
with no 30-day clock. The separate 30-day free trial has no export limit but expires. Academics needing
more can get heavily discounted academic licenses (email sales@coreform.com to lift the export cap).
Where it fits a free student FEA pipeline: choose a solver (application-focused like MOOSE, CalculiX,
OpenFOAM, FEBio, code_aster vs extensible frameworks like deal.II, FreeFEM, MFEM), model-prep/optimization
tools (Octave/Python/Julia/spreadsheets, and Dakota for design-of-experiments/optimization via Cubit's
Python "black-box" interface), meshing (Coreform Cubit Learn), visualization (ParaView for ~99% of users,
plus GLVis/VisIt for specific solvers), and a file format -- the recommendation is Exodus/Genesis, the
open Sandia format Cubit was built around, part of the SEACAS toolkit with Python readers and converters.
A real BYU example chained Cubit Learn -> MOOSE (tensor mechanics) -> ParaView, all free. Cubit also has
0D-to-3D meshing, conforming hex/tet/pyramid transitions, and smooth U-spline elements for IGA (build/fit
a U-spline basis on a volume) feeding Coreform's IGA solver. A recurring complex-geometry tip: Cubit's
ACIS kernel uses a 1e-6 absolute tolerance, so for tiny or highly curved parts whose smallest feature
approaches that tolerance, scale the model up (e.g. by ~1000 so the smallest curve is order 1-3), do all
decomposition/meshing, then scale back on export.

### Sources
- Introduction to Coreform Cubit (advanced meshing overview)
- 2025 update: Your first 15 minutes with Coreform Cubit (webinar)
- 2025 update: Your first hex mesh with the ITEM wizard in Power Tools
- 2025 update: Your second mesh with the ITEM wizard (includes decomposition)
- Introduction to Coreform Cubit: hex meshing for beginners with the ITEM wizard (webinar)
- 2025 update: Your first simple mesh with the Command Panel
- 2025 update: Your second hex mesh using the Command Panel
- 2025 update: Your first hex mesh using the command line
- Journal files and journaling in Cubit
- Your first 15 minutes using Coreform Cubit (webinar)
- 2025 update: Creating a mesh in Coreform Cubit with Python
- Importing and exporting custom toolbars in Coreform Cubit
- Coreform Cubit overview: Graphical User Interface
- Extended parsing in Coreform Cubit
- Mesh pillowing in Coreform Cubit
- Extended selection tool using Python
- Cubit Boolean operations
- Cubit nodeset creation with extended entity selection
- Cubit orthogonal smoothing
- Localize Cubit for use in non-English-speaking countries
- Setting vertex types in Cubit
- Using free surfaces to generate a mesh on a volume
- Using the X-ray selection tool in Cubit
- Why do some videos show Trelis?
- Coreform Cubit Associate: a valuable part of a free, efficient student FEA software pipeline (webinar)
- Introduction to Python scripting in Coreform Cubit (webinar)
- 2025 update: How to download and install Coreform Cubit
"""

FLEX_IGA = r"""
## Coreform IGA / Flex: isogeometric analysis & U-splines

This topic distills the public Coreform webinar series on isogeometric analysis (IGA), the U-spline technology, and the Coreform IGA / Coreform Flex solver. Coreform Flex is a commercial product (a license is required to run it); the conceptual and workflow knowledge captured here is drawn from Coreform's own public talks. Where the talks make business claims (cost/time savings, speedups), they are attributed as Coreform's stated claims rather than independently verified fact.

### What IGA is and why it exists

Isogeometric analysis is, at its core, just the finite element method (FEM) with a different choice of basis. The motivation is a well-known pain point in simulation workflows, traced repeatedly in the talks to an early-2000s Sandia National Labs study of analyst time: across a typical simulation, the actual solve is a small fraction of the effort, while the dominant cost is non-value-added model preparation: importing CAD, de-featuring (removing fillets, chamfers, small holes, embossings), healing/decomposing geometry so it can be meshed, building a hexahedral mesh, iterating on it, and cleaning up element quality. The Sandia data attributed on the order of 70-80% of analyst time to mesh-related work, and Coreform's presenters argue that in advanced engineering organizations with sophisticated assemblies the figure now often exceeds 90%. Two further problems compound this: (1) the customer derives no engineering value from the mesh itself, and (2) the de-featured "analysis solid model" is no longer the part that will actually be manufactured, so the simulation runs on an approximation, both geometrically (faceted/simplified surfaces) and topologically.

IGA's promise is to drive the mesh-related, non-value-added steps toward zero by analyzing fully-featured CAD directly, using the SAME smooth spline basis that represents the geometry. The slogan throughout the series is getting the analyst "as close to CAD" as possible. Coreform was founded in 2014 (presenters variously cite 2014 and 2016) with this goal; its founders had been investigating IGA since the early 2000s, and the technology rests on roughly 15+ years of IGA research. Coreform positions itself as having shipped the first commercial native IGA solver (Coreform IGA, first released December 2020), later evolving into the Coreform Flex product.

### The mathematical idea

FEM has two ingredients: a function space (where you place little basis functions in space, normally defined by a mesh) and a weak/energy form describing equilibrium. You discretize the weak form against the function space to assemble stiffness matrices, residual/forcing vectors, and mass matrices; the entries of those matrices are determined entirely by the choice of functions. Traditional FEA uses low-order, C0-continuous Lagrange "hat" functions. IGA instead uses smooth, higher-order splines (B-splines, NURBS, Bezier, and Coreform's U-splines) as the shape functions. A key framing in the talks: every mesh you have ever built IS already a spline, because a spline is just a piecewise polynomial; IGA leans into that fact and unlocks the rich mathematics of spline theory. Coreform's presenters state the strong mathematical claim that smooth higher-order splines are the best (most efficient/robust) approximation space one can choose; e.g. a degree-5, C4-continuous function (continuity = degree-1) is far more powerful per degree of freedom than a C0 element.

Consequences of the smooth, higher-order basis, as explained in the talks:
- Accuracy per DOF: solution fields (e.g. displacement) converge at rate p+1 in mesh refinement; derived quantities (strain, stress) converge one order slower, at rate p. So higher p means displacements AND stresses converge much faster. In an L2 sense, a cubic basis gives ~order-4 convergence for displacement and ~order-3 for stress, vs order-2 / order-1 for linear FEA. This matters for fatigue, work-hardening, and any stress-driven analysis. Coreform repeatedly shows convergence plots (e.g. a C-frame) where cubic U-splines essentially converge in well under ~10,000 elements while linear hexes have not converged at a million; for the Kansas City flat-flex-cable problem they cite reaching converged answers with under ~11,000 elements vs traditional meshes of ~10 million (~1000x fewer).
- Field smoothness: stresses and contact pressures come out smooth across element boundaries instead of C0-kinked; the nuclear fuel-rod contact example contrasts smooth U-spline contact pressure against the spurious checkerboard pattern from low-order FEM.
- Efficient quadrature: Gauss quadrature over-integrates smooth splines and proliferates points with degree. Coreform uses "function-maxima" quadrature: place one quadrature point at the maximum of each global (patch-level) basis function, with weights computed by moment fitting (and folded quadrature for efficient weight computation). Each degree increase adds only one basis function and thus one quadrature point, so an 8th-order patch can be exactly integrated with ~18 points vs ~50 for Gauss. This supersedes their earlier Greville-quadrature work and reduces over-integration (a cause of locking).
- Conditioning: with high-order C0 elements the solution spectrum gets polluted by non-physical "admissible" modes that add high-frequency noise and hurt iterative solvers. Adding smoothness cleans this up, yielding better-conditioned systems. As degree increases, traditional methods' condition numbers grow rapidly while Coreform reports their condition numbers stay flat or even decrease, enabling iterative solvers and (in explicit dynamics) larger stable time steps. They claim problems that customers could not solve on HPC clusters with thousands of cores were solved on a laptop because of this conditioning.
- Higher-order PDEs / thin structures: smooth higher-order bases capture bending (a higher-order effect) well, so shell-like behavior can be resolved with even less than one element through the thickness (demonstrated on the Scordelis-Lo roof and on thin pipes/clamps).
- Bigger, denser-but-smaller matrices: higher-order elements give denser per-element stiffness blocks, but the global system is dramatically smaller (their example: 18x18 vs 81x81), so net efficiency improves.

Bezier extraction is the linchpin technology. It is a linear transformation (encoded per-element as local matrices) from a standard C0 Bezier/Bernstein basis to the spline (U-spline/NURBS/T-spline/B-spline) basis. It generalizes the traditional FEM assembly process, letting a conventional FEA code consume spline elements, and lets Coreform switch between a global (patch) view of the spline and a local (per-element) view depending on the operation. The presenters repeatedly advise newcomers that the single best key to understanding IGA is to read the open literature on Bezier extraction. Coreform exports Bezier-extracted U-splines in a BEXT file (the name derives from Bezier extraction; LS-Dyna consumes this as IGA parts, e.g. element formulation 201) and is developing a newer JSON-based extraction format and an Exodus-based format (with Sandia / Idaho National Labs for the open-source MOOSE code).

### U-splines: the key Coreform technology

U-splines ("unstructured splines") are Coreform's patent-pending spline construction. They provide a smooth spline basis over UNstructured, mixed-degree, multi-patch / arbitrary-topology meshes that classic tensor-product NURBS and B-splines cannot represent with smoothness. The limitation they overcome: B-splines are built via Cox-de-Boor recursion (uniform degree, univariate) and extended to 2D/3D only by tensor products, which are inherently structured (u/v/w principal directions) and confined to a single rectangular patch. On any unstructured or multi-patch layout, smoothness is lost (drops to C0) at every patch interface. U-splines instead take, as input, just (1) the mesh layout and (2) the desired local spline properties (per-region polynomial degree and inter-element smoothness), and produce as output a spline basis plus the control points that fit the geometry, with maximal smoothness wherever the topology allows.

What U-splines enable, per the talks:
- Smooth bases on unstructured, mixed-topology meshes (including triangle/pyramid interface elements), with non-uniform degree, non-uniform smoothness, and local refinement, all in one algorithm that works identically across 1D/2D/3D.
- True local adaptivity via hierarchical refinement: the U-spline space is built with several nested levels (coarse/medium/fine, and more), and refinement just turns basis functions on/off locally. Unlike T-splines (whose local refinement can propagate dependencies through hanging nodes) and NURBS (whose local changes propagate across the whole patch, proliferating DOFs), U-splines aim for genuine local additivity.
- "Super-smooth intersections" instead of true hanging nodes: T-junctions that look like hanging nodes are actually fully supported in the basis (no tie constraint needed) and retain smoothness.
- Mixed continuity: C0 can be inserted deliberately (e.g. at material interfaces or to model cracks/discontinuities), exactly as a C0 FEM basis would, while the rest of the model stays smooth. At extraordinary points the basis is currently only continuous (C0), not smooth; Coreform notes approaches to make these at least C1.
- NURBS compatibility: U-splines are a generalization of NURBS/B-splines/T-splines, can be rational, and can be exported losslessly to a NURBS layout / STEP (converted to NURBS patches on export). Coreform was the company behind T-splines (sold to Autodesk ~a decade prior to these talks); U-splines were explicitly motivated by lessons from T-splines and built to be suitable for BOTH CAD surfacing and analysis. The "B" in B-spline stood for "basis," and Coreform frames U-splines as finally fulfilling that original vision.

In FEM terminology, U-splines support full h- (subdivision), p- (degree), and k- (smoothness) refinement, locally. Control points are the corollary of FEM nodes (one basis function per control point).

### Model preparation in Coreform Cubit

Coreform Cubit is Coreform's preprocessor, co-developed with Sandia (it descends from Sandia Cubit and the former commercial "Trellis"). It is a state-of-the-art structured/unstructured hex mesher, also does tet and hex-to-tet (pyramid transition) meshing, ships ACIS-based geometry cleanup, reads STEP / SAT natively (plus optional translators for SolidWorks, NX, CATIA, etc., and STL/PLY/OBJ), and is fully scriptable through a Python API where every GUI action echoes a journaled command. A free, non-commercial "Cubit Learn" license is offered, fully featured but capped at 50,000 elements on export.

The general U-spline build workflow (the bolded verbs are Cubit commands): MESH -> SET -> BUILD -> FIT -> EXPORT.
1. MESH the geometry with hexahedra (volumes) or quadrilaterals (surfaces). For a body-fit U-spline this means traditional prep first: de-feature (e.g. `surface ... remove`, with "select similar surfaces"), decompose via web-cutting (`webcut`) to build structure, composite away sliver/tangent surfaces with virtual topology (applied last, since virtual ops reduce downstream robustness), then imprint+merge for a conforming mesh.
2. SET the U-spline properties on the geometry: `set u_spline volume all degree 2 continuity 1` (default is degree 2 / continuity 1; max continuity = degree-1).
3. BUILD the basis: `build u_spline volume all as <id>`.
4. FIT the U-spline to the CAD via Bezier projection: `fit u_spline <id>`. This step diverges from low-order FEM (whose elements interpolate the geometry, needing no fit); Bezier projection fits the spline's Bezier components as closely as possible to the input CAD.
5. EXPORT: `export u_spline <id> ...` to BEXT (for LS-Dyna; add the `dyna_cards` option to emit node/side-set .k files), to VTK (surface, for ParaView), to the JSON Bezier-extraction format, or to STEP (lossless, but it becomes NURBS, no longer a U-spline). Assign boundary-condition sets (blocks/side sets/node sets, with globally unique IDs for the BEXT format) BEFORE building the U-spline so the algorithm captures that topology.

Gotchas the webinars call out (in the 2021.11-era releases): U-splines were initially supported only on swept hex meshes (map/submap/tet-primitive/sphere/polyhedron unsupported, since added); officially tested for quadratic and cubic (degree 2-3), with higher orders buildable but untested; a practical limit around ~10,000 (up to ~80,000) hex elements pending performance work; per-element local degree variation not yet exposed (different degrees allowed only on disconnected volumes); T-junctions and extraordinary-point smoothness still in progress. Sharp edges/creases can be "filleted over" by a smooth fit; workarounds shown include (a) web-cutting along the crease (giving C0 there), (b) using a polyhedron-scheme quad whose unshared corner sits on the crease vertex, or (c) for surfaces the semi-automatic `build u_spline crease group`. Over-partitioning is discouraged in those releases because volume boundaries were only continuous, not smooth. Meshes must currently be created in Cubit (importing arbitrary external quad/hex meshes was not officially supported) because the U-spline must be fit to the CAD. Coreform Lattice GC is a related Cubit module that inserts unit cells (jack, octet truss, truncated sphere, PVB self-supporting, or custom) into a U-spline background mesh, inheriting smooth, geometry-conforming, gap-free lattices that can be sliced in parallel for additive manufacturing and represented in a model-based-enterprise without rendering hundreds of thousands of STL cells.

### The Coreform IGA / Flex solver and workflow

Coreform IGA (later Coreform Flex) is built from the ground up as a native IGA solver, not an FEA code retrofitted to IGA. Architecture and capabilities described across the series:
- Solvers built on PETSc (parallel direct and iterative/Krylov linear solves); they also investigated MFEM (HPC libraries from national labs / exascale efforts) to pair with function-maxima quadrature for rapid assembly. Early demos ran single-core direct solves (PETSc parallel support was still being ported to Windows); the trimming process is itself parallelized.
- Physics/analysis types shown: linear elastostatics (the most mature/"released" capability); nonlinear structural mechanics including large/finite deformation, isotropic (J2) plasticity with piecewise-linear hardening, hyperelasticity (Neo-Hookean, Mooney-Rivlin) with near-incompressible formulations (pressure stabilization), and contact (self-contact, surface-to-surface, and a "general/base contact" needing only a single mechanical-contact interaction with optional Coulomb friction). Constraint enforcement starts from an augmented-Lagrangian framework specialized to penalty / Nitsche / mixed schemes. Time integration: implicit statics (Newton-Raphson, adaptive step sizes), a nonlinear-static "continuation" method, implicit dynamics (generalized-alpha / implicit midpoint), and explicit dynamics (central difference). Also: tie constraints, element death, multiple parts/assemblies, isotropic heat conduction and modal analysis (beta), with thermomechanical/multiphysics, RBE2/RBE3-style couplings, shells, and topology-optimization integration on the roadmap. User subroutines (custom boundary conditions, loads, and eventually user materials, plus paths like a welding torch trajectory) are written in Julia, chosen for being dynamically typed and just-in-time compiled (no separate compiler per platform) while near-C/Fortran performance.
- Input and I/O: human-readable JSON5 input decks (text-editable, not binary), VTK output for ParaView, HDF5 (.h5) probe output. Coreform Flex runs in the browser (built on web technologies to support cloud, internal cloud, or air-gapped/node-locked compute) and also as a standalone client; the native session format is .CF, and Cubit can export a Flex model to .CF.
- Model tree and probe output: Cubit's CAE model tree gives one-to-one correspondence between tree entries and the solver input deck (description, version, material model, function definitions, loads). The probe output evaluates any field (stress, displacement, etc.) at an arbitrary CAD location, NOT tied to a node or integration point, independent of the mesh discretization; probes can sample a single point, a line of N points, or report a field extremum (e.g. max von Mises and the displacement at that location). In a linear-statics demo a probe gave ~-2.3 ksi vs a ParaView peak ~-2.5 ksi.

The Flex Representation Method (FRM) and Flex meshing. FRM is Coreform's core innovation, described as the marriage of U-splines with the finite cell method (an immersed/embedded-domain, cut-cell technique using high-order bases, hierarchical refinement, adaptive integration at boundaries, and weak enforcement of immersed BCs). The key move: relax the requirement that the mesh conform to the geometry. The geometry is captured in the U-spline basis; a CAD body is immersed in a higher-order smooth spline background grid (e.g. quadratic, C1, rectilinear) and a fully-automated, parallel "volumetric trimming" / "flex meshing" step trims the CAD out of the spline grid to produce a volumetric simulation model that exactly represents the original CAD interior, not just its bounding surfaces. Because U-splines are unstructured, FRM spans a CONTINUUM of options: at one end, traditional de-feature-and-body-fit hex meshing; at the other, fully-immersed bounding-box meshing with zero manual prep; and anywhere in between ("locally immersed" / "partially body-fit"), e.g. body-fit a critical region (the outer radius of a piston, the tread of a tire, symmetry faces) for efficiency/accuracy-per-DOF while immersing the hard-to-mesh features (shoulder fillets, embossings). The analyst thereby calibrates the speed/accuracy trade-off per component. The minimum CAD information FRM needs is an inside/outside test, so it accepts not only B-rep CAD but also STL/faceted data, scan/point-cloud data (with a triangulation only for BC/contact imposition), and emerging implicit/additive/generative CAD kernels, all without changing the solver. Trimming is robust to "dirty" CAD (small gaps, lost quadrature points are a practical non-issue) provided the geometry is topologically correct (gaps beyond the model tolerance are respected as real disconnections). Boundary conditions and loads are applied directly to the CAD object, never to node/element sets, so changing a BC region or a tread design within the immersed envelope requires NO remesh; trimmed cells use adaptive (octree + function-maxima) quadrature to integrate accurately up to the cut boundary, and elements are treated merely as scaffolding for the basis functions and quadrature.

Rapid simulation across design iterations (no per-iteration remesh). Because one immersed/locally-immersed mesh covers a whole design envelope, many geometric variants that fit the same envelope reuse the SAME mesh. Demonstrated on a wheel sizing-optimization (spoke-width changes that alter CAD topology and even an invalid regenerated geometry leave the analysis unaffected) and especially on the GE jet-engine bracket grand-challenge (GrabCAD, 2013, ~630 entries). All entries fill the same design envelope, so the prep is identical: select the bolt/load surfaces (using "select similar surfaces"), immerse, and a simple Python script over Cubit's API + Coreform IGA can verify all of them "set-and-forget" with no manual intervention, no per-model tet/hex remeshing. Those models ran ~200,000-600,000 DOF, ~15 min to ~1 hr per model on a single core with a direct solver. This directly answers the original GE finding that evaluating the entries was itself the bottleneck (and some STL files were too large to mesh).

Predicting onset of failure with IGA. A dedicated webinar ("Predicting the onset of failure with IGA," May 2021) worked an undergraduate-textbook 1.5-ton hydraulic-press C-frame in gray cast iron (brittle; max-principal-stress failure criterion). The hand calculation (idealized curved-beam bending) predicted ~127.2 ksi at the inner radius, exceeding the ~50 ksi ultimate, hence failure. Coreform IGA, immersing the fillet-laden CAD in a "horseshoe" partially-immersed domain (and even a tight bounding box), recovered ~127-130 ksi maximum principal stress at the inner radius, in close agreement with both the hand calc and a near-million-element traditional FEA run, but with roughly 2-3 orders of magnitude fewer DOF (cited as ~1000x fewer DOF / ~1000x faster for that accuracy). A ductile (304/316-style steel) variant with perfect plasticity showed FRM handling very large plastic deformation robustly on an immersed mesh. The talk's framing: the textbook problem is for verification/intuition; the real payoff is on complex assemblies (e.g. a helicopter rotor) where load paths and stresses are NOT obvious and hand calcs cannot be trusted.

The "next-generation FEA solver on four challenging problems" demo. This webinar exercised Flex on a spread of physics: (1) a linear-static, additively-manufactured fin-stock heat-exchanger chamber with internal angled fins and fillets, immersed in a rectilinear grid, refined just by typing a smaller mesh size; (2) a topology-optimized radar dish + support (organic surfaces) as a tie-constrained assembly for natural-frequency extraction (linear statics); (3) the flat-flex-cable bending / delamination problem (large deformation, where body-fitting would require projecting every trace through the whole geometry, badly distorting the mesh; immersing each component eliminates that); and (4) the direct-ink-write (DIW) pad live demo: an implicit-dynamic, nearly-incompressible (Neo-Hookean, Poisson up to ~0.499), self-contacting compression of an additively-manufactured silicone-thread engineered-foam pad. The DIW pad is the showcase of full automation: traditionally it needs ~thousands of sub-cells and ~100 web-cuts to hand-build a quality hex mesh, whereas Flex just "drops in a box" and immerses (with body-fit platens where meshing is trivial), then a Python script over the Flex + Cubit APIs sweeps thread diameters / spacings / pad-volume ratios (examples on Coreform's public GitHub, "coreform-llc/flex-python-examples"), with reaction forces probed and fit. A separate "introducing Coreform Flex" webinar adds a NAFEMS trunnion-supported-pipe linear-static benchmark (verified against the published multi-code solution), a bolted pipe-repair-clamp implicit problem (large deformation, friction contact, a nearly-incompressible gasket, metal plasticity, only the bolt threads removed), and a ball-drop explicit-dynamic test (nonlinear, contact, plasticity, fillets immersed). Coreform also referenced the Sandia fracture-challenge problem to show J2 plasticity + friction contact matching experimental reaction-force data.

### Positioning and the March 2025 advancements

Coreform Flex is presented as Coreform's flagship product: a single tool bundling the IGA preprocessor, the native solver, and post-processing, aimed at simulating fully-featured CAD with no required de-featuring, while remaining a true generalization of FEM (standard Galerkin, supports in principle any FEM capability: material nonlinearity, plasticity, contact, fracture, element death, multiphysics). The product family: Coreform Cubit (traditional meshing/geometry prep, competing with HyperMesh / ANSA), Coreform Flex (the IGA solver, competing with general FEA solvers such as Abaqus / LS-Dyna), and Coreform Suite (both together). Coreform claims compatibility/interoperability via Bezier extraction (consumable in LS-Dyna; an Abaqus/CAE integration was disclosed as in active development, letting existing FEA workflows incorporate IGA and connect trimmed IGA regions to traditional beams/shells/connectors "for free"). Coreform repeatedly notes it implements features only against concrete customer requests; primary markets cited are automotive and defense (also aerospace, nuclear/energy), where most problems are large nonlinear/dynamic assemblies with contact.

Stated benefits (attributed to Coreform): drastically reduced model-prep time (some customers' month-long model builds reportedly compressed to a day, two days, or a week); elimination of geometric error from de-featuring; superior accuracy per DOF and superior robustness for nonlinear/dynamic/large-deformation problems; better-conditioned linear systems and larger explicit time steps (for explicit automotive crash work, they claim a roughly order-of-magnitude larger critical time step versus linear FEA, from both coarser high-order meshes and U-spline spectral tailoring that keeps the critical time step roughly constant across polynomial degree); and unlocking design-space exploration / automated optimization because robust, user-intervention-free model building can be scripted (Flex importable as a Python module alongside, e.g., PyTorch).

The March 2025 webinar ("the latest advancements in Coreform Flex") reiterated that Coreform Flex is now available for testing and purchase as the first solver built natively on IGA, with linear elastostatics as the most robust released capability and full nonlinear structural mechanics (including contact) close to leaving beta. New/updated points from that session: a refreshed GUI with mesh-manager and procedure/breadcrumb model-tree wizards; "primitive" background-mesh types beyond the rectilinear grid (cylinder/sphere/annulus in progress, plus the ability to author any background mesh, body-fit or partial, in Cubit and import it); a published, CI-generated verification manual at docs.coreform.com (cantilever beam, large deformation, near-incompressible elasticity, plasticity, mechanical contact) with results regenerated from live CI runs; default maximally-smooth splines with selectable degree up to quartic C3 and selectable custom continuity (e.g. P2C0); element death (not yet activation) and basic element-failure with sub-element resolution down to a quadrature point; the active Abaqus/CAE integration; the worked door-lock and C-frame example models; and a candid estimate that hand-meshing the showcased complex cast part traditionally would take a first-time expert ~20 hours (and that some "dirty" parts shown could not be tet-meshed by any available mesher), versus an automated trim in seconds. Local adaptivity / a-posteriori error estimation and mesh adaptation remain on the roadmap (today, mesh convergence is studied by running several uniform mesh sizes, with per-part size control as the only built-in adaptivity).

### Sources

- 001_What is IGA?
- 002_Explaining IGA: a brief technical introduction to IGA, U-splines, and Flex IGA
- 003_Introducing Coreform IGA
- 004_Rapid simulation of multiple design iterations with Coreform IGA
- 005_Coreform IGA update webinar: model tree, probe output, linear statics
- 006_Introduction to Coreform IGA
- 007_U-splines for IGA model prep in Coreform Cubit
- 008_IGA model preparation in Coreform Cubit Webinar
- 009_Coreform IGA: Predicting the onset of failure with IGA
- 010_Clip from Tire Society presentation
- 011_The mathematical idea underlying isogeometric analysis (IGA)
- 012_Introducing Coreform Flex for rapid simulation of fully-featured CAD models
- 013_Why everyone is so excited about isogeometric analysis (IGA): a brief explanation
- 014_A demonstration of next-generation FEA solver Coreform Flex on four challenging problems
- 015 / 016_The latest advancements in Coreform Flex (March 2025)
"""

CLIPS_HIGHLIGHTS = r"""
## Clips: short-form highlights (geodynamics/planetary + value framing)

These are Coreform's short-form "Clips" (46 shorts). MOST are brief excerpts of
material already captured at length in other topics -- see those for depth:
getting_started (wizard, first hex mesh); advanced_meshing (hybrid hex/pyramid/tet
mesh, fluid region from bounding surfaces, 4000-part imprint&merge, large-assembly
selection, toolbar/PySide6 customization); solver_workflows + python_automation
(MOOSE meshing, parameter-sweep automation, command vs read-only API); neutronics_fusion
(Cubit as preprocessor for Cardinal/OpenMC, CSG import, DAGMC surface mesh, conformal
volume mesh); and the Flex/IGA topic ("why FEA is painful", FEA built on CAD math /
isogeometric analysis, ~50% project cost savings, single-source-of-truth geometry).
This topic only synthesizes the genuinely DISTINCT content below.

### Geodynamics / planetary-science workflows (new domain)

A planetary-science group uses Cubit (driven by Python + journal files in batch mode)
to mesh icy-moon crusts -- the worked example is Saturn's moon Enceladus.

- Base body: in Cubit, define a spherical shell by outer radius (Enceladus ~252 km)
  and a shell thickness, which fixes the inner radius. Export to Exodus (EXO).
- Crustal-thickness variation: a Python script loads the EXO node coordinates,
  perturbs the outer and inner surface node positions to impose lateral thickness
  variations, then writes the EXO back out. So Cubit makes the clean shell; Python
  edits the geometry as nodal data.
- Geodynamic faults (the "tiger stripes"): create a finite volume from latitude/
  longitude points on the outer surface, project it into the shell, then use Cubit's
  cut feature to merge the two geometries and form fault interfaces. The four main
  stripes plus small minor "splays" (folds at the stripe tips) are all built this way.
  The whole sequence runs unattended from a journal file in batch mode; the body is
  Tet-meshed with refinement concentrated around the fault interfaces.
- Quality repair: shells with strong lateral thickness variation risk overlapping
  tetrahedra / coincident nodes. Cubit is asked to find cells with negative Jacobians
  and untangle them, which is what makes a clean mesh achievable on these bodies.
- Strain-driven mesh refinement: run an initial sim to get strain over the outer
  surface, write that as a sizing variable into the EXO, re-import to Cubit, delete the
  old mesh and remesh -> finer resolution where strain is high, for accurate strain.
- Why Cubit here: the methodology is a long flowchart of chained steps; any slow stage
  bottlenecks the whole pipeline of repeated simulations. Cubit's batch-mode speed and
  scripting efficiency keep the geometry/mesh steps from being the bottleneck.

### Value framing and a concrete tip

Recurring reasons an engineer cites for Cubit (positioning, distilled, no fluff):
extensive online documentation and deliberate extensibility; an active forum with
developer participation; a strong Python scripting + automatic-journaling workflow,
where anything the GUI does is also reachable from Python ("first-class" command
access, nothing hidden) and verbose output makes debugging/sharing easy; GUI
customization so a recurring workflow becomes a personal toolbar; and auto-journaling
that reconstructs "how did I do that" and converts straight into automation. The
recommended hand-meshing order is import -> heal -> web cut -> composite -> imprint &
merge -> mesh, ideally walked via the item wizard for discipline; similarity selection
speeds decomposing arrayed/assembly parts, and order of operations matters (wrong order
forces backtracking or bad meshes).

Concrete tip not in other topics -- large-assembly selection qualifiers: build named
groups (add volumes/surfaces/nodes/vertices to a group); "select similar volumes" for
repeated parts; and compose selections with "except", "in", and "with". Examples:
select every volume in group 5 except those in block 2 (exclude a different material);
"select surface in volume 4" for that volume's faces; and "... with is_merged = false"
to grab only the un-imprinted/-merged (e.g. air-exposed) surfaces. If an imprint goes
wrong, "regularize" reverts a volume back to its un-imprinted CAD.

### Sources

- Using python scripting to generate geodynamic faults with Coreform Cubit
- Generating planetary satellite crustal thickness variations with Coreform Cubit
- Using Coreform Cubit to identify malformed elements
- A geodynamics workflow relies on the efficiency of Coreform Cubit
- Mesh refinement with Coreform Cubit in a geodynamics workflow
- An Engineer's Top Reasons for Using Coreform Cubit
- An engineer's Coreform Cubit workflow and tips
- Useful commands when meshing a large assembly with Coreform Cubit
"""

_TOPICS = {
    "python_automation": PYTHON_AUTOMATION,
    "meshing_strategy": MESHING_STRATEGY,
    "solver_workflows": SOLVER_WORKFLOWS,
    "getting_started": GETTING_STARTED,
    "advanced_meshing": ADVANCED_MESHING,
    "ml_and_gui": ML_AND_GUI,
    "neutronics_fusion": NEUTRONICS_FUSION,
    "domain_applications": DOMAIN_APPLICATIONS,
    "getting_started_2025": GETTING_STARTED_2025,
    "flex_iga": FLEX_IGA,
    "clips_highlights": CLIPS_HIGHLIGHTS,
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
    "getting_started": "getting_started", "intro": "getting_started",
    "beginner": "getting_started", "item": "getting_started",
    "item_wizard": "getting_started", "command_panel": "getting_started",
    "first_mesh": "getting_started", "gui": "getting_started",
    "advanced_meshing": "advanced_meshing", "boundary_layer": "advanced_meshing",
    "boundary_layers": "advanced_meshing", "hybrid": "advanced_meshing",
    "fluid": "advanced_meshing", "fluid_region": "advanced_meshing",
    "assembly": "advanced_meshing", "large_assembly": "advanced_meshing",
    "imprint": "advanced_meshing", "tolerant_imprint": "advanced_meshing",
    "ml": "ml_and_gui", "machine_learning": "ml_and_gui",
    "classification": "ml_and_gui", "defeaturing": "ml_and_gui",
    "defeature": "ml_and_gui", "custom_gui": "ml_and_gui",
    "toolbar": "ml_and_gui", "pyside6": "ml_and_gui",
    "neutronics_fusion": "neutronics_fusion", "neutronics": "neutronics_fusion",
    "dagmc": "neutronics_fusion", "openmc": "neutronics_fusion",
    "mcnp": "neutronics_fusion", "cardinal": "neutronics_fusion",
    "fusion": "neutronics_fusion", "nuclear": "neutronics_fusion",
    "domain_applications": "domain_applications", "domain": "domain_applications",
    "lattice": "domain_applications", "tire": "domain_applications",
    "geomechanics": "domain_applications", "irazu": "domain_applications",
    "level_set": "domain_applications", "bio": "domain_applications",
    "associate": "domain_applications", "learn": "domain_applications",
    "student": "domain_applications", "3d_printing": "domain_applications",
    "getting_started_2025": "getting_started_2025", "gs2025": "getting_started_2025",
    "power_tools": "getting_started_2025", "powertools": "getting_started_2025",
    "2025_update": "getting_started_2025", "journaling": "getting_started_2025",
    "journal": "getting_started_2025", "journal_editor": "getting_started_2025",
    "pillowing": "getting_started_2025", "pillow": "getting_started_2025",
    "vertex_type": "getting_started_2025", "vertex_types": "getting_started_2025",
    "smoothing": "getting_started_2025", "orthogonal_smoothing": "getting_started_2025",
    "xray": "getting_started_2025", "x_ray": "getting_started_2025",
    "extended_selection": "getting_started_2025", "extended_parsing": "getting_started_2025",
    "free_surfaces": "getting_started_2025", "localization": "getting_started_2025",
    "trelis": "getting_started_2025", "install": "getting_started_2025",
    "download": "getting_started_2025", "cubit_learn": "getting_started_2025",
    "associate_edition": "getting_started_2025",
    "flex_iga": "flex_iga", "flex": "flex_iga", "iga": "flex_iga",
    "isogeometric": "flex_iga", "isogeometric_analysis": "flex_iga",
    "u_spline": "flex_iga", "u_splines": "flex_iga", "usplines": "flex_iga",
    "uspline": "flex_iga", "splines": "flex_iga", "spline": "flex_iga",
    "coreform_iga": "flex_iga", "bezier": "flex_iga", "bezier_extraction": "flex_iga",
    "frm": "flex_iga", "nurbs": "flex_iga", "bext": "flex_iga",
    "clips_highlights": "clips_highlights", "clips": "clips_highlights",
    "clip": "clips_highlights", "shorts": "clips_highlights",
    "geodynamics": "clips_highlights", "geodynamic": "clips_highlights",
    "planetary": "clips_highlights", "enceladus": "clips_highlights",
    "geology": "clips_highlights",
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

Batch 2 topics:
  getting_started    - GUI layout, tet-vs-hex recipes, ITEM wizard pipeline,
                       Command Panel first-mesh, blocks/sidesets/nodesets,
                       beginner decision points.
  advanced_meshing   - boundary-layer tool (intersections, y+, internal
                       continuity), hybrid hex-pyramid-tet, fluid-region
                       extraction, 4000-part large-assembly imprint/merge + tolerant
                       imprint + group strategy + HDF5 Exodus.
  ml_and_gui         - ML defeaturing (predict tet quality), ML part
                       classification & reduction (random forest; bolt recipes;
                       custom training), custom-GUI intro (PySide6 toolbars).

Aliases: python/automation/dakota -> python_automation; hex/decomposition/sculpt -> meshing_strategy;
moose/calculix/openfoam/exodus -> solver_workflows; item/beginner/gui -> getting_started;
boundary_layer/hybrid/fluid/assembly -> advanced_meshing; ml/classification/defeature/toolbar -> ml_and_gui.

Batch 3 topics:
  neutronics_fusion  - DAGMC watertight export + Exodus coupling (one Cubit model
                       -> OpenMC/MCNP transport + MOOSE/Cardinal FE), faceting
                       (deviation angle), block/sideset metadata, Open FUSION cuts.
  domain_applications- sculpt overlay hex (organic/CT/STL) + caveats, evolving
                       level-set remesh, Lattice GC (U-spline lattices for 3D
                       printing), tire 2D->3D cross-section, geomechanics/Irazu
                       FDEM, and Cubit Learn/Associate (free edition, 50k export cap).
Aliases: dagmc/openmc/mcnp/fusion -> neutronics_fusion; sculpt(organic)/lattice/tire/
irazu/learn/student -> domain_applications. ALL 54 Tutorials-playlist videos ingested
Additional playlists (3 more, ingested separately):
  getting_started_2025 - Getting-started playlist (2025.8 UI refresh): three first-mesh
                       paths (Power Tools/ITEM, Command Panel, command line), journaling +
                       Journal Editor, Python (cubit.cmd / extended selection), extended
                       parsing, pillowing, vertex types, orthogonal smoothing, X-ray select,
                       localization, Trelis history, Coreform Associate/Learn free edition.
  flex_iga           - Coreform IGA / Flex playlist: isogeometric analysis, U-splines (Bezier
                       extraction, BEXT), the FRM immersed/locally-immersed workflow,
                       MESH->SET->BUILD->FIT->EXPORT in Cubit, rapid design iteration (GE
                       bracket), C-frame failure prediction, function-maxima quadrature.
                       (Flex = licensed product; knowledge from public talks, claims attributed.)
  clips_highlights   - Clips playlist (short excerpts, mostly cross-ref other topics): the new
                       geodynamics/planetary domain (Enceladus crust, Python faults), value
                       framing, large-assembly selection qualifiers (except/in/with is_merged).
Aliases: power_tools/journaling/pillowing/associate_edition/trelis -> getting_started_2025;
iga/isogeometric/u_splines/flex/bezier/frm/nurbs -> flex_iga; clips/geodynamics/planetary -> clips_highlights.

(C++ SDK lives in cpp_sdk). "coreform_all" concatenates every topic.
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
                           ("python_automation", "meshing_strategy", "solver_workflows",
                            "getting_started", "advanced_meshing", "ml_and_gui",
                            "neutronics_fusion", "domain_applications",
                            "getting_started_2025", "flex_iga", "clips_highlights"))
    resolved = _ALIASES.get(topic, topic)
    if resolved in _TOPICS:
        return _TOPICS[resolved]
    return (f"Unknown coreform topic {topic!r}. Available: index, all, "
            f"{sorted(_TOPICS)}; aliases: {sorted(_ALIASES)}.")
