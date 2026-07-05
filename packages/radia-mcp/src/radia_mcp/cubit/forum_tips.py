"""
Cubit practical techniques sourced from the Coreform Cubit forum.

Covers mesh quality improvement, biased/graded meshing, geometry
decomposition for hex meshing, sweep tips, parallel meshing,
Python API patterns, geometry healing/import, boundary layer meshing,
Sculpt workflows, scripting performance, quality diagnostics,
solver integration workflows, and advanced meshing techniques
collected from real user problems and solutions.

Source: https://forum.coreform.com/c/cubit/9
"""

FORUM_MESH_QUALITY = """
# Mesh Quality Improvement Techniques

## Condition Number Smoothing (Fix Negative Jacobians)

After meshing, if elements have bad quality or negative Jacobians, apply
condition number smoothing:

```python
# Standard smoothing
cubit.cmd("volume all smooth scheme condition number beta 2.0 cpu 0.5")
cubit.cmd("smooth volume all")
```

| Parameter | Description |
|-----------|-------------|
| `beta` | Aggressiveness (2.0 is typical) |
| `cpu` | Max CPU time in seconds per volume |

This is **the most common fix** for bad elements after meshing complex
Boolean/webcut geometry with polyhedral or tetmesh schemes.

Source: forum.coreform.com/t/2636 (Bad Quality Shear Hex Generated)

## Quality Diagnostics

```python
# Check shape quality (0=bad, 1=ideal)
cubit.cmd("quality volume all shape")

# Check aspect ratio
cubit.cmd("quality volume all aspect ratio")

# Check Jacobian (negative = inverted element)
cubit.cmd("quality volume all jacobian")

# Check scaled Jacobian (normalized, -1 to 1)
cubit.cmd("quality volume all scaled jacobian")
```

## Smoothing After Hex Meshing

For hex meshes on complex decomposed geometry:
```python
cubit.cmd("volume all scheme auto")
cubit.cmd("mesh volume all")

# Post-mesh smoothing
cubit.cmd("volume all smooth scheme condition number beta 2.0 cpu 0.5")
cubit.cmd("smooth volume all")
```

## Node Constraint for 2nd Order Curved Elements

To force mid-edge nodes to conform to curved geometry:
```python
cubit.cmd("set node constraint on")
```

Without this, 2nd order mid-nodes form straight lines between corner nodes
rather than following the surface curvature.

Source: forum.coreform.com/t/2461 (Smooth curved 2nd-order elements)
"""

FORUM_BIASED_GRADED = """
# Biased and Graded Meshing Techniques

## Curve Bias (Size Gradation Along a Direction)

Create mesh that transitions from fine to coarse along a curve:

```python
# Bias from vertex 1 (fine) to vertex 2 (coarse)
cubit.cmd("curve 4 scheme bias fine size 0.005 coarse size 0.05 start vertex 1")

# Propagate the bias through the volume
cubit.cmd("propagate curve bias volume 1")
```

Source: forum.coreform.com/t/2700 (Cylindrical geometry refined near top)

## Circle Scheme (Radial Gradation)

Control radial mesh density on circular surfaces:

```python
# 16 angular divisions, inner core at 5% of radius
cubit.cmd("surface 2 scheme circle interval 16 fraction 0.05")
```

| Parameter | Description |
|-----------|-------------|
| `interval` | Number of elements around circumference |
| `fraction` | Fraction of radius for inner mesh core (0-1) |

Source: forum.coreform.com/t/2700

## Size Auto Factor (Automatic Size Control)

Let Cubit auto-determine element size with a scale factor:

```python
# Coarser (factor 10 = larger elements)
cubit.cmd("volume 1 size auto factor 10")

# Finer (factor 3 = smaller elements)
cubit.cmd("volume 2 size auto factor 3")
```

The factor is relative: higher = coarser, lower = finer. Useful when
different regions need different mesh densities:

```python
# Fine mesh near feature, coarse elsewhere
cubit.cmd("volume 5 6 7 8 size auto factor 5")    # Near feature
cubit.cmd("volume 1 2 3 4 size auto factor 10")   # Far from feature
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
```

Source: forum.coreform.com/t/2503 (Crack model with graded mesh)

## Surface Refinement by Depth

Refine near a specific surface, with decreasing density away from it:

```python
cubit.cmd("mesh volume all")
cubit.cmd("refine volume all using surface 7 depth 2")
```

The `depth` parameter controls how many layers of elements are refined.

Source: forum.coreform.com/t/395 (6822 views - Mesh Refinement in one direction)
"""

FORUM_HEX_DECOMPOSITION = """
# Geometry Decomposition for Hex Meshing

## Webcut + Imprint + Merge Pattern

The standard pattern for decomposing geometry into hex-meshable sub-volumes:

```python
# 1. Create or import geometry
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("create sphere radius 0.5")

# 2. Boolean operations
cubit.cmd("subtract volume 2 from volume 1 keep_tool")

# 3. Webcut along symmetry/feature planes
cubit.cmd("webcut volume all with plane xplane offset 0")
cubit.cmd("webcut volume all with plane yplane offset 0")
cubit.cmd("webcut volume all with plane zplane offset 0")

# 4. Imprint and merge (CRITICAL for conformal mesh)
cubit.cmd("imprint volume all")
cubit.cmd("merge volume all")

# 5. Mesh
cubit.cmd("volume all scheme auto")
cubit.cmd("mesh volume all")
```

Source: forum.coreform.com/t/2636, /t/2398

## Torus: Webcut to Enable Sweep

A torus cannot be sweep-meshed as-is (no cross-section). Cut it first:

```python
cubit.cmd("create torus major radius 1 minor radius 0.1")
cubit.cmd("webcut volume all with plane yplane offset 0")
cubit.cmd("imprint volume all")
cubit.cmd("merge volume all")
cubit.cmd("mesh volume all")
```

Source: forum.coreform.com/t/2632 (Meshing a torus)

## Brick with Embedded Sphere

Use webcuts + separate sizing per volume:

```python
cubit.cmd("brick x 1")
cubit.cmd("sphere radius 0.05")

# Webcut for symmetric decomposition
cubit.cmd("webcut volume all with plane xplane offset 0")
cubit.cmd("webcut volume all with plane yplane offset 0")
cubit.cmd("webcut volume all with plane zplane offset 0")

# Remove unwanted octants (for quarter/eighth symmetry)
cubit.cmd("delete volume with x_coord < 0")
cubit.cmd("delete volume with y_coord < 0")
cubit.cmd("delete volume with z_coord < 0")

# Boolean
cubit.cmd("subtract body 6 from body 5 keep_tool")
cubit.cmd("merge volume all")

# Different mesh density per region
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume 6 size auto factor 9")    # Sphere: fine
cubit.cmd("volume 5 size auto factor 10")   # Brick: coarse
cubit.cmd("mesh volume all")
```

Source: forum.coreform.com/t/2398 (Brick with Embedded Sphere)

## Angle-Based Compositing (Fillet/Sliver Cleanup)

Simplify geometry by compositing faces connected at small dihedral angles:

```
Workflow:
1. Select faces by dihedral angle (~1 degree max)
2. Keep anchor surfaces at boundaries
3. Composite the selection
4. Mesh coarse first, evaluate, then refine
```

For solid-mechanics workflows, higher-order tets can represent curvature
without needing the small CAD features resolved in the mesh.

Source: forum.coreform.com/t/2692 (Angle-based compositing)

## When Hex Meshing Fails: "Not mappable, submappable or sweepable"

```
WARNING: Volume 5 must have its meshing scheme explicitly specified;
it is not automatically mappable, submappable or sweepable.
```

**Solutions** (try in order):
1. Further webcut decomposition
2. Set `scheme tetmesh` for problematic volumes
3. Use `scheme polyhedron` (poly-hex) for complex volumes
4. Give up on hex, use `volume all scheme tetmesh`

```python
# Try polyhedron for complex volumes
cubit.cmd("volume 5 scheme polyhedron")
cubit.cmd("mesh volume 5")

# Post-mesh smoothing often needed
cubit.cmd("volume 5 smooth scheme condition number beta 2.0 cpu 0.5")
cubit.cmd("smooth volume 5")
```

Source: forum.coreform.com/t/2595
"""

FORUM_SWEEP_TIPS = """
# Sweep Meshing Tips

## Exact Sweep with Transform Translate

When sweep meshing deforms the target surface mesh (source and target don't
match), use `sweep transform translate`:

```python
cubit.cmd("volume 1 redistribute nodes off")
cubit.cmd("volume 1 scheme sweep source surface 1 target surface 2 "
          "sweep transform translate propagate bias")
cubit.cmd("volume 1 autosmooth target off")
cubit.cmd("mesh volume 1")
```

| Setting | Description |
|---------|-------------|
| `redistribute nodes off` | Don't redistribute nodes during sweep |
| `sweep transform translate` | Project source mesh exactly to target |
| `autosmooth target off` | Don't smooth the target surface |

Source: forum.coreform.com/t/469 (3067 views - Sweep source/target matching)

## Copy Mesh for Periodic/Cyclic Boundary Conditions

For meshes requiring identical faces (e.g., cyclic symmetry in OpenFOAM):

```python
# 1. Create matching topology by copying/imprinting
cubit.cmd("volume 1 copy rotate 60 about z")
cubit.cmd("volume 1 copy rotate -60 about z")
cubit.cmd("imprint volume all")
cubit.cmd("delete volume 2 3")

# 2. Set mesh size and mesh first surface
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size auto factor 5")
cubit.cmd("mesh surface 1")   # Mesh one side first

# 3. Copy mesh to opposite surface
cubit.cmd("copy mesh surface 1 onto surface 4")

# 4. Mesh the volume
cubit.cmd("mesh volume 1")
```

The key: imprinting ensures matching topology on both faces.

Source: forum.coreform.com/t/2507 (Copy surface mesh for cyclic condition)

## Copy Mesh for Kelvin Transform (Periodic Sphere BC)

For NGSolve Periodic Kelvin transform, two spheres need identical surface meshes.
ACIS spheres have no curves/vertices, so `copy mesh surface` fails without webcut.

```python
R = 0.060       # sphere radius
offset = 0.300  # exterior sphere offset

# 1. Create sphere + webcut to add curves at equator
cubit.cmd(f"create sphere radius {R}")
cubit.cmd("webcut volume 1 with plane zplane")
# -> volume 1 (upper half), volume 2 (lower half)
# Webcut creates a curve at the equator for mesh copy mapping

# 2. Copy halves to create exterior sphere (no mesh)
cubit.cmd(f"volume 1 2 copy move x {offset} nomesh")
# -> volumes 3, 4 (exterior copies with identical topology)

# 3. Mesh interior halves
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume 1 2 size 0.010")
cubit.cmd("mesh volume 1 2")

# 4. Copy mesh from interior to exterior hemisphere surfaces
# Must specify source/target curve + vertex for mapping
for src_vid, dst_vid in [(1, 3), (2, 4)]:
    src_sids = cubit.get_relatives("volume", src_vid, "surface")
    dst_sids = cubit.get_relatives("volume", dst_vid, "surface")
    # Hemisphere = largest area surface (not the flat equator cut)
    src_hemi = max(src_sids, key=lambda s: cubit.surface(s).area())
    dst_hemi = max(dst_sids, key=lambda s: cubit.surface(s).area())
    src_c = cubit.get_relatives("surface", src_hemi, "curve")[0]
    dst_c = cubit.get_relatives("surface", dst_hemi, "curve")[0]
    src_v = int(cubit.curve(src_c).vertices()[0].entity_name().split()[1])
    dst_v = int(cubit.curve(dst_c).vertices()[0].entity_name().split()[1])
    cubit.cmd(
        f"copy mesh surface {src_hemi} onto surface {dst_hemi} "
        f"source curve {src_c} source vertex {src_v} "
        f"target curve {dst_c} target vertex {dst_v}")

# 5. Mesh exterior volumes (boundary mesh already copied)
cubit.cmd("volume 3 4 size 0.015")
cubit.cmd("mesh volume 3 4")

# 6. Create surface blocks for NGSolve periodic identification
cubit.cmd("set duplicate block elements on")
# kelvin_int: hemisphere surfaces of interior
# kelvin_ext: hemisphere surfaces of exterior
```

Key points:
- **Webcut is required**: ACIS spheres have 0 curves; copy mesh needs curve/vertex
- **Volume copy nomesh**: ensures identical topology between interior/exterior
- **Node positions match exactly** (max dist ~1e-17) after copy mesh
- NGSolve `IdentifyPeriodicBoundaries("kelvin", "kelvin_int", Trafo(...), point_tolerance=1e-3)`
  then `Mesh(ngmesh)` re-wrap + `Periodic(HCurl)` constrains matching DOFs

Source: public-safe curated corpus (verified 2026-03-29)

## Interval Matching: Over-Constrained Systems

When sweep meshing fails with "Matching intervals UNSUCCESSFUL":

```
ERROR: Matching intervals <<< UNSUCCESSFUL >>>
```

**Cause**: Linking surfaces between already-meshed volumes are over-constrained.
Setting additional hard intervals on curves that are already determined by
adjacent meshes leads to conflicting constraints.

**Fix**: Don't set intervals on curves that belong to linking surfaces between
meshed and unmeshed volumes. Let Cubit resolve them automatically.

Source: forum.coreform.com/t/2709
"""

FORUM_PARALLEL_MESHING = """
# Parallel Mesh Generation

## Tetmesh: Built-in Multi-Threading

The tetmesh scheme defaults to multi-threaded (HPC) meshing:

```python
# Enable (default)
cubit.cmd("set tetmesher hpc on")

# Disable (single-threaded, for debugging)
cubit.cmd("set tetmesher hpc off")
```

Source: forum.coreform.com/t/2368, /t/2390

## Process-Based Parallelization (Python)

For meshing assemblies with independent parts, use Python's multiprocessing:

```python
import multiprocessing
from pathlib import Path

def mesh_part(cub5_file):
    \"\"\"Mesh a single .cub5 file in a subprocess.\"\"\"
    import cubit
    cubit.init(['cubit', '-nojournal', '-batch'])
    cubit.cmd(f'open "{cub5_file}"')
    cubit.cmd("volume all scheme tetmesh")
    cubit.cmd("volume all size auto factor 5")
    cubit.cmd("mesh volume all")
    output = str(cub5_file).replace('.cub5', '.exo')
    cubit.cmd(f'export mesh "{output}" overwrite')

if __name__ == '__main__':
    parts = list(Path("parts/").glob("*.cub5"))
    with multiprocessing.Pool(4) as pool:
        pool.map(mesh_part, parts)
```

Each process gets its own Cubit instance with a separate license checkout.
Works best when assembly parts are independent (no shared surfaces).

Source: forum.coreform.com/t/2391 (Parallel Execution of Coreform Cubit)
"""

FORUM_PYTHON_API = """
# Python API Patterns from Forum

## Get Volume/Surface Names

```python
vol_ids = cubit.parse_cubit_list('volume', 'all')
for vid in vol_ids:
    vol = cubit.volume(vid)
    name = vol.entity_name()
    print(f"Volume {vid}: '{name}'")
```

Also works for surfaces, curves, vertices:
```python
surf_ids = cubit.parse_cubit_list('surface', 'all')
for sid in surf_ids:
    name = cubit.get_entity_name('surface', sid)
    print(f"Surface {sid}: '{name}'")
```

Source: forum.coreform.com/t/2704

## get_expanded_connectivity for Solver Mapping

The `get_expanded_connectivity()` returns nodes in Exodus element order.
Different solvers use different node orderings:

```python
# Get Exodus-ordered connectivity
nodes = cubit.get_expanded_connectivity("hex", hex_id)

# Map from Exodus HEX20 to CalculiX C3D20 ordering:
# exodus_to_calculix = [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]
# mapped = [nodes[i] for i in exodus_to_calculix]
```

Source: forum.coreform.com/t/2446 (Numbering After Webcutting)

## Batch Mode with Python Script

```bash
# Run a .jou file in batch mode (use CUBIT_PATH or full path)
"%CUBIT_PATH%\\coreform_cubit.exe" -batch -nographics -nojournal workflow.jou
# Or: "C:\\Program Files\\Coreform Cubit 2025.12\\bin\\coreform_cubit.exe" -batch -nographics -nojournal workflow.jou

# Run a Python script directly
"%CUBIT_PATH%\\coreform_cubit.exe" -nojournal script.py

# Alternative: Use system Python with CUBIT_PATH for Cubit + NGSolve access
set CUBIT_PATH="C:\\Program Files\\Coreform Cubit 2025.12\\bin"
python script.py
```

Add `-nojournal` to prevent Cubit from creating cubitXXX.jou files.

Source: forum.coreform.com/t/2339

## STL/Faceted Geometry: Set Geometry Engine

When mixing STL imports with CSG operations (brick, cylinder, etc.):

```python
# Import STL files first
cubit.cmd('import stl "part1.stl" feature_angle 135.0 merge stitch')
cubit.cmd('import stl "part2.stl" feature_angle 135.0 merge stitch')

# Switch to facet engine BEFORE creating CSG primitives
cubit.cmd("set geometry engine facet")
cubit.cmd("create brick bounding box volume all extended percentage 50")

# Boolean operations now work
cubit.cmd("subtract volume 1 from volume 3 keep")
```

Without `set geometry engine facet`, you get: "Objects for Booleans must
first be classified" error.

Source: forum.coreform.com/t/2508 (Boolean operations with STL files)

## Element Numbering Consistency

After webcuts, element and node IDs may change. To get consistent numbering:

```python
# Compress IDs to remove gaps
cubit.cmd("compress ids")
```

For Exodus export, node/element numbering follows the order:
vertices → curves → surfaces → volumes.

Source: forum.coreform.com/t/2446, /t/2711

## Freeze Curve Mesh Before Surface Mesh

To preserve curve mesh during surface meshing (e.g., for periodic meshes):

```python
# Mesh curves first
cubit.cmd("curve 1 2 3 4 interval 10")
cubit.cmd("mesh curve 1 2 3 4")

# Freeze the curve mesh
cubit.cmd("set mesh freeze on")

# Now mesh the surface - curve mesh is preserved
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")
```

Source: forum.coreform.com/t/2100 (2847 views - Freeze curve mesh)
"""

FORUM_EXPORT_TIPS = """
# Export Tips from Forum

## Abaqus Export with Parts

```python
# Flat format (default)
cubit.cmd('export abaqus "mesh.inp" mesh_only dimension 3 overwrite everything')

# With separate parts per block
cubit.cmd('export abaqus "mesh.inp" mesh_only dimension 3 overwrite partial')
```

Source: forum.coreform.com/t/2622 (Experience with Abaqus mesh export)

## Exodus Sizing Function from Python

For adaptive mesh refinement, create sizing functions via Exodus:

```python
# 1. Export initial mesh as Exodus
cubit.cmd('export mesh "initial.exo" overwrite')

# 2. Add sizing data using netCDF4 (Python)
import netCDF4
ds = netCDF4.Dataset("initial.exo", "r+")
# Add time-based nodal variable for sizing
# See SEACAS documentation for Exodus format
ds.close()

# 3. Import sizing function back into Cubit
cubit.cmd('sizing function type exodus "initial.exo" variable "sizing" time 1')
```

Source: forum.coreform.com/t/1035 (2624 views - Exodus sizing function)

## LS-DYNA .k Export

Standard Cubit LS-DYNA export may not include all boundary condition sets.
Verify node sets and side sets are properly defined before export:

```python
cubit.cmd("nodeset 1 surface with x_coord 0")
cubit.cmd("nodeset 1 name 'fixed_nodes'")
cubit.cmd('export lsdyna "mesh.k" overwrite')
```

Source: forum.coreform.com/t/2386 (Exporting .k files for LS-DYNA)
"""


FORUM_GEOMETRY_HEALING = """
# Geometry Import, Healing, and Repair

## STEP Import Workflow

When importing STEP files, tolerances matter. STEP files specify
`UNCERTAINTY_MEASURE_WITH_UNIT` (e.g., 0.01 mm). Features smaller
than this tolerance may be ignored.

```python
# Basic STEP import
cubit.cmd('import step "model.step" heal')

# If surfaces are disconnected after import, use autoheal
cubit.cmd('import step "model.step"')
cubit.cmd("healer autoheal volume all rebuild")
```

Source: forum.coreform.com/t/1254 (2945 views - .STEP file Processing)

## STL/OBJ Import and Faceted Geometry

For faceted geometry (STL, OBJ), use `feature_angle` to control
surface segmentation:

```python
# Import STL with feature angle control
cubit.cmd('import stl "model.stl" feature_angle 135.0 merge stitch')

# For multiple STL files representing different regions
cubit.cmd('import stl "region1.stl" feature_angle 135.0 merge stitch')
cubit.cmd('import stl "region2.stl" feature_angle 135.0 merge stitch')

# Switch to facet geometry engine for CSG + STL combinations
cubit.cmd("set geometry engine facet")
cubit.cmd("create brick bounding box volume all extended percentage 50")
```

Source: forum.coreform.com/t/2508, /t/2613

## Bounding Imported CAD Geometry (Void/Cavity Detection)

When importing CAD with internal cavities, Cubit may not recognize
internal voids. Fix with healing:

```python
# Method 1: Simple heal
cubit.cmd('import step "model.step" heal')

# Method 2: Autoheal with rebuild (more thorough)
cubit.cmd('import step "model.step"')
cubit.cmd("healer autoheal volume all rebuild")

# Then proceed with bounding/subtraction
cubit.cmd("create brick bounding box volume 1 tight")
cubit.cmd("subtract body 1 from body 2 keep_tool")
```

Without healing, `subtract` may produce only one volume instead of
separating the void region.

Source: forum.coreform.com/t/1407 (3715 views - Bounding imported CAD)

## Stitch for Faceted Surface Consolidation

When importing ACIS files with many disconnected spline surfaces:

```python
# Stitch surfaces with tolerance
cubit.cmd("stitch body all tolerance 0.006")
```

Better than `composite create surface all` for performance.

Source: forum.coreform.com/t/1293 (3088 views - ACIS import surfaces)

## Imprint + Merge After Boolean Operations

**Critical step** for conformal meshing between adjacent volumes:

```python
cubit.cmd("imprint volume all")
cubit.cmd("merge volume all")
```

Without merge, adjacent volumes will have duplicate nodes at shared
interfaces, creating non-conformal meshes.

Source: multiple forum topics
"""

FORUM_SCULPT_WORKFLOW = """
# Sculpt (Overlay Grid) Meshing Workflow

## Recommended Sculpt Workflow

1. Import geometry (STL/STEP) into Coreform Cubit
2. Set working directory to a scratch folder
3. Create initial Sculpt input file from Cubit GUI
4. Run Sculpt from the command line (not from Cubit)
5. Recombine partitioned Exodus output

```python
# Step 1: Import and prepare geometry
cubit.cmd('import step "assembly.step" heal')

# Step 2: Set scratch directory
cubit.cmd('cd "/path/to/scratch"')

# Step 3: Export Sculpt input via Cubit
# (Or create manually as a .diatom file)

# Step 4: Run Sculpt from command line (outside Cubit)
# sculpt.exe --num_procs 8 -i input_file.diatom

# Step 5: Import result back
cubit.cmd('import mesh "sculpt_output.e" no_geom')
```

### Sculpt Command-Line Help

```bash
# Show all commands
sculpt.exe -h

# Get help on specific parameter
sculpt.exe -h <parameter>
```

Source: forum.coreform.com/t/1114 (3231 views - Sculpt Workflow Suggestion)

## Sculpt with "capture" Option (Negative Jacobians)

When Sculpt's `capture` option creates negative Jacobians, apply
post-mesh smoothing:

```python
cubit.cmd("volume all smooth scheme condition number beta 2.0 cpu 0.5")
cubit.cmd("smooth volume all")
```

Source: forum.coreform.com/t/2689

## Sculpt HEX27 Support

Sculpt supports HEX27 elements. After sculpt meshing:

```python
cubit.cmd("block all element type hex27")
```

Source: forum.coreform.com/t/2454 (700 views)
"""

FORUM_BOUNDARY_LAYER = """
# Boundary Layer Meshing

## Basic Boundary Layer Setup

```python
# Create boundary layer on surfaces for CFD
cubit.cmd("create boundary_layer 1")
cubit.cmd("boundary_layer 1 first_row_height 0.001")
cubit.cmd("boundary_layer 1 growth_factor 1.2")
cubit.cmd("boundary_layer 1 num_rows 10")
cubit.cmd("boundary_layer 1 add surface 3 4 5 volume 1")

# Mesh the volume
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
```

## Fix: "Reversal Nodes Found" Error

When boundary layer fails with "Reversal Nodes Found" and
"Boundary layer on surface is invalid", two root causes exist:

1. **Merged surfaces between volumes** create interface complications
2. **Multiple individual surfaces** in boundary layer definition

**Solution**: Delete merged volume early + use composite surface:

```python
# Step 1: Create geometry with Boolean merge
cubit.cmd("create brick 10 5 5")
cubit.cmd("create brick 10 5 5")
cubit.cmd("move vol 2 x 10")
cubit.cmd("merge vol all")

# Step 2: Delete the merged volume IMMEDIATELY (key step!)
cubit.cmd("delete vol 2")

# Step 3: Create composite from multiple BL surfaces
# Instead of referencing surfaces individually (causes reversal nodes),
# consolidate them into a single composite surface
cubit.cmd("composite create surface 4 5 10 11 12 15 16 17")

# Step 4: Set up boundary layer on the composite surface
cubit.cmd("create boundary_layer 1")
cubit.cmd("boundary_layer 1 first_row_height 0.001")
cubit.cmd("boundary_layer 1 growth_factor 1.2")
cubit.cmd("boundary_layer 1 num_rows 10")
# Reference ONLY the composite surface (e.g., surface 24)
cubit.cmd("modify boundary_layer 1 add surface 24 volume 1")

# Step 5: Mesh, then reflect and merge afterward
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("volume 1 copy reflect z")
cubit.cmd("merge vol all")
```

**Key insights**:
- `composite create surface` groups multiple surfaces into one entity,
  preventing BL conflicts at surface intersections
- Delete merged volumes BEFORE meshing, not after
- Reflect and merge AFTER meshing (not before)

Source: forum.coreform.com/t/2579 (Reversal Nodes Found)

## Boundary Layer + Sweep

When boundary layer on surface cannot be swept:

```python
# Mesh the boundary layer surface first
cubit.cmd("mesh surface 1")

# Then sweep the volume
cubit.cmd("volume 1 scheme sweep source surface 1 target surface 2")
cubit.cmd("mesh volume 1")
```

Source: forum.coreform.com/t/1532 (1606 views)
"""

FORUM_PYTHON_PERFORMANCE = """
# Python Scripting Performance Tips

## Performance Settings

```python
# Turn off messages
cubit.cmd("info off")
cubit.cmd("warning off")

# Turn off undo tracking
cubit.cmd("undo off")

# Pause graphics (fastest)
cubit.cmd("graphics off")
# ... do work ...
cubit.cmd("graphics on")

# Or: pause + flush (for GUI use)
cubit.cmd("graphics pause")
# ... do work ...
cubit.cmd("graphics flush")

# Reduce graphics faceting (if graphics needed)
cubit.cmd("graphics tolerance angle 15")

# Don't use transparent mode during computation
```

Source: forum.coreform.com/t/1949 (Tips for maximizing scripting performance)

## Precompute and Reuse Data

```python
import numpy as np
import time

# BAD: calling API in inner loop
node_ids = cubit.get_volume_nodes(1)
for nid in node_ids:
    coord = cubit.get_nodal_coordinates(nid)  # API call each time
    # do something with coord

# GOOD: batch fetch, then process
node_ids = cubit.get_volume_nodes(1)
all_coords = np.array([cubit.get_nodal_coordinates(nid) for nid in node_ids])
# Now process all_coords with NumPy vectorized operations
```

Source: forum.coreform.com/t/1949

## Capturing cubit.cmd() Output

`cubit.cmd()` returns only `True`/`False`. To capture console output,
use Python API methods directly instead:

```python
# This only returns True:
result = cubit.cmd("list vol 1")  # result = True

# Use API methods instead:
vol_ids = cubit.parse_cubit_list("volume", "all")
for vid in vol_ids:
    name = cubit.get_entity_name("volume", vid)
    center = cubit.get_center_point("volume", vid)
    print(f"Volume {vid}: '{name}' at {center}")
```

Source: forum.coreform.com/t/359 (5692 views - Capturing cubit.cmd output)

## Tracking Entities After Boolean Operations

After boolean operations, entity IDs change. Strategies:

```python
# Strategy 1: Use names (best)
cubit.cmd("volume 1 name 'matrix'")
cubit.cmd("volume 2 name 'inclusion'")
cubit.cmd("subtract volume 2 from volume 1 keep_tool")

# Find by name after boolean
vol_ids = cubit.parse_cubit_list("volume", "all")
for vid in vol_ids:
    name = cubit.get_entity_name("volume", vid)
    print(f"Volume {vid}: {name}")

# Strategy 2: Track last ID
last_id_before = cubit.get_last_id("volume")
cubit.cmd("subtract volume 2 from volume 1")
last_id_after = cubit.get_last_id("volume")
new_volumes = range(last_id_before + 1, last_id_after + 1)
```

Source: forum.coreform.com/t/1356 (1811 views - Tracking Entities)

## Vertex/Surface Coordinate Queries

```python
# Get vertex position
v = cubit.vertex(vertex_id)
coord = v.coordinates()
# or: coord = cubit.get_center_point("vertex", vertex_id)

# Get surface normal at a point
s = cubit.surface(surface_id)
coord = cubit.get_nodal_coordinates(node_id)
normal = s.normal_at(coord)

# Get entity center point (works for any entity type)
center = cubit.get_center_point("volume", vol_id)
center = cubit.get_center_point("surface", surf_id)
```

Source: forum.coreform.com/t/1380 (1743 views), /t/1284 (2136 views)
"""

FORUM_QUALITY_DIAGNOSTICS = """
# Quality Diagnostics and Element Selection

## Extract Quality Values for All Elements

```python
# Get all hex IDs in a volume
hex_ids = cubit.parse_cubit_list("hex", "all in volume 1")

# Query quality for each element
quality_values = []
for hid in hex_ids:
    q = cubit.get_quality_value("hex", hid, "scaled jacobian")
    quality_values.append(q)

# Find bad elements
bad_ids = [hex_ids[i] for i, q in enumerate(quality_values)
           if q < 0]
print(f"Elements with negative Jacobian: {bad_ids}")
```

Source: forum.coreform.com/t/2456, /t/1113

## Available Quality Metrics

| Metric | Description | Ideal |
|--------|-------------|-------|
| `scaled jacobian` | Normalized Jacobian (-1 to 1) | 1.0 |
| `shape` | Element shape quality (0 to 1) | 1.0 |
| `aspect ratio` | Edge length ratio (1 to inf) | 1.0 |
| `jacobian` | Raw Jacobian (negative = inverted) | > 0 |
| `condition` | Condition number (1 to inf) | 1.0 |

```python
# Visualize bad elements
cubit.cmd("quality volume all scaled jacobian high 0 draw mesh")

# Select elements by coordinate range
cubit.cmd("draw tet all with x_coord > -0.1 And x_coord < 0.1")

# Refine elements near a plane
cubit.cmd("refine tet all with x_coord > -0.1 And x_coord < 0.1 numsplit 1")
```

Source: forum.coreform.com/t/1113 (1436 views), /t/2046

## Global Element IDs vs Cubit IDs

```python
# Get the Cubit hex ID
hex_id = cubit.get_volume_hexes(1)

# Convert to global (Exodus/Abaqus) element ID
global_id = cubit.get_global_element_id("hex", hex_id[0])

# Compress IDs to remove gaps after webcut
cubit.cmd("compress ids")
```

Source: forum.coreform.com/t/1113, /t/2446

## Reverse Element Normals (Right-Hand Rule)

When a solver complains about left-handed elements:

```python
# Check normals
cubit.cmd("draw face all normal")

# Reverse normals on a surface
cubit.cmd("reverse face in surface 1 all")
```

Source: forum.coreform.com/t/1901 (Mirroring elements)
"""

FORUM_SOLVER_WORKFLOWS = """
# Solver Integration Workflows

## FEniCS / FEniCSx (dolfinx)

Export as Exodus, convert to XDMF using `mesh_converter`:

```python
# In Cubit: export Exodus mesh
cubit.cmd("block 1 add volume 1")
cubit.cmd("sideset 1 add surface 3 4")
cubit.cmd('export mesh "mesh.exo" overwrite')
```

```python
# Outside Cubit: convert Exodus → XDMF
# pip install mesh_converter  (or use jorgensd/mesh_converter)
# mesh_converter mesh.exo mesh.xdmf

# In dolfinx:
from dolfinx import io
from mpi4py import MPI
with io.XDMFFile(MPI.COMM_WORLD, "mesh.xdmf", "r") as f:
    msh = f.read_mesh(name="Mesh")
    block_data = f.read_meshtags(msh, name="Mesh")
```

Source: forum.coreform.com/t/2172 (2082 views - Workflow for FEniCS)

## VTK / VTU Export

Cubit does **not** have a built-in VTK writer. Options:

1. **MeshIO** conversion:
   ```python
   # Export as Exodus from Cubit, then convert
   cubit.cmd('export mesh "mesh.exo" overwrite')
   # meshio convert mesh.exo mesh.vtu
   ```

3. **ParaView** can read Exodus directly.

Source: forum.coreform.com/t/2025 (1658 views - VTK/VTU export)

## CalculiX (C3D20 Node Ordering)

Cubit uses Exodus node ordering. CalculiX uses different ordering
for HEX20 (C3D20). Map mid-face nodes:

```python
# Exodus HEX20 → CalculiX C3D20 reordering
# Nodes 12-15 and 16-19 are swapped
exodus_to_ccx = [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]

nodes = cubit.get_expanded_connectivity("hex", hex_id)
ccx_nodes = [nodes[i] for i in exodus_to_ccx]
```

Source: forum.coreform.com/t/2446, /t/2633

## OpenFOAM Export

```python
cubit.cmd('export openfoam "/path/to/case" overwrite')
```

**Limitations**: Materials/block data are not exported to OpenFOAM format.
Only boundary, faces, neighbour, owner, and points files are created.

Source: forum.coreform.com/t/379 (6064 views - OpenFOAM export)

## Group Names and Material Properties (Custom Export)

```python
# Get group names and IDs
names_ids = cubit.group_names_ids()
for name, gid in names_ids:
    print(f"Group '{name}' = ID {gid}")

# Get material from block
material_id = cubit.get_block_material(block_id)
material_name = cubit.get_material_name(material_id)
```

Source: forum.coreform.com/t/2096 (897 views)
"""

FORUM_ADVANCED_MESHING = """
# Advanced Meshing Techniques

## Volume with Internal Void (Cylinder Inside Brick)

```python
# Method 1: Decompose for hex meshing
cubit.cmd("brick x 10")
cubit.cmd("cylinder height 3 radius 2")
cubit.cmd("subtract body 2 from body 1 keep")
cubit.cmd("delete volume 1")
# Webcut to create sweepable volumes
cubit.cmd("webcut volume 3 with sheet extended from surface 7")
cubit.cmd("merge volume all")
cubit.cmd("volume all size auto factor 2")
cubit.cmd("mesh volume all")

# Method 2: Tet + THex (simpler but lower quality)
cubit.cmd("merge volume all")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
cubit.cmd("thex volume all")
```

Source: forum.coreform.com/t/428 (5176 views)

## Crack/Fracture Geometry (Embedded Surface)

Create a surface inside a volume for crack simulation:

```python
# Method: Non-manifold body via subtract
cubit.cmd("brick x 10 y 10 z 10")
cubit.cmd("surface 6 copy scale x 0.0001 y 0.5 z 0.5")
# Subtract sheet to create non-manifold body
cubit.cmd("subtract body 2 from body 1")
# The tet-mesher handles non-manifold correctly
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
```

The key: subtracting a zero-thickness sheet creates a non-manifold
body where the tet-mesher respects the internal surface.

Source: forum.coreform.com/t/312 (31087 views!), /t/1418 (1489 views)

## Hemisphere Hex Meshing

```python
cubit.cmd("sphere radius 1")
cubit.cmd("brick x 8")
# Quarter symmetry
cubit.cmd("webcut volume all with plane xplane offset 0")
cubit.cmd("delete volume 3 4")
cubit.cmd("subtract volume 1 from volume 2 keep_tool")
cubit.cmd("compress ids")
cubit.cmd("imprint volume all")
cubit.cmd("merge volume all")
# Set scheme
cubit.cmd("volume all scheme auto")
cubit.cmd("volume all size auto factor 5")
cubit.cmd("mesh volume all")
```

For hex meshing curved geometry: the decomposition-sweep strategy
is essential. Watch the Coreform hex meshing webinar for the
fundamental sweep-based approach.

Source: forum.coreform.com/t/2200 (1665 views)

## Dome Volume Meshing

For dome shapes created from swept cross-sections, decompose
the dome into multiple sweep-meshable sections. Key insight:
do not try to sweep through a singularity at the top. Instead:

1. Leave an opening at the dome apex
2. Create a separate cylinder/plug for the opening
3. Webcut the dome into half and mesh each section
4. Use imprint+merge for conformal interface

Source: forum.coreform.com/t/1044 (4055 views)

## Plate with Circular Hole (Quad Mesh)

For high-quality quad meshes around circular holes, use `polyhedron`
scheme instead of default `paving`:

```python
cubit.cmd("surface all scheme polyhedron")
cubit.cmd("surface all size 2")
cubit.cmd("mesh surface all")
```

Or decompose with O-grid for perfect radial meshes:

```python
# Create geometry
cubit.cmd("create surface rectangle width 100 height 50 zplane")
cubit.cmd("create surface circle radius 10 zplane")
cubit.cmd("subtract surface 2 from surface 1")

# For structured mesh, webcut + circle scheme:
cubit.cmd("surface all scheme polyhedron")
cubit.cmd("mesh surface all")
```

Source: forum.coreform.com/t/2428 (1505 views - Mesh Plate with Hole)

## Topography Surface from Structured Data

Create a CAD surface from (x, y, z) grid data:

```python
import numpy as np

# Create vertices from data
for row in data:
    for x, y, z in row:
        cubit.cmd(f"create vertex {x} {y} {z}")

# Create spline curves through rows of vertices
for row_vertices in vertex_rows:
    ids = " ".join(str(v) for v in row_vertices)
    cubit.cmd(f"create curve spline location vertex {ids} delete")

# Skin curves into surface
curve_ids = " ".join(str(c) for c in all_curves)
cubit.cmd(f"create surface skin curve {curve_ids}")

# Mesh the surface
cubit.cmd(f"surface {surf_id} size auto factor 5")
cubit.cmd(f"mesh surface {surf_id}")
```

Source: forum.coreform.com/t/2568 (741 views - Topography surface)

## External Faces Sideset (All Unmerged Surfaces)

```python
# All unmerged surfaces are external
sids = cubit.parse_cubit_list("surface", "all with is_merged=false")
face_ids = cubit.parse_cubit_list("face",
    f"in surface {' '.join(str(s) for s in sids)}")
cubit.cmd(f"sideset 1 add surface {' '.join(str(s) for s in sids)}")
cubit.cmd(f"nodeset 1 add surface {' '.join(str(s) for s in sids)}")
```

Source: forum.coreform.com/t/2581 (464 views)

## THex: Convert Tet Mesh to Hex

```python
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
cubit.cmd("thex volume all")

# Post-smoothing usually needed
cubit.cmd("volume all smooth scheme condition number beta 2.0 cpu 0.5")
cubit.cmd("smooth volume all")
```

**Note**: THex creates lower quality hexes than direct hex meshing.
Quality depends heavily on the input tet mesh quality.

Source: forum.coreform.com/t/997 (4075 views - Negative Jacobians)

## High-Order Element Conversion via Exodus Reimport

To convert imported meshes to high-order while preserving geometry:

```python
# 1. Export current mesh as Exodus
cubit.cmd('export mesh "mesh.exo" overwrite')

# 2. Reimport with geometry interpolation
cubit.cmd('import mesh geometry "mesh.exo" feature_angle 135')
# Use "linear", "spline", or "gradient" interpolation

# 3. Convert to high order
cubit.cmd("block all element type hex20")

# 4. Smooth to fix any inversions
cubit.cmd("volume all smooth scheme condition number beta 2.0 cpu 0.5")
cubit.cmd("smooth volume all")
```

Source: forum.coreform.com/t/863 (2143 views)
"""

FORUM_VENDOR_TIPS = """
# Vendor Tips (M&T / Coreform Support)

Knowledge from official Coreform support interactions and vendor guidance.

## Pressure Boundary Condition Bug

**CRITICAL**: Pressure boundary conditions on Curves do NOT export correctly
to Nastran format. This is a known bug (reported to Coreform developer as of 2021).

```python
# CORRECT: Pressure on Surface
cubit.cmd('create pressure name "load" on surface 1 magnitude 998')
# -> Exports as PLOAD4 card in Nastran

# WRONG: Pressure on Curve (bug - does not export correctly)
cubit.cmd('create pressure name "load" on curve 5 magnitude 998')
```

**Workaround**: Always apply pressure BCs to surfaces, not curves.

## Node Numbering in Nastran Export

Node ordering in Nastran export is NOT ascending by default:
- Elements are written in ascending order
- Each element's node IDs are assigned starting from the block's start_id
- Use `renumber` to get clean sequential numbering:

```
renumber node all in Volume all start_id 1 uniqueids
```

## Offset Webcut (Advanced Decomposition)

The `offset_tool` webcut variant cuts inward from surface boundaries.
Especially useful with the `loft` option for complex surfaces:

```
webcut Volume 7 offset_tool Surface 107 offset 1 inward merge
webcut Volume 20 offset_tool Surface 37 104 103 offset 1 loft inward merge
```

**Tip**: The `loft` keyword handles surfaces that are not simply extrudable.

Source: M&T vendor support (April 2024)
"""

FORUM_RENUMBERING = """
# Node and Element Renumbering

## Basic Renumbering

```
renumber node all start_id 300 uniqueids
renumber node all in Volume all start_id 1 uniqueids
```

## Compress (After Deletions)

After deleting entities, IDs may have gaps. Use `compress` to renumber:

```
compress       # renumber all entity IDs to be sequential
```

## Notes

- `uniqueids` ensures no duplicate IDs
- Element ordering follows block assignment order, not spatial order
- For clean exports, renumber before `export nastran`
"""

FORUM_MESH_COPY = """
# Mesh Copy Between Surfaces

Copy mesh from one surface to another (useful for conformal meshes
across interfaces, periodic boundaries, or symmetric models).

## Basic Mesh Copy

```
copy mesh surface 1 onto surface 4 source curve 2 source vertex 2 target curve 14 target vertex 14 smooth
```

**Parameters**:
- `source curve` / `source vertex`: reference edge and vertex on source surface
- `target curve` / `target vertex`: corresponding edge and vertex on target surface
- `smooth`: apply smoothing after copy (recommended)

## Multi-Surface Copy

```
copy mesh surface 3 onto surface 12 source curve 8 source vertex 1 target curve 10 target vertex 4
```

## Use Cases

- Ensure conformal mesh at material interfaces after imprint/merge
- Create periodic mesh pairs (same topology on both sides)
- Copy from a well-meshed surface to a similar surface

## Tips

- Source and target surfaces must have the same number of sides
- Curve/vertex pairs define the orientation mapping
- Use `list surface N` to check surface topology before copying
"""

FORUM_ENTITY_SELECTION = """
# Entity Selection and Filtering Tips

## "except" Keyword

Exclude specific entities from a selection:

```
surface all except 81 66 75          # all surfaces except 81, 66, 75
draw surface all except 5            # draw all except surface 5
```

## Spatial Filtering

```
draw surface with radius=2           # surfaces with specific radius
draw hex with x_min > 0              # elements in positive x half
draw tet with z_coord < 0.5          # tets below z=0.5
draw block 3 with x_min < -1500 color red
```

## Name-Based Selection (Wildcards)

```
mesh volume with name "coil_1*"      # match coil_1, coil_10, coil_1a, ...
mesh volume with name "coil_*"       # match all coil_ prefixes
body in surface with name "Tear_Fault_13@A"
```

## Python API Equivalent

```python
body_id = cubit.parse_cubit_list("body", "in surface with name 'name'")
sids = cubit.parse_cubit_list("surface", "all with is_merged=false")
```

## Groups for Temporary Collections

```python
tmpGroup = cubit.create_new_group()
cubit.cmd(f"group {tmpGroup} add curve all")
curves = cubit.get_group_curves(tmpGroup)
cubit.delete_group(tmpGroup)
```

## HARD Curve Detection

After manual interval assignments, curves become "HARD" locked:

```python
hard_curves = []
for c_id in cubit.parse_cubit_list("curve", "all"):
    if cubit.get_mesh_interval_firmness("curve", c_id) == "HARD":
        hard_curves.append(c_id)
cubit.cmd("draw curve " + " ".join(map(str, hard_curves)))
```
"""


def get_forum_tips(topic: str = "all") -> str:
	"""Return forum-sourced Cubit tips by topic."""
	topics = {
		"mesh_quality": FORUM_MESH_QUALITY,
		"biased_graded": FORUM_BIASED_GRADED,
		"hex_decomposition": FORUM_HEX_DECOMPOSITION,
		"sweep_tips": FORUM_SWEEP_TIPS,
		"parallel_meshing": FORUM_PARALLEL_MESHING,
		"python_api": FORUM_PYTHON_API,
		"export_tips": FORUM_EXPORT_TIPS,
		"geometry_healing": FORUM_GEOMETRY_HEALING,
		"sculpt_workflow": FORUM_SCULPT_WORKFLOW,
		"boundary_layer": FORUM_BOUNDARY_LAYER,
		"python_performance": FORUM_PYTHON_PERFORMANCE,
		"quality_diagnostics": FORUM_QUALITY_DIAGNOSTICS,
		"solver_workflows": FORUM_SOLVER_WORKFLOWS,
		"advanced_meshing": FORUM_ADVANCED_MESHING,
		"vendor_tips": FORUM_VENDOR_TIPS,
		"vendor": FORUM_VENDOR_TIPS,  # alias
		"renumbering": FORUM_RENUMBERING,
		"renumber": FORUM_RENUMBERING,  # alias
		"mesh_copy": FORUM_MESH_COPY,
		"entity_selection": FORUM_ENTITY_SELECTION,
		"selection": FORUM_ENTITY_SELECTION,  # alias
		"filtering": FORUM_ENTITY_SELECTION,  # alias
	}

	topic = topic.lower().strip()
	if topic == "all":
		return "\n\n".join(topics.values())
	elif topic in topics:
		return topics[topic]
	else:
		return (
			f"Unknown topic: '{topic}'. "
			f"Available: all, {', '.join(topics.keys())}"
		)
