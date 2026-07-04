"""
GMSH Post-Processing Specification for Radia Panels.

This is the SINGLE SOURCE OF TRUTH for what the GMSH output from
Radia panels must look like. Every calc_*.py that writes GMSH
output MUST satisfy ALL items in GMSH_POST_SPEC below. The
gmsh-verify skill checks a subset of these automatically; the
rest are verified by visual inspection against the reference
screenshots.

Call via MCP: ``gmsh_post_spec()``
"""

GMSH_POST_SPEC = """\
# GMSH Post-Processing Specification (Radia Panels)

ALL items are mandatory. No exceptions, no "future work".

## 1. File Format

- **v4.1 ONLY** (CLAUDE.md policy: "Radia-NGSolve panels output v4.1 only")
- Single .msh file (NO two-file Merge — element ID collision kills GMSH)
- Node IDs offset between physical groups (volume 1..N, surface N+1..M)
- Companion .geo file (Open GMSH button opens this, not the raw .msh)

## 2. Physical Groups (BEM post with air field)

Three groups in one file:

| Group | Dim | Content | NodeData |
|-------|-----|---------|----------|
| "air" | 3D (tets) | Volume mesh (Biot-Savart air box) | **B** vector |
| "coil_surface" | 2D (tris) | Coil surface mesh | **J** vector |
| "coil_wire" | 1D (lines) | Coil wireframe edges | (none) |

If workpiece coupling is on, add:

| Group | Dim | Content | NodeData |
|-------|-----|---------|----------|
| "workpiece_surface" | 2D (tris) | Workpiece surface mesh | **q** scalar (sigma*|J|^2/2) |

## 3. Mesh Curving

- Coil surface mesh MUST use **high-order elements** matching the
  .vol's curve order (mesh.Curve(N) where N = .vol's curvedelements
  order, typically 2).
- This means Tri6 (order 2), NOT Tri3 (order 1).
- The coil surface mesh is extracted from mesh_full (which has
  Curve(N) activated at load time), so the GmshPostExport high-order
  path (via _compute_ho_bnd_nodes or GetTrafo) must be exercised.
- Volume mesh (air box) stays order 1 (no geometry to curve).

## 4. Display Options (.geo launch + exact .opt sidecars)

Write a ``case.geo`` launch recipe next to the .msh file. GMSH auto-loads
``case.geo.opt`` when ``case.geo`` is opened. Also write ``case.msh.opt``
for raw mesh/data inspection when a user intentionally opens the .msh
directly. A plain ``case.opt`` is not the Explorer auto-load contract.

The panel's Open GMSH button should open ``case.geo``. If an older panel
passes ``case.msh``, the launcher should prefer the sibling ``case.geo``.

The .opt MUST contain ALL of the following:

```
Mesh.NumSubEdges = 4;       // render curved elements smoothly
Mesh.VolumeEdges = 0;       // hide air tet edges (black cube)
Mesh.VolumeFaces = 0;       // hide air tet faces
Mesh.SurfaceFaces = 1;      // show coil surface as smooth body

// B view
View[0].VectorType = 4;     // 3D arrow
View[0].ArrowSizeMin = 20;  // FIXED 20 px (not auto-scale!)
View[0].ArrowSizeMax = 20;

// J view
View[1].VectorType = 4;
View[1].ArrowSizeMin = 20;
View[1].ArrowSizeMax = 20;
```

**Mesh.Volumes does NOT exist** in GMSH 4.x. Never write it.

ArrowSizeMin = ArrowSizeMax = 20 is THE critical setting. Without
it the arrows are ~3 px and invisible. This was the "kirei" knob
from v3.6.1 (commit 5e0b69c).

The Open GMSH button handler should call ``gmsh.open(geo_path)``.
GMSH auto-discovers ``case.geo.opt`` and applies the options. Do not
rely on ``case.opt`` or Windows ``UserChoice`` associations.

## 5. NodeData Requirements

- B: 3-component vector on ALL volume nodes. Real part of
  BiotSavartCF for AC problems.
- J: 3-component vector on ALL coil surface nodes. Per-vertex
  averaged from per-element HDivSurface integration.
- q (workpiece): 1-component scalar on workpiece surface nodes.
  q = sigma * (|J_re|^2 + |J_im|^2) / 2 [W/m^2].
- NO scalar |B| or |J| views (redundant — GMSH vector view shows
  magnitude via arrow color).

## 6. What "kirei" Looks Like (v3.6.1 Reference)

The target visualization has these visual properties:
- Coil is a smooth shaded torus (Tri6 curved elements)
- J arrows sit on the coil surface, tangential, colored by magnitude
- B arrows float in the air around the coil, 3D arrows, colored by |B|
- Air mesh is INVISIBLE (no black edges, no tet faces)
- Coil wireframe visible as thin lines showing the mesh topology
- Workpiece (if present): smooth sphere with heating density color map

## 7. FEM Post (calc_fem_kelvin.py)

Same requirements except:
- Volume mesh is mesh_full (not a separate air box)
- B is curl(A) evaluated per-vertex via CoefficientFunction
- J is no longer exported (PEEC filament source, not a volume J)
- No separate coil wireframe (mesh_full already has the coil)
- Companion .geo hides Mesh.VolumeEdges + Mesh.VolumeFaces +
  Mesh.SurfaceEdges (keep SurfaceFaces for the material boundary)

## 8. Invariants (Never Break These)

- v4.1 format
- Single file (no Merge)
- ArrowSizeMin = ArrowSizeMax = 20
- No Mesh.Volumes in .geo/.opt
- No |B| or |J| scalar views
- Coil surface MUST be curved (Tri6+)
- B and J arrows MUST be visible at default zoom
- gmsh-verify must PASS before deploying
"""


def get_gmsh_post_spec():
    return GMSH_POST_SPEC
