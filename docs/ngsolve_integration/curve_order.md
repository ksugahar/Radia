# Geometry order (`Mesh.Curve`) in Radia workflows

This page states how to obtain and keep a curved (high-order) geometry map for
NGSolve meshes used by Radia, and why `Mesh.Curve` must not be called on a
mesh loaded from a `.vol`. Behaviour was measured on NGSolve/Netgen 6.2.2607
and is pinned by `tests/test_curve_order_contract.py`.

Three different things are often called "order":

| Term | Meaning | Set by |
| --- | --- | --- |
| geometry order | polynomial degree of the element map (curved elements) | `Mesh.Curve(p)` on live CAD, or the order baked into a `.vol` export |
| FE order | degree of the finite-element space | `H1/HCurl/HDiv(mesh, order=k)` |
| `curvaturesafety` | mesh *size* refinement near curved CAD faces (a meshing parameter) | `GenerateMesh(..., curvaturesafety=c)` |

Raising the FE order does not curve the geometry, and `curvaturesafety` does
not curve anything; it only makes flat elements smaller near curvature.

## What a loaded `.vol` does

`ngsolve.Mesh("x.vol")` applies the `curvedelements` section stored in the file.
The saved geometry order is in effect immediately; `GetCurveOrder()` reports it.
A Netgen/OCC `.vol` also embeds its OCC geometry after the `endmesh` line; a
Cubit `export netgen` file does not.

Calling `Curve(p)` on a loaded mesh rebuilds the curved nodes from the geometry
Netgen associates with the mesh. Measured outcomes (sphere, `R = 10 mm`,
`maxh = R/2`, boundary-area error):

| `.vol` content | as loaded | `Curve(p)` after loading |
| --- | --- | --- |
| Netgen/OCC, saved at order 3, CAD in the same units | +0.012% (order 3) | correct for `p >= 2` (order 2: −0.12%); `p = 1` straightens it (−5.3%) |
| Netgen/OCC, saved flat, CAD in the same units | −5.3% (order 1) | correct (order 2: −0.12%) |
| meshed in mm, `ngmesh.Scale(1e-3)`, saved | −5.3% (order 1) | **+7.9e8 %** (nodes projected onto the mm CAD) |
| CAD-less (Cubit-style) curved order 3, fresh process | +0.012% (order 3) | **−5.3% for every `p`** while `GetCurveOrder()` reports `p` |
| same CAD-less file after another `.vol` with CAD was loaded in the process | +0.012% | **curved onto that other geometry** (+7.9e8 % when it was the mm sphere) |

A `.vol` with embedded CAD keeps that CAD attached to its mesh even after other
files are loaded. A CAD-less mesh has none, and `GetGeometry()` then returns the
most recently loaded geometry, so the result of a post-load `Curve` depends on the
history of the process. This is why earlier notes in this repository
disagreed ("flattens", "works", "corrupts the mapping"): each described one
row of the table.

The 2026 induction-heating thermal mesh that reported a 6.9e6 m² boundary
(2.97e8 × the true 0.0233 m²) is the third row: it was meshed in millimetres
and rescaled with `ngmesh.Scale(1e-3)`.

## Rules

1. **Never call `Mesh.Curve` on a loaded `.vol`.** The stored order is already
   used. If a higher order is required, call
   `radia.mesh_curve.ensure_curve_order(mesh, order, vol_path=path)` inside the
   caller's `TaskManager`. It keeps a sufficient stored order (`"kept"` does not
   re-examine the geometry), takes the CAD provenance from the file because
   `GetGeometry()` cannot tell embedded from borrowed geometry (an OCC archive
   after `endmesh`: `TextOutArchive`, `netgen::OCCGeometry`, `CASCADE
   Topology`; any other trailing text is refused), requires the current
   vertices to equal the points stored in that file (a mesh moved or rescaled
   after loading is refused), rejects a CAD box that differs from the mesh box
   by more than a factor of 2 or is shifted (mesh rescaled before saving), and
   rejects a curving that changes boundary or domain measure by
   more than a factor of 3. These are sanity checks for gross errors, not a
   proof of a valid map. After a rejection with `mesh_modified=True` the mesh
   must be reloaded.
2. **Curve where the CAD lives.** Call `Curve(p)` right after
   `OCCGeometry(shape).GenerateMesh(...)` in the same process, inside
   `TaskManager`, then save. For Cubit, export the order:
   `export netgen "<f>.vol" order <p>`.
3. **Mesh in final units.** Scale the shape, not the mesh:
   `shape = shape.Scale(Pnt(0, 0, 0), 1e-3)` before `GenerateMesh`.
   `ngmesh.Scale` after meshing leaves the embedded CAD in the old units
   (lint rule `netgen-scale-after-generate`).
4. **A saved `GridFunction` belongs to the geometry map it was computed on.**
   Load the same `.vol` and do not curve it on the reader side; record the
   geometry order next to the field (e.g. `geometry_curve_order` in the result
   JSON) and check `GetCurveOrder()` against it when reading.
5. **`GetCurveOrder()` is not evidence.** After any change of geometry order,
   compare boundary area and volume before and after (`mesh_measures`).
6. **Use the geometry order the method needs.** FE order and geometry order
   are independent. Curved boundaries that the method relies on need a curved
   geometry (the Kelvin shell agreed with a 2D axisymmetric reference to about
   1.15% at geometry order 2); a flat geometry at FE order `p` keeps an `O(h²)`
   geometric error that a higher `p` cannot remove. Flat models (boxes) need no
   curving at any FE order. The accelerator-magnet panel therefore requires the
   order only when a Kelvin shell is present and otherwise keeps the stored
   order; the Kelvin benchmark panel always requires it.

## Workpiece boundary-element surfaces are flat today

The P1 workpiece BEM routes (`calc_inductance` panels, `bem_coupled_solver`,
`peec_coupled_bem_solver`, the loop/cohomology extension) extract the boundary
by copying vertices into a new surface mesh (`_extract_bnd_only_inline`,
`_extract_surface_mesh_filtered`). The curved-element data and the CAD
association are dropped, so these BEM solves always use flat P1 panels,
whatever order the `.vol` carries. `extract_surface_curved` accepts
`geom_order` but evaluates only `GetTrafo` at six reference points, so the
argument is nominal: the P2 panel geometry comes from the parent mesh's own
element maps, so a flat parent gives flat panels.
`calc_inductance` records `wp_vol_curve_order` but does not use it.

Consequences: the BEM surface area and normals carry an `O(h²/R²)` error at
small radii (bores, fillets), which shows up in mesh-convergence studies.
Curved workpiece panels require a P2 (Tri6) surface built from the parent
mesh's element maps; this is not implemented for the P1 loop extension.

## Lint rules

| Rule | Severity | Flags |
| --- | --- | --- |
| `ngsolve-curve-after-vol-import` | HIGH | `.Curve()` on a mesh created from a `.vol` literal, path variable or `.Load` |
| `netgen-scale-after-generate` | HIGH | `.Scale()` anywhere in a scope on a name that is ever bound to a `GenerateMesh` result (or its alias); flow-insensitive, so a name reused for a shape is reported too |
| `ngsbem-missing-curvaturesafety` | MODERATE | ngsbem meshes without `curvaturesafety` |
| `curve-without-export-curved` (cubit lint) | MODERATE | `.Curve(p>=2)` without a curved Cubit export |

The `.vol` rule cannot see `Mesh(vol_file)` with a variable whose value comes
from a function argument; route such code through `ensure_curve_order`.

## Known places still to review

These call `Curve` on meshes that may have been loaded; they were not changed
together with this page and should move to `ensure_curve_order`:

- `src/radia/vim/_vim.py` `_curve_mesh` (no-op when the stored order suffices,
  otherwise curves a possibly loaded mesh),
- `packages/radia-mcp/.../vol2d_scalar.py`, `vol2d_thermal.py`,
  `vol2d_circuit.py` (`mesh.Curve(order)` after loading a `.vol`),
- `gmsh_post_export.vol2msh(..., mesh_curve_order>=2)`.

`curvaturesafety` advice also differs between `ngsbem_inductance` (≥ 1,
degenerate elements otherwise) and `docs/bem_extractor` (0.5); both statements
concern mesh size and BEM conditioning, not geometry order.

## Measurement scripts

`tests/test_curve_order_contract.py` runs each scenario in a fresh interpreter
so that the process-global geometry cannot leak between cases.
