# Cubit Mesh Export Examples

Example scripts for the radia Cubit C++ plugin (`cubit_mesh_export.ccm`), organized by export format.

Presentation notebook:
[`docs/cubit_mesh_export/cubit_mesh_export_showcase.ipynb`](../../docs/cubit_mesh_export/cubit_mesh_export_showcase.ipynb)
is the result-saved talk layer. It combines the validation summary JSONs,
high-order curved-hex NGSolve results, and source excerpts for Q&A.

Mesh-evaluation demo:
[`docs/cubit_mesh_export/netgen/p_convergence_demo.ipynb`](netgen/p_convergence_demo.ipynb)
is the result-saved p-convergence demonstration. It records a Cubit batch run
that calls APREPRO export commands through `cubit.cmd(...)` and evaluates the
exported `.vol` files with `src/radia/panels/calc_mesh_eval.py`. This is a
documentation demo, not a Cubit menu action or engineering design panel.

## Folder Structure

| Folder | Format | Description |
|--------|--------|-------------|
| [gmsh/](gmsh/) | Gmsh (.msh) | v4.1, order 1-3 |
| [nastran/](nastran/) | Nastran BDF (.bdf) | CTETRA/CHEXA/CTRIA, order 1-2 |
| [other_formats/](other_formats/) | Custom formats | FreeFEM, ANSYS CDB, Lukas 2D |

## Validation examples

| Example | Shows |
|---|---|
| [`validation_vol_surface_closure.py`](validation_vol_surface_closure.py) | Netgen `.vol` boundary orientation, vector-area closure, and boundary/tet volume agreement |
| [`validation_vol_fem_bem_topology.py`](validation_vol_fem_bem_topology.py) | First-order FEM/BEM topology view: closed surface edge manifold, Euler characteristic, compact scalar-BEM nodes, RWG-to-HCurl edge trace |
| [`validation_vol_p1_trace_matrix.py`](validation_vol_p1_trace_matrix.py) | First-order H1 FEM to scalar BEM trace matrix as one-based sparse COO rows/cols/values |
| [`validation_p1_surface_triangle_element_summary.py`](validation_p1_surface_triangle_element_summary.py) | P1 boundary-triangle teaching block: geometry, mass, stiffness, load vector, and one-based sparse triplets |
| [`validation_vol_boundary_oriented_edges.py`](validation_vol_boundary_oriented_edges.py) | First-order boundary triangle oriented-edge rows for readable RWG-style trace assembly and HCurl edge matching |
| [`validation_vol_mesh_health.py`](validation_vol_mesh_health.py) | Netgen `.vol` first-order FEM/BEM mesh health: shape quality, surface closure, boundary-to-tet face consistency, and worst-element rows |
| [`validation_vol_tet_quality.py`](validation_vol_tet_quality.py) | Netgen `.vol` tetrahedron quality: edge ratio, inradius, circumradius, radius-ratio quality, corner-normalized Jacobian quality, and optional real Cubit export evaluation |
| [`validation_vol_surface_triangle_quality.py`](validation_vol_surface_triangle_quality.py) | Netgen `.vol` boundary-triangle quality: area, edge ratio, inradius/circumradius quality, angle range, sliver detection, and optional real Cubit export evaluation |
| [`validation_vol_boundary_normal_vectors.py`](validation_vol_boundary_normal_vectors.py) | Netgen `.vol` boundary normal/vector-area rows for Maxwell-stress force integration over named sidesets |
| [`validation_vol_boundary_pressure_force.py`](validation_vol_boundary_pressure_force.py) | Netgen `.vol` boundary pressure-force rows: scalar pressure times oriented vector area over named sidesets |
| [`validation_vol_boundary_pressure_moment.py`](validation_vol_boundary_pressure_moment.py) | Netgen `.vol` boundary pressure-force/moment rows: triangle-centroid moment integration and generic resultant reduction |
| [`validation_vol_boundary_pressure_resultant.py`](validation_vol_boundary_pressure_resultant.py) | Netgen `.vol` boundary pressure resultant summary: closed-surface pressure cancellation and one-sided pressure force/moment |
| [`validation_vol_boundary_traction_moment.py`](validation_vol_boundary_traction_moment.py) | Netgen `.vol` boundary vector-traction force/moment rows: constant global traction over named sidesets |
| [`validation_vol_boundary_inventory.py`](validation_vol_boundary_inventory.py) | Netgen `.vol` named-boundary inventory for Cubit/Coreform sidesets: per-boundary area, triangle count, and trace-node ids |
| [`validation_vol_boundary_edge_inventory.py`](validation_vol_boundary_edge_inventory.py) | Netgen `.vol` boundary-local edge inventory: separates sideset perimeter edges from triangulation diagonals |
| [`validation_vol_boundary_condition_assignment.py`](validation_vol_boundary_condition_assignment.py) | Netgen `.vol` boundary-condition assignment audit: map condition labels by boundary number/name and catch missing or unknown keys |
| [`validation_vol_material_interface.py`](validation_vol_material_interface.py) | Netgen `.vol` material/interface inventory: material volumes, exterior/interface areas, and `domin/domout` boundary incidence |
| [`validation_vol_boundary_tet_face_incidence.py`](validation_vol_boundary_tet_face_incidence.py) | Netgen `.vol` boundary triangle to tetrahedron face incidence: exterior/interface adjacency, orphan detection, and `domin/domout` material consistency |

```powershell
python validation_vol_surface_closure.py
python validation_vol_fem_bem_topology.py
python validation_vol_p1_trace_matrix.py
python validation_p1_surface_triangle_element_summary.py
python validation_vol_boundary_oriented_edges.py
python validation_vol_mesh_health.py
python validation_vol_tet_quality.py
python validation_vol_surface_triangle_quality.py
python validation_vol_boundary_normal_vectors.py
python validation_vol_boundary_pressure_force.py
python validation_vol_boundary_pressure_moment.py
python validation_vol_boundary_pressure_resultant.py
python validation_vol_boundary_traction_moment.py
python validation_vol_boundary_inventory.py
python validation_vol_boundary_edge_inventory.py
python validation_vol_boundary_condition_assignment.py
python validation_vol_material_interface.py
python validation_vol_boundary_tet_face_incidence.py
```

Each folder contains:
- Example Python scripts (using `cubit.cmd('export ...')`)
- Pre-generated sample output files
- README.md with format-specific documentation

## Running Examples

```powershell
python <folder>/<script>.py
```

Requires Cubit installed with the radia plugin (`radia-setup`).

## Export Commands (C++ .ccm plugin)

```
export gmsh "file.msh" [order {1|2}] [dimension {2|3}] [overwrite]
export nastran_bdf "file.bdf" [order {1|2}] [dimension {2|3}] [nopyramid] [overwrite]
export vtk "file.vtk" [order {1|2}] [dimension {2|3}] [overwrite]
```

## See Also

- [src/cubit_plugin/](../../src/cubit_plugin/) - C++ plugin source
- [tests/cubit/](../../tests/cubit/) - Automated tests (42/42 PASS)
