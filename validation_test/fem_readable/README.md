# Readable FEM validation corpus

Small validation-class FEM checks written to be easy to translate into
first-order teaching scripts.  The focus is readability and explicit assembly
steps rather than production solver performance.

Docs surface: the result-saved notebook
[`docs/fem_readable/fem_readable_validation_archive.ipynb`](../../docs/fem_readable/fem_readable_validation_archive.ipynb)
collects the 12 scripts, their existing summary JSON payloads, full source, and
SHA-256 hashes. The executable lane is this `validation_test/fem_readable`
directory; the notebook is the human-facing rendered archive.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_p1_tet_patch_test.py`](validation_p1_tet_patch_test.py) | P1 tetrahedron stiffness assembly on a 12-tet cube star patch; affine Dirichlet patch tests on regular and distorted interior nodes | `p1_tetrahedron_stiffness`, `p1_tetrahedron_geometry` |
| [`validation_p1_tet_poisson_convergence.py`](validation_p1_tet_poisson_convergence.py) | Validation-class 3D Poisson manufactured-solution convergence on structured P1 tetra meshes | `p1_tetrahedron_stiffness`, `p1_tetrahedron_geometry` |
| [`validation_p1_surface_trace_coupling.py`](validation_p1_surface_trace_coupling.py) | Validation-class `.vol` FEM/BEM trace example: boundary P1 surface mass, surface stiffness, load vector, closed-surface normal checks, and shared H1/BEM node ids | `parse_netgen_tri_tet_vol`, `p1_surface_triangle_*`, `first_order_fem_bem_topology` |
| [`validation_p1_surface_single_layer_moments.py`](validation_p1_surface_single_layer_moments.py) | Validation-class P1 surface-density moments for a readable Laplace single-layer BEM far-field gate | `p1_surface_triangle_density_moments`, `laplace_single_layer_far_potential` |
| [`validation_p1_tet_robin_trace_system.py`](validation_p1_tet_robin_trace_system.py) | Validation-class `.vol` P1 volume/boundary trace system: dense tet stiffness plus boundary Robin mass/load, with constant-solution and pure-Neumann gates | `parse_netgen_tri_tet_vol`, `assemble_p1_tet_robin_system`, `first_order_fem_bem_topology` |
| [`validation_p1_tet_flux_trace.py`](validation_p1_tet_flux_trace.py) | Validation-class P1 tetrahedron gradient/flux trace: affine field gradient, physical flux, outward face Neumann rows, and stiffness-energy identity | `p1_tetrahedron_gradient`, `p1_tetrahedron_flux`, `p1_tetrahedron_boundary_fluxes` |
| [`validation_p1_tet_face_trace_projection.py`](validation_p1_tet_face_trace_projection.py) | Validation-class P1 tetrahedron face trace projection: shared `.vol` node ids, surface mass projection `int u_h N_j dS`, and trace `L2` energy | `parse_netgen_tri_tet_vol`, `p1_tetrahedron_face_trace_summary` |
| [`validation_surface_triangle_constant_traction_load.py`](validation_surface_triangle_constant_traction_load.py) | Validation-class P1 surface-triangle constant vector traction load: integrated force, equivalent nodal loads, and force/moment preservation | `surface_triangle_constant_traction_load_summary`, `force_moment_resultant_summary` |
| [`validation_surface_maxwell_force_trace.py`](validation_surface_maxwell_force_trace.py) | Validation-class `.vol` surface Maxwell traction trace: oriented boundary triangles, integrated force, and P1 equivalent nodal force loads | `parse_netgen_tri_tet_vol`, `surface_triangle_maxwell_traction_summary`, `air_gap_maxwell_pressure` |
| [`validation_tetrahedron_lorentz_force_load.py`](validation_tetrahedron_lorentz_force_load.py) | Validation-class P1 tetrahedron Lorentz body-force load: constant `J x B`, integrated force, equivalent nodal loads, and force/moment preservation | `tetrahedron_lorentz_force_summary`, `force_moment_resultant_summary` |
| [`validation_hcurl_rwg_trace_orientation.py`](validation_hcurl_rwg_trace_orientation.py) | Validation-class `.vol` first-order HCurl/RWG trace example on a four-tet star mesh with one interior node; checks RWG edge orientation, RWG-to-HCurl ids, and interior-edge exclusion | `parse_netgen_tri_tet_vol`, `first_order_fem_bem_topology`, `surface_edge_manifold_summary` |
| [`validation_vol_boundary_components.py`](validation_vol_boundary_components.py) | Validation-class `.vol` FEM/BEM block example: disconnected boundary components, trace node blocks, Euler checks, and RWG-to-HCurl edge trace mapping | `parse_netgen_tri_tet_vol`, `surface_connected_components`, `first_order_fem_bem_topology` |

```powershell
python validation_p1_tet_patch_test.py
python validation_p1_tet_poisson_convergence.py
python validation_p1_surface_trace_coupling.py
python validation_p1_surface_single_layer_moments.py
python validation_p1_tet_robin_trace_system.py
python validation_p1_tet_flux_trace.py
python validation_p1_tet_face_trace_projection.py
python validation_surface_triangle_constant_traction_load.py
python validation_surface_maxwell_force_trace.py
python validation_tetrahedron_lorentz_force_load.py
python validation_hcurl_rwg_trace_orientation.py
python validation_vol_boundary_components.py
```

The checks are public-safe clean-room formulas. Private teaching prototypes can
mirror the scripts, but no private repository code is included here.
