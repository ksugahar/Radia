# Readable FEM examples

Small validation-class FEM examples written to be easy to translate into
MATLAB/Gypsilab-style teaching scripts.  The focus is readability and explicit
assembly steps rather than production solver performance.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_p1_tet_patch_test.py`](validation_p1_tet_patch_test.py) | P1 tetrahedron stiffness assembly on a 12-tet cube star patch; affine Dirichlet patch tests on regular and distorted interior nodes | `p1_tetrahedron_stiffness`, `p1_tetrahedron_geometry` |
| [`validation_p1_tet_poisson_convergence.py`](validation_p1_tet_poisson_convergence.py) | Validation-class 3D Poisson manufactured-solution convergence on structured P1 tetra meshes | `p1_tetrahedron_stiffness`, `p1_tetrahedron_geometry` |
| [`validation_p1_surface_trace_coupling.py`](validation_p1_surface_trace_coupling.py) | Validation-class `.vol` FEM/BEM trace example: boundary P1 surface mass, surface stiffness, load vector, closed-surface normal checks, and shared H1/BEM node ids | `parse_netgen_tri_tet_vol`, `p1_surface_triangle_*`, `first_order_fem_bem_topology` |
| [`validation_p1_tet_robin_trace_system.py`](validation_p1_tet_robin_trace_system.py) | Validation-class `.vol` P1 volume/boundary trace system: dense tet stiffness plus boundary Robin mass/load, with constant-solution and pure-Neumann gates | `parse_netgen_tri_tet_vol`, `assemble_p1_tet_robin_system`, `first_order_fem_bem_topology` |
| [`validation_p1_tet_flux_trace.py`](validation_p1_tet_flux_trace.py) | Validation-class P1 tetrahedron gradient/flux trace: affine field gradient, physical flux, outward face Neumann rows, and stiffness-energy identity | `p1_tetrahedron_gradient`, `p1_tetrahedron_flux`, `p1_tetrahedron_boundary_fluxes` |
| [`validation_surface_maxwell_force_trace.py`](validation_surface_maxwell_force_trace.py) | Validation-class `.vol` surface Maxwell traction trace: oriented boundary triangles, integrated force, and P1 equivalent nodal force loads | `parse_netgen_tri_tet_vol`, `surface_triangle_maxwell_traction_summary`, `air_gap_maxwell_pressure` |
| [`validation_tetrahedron_lorentz_force_load.py`](validation_tetrahedron_lorentz_force_load.py) | Validation-class P1 tetrahedron Lorentz body-force load: constant `J x B`, integrated force, and equivalent nodal loads | `tetrahedron_lorentz_force_summary` |
| [`validation_hcurl_rwg_trace_orientation.py`](validation_hcurl_rwg_trace_orientation.py) | Validation-class `.vol` first-order HCurl/RWG trace example on a four-tet star mesh with one interior node; checks RWG edge orientation, RWG-to-HCurl ids, and interior-edge exclusion | `parse_netgen_tri_tet_vol`, `first_order_fem_bem_topology`, `surface_edge_manifold_summary` |
| [`validation_vol_boundary_components.py`](validation_vol_boundary_components.py) | Validation-class `.vol` FEM/BEM block example: disconnected boundary components, trace node blocks, Euler checks, and RWG-to-HCurl edge trace mapping | `parse_netgen_tri_tet_vol`, `surface_connected_components`, `first_order_fem_bem_topology` |

```powershell
python validation_p1_tet_patch_test.py
python validation_p1_tet_poisson_convergence.py
python validation_p1_surface_trace_coupling.py
python validation_p1_tet_robin_trace_system.py
python validation_p1_tet_flux_trace.py
python validation_surface_maxwell_force_trace.py
python validation_tetrahedron_lorentz_force_load.py
python validation_hcurl_rwg_trace_orientation.py
python validation_vol_boundary_components.py
```

The examples are public-safe clean-room formulas. MATLAB/Lukas/Gypsilab
prototypes can mirror the scripts, but no private MATLAB repository code is
included here.
