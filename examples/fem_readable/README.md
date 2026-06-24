# Readable FEM examples

Small validation-class FEM examples written to be easy to translate into
MATLAB/Gypsilab-style teaching scripts.  The focus is readability and explicit
assembly steps rather than production solver performance.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_p1_tet_patch_test.py`](validation_p1_tet_patch_test.py) | P1 tetrahedron stiffness assembly on a 12-tet cube star patch; affine Dirichlet patch tests on regular and distorted interior nodes | `p1_tetrahedron_stiffness`, `p1_tetrahedron_geometry` |
| [`validation_p1_tet_poisson_convergence.py`](validation_p1_tet_poisson_convergence.py) | Validation-class 3D Poisson manufactured-solution convergence on structured P1 tetra meshes | `p1_tetrahedron_stiffness`, `p1_tetrahedron_geometry` |
| [`validation_p1_surface_trace_coupling.py`](validation_p1_surface_trace_coupling.py) | Validation-class `.vol` FEM/BEM trace example: boundary P1 surface mass, surface stiffness, load vector, closed-surface normal checks, and shared H1/BEM node ids | `parse_netgen_tri_tet_vol`, `p1_surface_triangle_*`, `first_order_fem_bem_topology` |
| [`validation_hcurl_rwg_trace_orientation.py`](validation_hcurl_rwg_trace_orientation.py) | Validation-class `.vol` first-order HCurl/RWG trace example on a four-tet star mesh with one interior node; checks RWG edge orientation, RWG-to-HCurl ids, and interior-edge exclusion | `parse_netgen_tri_tet_vol`, `first_order_fem_bem_topology`, `surface_edge_manifold_summary` |

```powershell
python validation_p1_tet_patch_test.py
python validation_p1_tet_poisson_convergence.py
python validation_p1_surface_trace_coupling.py
python validation_hcurl_rwg_trace_orientation.py
```

The examples are public-safe clean-room formulas. MATLAB/Lukas/Gypsilab
prototypes can mirror the scripts, but no private MATLAB repository code is
included here.
