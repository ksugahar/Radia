# Radia MCP Server

Model Context Protocol server for Radia electromagnetic simulation (standalone).

## Overview

This server provides Radia magnetostatics simulation capabilities through MCP, **without NGSolve dependency**. For NGSolve integration, use the companion `mcp_server_ngsolve` server.

## Features

- **Geometry**: Magnets, coils, containers
- **Materials**: Linear, nonlinear, permanent magnets
- **Solver**: LU, BiCGSTAB, HACApK
- **Fields**: B, H, A, Phi calculations
- **Workspace**: Export to shared workspace for NGSolve coupling

## Installation

```bash
# Ensure Radia is built
cd S:\Radia\01_Github
powershell -ExecutionPolicy Bypass -File BuildMSVC.ps1

# Install MCP SDK
pip install mcp
```

## Tools (35 total)

### Geometry (6)
- `radia_geometry_create_recmag` - Rectangular magnet
- `radia_geometry_create_hexahedron` - Hexahedral element
- `radia_geometry_create_tetrahedron` - Tetrahedral element
- `radia_geometry_create_racetrack_coil` - Race-track coil
- `radia_geometry_create_container` - Group objects
- `radia_geometry_set_units` - Set units (m/mm)

### Material (6)
- `radia_material_create_linear` - Linear material
- `radia_material_create_linear_anisotropic` - Anisotropic material
- `radia_material_create_nonlinear` - B-H curve material
- `radia_material_create_permanent_magnet_fixed` - Fixed PM
- `radia_material_create_permanent_magnet_linear` - Linear PM
- `radia_material_apply` - Apply material to object

### Solver (5)
- `radia_solver_solve` - Run solver
- `radia_solver_set_hacapk_params` - H-matrix parameters
- `radia_solver_set_relax_param` - Relaxation parameter
- `radia_solver_create_background_field` - Background field
- `radia_solver_get_hacapk_stats` - Solver statistics

### Field (5)
- `radia_field_compute` - Field at point
- `radia_field_compute_line` - Field along line
- `radia_field_compute_grid` - Field on grid
- `radia_field_force` - Force calculation
- `radia_field_torque` - Torque calculation

### Export & Workspace (4)
- `radia_export_vtk_legacy` - VTK export
- `radia_workspace_create_session` - New session
- `radia_workspace_export_object` - Export to workspace
- `radia_workspace_list_sessions` - List sessions

### NGSolve Integration (5)
- `radia_ngsolve_create_field` - Create field data for NGSolve CoefficientFunction
- `radia_ngsolve_batch_evaluate` - Batch field evaluation with H-matrix (10-100x speedup)
- `radia_ngsolve_enable_hmatrix` - Enable H-matrix acceleration
- `radia_ngsolve_export_field_grid` - Export structured grid to NPZ format
- `radia_ngsolve_mesh_field_sample` - Sample field at mesh vertices

### Diagnostic & Debugging (4)
- `radia_server_info` - Get server version and status
- `radia_list_objects` - List all objects in server state
- `radia_get_object_info` - Get detailed information about an object
- `radia_clear_state` - Clear all objects from server state (reset)

## Usage

### Claude Desktop Configuration

```json
{
  "mcpServers": {
    "radia": {
      "command": "python",
      "args": ["-m", "mcp_server_radia.server"],
      "cwd": "S:\\Radia\\01_Github",
      "env": {
        "PYTHONPATH": "S:\\Radia\\01_Github\\build\\Release;S:\\Radia\\01_Github\\src\\radia"
      }
    }
  }
}
```

### Standalone Example

```
User: "Create a 10cm cubic NdFeB magnet"

Claude calls:
  radia_geometry_set_units(units="m")
  radia_geometry_create_recmag(
    center=[0, 0, 0],
    dimensions=[0.1, 0.1, 0.1],
    magnetization=[0, 0, 954930]
  )
```

### Export for NGSolve

```
User: "Export the magnet to workspace for NGSolve analysis"

Claude calls:
  radia_workspace_export_object(
    object_name="magnet",
    export_geometry=True,
    export_fields=True
  )

Result: Session ID that NGSolve server can import
```

### NGSolve Integration with H-Matrix

```
User: "Evaluate field at 10000 points using H-matrix acceleration"

Claude calls:
  radia_ngsolve_enable_hmatrix(enable=True, precision=1e-6)
  radia_ngsolve_batch_evaluate(
    radia_object_name="magnet",
    field_type="b",
    points=[[x1,y1,z1], ...],  # 10k points
    use_hmatrix=True
  )

Result: Field values with 10-100x speedup vs standard evaluation
```

### Export Field Grid for NGSolve Import

```
User: "Export B field on 20x20x20 grid from -0.1 to 0.1m cube"

Claude calls:
  radia_ngsolve_export_field_grid(
    radia_object_name="magnet",
    field_type="b",
    bbox_min=[-0.1, -0.1, -0.1],
    bbox_max=[0.1, 0.1, 0.1],
    grid_size=[20, 20, 20],
    output_file="field_grid.npz"
  )

Result: NPZ file compatible with scipy.interpolate for NGSolve
```

## Validation Examples

### Rotating Magnet Eddy Current Analysis

Complete validation examples demonstrating Radia-NGSolve coupling for transient eddy current analysis are available in [`examples/NGSolve_Integration/rotating_magnets/`](../examples/NGSolve_Integration/rotating_magnets/).

**Physical Model:**
- Rotating 1mm³ permanent magnet (Br = 0.2 T) moving over 0.5mm copper plate
- 180 timesteps with 4°/step rotation (2 full rotations)
- Movement range: X = -6mm to 4mm, fixed height Y = 2mm

**Two Formulation Comparison:**

1. **A-Φ Method** (Vector-Scalar Potential)
   - Vector potential: A = A_ext + A_r
   - Radia provides A_ext via `'a'` field type
   - HCurl(nograds=True) + H1 mixed formulation
   - Direct computation of B = curl(A), E = -∂A/∂t - grad(Φ)
   - File: [`comparison_A_Phi_method.py`](../examples/NGSolve_Integration/rotating_magnets/comparison_A_Phi_method.py)

2. **T-Ω Method** (Current-Magnetic Scalar Potential)
   - Current potential: J = curl(T)
   - Radia provides H_ext via `'h'` field type
   - HCurl(nograds=True) + H1 mixed formulation
   - Conductor-only T field (DOF reduction)
   - File: [`comparison_T_Omega_method.py`](../examples/NGSolve_Integration/rotating_magnets/comparison_T_Omega_method.py)

**Validation Results:**
- Maxwell relation verification: curl(A_ext) ≈ B_ext/μ₀ (< 0.1% error)
- Both formulations reproduce similar eddy current patterns
- Magnetic energy and Joule loss calculations
- CSV outputs: field data, statistics, energy metrics

**Radia Integration Features Demonstrated:**
- External field provision: `radia_ngsolve_create_field()` for A_ext, B_ext, H_ext
- Time-dependent analysis: 180-step transient simulation
- Batch field evaluation: efficient multi-point field computation
- Workspace coupling: Radia geometry → NGSolve FEM solver

For detailed documentation, see the [README.md](../examples/NGSolve_Integration/rotating_magnets/README.md) in the validation directory.

## Troubleshooting

### Solver Issues

**Problem: Solver does not converge (rad.Solve() returns error)**

**原因 (Causes):**
- Material permeability too high or nonlinear B-H curve poorly defined
- Geometry has very small or degenerate elements
- Insufficient relaxation for nonlinear materials

**解決策 (Solutions):**
```python
# 1. Increase relaxation parameter for nonlinear materials
radia_solver_set_relax_param(relaxation_param=0.5)  # Default 1.0 → reduce to 0.5

# 2. Check geometry for degenerate elements
radia_list_objects()  # List all objects
radia_get_object_info(object_name="magnet")  # Check specific object

# 3. Use BiCGSTAB for large problems
radia_solver_solve(
    object_name="magnet",
    method="bicgstab",  # Not "lu"
    precision=1e-6,
    max_iter=1000
)

# 4. For H-matrix solver
radia_solver_set_hacapk_params(
    precision=1e-6,
    max_iter=1000,
    use_hmatrix=True
)
```

**Problem: Solver is too slow**

**原因 (Causes):**
- Too many elements (mesh too fine)
- Using LU solver for large problems
- H-matrix not enabled for large systems

**解決策 (Solutions):**
```python
# 1. Enable H-matrix acceleration (O(N log N) vs O(N²))
radia_solver_set_hacapk_params(
    precision=1e-6,
    use_hmatrix=True
)

# 2. Use iterative solver for large problems
radia_solver_solve(
    object_name="magnet",
    method="bicgstab",
    precision=1e-6
)

# 3. Check element count
radia_list_objects()
# If > 500 elements, consider coarser subdivision
```

### Field Evaluation Issues

**Problem: Field values are zero or incorrect**

**原因 (Causes):**
- Solver not run before field evaluation
- Evaluation point inside magnet geometry (singularity)
- Units mismatch (Radia in mm, NGSolve in m)

**解決策 (Solutions):**
```python
# 1. ALWAYS set units before creating geometry
radia_geometry_set_units(units="m")  # Required for NGSolve coupling

# 2. Run solver before evaluating fields
radia_solver_solve(object_name="magnet", method="bicgstab")

# 3. Evaluate fields OUTSIDE magnet geometry
# Rule: evaluation point > 0.001m from magnet surfaces
radia_field_compute(
    object_name="magnet",
    point=[0.12, 0, 0],  # Outside 0.1m magnet
    field_type="b"
)
```

**Problem: H-matrix batch evaluation slower than expected**

**原因 (Causes):**
- H-matrix not properly enabled
- Batch size too small (< 100 points)
- Precision too tight (< 1e-8)

**解決策 (Solutions):**
```python
# 1. Explicitly enable H-matrix
radia_ngsolve_enable_hmatrix(
    enable=True,
    precision=1e-6  # Not 1e-10
)

# 2. Use large batches (> 1000 points for best speedup)
points = [[x, y, z] for ... in range(10000)]  # Large batch
radia_ngsolve_batch_evaluate(
    radia_object_name="magnet",
    field_type="b",
    points=points,
    use_hmatrix=True
)

# 3. Check speedup
# Expected: 10-100× faster for N > 1000 points
```

### NGSolve Integration Issues

**Problem: Field grid export produces NaN values**

**原因 (Causes):**
- Evaluation grid extends inside magnet geometry
- Grid resolution too fine near magnet surfaces
- Solver precision insufficient

**解決策 (Solutions):**
```python
# 1. Keep grid outside magnet geometry
# For 0.1m magnet, use bbox_min/max > ±0.11m
radia_ngsolve_export_field_grid(
    radia_object_name="magnet",
    field_type="b",
    bbox_min=[-0.15, -0.15, -0.15],  # Outside magnet
    bbox_max=[0.15, 0.15, 0.15],
    grid_size=[20, 20, 20],
    output_file="field_grid.npz"
)

# 2. Reduce grid resolution near boundaries
# Use grid_size=[15, 15, 15] instead of [30, 30, 30]

# 3. Increase solver precision
radia_solver_solve(
    object_name="magnet",
    method="bicgstab",
    precision=1e-8  # Was 1e-6 → tighten
)
```

**Problem: Workspace export fails or session not found**

**原因 (Causes):**
- Object not created or named incorrectly
- Workspace directory not accessible
- Insufficient permissions

**解決策 (Solutions):**
```python
# 1. Verify object exists in state
radia_list_objects()  # Check object is listed

# 2. Use exact object name
radia_workspace_export_object(
    object_name="magnet",  # Must match creation name
    export_geometry=True,
    export_fields=True
)

# 3. List created sessions
radia_workspace_list_sessions()
# Copy session_id for NGSolve import

# 4. Check workspace directory
# Verify: S:\Radia\01_Github\mcp_shared\ is writable
```

### Geometry Creation Issues

**Problem: Object creation fails with invalid geometry**

**原因 (Causes):**
- Vertex ordering incorrect for hexahedron/tetrahedron
- Dimensions too small or negative
- Overlapping geometry (container issues)

**解決策 (Solutions):**
```python
# 1. Check hexahedron vertex ordering
# Must follow Radia convention:
# vertices[0-3]: bottom face (z=z_min), counterclockwise
# vertices[4-7]: top face (z=z_max), counterclockwise
vertices = [
    [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],  # Bottom
    [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1],  # Top
]

# 2. Verify dimensions are positive
radia_geometry_create_recmag(
    center=[0, 0, 0],
    dimensions=[0.1, 0.1, 0.1],  # All > 0
    magnetization=[0, 0, 954930]
)

# 3. Check container objects
radia_list_objects()  # Verify all sub-objects exist
```

### Diagnostic Tools Usage

**サーバ状態確認 (Check server state):**
```python
# 1. Get server information
radia_server_info()
# Returns: version, Radia availability, object count

# 2. List all objects in memory
radia_list_objects()
# Shows all Radia objects with IDs

# 3. Get detailed object information
radia_get_object_info(object_name="magnet")
# Returns: type, Radia ID, properties

# 4. Reset server state (clear all objects)
radia_clear_state(confirm=True)
# WARNING: Deletes all objects from memory
```

## Development Policy

### Branch Management
- **Feature branches**: Create feature branches for Pull Request purposes
  - Feature branches are temporary and should be deleted after PR approval and merge
  - Naming convention: `feature/mcp-tools-enhancement`, `fix/hmatrix-bug`
- **Master branch**: Always kept in sync with the latest version
  - Research lab internal policy: master branch is continuously updated with latest stable code
  - Direct commits to master are allowed for internal development
  - Master branch always reflects the current working state

### Pull Request Workflow (for upstream contributions)
1. Create a feature branch from master
2. Make your changes and commit with descriptive messages
3. Push to your fork and create a Pull Request
4. After PR approval and merge, delete the feature branch
5. Keep master branch clean and stable

### Internal Development (Research Lab)
- Master branch is the primary development branch
- Always keep master in sync with the latest version
- Feature branches are used only for PR submissions to upstream
- Internal changes can be committed directly to master

## See Also

- [NGSolve MCP Server](https://github.com/ksugahar/ngsolve) - Companion server for NGSolve FEM
- [Dual Server Deployment](../docs/DUAL_SERVER_DEPLOYMENT.md)
- [MCP Server Architecture](../docs/MCP_SERVER_ARCHITECTURE.md)
