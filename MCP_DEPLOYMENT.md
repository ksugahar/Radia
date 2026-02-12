# MCP Server Deployment Summary

Two independent MCP servers have been deployed to support Radia and NGSolve electromagnetic simulation.

## Deployment Locations

### Radia MCP Server
**Location:** `S:\Radia\01_Github\mcp_server_radia\`

**Features:**
- 26 tools across 5 categories
- Standalone Radia magnetostatics simulation
- Workspace export for NGSolve coupling

**Tools:**
- Geometry (6): RecMag, Hexahedron, Tetrahedron, Racetrack Coil, Container, Units
- Material (6): Linear, Anisotropic, Nonlinear, PM Fixed, PM Linear, Apply
- Solver (5): Solve, HACApK params, Relaxation, Background field, Stats
- Field (5): Compute, Line, Grid, Force, Torque
- Workspace (4): Create session, Export object, Export fields, List sessions

### NGSolve MCP Server
**Location:** `S:\NGSolve\01_Github\mcp_server_ngsolve\`

**Features:**
- 15 tools across 3 categories
- Mesh generation and import
- Radia field coupling via shared workspace
- Kelvin transformation for infinite domains

**Tools:**
- Mesh (4): Box, Cylinder, Import, Get info
- Radia Coupling (4): Import object, Get field data, Create interpolated field, List objects
- Kelvin Transformation (7): Create mesh, Solve, Compute energy, Export VTK, Compare analytical, Adaptive (planned), Check availability

### Shared Workspace
**Location:** `S:\Radia\01_Github\mcp_shared\`
**Symbolic Link:** `S:\NGSolve\01_Github\mcp_shared -> S:\Radia\01_Github\mcp_shared`

**Purpose:**
- Inter-server communication
- Session-based data exchange
- VTK geometry export
- NPZ field data export

## Claude Desktop Configuration

### Dual Server Configuration

Add both servers to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "radia": {
      "command": "python",
      "args": ["-m", "mcp_server_radia.server"],
      "cwd": "S:\Radia\01_Github",
      "env": {
        "PYTHONPATH": "S:\Radia\01_Github\build\Release;S:\Radia\01_Github\src\radia"
      }
    },
    "ngsolve": {
      "command": "python",
      "args": ["-m", "mcp_server_ngsolve.server"],
      "cwd": "S:\NGSolve\01_Github",
      "env": {
        "PYTHONPATH": "S:\NGSolve\01_Github;S:\Radia\01_Github"
      }
    }
  }
}
```

## Usage Workflows

### Workflow 1: Standalone Radia Analysis

```
1. Create Radia geometry
   radia_geometry_create_recmag(...)

2. Set material properties
   radia_material_create_linear(...)
   radia_material_apply(...)

3. Solve magnetostatic problem
   radia_solver_solve(...)

4. Compute fields
   radia_field_compute_line(...)

5. Export results
   radia_export_vtk_legacy(...)
```

### Workflow 2: Radia + NGSolve Coupling

```
Radia Server:
1. Create and solve Radia model
   radia_geometry_create_recmag(...)
   radia_solver_solve(...)

2. Export to workspace
   radia_workspace_export_object(
     object_name="magnet",
     export_geometry=True,
     export_fields=True
   )
   → Returns session_id

NGSolve Server:
3. Create NGSolve mesh
   ngsolve_mesh_create_box(...)

4. Import Radia object
   ngsolve_radia_import_object(
     session_id=<from_radia>,
     radia_object_name="magnet"
   )

5. Create interpolated field
   ngsolve_radia_create_interpolated_field(
     session_id=<from_radia>,
     radia_object_name="magnet",
     mesh_name="fem_mesh"
   )

6. Use in NGSolve FEM analysis
```

### Workflow 3: Kelvin Transformation for Infinite Domains

```
1. Create mesh with Kelvin transformation
   kelvin_create_mesh_with_transform(
     geometry_type="sphere",
     magnetic_region_size=0.5,
     inner_air_radius=1.0,
     kelvin_radius=1.5,
     maxh=0.05
   )

2. Solve using Ω-Reduced Ω method
   kelvin_omega_reduced_omega_solve(
     mesh_name="kelvin_mesh",
     source_field_type="uniform",
     source_field_params={"H0": 1.0, "direction": [0, 0, 1]},
     permeability=100.0,
     use_kelvin=True
   )

3. Compute perturbation energy
   kelvin_compute_perturbation_energy(
     solution_name="kelvin_mesh_omega_solution"
   )

4. Compare with analytical solution
   kelvin_compare_analytical(
     solution_name="kelvin_mesh_omega_solution",
     geometry_params={
       "sphere_radius": 0.5,
       "mu_r": 100.0,
       "H0": 1.0
     }
   )

5. Export to VTK
   kelvin_export_vtk(
     solution_name="kelvin_mesh_omega_solution",
     output_file="results/sphere_kelvin"
   )
```

## Architecture

### Server Independence
- Each server runs in its own process
- No direct inter-process communication
- Communication through file-based shared workspace

### Shared Workspace Structure
```
~/.radia-ngsolve-workspace/
└── sessions/
    └── <session_id>/
        ├── metadata.json
        ├── geometry.vtk
        └── field_data.npz
```

### Session Management
- Radia creates sessions via `radia_workspace_create_session()`
- NGSolve accesses sessions via session_id
- Multiple objects can be stored in one session

## Development and Testing

### Testing Radia Server Standalone
```bash
cd S:\Radia\01_Github
python -m mcp_server_radia.server
```

### Testing NGSolve Server Standalone
```bash
cd S:\NGSolve\01_Github
set PYTHONPATH=S:\NGSolve\01_Github;S:\Radia\01_Github
python -m mcp_server_ngsolve.server
```

### Testing Workspace Communication
```python
# Create test session from Radia
from mcp_shared.workspace import SharedWorkspace
ws = SharedWorkspace()
session_id = ws.create_session("test_session")
print(f"Session ID: {session_id}")

# Access from NGSolve
ws2 = SharedWorkspace()
info = ws2.get_session_info(session_id)
print(f"Session info: {info}")
```

## Repository Structure

### Radia Repository (S:\Radia\01_Github)
```
├── mcp_server_radia/
│   ├── __init__.py
│   ├── server.py
│   ├── README.md
│   └── tools/
│       ├── __init__.py
│       ├── geometry_tools.py
│       ├── material_tools.py
│       ├── solver_tools.py
│       ├── field_tools.py
│       └── export_tools.py
├── mcp_shared/
│   ├── __init__.py
│   └── workspace.py
└── MCP_DEPLOYMENT.md
```

### NGSolve Repository (S:\NGSolve\01_Github)
```
├── mcp_server_ngsolve/
│   ├── __init__.py
│   ├── server.py
│   ├── README.md
│   └── tools/
│       ├── __init__.py
│       ├── mesh_tools.py
│       ├── radia_coupling_tools.py
│       └── kelvin_transform_tools.py
└── mcp_shared -> S:\Radia\01_Github\mcp_shared (symlink)
```

## Known Limitations

### Current
- Adaptive mesh refinement not yet implemented (planned)
- VTK export is legacy format (ASCII)
- Field interpolation uses linear interpolation only

### Future Enhancements
- A-formulation implementation for NGSolve
- Advanced adaptive mesh refinement with error estimation
- Halbach array geometry tools for Radia
- Optimization and parameter sweep tools
- Parallel solver support

## Support and Documentation

### Radia MCP Server
- README: `S:\Radia\01_Github\mcp_server_radia\README.md`
- Tool documentation: See individual tool files in `tools/`

### NGSolve MCP Server
- README: `S:\NGSolve\01_Github\mcp_server_ngsolve\README.md`
- Kelvin transformation examples: `S:\NGSolve\NGSolve\2025_12_14_Kelvin変換\`

### Shared Workspace
- Module: `S:\Radia\01_Github\mcp_shared\workspace.py`

## Version History

- **2026-02-12**: Initial deployment
  - Radia MCP Server: 26 tools
  - NGSolve MCP Server: 15 tools (including enhanced Kelvin transformation)
  - Shared workspace for server communication

## Authors

Based on existing Radia and NGSolve implementations, wrapped for MCP by Claude Code.
