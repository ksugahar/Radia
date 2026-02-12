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

## Tools (26 total)

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

## See Also

- [Dual Server Deployment](../docs/DUAL_SERVER_DEPLOYMENT.md)
- [MCP Server Architecture](../docs/MCP_SERVER_ARCHITECTURE.md)
