# Tests for Radia Cubit Plugin (cubit_mesh_curver)

## Running Tests

Tests should be run using **system Python** with the `CUBIT_PATH` environment variable set to point to Cubit's bin directory.

### Setting Up Cubit Access

Either add Cubit's `bin` directory to your system PATH, or set the `CUBIT_PATH` environment variable:

```bash
# Windows (PowerShell) -- Coreform Cubit 2025.12+ (PySide6 plugin)
$env:CUBIT_PATH = "C:/Program Files/Coreform Cubit 2025.12/bin"

# Windows (cmd)
set CUBIT_PATH=C:/Program Files/Coreform Cubit 2025.12/bin

# Linux/Mac
export CUBIT_PATH=/path/to/cubit/bin
```

### Using pytest (Recommended)

```bash
python -m pytest tests/
```

### Running Individual Tests

```bash
python tests/cubit/test_gmsh_export.py
python tests/cubit/test_vtk_auto_order.py
```

## Test Files

| File | Description |
|------|-------------|
| `test_basic.py` | Module import and function signature tests |
| `test_gmsh_export.py` | Gmsh v4.1 export tests (format, 1st/2nd order, mixed, $Entities, $PhysicalNames) |
| `test_nastran_export.py` | Nastran BDF export tests |
| `test_meg_export.py` | MEG format export tests |
| `test_vtk_auto_order.py` | VTK auto element order detection tests |
| `test_vtk_node_ordering.py` | VTK node ordering validation tests |
| `test_vtu_export.py` | VTU XML format export tests |
| `test_netgen_export.py` | Netgen mesh export tests |
| `test_netgen_first_order.py` | First-order Netgen mesh tests |
| `test_netgen_mixed_elements.py` | Mixed element type tests |
| `test_netgen_with_ngsolve.py` | NGSolve integration tests (requires ngsolve) |
| `test_geometry_blocks.py` | Geometry-based block export tests |
| `test_block_geometry_api.py` | Cubit API investigation for geometry blocks |
| `test_mixed_element_warning.py` | Mixed element type warning tests |
| `test_setgeominfo.py` | SetGeomInfo API tests |
| `test_setgeominfo_uv.py` | SetGeomInfo UV parameter tests |
| `test_curve_workflow.py` | Cubit-to-NGSolve high-order curving workflow |

## Notes

- Most tests require Coreform Cubit installation (set `CUBIT_PATH` env var)
- NGSolve tests require separate ngsolve installation (`pip install ngsolve`)
- **Important**: When both NGSolve and Cubit are used in the same script, NGSolve must be imported BEFORE adding Cubit to `sys.path` and importing `cubit`. This avoids DLL conflicts on Windows.
- Tests create temporary mesh files that are cleaned up after execution

## Test Structure

Each test file follows this pattern:

```python
import sys
import os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

import cubit_mesh_curver

# Test code here...
```
