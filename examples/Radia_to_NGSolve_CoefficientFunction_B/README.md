# Radia-NGSolve Integration Examples

This directory contains examples demonstrating the integration of Radia magnetic field calculations with NGSolve finite element analysis.

## Overview

The `rad_ngsolve` module provides a seamless interface between Radia and NGSolve, allowing Radia magnetic fields to be used as NGSolve CoefficientFunctions. This enables:
- Direct evaluation of Radia fields on NGSolve meshes
- Integration of Radia fields into FEM simulations
- Automatic unit conversion (Radia: mm, NGSolve: m)

## Files

### Main Examples

1. **demo_field_types.py**
   - Demonstrates all supported field types
   - Shows how to use the unified `RadiaField` interface
   - Compares NGSolve evaluation with direct Radia calculation
   - **Field types:** 'b' (flux density), 'h' (magnetic field), 'a' (vector potential), 'm' (magnetization)

2. **visualize_field.py**
   - Field visualization and comparison tool
   - Compares CoefficientFunction vs GridFunction evaluation
   - Exports fields to VTK format for Paraview visualization
   - Includes mesh refinement study

2. **export_radia_geometry.py**
   - Exports Radia geometry to VTK format
   - Uses existing `radia_vtk_export` module
   - Automatic mm to m unit conversion

### Output Files

3. **radia_components.vtk**
   - Example VTK output from export_radia_geometry.py
   - Can be opened in Paraview for visualization

4. **radia_field.pvsm**
   - ParaView state file for field visualization
   - Opens pre-configured visualization setup

### Utility Scripts

4. **export_radia_geometry.py**
   - Exports Radia geometry to VTK format
   - Uses existing `radia_vtk_export` module
   - Automatic mm to m unit conversion

### Output Files

5. **radia_components.vtk**
   - Example VTK output from export_radia_geometry.py
   - Can be opened in Paraview for visualization

## Quick Start

### Basic Usage

```python
import radia as rad
import rad_ngsolve
from ngsolve import *

# Create Radia magnet
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

# Create field CoefficientFunctions
B_cf = rad_ngsolve.RadiaField(magnet, 'b')  # Magnetic flux density
H_cf = rad_ngsolve.RadiaField(magnet, 'h')  # Magnetic field
A_cf = rad_ngsolve.RadiaField(magnet, 'a')  # Vector potential
M_cf = rad_ngsolve.RadiaField(magnet, 'm')  # Magnetization

# Use in NGSolve
mesh = Mesh(...)
fes = VectorH1(mesh, order=2)
gf = GridFunction(fes)
gf.Set(B_cf)
```

### Field Types

| Type | Description | Units | Radia Command |
|------|-------------|-------|---------------|
| 'b' | Magnetic flux density | Tesla (T) | rad.Fld(obj, 'b', pt) |
| 'h' | Magnetic field | A/m | rad.Fld(obj, 'h', pt) |
| 'a' | Vector potential | T*m | rad.Fld(obj, 'a', pt) |
| 'm' | Magnetization | A/m | rad.Fld(obj, 'm', pt) |

## Running Examples

```bash
# Field types demonstration
python demo_field_types.py

# Field visualization and comparison
python visualize_field.py --method both --maxh 0.005

# Mesh convergence study
python test_sphere_in_quadrupole.py

# Export geometry to VTK
python export_radia_geometry.py
```

## Unit Conversion

**Important:** Radia uses millimeters, NGSolve uses meters.
- The `rad_ngsolve` module automatically converts coordinates: m -> mm (x1000)
- The `radia_vtk_export` module automatically converts coordinates: mm -> m (/1000)
- Field values remain in SI units (Tesla, A/m, etc.)

## Requirements

- Python 3.8+
- NGSolve 6.2+
- Radia
- NumPy

## Notes

- All examples use the unified `RadiaField` interface
- Legacy interfaces (RadBfield, RadHfield, etc.) have been removed
- For performance-critical applications, use GridFunction interpolation
- For exact values at specific points, use CoefficientFunction direct evaluation

## Visualization

VTK output files can be visualized in Paraview:
1. Open Paraview
2. File -> Open -> select .vtk or .vtu file
3. Click 'Apply' in Properties panel
4. Select field to visualize from dropdown menu

## Author

Claude Code, 2025
