# Radia Background Field Docs

This directory contains the result-saved docs layer and notebook-coupled helper
scripts demonstrating how to use Python callback functions as background
magnetic fields in Radia simulations using `rad.ObjBckg()`.

## Overview

Radia's `ObjBckg` function allows you to define arbitrary background magnetic fields using Python callback functions. This enables:
- Integration of analytically defined fields (quadrupole, sextupole, etc.)
- Coupling with external field solvers
- Custom field distributions for specific applications

## Files

### Mesh Generation

#### **cubit_to_nastran.py**
   - Generates high-quality tetrahedral mesh of sphere using Cubit
   - Exports to Nastran .bdf format
   - Sphere radius: 10 mm, element size: ~1 mm
   - Uses tetrahedral mesh (always convex, required for Radia)
   - Requires Coreform Cubit installation

### Example Scripts

#### 1. **quadrupole_analytical.py**
   - Simple quadrupole background field example with nonlinear material (MatSatIsoFrm)
   - Tests B→H conversion and solver convergence with background fields
   - Verifies B/H = μ₀ at multiple far-field points

#### 2. **sphere_in_quadrupole.py**
   - Analytical solution comparison for magnetizable cube in quadrupole field
   - Uses nonlinear material (MatSatIsoFrm)
   - Evaluates 11 test points at distances 20-100 mm
   - Includes error statistics grouped by distance

#### 3. **permeability_comparison.py**
   - Compares accuracy across different permeability values
   - Tests with μᵣ = 10, 100, 1000 using linear material (MatLin)
   - 11 test points per permeability value
   - Demonstrates accuracy across a wide permeability range


## Quick Start

### Using Callback Function for Background Field

```python
import radia as rd
import numpy as np

# Radia always uses meters

# Define background field function
def quadrupole_field(pos):
	"""
	Quadrupole field: B = g*(y*ex + x*ey)

	Args:
		pos: [x, y, z] in meters (Radia always uses meters)

	Returns:
		[Bx, By, Bz] in Tesla
	"""
	x, y, z = pos
	g = 10.0  # Gradient in T/m
	Bx = g * y
	By = g * x
	Bz = 0.0
	return [Bx, By, Bz]

# Create background field source
background = rd.ObjBckg(quadrupole_field)

# Create magnetizable object using ObjHexahedron
mm = 1e-3
half = 5 * mm
vertices = [[-half,-half,-half], [half,-half,-half], [half,half,-half], [-half,half,-half],
            [-half,-half,half], [half,-half,half], [half,half,half], [-half,half,half]]
cube = rd.ObjHexahedron(vertices, [0, 0, 0])

# Apply linear isotropic material (mu_r = 1000)
mat = rd.MatLin(1000)
rd.MatApl(cube, mat)

# Combine with background field
system = rd.ObjCnt([cube, background])

# Solve
rd.Solve(system, 0.0001, 10000)

# Evaluate total field (object + background)
B_total = rd.Fld(system, 'b', [20*mm, 0, 0])
```

## Background Field Function Requirements

### Function Signature

```python
def my_field(pos):
	"""
	Args:
		pos: [x, y, z] in meters (Radia always uses meters)

	Returns:
		[Bx, By, Bz] in Tesla
	"""
	x, y, z = pos
	# ... compute field ...
	return [Bx, By, Bz]
```

### Important Notes

1. **Units**:
   - Input position units are always in meters (Radia always uses meters)
   - Output: **Magnetic flux density B in Tesla**
   - Internal conversion: Radia automatically converts B→H using H = B/μ₀

2. **Return Type**:
   - Must return a list or tuple of 3 numbers: `[Bx, By, Bz]`
   - Alternative: Return dict `{'B': [Bx, By, Bz], 'A': [Ax, Ay, Az]}` for both B and vector potential A

3. **Physical Quantities**:
   - **B field (Tesla)**: Magnetic flux density - what the callback returns
   - **H field (A/m)**: Magnetic field intensity - automatically computed as H = B/μ₀
   - **μ₀ = 1.25663706212×10⁻⁶ T/(A/m)**: Permeability of free space
   - Conversion factor: 1/μ₀ = 795774.715459 (A/m)/T

4. **Thread Safety**:
   - Function will be called multiple times during field computation
   - Should be stateless or thread-safe

## Common Background Field Types

All examples below use positions in meters (Radia always uses meters).

### Uniform Field

```python
def uniform_field(pos):
	return [0.0, 1.0, 0.0]  # 1 T in Y direction
```

### Gradient Field (Dipole-like)

```python
def gradient_field(pos):
	x, y, z = pos  # meters
	g = 10.0  # T/m
	# div(B) = ∂(gz)/∂x + ∂(0)/∂y + ∂(-gx)/∂z = 0 + 0 + 0 = 0
	return [g * z, 0.0, -g * x]
```

### Quadrupole Field

```python
def quadrupole_field(pos):
	x, y, z = pos  # meters
	g = 10.0  # T/m
	return [g * y, g * x, 0.0]
```

### Sextupole Field

```python
def sextupole_field(pos):
	x, y, z = pos  # meters
	k = 100.0  # T/m^2
	Bx = k * x * y
	By = k / 2.0 * (x**2 - y**2)
	return [Bx, By, 0.0]
```

## Limitations and Notes

1. **Binary Serialization**: `rd.DumpBin`/`rd.Parse` not supported for callback functions
2. **B→H Conversion**: Callback returns B (Tesla), Radia automatically converts to H = B/μ₀
   - Verified working for standalone background field sources
   - Test scripts validate this conversion
3. **Stray Fields**: A magnetized object produces stray fields that extend beyond its surface.
   When comparing with analytical background fields, test points should be placed far from
   the object to minimize the stray field contribution.

## Comparison with NGSolve Integration

| Feature | Background Field (this folder) | NGSolve Integration |
|---------|-------------------------------|---------------------|
| Direction | Python → Radia | Radia → NGSolve |
| Use Case | Add external fields to Radia | Use Radia fields in NGSolve FEM |
| Function | `rd.ObjBckg()` | `rad.RadiaField()` |
| Input | Python callback | Radia object |
| Output | Radia field source | NGSolve CoefficientFunction |
| Location | `docs/background_fields/` | `examples/ngsolve_integration/` |

## Requirements

- Python 3.8+
- Radia with CoefficientFunction support
- NumPy
- Cubit (optional, for mesh generation)

## References

- Result-saved notebook: `docs/background_fields/background_fields.ipynb`
- Synchronized result JSON: `docs/background_fields/background_fields_results.json`
- Radia to NGSolve examples: `examples/ngsolve_integration/`
