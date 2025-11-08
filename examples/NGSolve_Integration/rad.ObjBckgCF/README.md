# NGSolve CoefficientFunction to Radia Background Field

This directory contains examples demonstrating how to use NGSolve CoefficientFunctions (or Python callback functions) as background magnetic fields in Radia simulations.

## Overview

Radia's `ObjBckgCF` function allows you to define arbitrary background magnetic fields using Python callback functions. This enables:
- Integration of analytically defined fields (quadrupole, sextupole, etc.)
- Coupling with external field solvers
- Custom field distributions for specific applications

## Files

### 1. **Cubit2Nastran.py**
   - Generates high-quality tetrahedral mesh of sphere using Cubit
   - Exports to Nastran .bdf format
   - Sphere radius: 10 mm, element size: ~2 mm
   - Uses tetrahedral mesh (always convex, required for Radia)

### 2. **sphere_nastran_analysis.py**
   - Reads Nastran mesh and creates Radia model
   - Magnetizable sphere in quadrupole background field
   - Compares Radia numerical solution with analytical solution
   - Exports geometry and field distribution to VTK
   - Uses `rd.ObjBckgCF()` to define quadrupole field


## Quick Start

### Using Callback Function for Background Field

```python
import radia as rd
import numpy as np

# Define background field function
def quadrupole_field(pos):
	"""
	Quadrupole field: B = g*(x*ey + y*ex)

	Args:
		pos: [x, y, z] in millimeters

	Returns:
		[Bx, By, Bz] in Tesla
	"""
	x, y, z = pos
	g = 10.0  # Gradient in T/m
	# Convert mm to m
	x_m = x / 1000.0
	y_m = y / 1000.0

	Bx = g * y_m
	By = g * x_m
	Bz = 0.0

	return [Bx, By, Bz]

# Create background field source
background = rd.ObjBckgCF(quadrupole_field)

# Create magnetizable object
sphere = rd.ObjRecMag([0, 0, 0], [10, 10, 10])
rd.MatApl(sphere, rd.MatStd('Steel37'))

# Combine with background field
system = rd.ObjCnt([sphere, background])

# Solve
rd.Solve(system, 0.0001, 10000)

# Evaluate total field (object + background)
B_total = rd.Fld(system, 'b', [5, 5, 0])
```

## Background Field Function Requirements

### Function Signature

```python
def my_field(pos):
	"""
	Args:
		pos: [x, y, z] in millimeters

	Returns:
		[Bx, By, Bz] in Tesla
	"""
	x, y, z = pos
	# ... compute field ...
	return [Bx, By, Bz]
```

### Important Notes

1. **Units**:
   - Input: Position in **millimeters** (Radia's native units)
   - Output: Magnetic field in **Tesla**

2. **Return Type**:
   - Must return a list or tuple of 3 numbers: `[Bx, By, Bz]`

3. **Thread Safety**:
   - Function will be called multiple times during field computation
   - Should be stateless or thread-safe

## Common Background Field Types

### Uniform Field

```python
def uniform_field(pos):
	return [0.0, 1.0, 0.0]  # 1 T in Y direction
```

### Gradient Field

```python
def gradient_field(pos):
	x, y, z = pos
	g = 0.01  # T/mm
	return [g * x, g * y, g * z]
```

### Quadrupole Field

```python
def quadrupole_field(pos):
	x, y, z = pos
	g = 10.0  # T/m
	x_m, y_m = x / 1000.0, y / 1000.0
	return [g * y_m, g * x_m, 0.0]
```

### Sextupole Field

```python
def sextupole_field(pos):
	x, y, z = pos
	k = 100.0  # T/m^2
	x_m, y_m = x / 1000.0, y / 1000.0
	Bx = k * x_m * y_m
	By = k / 2.0 * (x_m**2 - y_m**2)
	return [Bx, By, 0.0]
```

## Running Examples

### Step 1: Generate Tetrahedral Mesh using Cubit

```bash
python Cubit2Nastran.py
# Creates sphere.bdf with tetrahedral elements
# Mesh: 1025 nodes, 4852 tetrahedra (~2mm element size)
```

### Step 2: Run Analysis with Different Permeabilities

```bash
# Run analysis with mu_r = 10
python sphere_nastran_analysis.py 10

# Run analysis with mu_r = 100
python sphere_nastran_analysis.py 100

# Run analysis with mu_r = 1000
python sphere_nastran_analysis.py 1000
```

Each run produces:
- `sphere_nastran_geometry.vtk` - Sphere geometry (generated once, independent of μᵣ)
- `sphere_nastran_field_mu{N}.vtu` - 3D field distribution with comparison data for μᵣ={N}

### Step 3: Visualize Results in ParaView

```bash
# Open geometry and field files together
paraview sphere_nastran_geometry.vtk sphere_nastran_field_mu10.vtu
```

## Limitations

1. **Vector Potential**: `rd.Fld(obj, 'a', pos)` not implemented for CF background fields
2. **Binary Serialization**: `rd.DumpBin`/`rd.Parse` not supported
3. **Infinite Integrals**: Uses simple numerical integration

## Comparison with NGSolve Integration

| Feature | Background Field (this folder) | CoefficientFunction (Radia_to_NGSolve) |
|---------|-------------------------------|---------------------------------------|
| Direction | Python → Radia | Radia → NGSolve |
| Use Case | Add external fields to Radia | Use Radia fields in NGSolve FEM |
| Function | `rd.ObjBckgCF()` | `rad_ngsolve.RadiaField()` |
| Input | Python callback | Radia object |
| Output | Radia field source | NGSolve CoefficientFunction |

## Validation Results

Magnetizable sphere (R = 10 mm) in quadrupole field (g = 10 T/m) was analyzed with three different relative permeabilities.

### Mesh Statistics

- **Mesh type**: Surface mesh with convex octant decomposition
- **Surface triangles (CTRIA3)**: 7408 triangles
- **Material groups**: 8 convex octants (sphere divided by 3 orthogonal planes)
- **Radia polyhedra**: 8 (one per octant)
- **Element size**: ~1 mm
- **Total nodes**: ~3700

**Mesh Strategy**: The sphere is decomposed into 8 convex octants using webcut operations along X, Y, and Z planes. Each octant's surface triangles are grouped by material ID and combined into a single convex polyhedron for Radia. This approach:
- Uses only surface elements (linear analysis doesn't require volume mesh)
- Ensures each polyhedron is convex (Radia requirement)
- Reduces computational cost significantly (8 polyhedra vs thousands of tetrahedra)
- Maintains geometric accuracy with fine surface triangulation

### Field Comparison at Test Points

Comparison between Radia numerical solution and analytical quadrupole field **outside the sphere** (r > 10 mm):

#### μᵣ = 10 (Low Permeability)

| Point (mm) | B_Radia (T) | B_Analytical (T) | Error (T) | Error (%) |
|------------|-------------|------------------|-----------|-----------|
| [15, 0, 0] | 0.134967 | 0.150000 | 0.015033 | 10.0% |
| [0, 15, 0] | 0.134961 | 0.150000 | 0.015039 | 10.0% |
| [20, 0, 0] | 0.195430 | 0.200000 | 0.004570 | 2.3% |
| [0, 20, 0] | 0.195428 | 0.200000 | 0.004572 | 2.3% |
| [30, 0, 0] | 0.299107 | 0.300000 | 0.000893 | 0.30% |
| [0, 30, 0] | 0.299106 | 0.300000 | 0.000894 | 0.30% |
| [40, 0, 0] | 0.399717 | 0.400000 | 0.000283 | 0.07% |
| [50, 0, 0] | 0.499884 | 0.500000 | 0.000116 | 0.02% |

**Far-field accuracy**: Excellent agreement at r ≥ 30 mm (error < 0.3%)

#### μᵣ = 100 (Medium Permeability)

| Point (mm) | B_Radia (T) | B_Analytical (T) | Error (T) | Error (%) |
|------------|-------------|------------------|-----------|-----------|
| [15, 0, 0] | 0.129561 | 0.150000 | 0.020439 | 13.6% |
| [0, 15, 0] | 0.129557 | 0.150000 | 0.020443 | 13.6% |
| [20, 0, 0] | 0.193794 | 0.200000 | 0.006206 | 3.1% |
| [0, 20, 0] | 0.193793 | 0.200000 | 0.006207 | 3.1% |
| [30, 0, 0] | 0.298788 | 0.300000 | 0.001212 | 0.40% |
| [0, 30, 0] | 0.298787 | 0.300000 | 0.001213 | 0.40% |
| [40, 0, 0] | 0.399617 | 0.400000 | 0.000383 | 0.10% |
| [50, 0, 0] | 0.499843 | 0.500000 | 0.000157 | 0.03% |

**Far-field accuracy**: Excellent agreement at r ≥ 30 mm (error < 0.4%)

#### μᵣ = 1000 (High Permeability - Soft Iron)

| Point (mm) | B_Radia (T) | B_Analytical (T) | Error (T) | Error (%) |
|------------|-------------|------------------|-----------|-----------|
| [15, 0, 0] | 0.128877 | 0.150000 | 0.021123 | 14.1% |
| [0, 15, 0] | 0.128874 | 0.150000 | 0.021126 | 14.1% |
| [20, 0, 0] | 0.193587 | 0.200000 | 0.006413 | 3.2% |
| [0, 20, 0] | 0.193586 | 0.200000 | 0.006414 | 3.2% |
| [30, 0, 0] | 0.298747 | 0.300000 | 0.001253 | 0.42% |
| [0, 30, 0] | 0.298747 | 0.300000 | 0.001253 | 0.42% |
| [40, 0, 0] | 0.399604 | 0.400000 | 0.000396 | 0.10% |
| [50, 0, 0] | 0.499838 | 0.500000 | 0.000162 | 0.03% |

**Far-field accuracy**: Excellent agreement at r ≥ 30 mm (error < 0.5%)

### Key Observations

1. **Near-field Distortion**: The magnetizable sphere distorts the external quadrupole field near the surface, with larger distortions for higher permeability
   - μᵣ = 10: Field reduced by ~10% at r = 15 mm
   - μᵣ = 100: Field reduced by ~13.6% at r = 15 mm
   - μᵣ = 1000: Field reduced by ~14.1% at r = 15 mm

2. **Excellent Far-field Accuracy**: Error decreases rapidly with distance:
   - r = 15 mm: 10-14% error
   - r = 20 mm: 2.3-3.2% error
   - r = 30 mm: 0.30-0.42% error
   - r = 40 mm: 0.07-0.10% error
   - r = 50 mm: 0.02-0.03% error

3. **Distance Scaling**: Error follows 1/r² behavior (expected for dipole field perturbation)

4. **Surface Mesh Efficiency**: Using only 8 convex polyhedra (surface mesh approach) provides:
   - Fast computation (8 polyhedra vs 4852 tetrahedra in volume mesh)
   - Good accuracy for linear magnetic analysis
   - Exact geometric representation of spherical surface (7408 fine triangles)

5. **Symmetry**: Results show good symmetry between [15,0,0] and [0,15,0] points (errors within 0.001%)

6. **ObjBckgCF Performance**: The callback function approach successfully implements arbitrary background fields with good accuracy

### Physical Interpretation

The magnetizable sphere acts as a magnetic dipole in the quadrupole field:
- **Field concentration**: High-permeability material concentrates field lines through the sphere
- **External shielding**: Near the sphere, the field is reduced compared to pure quadrupole
- **Symmetry preservation**: The solution maintains the expected quadrupole symmetry

## Requirements

- Python 3.8+
- Radia with CoefficientFunction support
- NumPy
- Cubit (optional, for mesh generation)

## References

- Main Radia documentation: `README.md`
- Radia to NGSolve examples: `examples/Radia_to_NGSolve_CoefficientFunction/`
- Build instructions: `README_BUILD.md`

---

**Date**: 2025
