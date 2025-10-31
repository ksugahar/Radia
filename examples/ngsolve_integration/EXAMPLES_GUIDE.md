# rad_ngsolve Examples Guide

Quick reference guide for using the rad_ngsolve C++ module examples.

## Setup

```powershell
# 1. Build radia.pyd
cd S:\radia\01_GitHub
.\Build.ps1

# 2. Build rad_ngsolve.pyd
.\build_ngsolve.ps1

# 3. Verify NGSolve installation
python -c "import ngsolve; print('NGSolve OK')"
```

## Example Files

### 1. simple_example.py ⭐ START HERE
**Purpose:** Simplest possible example - Quick verification

**Run:**
```bash
python simple_example.py
```

**What it does:**
- Step 1: Import ngsolve
- Step 2: Import radia  
- Step 3: Import rad_ngsolve
- Step 4: Create simple magnet
- Step 5: Create CoefficientFunction

**Expected output:**
```
SUCCESS: rad_ngsolve is working correctly!
```

**Runtime:** ~1 second

---

### 2. test_rad_ngsolve_full.py
**Purpose:** Verify rad_ngsolve installation and basic functionality

**Run:**
```bash
python test_rad_ngsolve_full.py
```

**What it tests:**
- Module import
- Radia magnet creation (10mm NdFeB cube)
- RadBfield, RadHfield, RadAfield creation
- Field evaluation on NGSolve mesh
- Direct field calculation

**Expected output:**
```
SUCCESS: Full integration test passed
C++ CoefficientFunction implementation working correctly!
```

---

### 3. example_rad_ngsolve_demo.py
**Purpose:** Complete demonstration of all capabilities

**Run:**
```bash
python example_rad_ngsolve_demo.py
```

**What it demonstrates:**
- Step 1: Radia magnetic structure (20mm NdFeB cube)
- Step 2: NGSolve mesh creation ([-40,40] x [-40,40] mm^2)
- Step 3: CoefficientFunction creation
- Step 4: Field evaluation on mesh
- Step 5: Field integrals (Int |B|^2, Int Bz)
- Step 6: Field sampling along z-axis

**Expected output:**
```
[Step 1] Creating Radia magnetic structure...
  [OK] Field at center: Bz = 0.7843 T

[Step 6] Sampling field along z-axis...
  z [mm]    Bz [T]
	 0.0     0.78431
	10.0     0.51283
	20.0     0.15857
	30.0     0.05336

SUCCESS: Demonstration completed
```

---

### 4. visualize_field.py
**Purpose:** Interactive 3D visualization of Radia fields using Netgen GUI

**Run:**
```bash
# Open Netgen GUI with Radia field (will be zero - known limitation)
python visualize_field.py

# Open Netgen GUI with TEST field (working visualization demo)
python visualize_field.py --test

# Skip GUI, export VTK file only (for automated testing)
python visualize_field.py --no-gui
```

**Important**: Radia fields evaluate to zero on 3D meshes (known limitation).
Use `--test` flag to see working visualization, or use 2D meshes for Radia fields.

**What it does:**
- Step 1: Creates Radia magnet (20mm NdFeB cube)
- Step 2: Creates NGSolve CoefficientFunction
- Step 3: Generates 3D Netgen mesh (100x100x100mm box)
- Step 4: Evaluates field on mesh
- Step 5: Computes field statistics
- Step 6: Opens Netgen GUI with `netgen.gui.StartGUI()`

**Drawing with ngsolve.Draw():**
```python
# Draw mesh
Draw(mesh)

# Draw GridFunction (no mesh parameter needed)
Draw(gf_B_mag)

# Draw CoefficientFunction (requires name parameter)
Draw(B_mag, mesh, "B_magnitude")

# Update display
Redraw()

# Start GUI event loop
import netgen.gui
netgen.gui.StartGUI()
```

**GUI Usage:**
- Netgen GUI window will open showing the mesh and field
- Use the View menu to select visualization options
- Rotate: Left mouse button drag
- Zoom: Mouse wheel
- Close the GUI window to exit the program

**Expected output:**
```
[Drawing] Setting up visualization...
  1. Mesh drawn
  2. Field magnitude (GridFunction) drawn
  3. Field magnitude (CoefficientFunction) drawn
  4. Display updated

======================================================================
Starting Netgen GUI...
(Program will wait until you close the GUI window)
======================================================================
```

**Fallback:**
- If GUI fails, automatically exports to VTK file (radia_field.vtk)
- Use `--no-gui` flag to skip GUI and export VTK only

**Note:**
- Field values on 3D mesh evaluate to zero (known limitation)
- Visualization still displays mesh geometry correctly
- For accurate field values, use 2D meshes (see example_rad_ngsolve_demo.py)

---

## Common Usage Pattern

```python
import sys
sys.path.insert(0, r'S:\radia\01_GitHub\build\Release')
sys.path.insert(0, r'S:\radia\01_GitHub\dist')

# STEP 1: Import ngsolve FIRST
import ngsolve
from ngsolve import Mesh, H1, GridFunction

# STEP 2: Import radia and rad_ngsolve
import radia as rad
import rad_ngsolve

# STEP 3: Create Radia geometry
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 20], [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

# STEP 4: Create CoefficientFunction
B_field = rad_ngsolve.RadBfield(magnet)

# STEP 5: Use in NGSolve
mesh = Mesh(...)
fes = H1(mesh, order=2, dim=3)
gf = GridFunction(fes)
gf.Set(B_field)  # Evaluate on mesh

# STEP 6: Compute integrals
from ngsolve import Integrate, dx
integral = Integrate(B_field[2], mesh)  # Integrate Bz
```

## Quick Reference

### Available Functions
- `rad_ngsolve.RadBfield(obj)` - Magnetic flux density [T]
- `rad_ngsolve.RadHfield(obj)` - Magnetic field intensity [A/m]
- `rad_ngsolve.RadAfield(obj)` - Vector potential [T·mm]

### Component Access
```python
B = rad_ngsolve.RadBfield(magnet)
Bx = B[0]  # x-component
By = B[1]  # y-component
Bz = B[2]  # z-component
B_mag_sq = Bx**2 + By**2 + Bz**2
```

### Integration
```python
from ngsolve import Integrate, dx

# Integrate field components
Bz_integral = Integrate(B[2], mesh)

# Integrate field magnitude
B_squared = B[0]**2 + B[1]**2 + B[2]**2
energy = Integrate(B_squared / (2*mu_0), mesh)
```

## Troubleshooting

### ImportError: DLL load failed
**Problem:** rad_ngsolve cannot load libngsolve.dll

**Solution:**
```python
# Import ngsolve FIRST
import ngsolve  # This loads libngsolve.dll
import rad_ngsolve  # Now this works
```

### ModuleNotFoundError: No module named 'rad_ngsolve'
**Problem:** rad_ngsolve.pyd not built

**Solution:**
```powershell
cd S:\radia\01_GitHub
.\build_ngsolve.ps1
```

### Field values are zero
**Problem:** Mesh or evaluation points in wrong location

**Check:**
1. Radia magnet is solved: `rad.Solve(magnet, 0.0001, 10000)`
2. Mesh overlaps with field region
3. Coordinate units match (Radia uses mm)

## Testing

Run comprehensive tests:
```bash
# In project root
python tests/test_rad_ngsolve.py

# Or with pytest
pytest tests/test_rad_ngsolve.py -v
```

## Performance Notes

The C++ implementation provides:
- **Zero Python overhead**: Direct C API calls
- **Native NGSolve integration**: True CoefficientFunction
- **Thread-safe**: Can be used in parallel assembly
- **Fast**: No Python callback per evaluation point

Typical performance:
- Field evaluation: ~1-10 μs per point
- Mesh evaluation (1000 DOFs): ~10-100 ms
- Integration: Native NGSolve speed

## Next Steps

1. Run `test_rad_ngsolve_full.py` to verify installation
2. Run `example_rad_ngsolve_demo.py` to see capabilities
3. Modify examples for your specific geometry
4. Check `../../tests/test_rad_ngsolve.py` for more examples
5. Read `../../RAD_NGSOLVE_BUILD_SUCCESS.md` for complete documentation

## Support

- **Build issues**: Check `../../RAD_NGSOLVE_BUILD_SUCCESS.md`
- **API reference**: See `../../src/python/rad_ngsolve.cpp`
- **Test examples**: See `../../tests/test_rad_ngsolve.py`
