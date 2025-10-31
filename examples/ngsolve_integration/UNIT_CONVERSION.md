# Unit Conversion in rad_ngsolve

## Overview

**NGSolve** and **Radia** use different length units:
- **NGSolve**: Meters (m)
- **Radia**: Millimeters (mm)

The `rad_ngsolve` module automatically handles this conversion.

## Automatic Conversion

When you evaluate a Radia field in NGSolve, the conversion happens automatically:

```python
# Create Radia magnet (mm coordinates)
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])  # 20x20x30 mm

# Create NGSolve mesh (meter coordinates)
geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05)))  # ±50mm in meters
mesh = Mesh(geo.GenerateMesh(maxh=0.03))

# Create coefficient function (handles conversion automatically)
B_cf = rad_ngsolve.RadBfield(magnet)

# Evaluate at NGSolve point (meters)
mesh_pt = mesh(0, 0, 0.02)  # 20mm in meters
B = B_cf(mesh_pt)  # Automatically converts 0.02m -> 20mm for Radia
```

## Implementation Details

The conversion is implemented in `rad_ngsolve.cpp`:

```cpp
virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
	                 FlatVector<> result) const override
{
	auto pnt = mip.GetPoint();
	int dim = pnt.Size();

	// Convert NGSolve coordinates (m) to Radia coordinates (mm)
	py::list coords;
	coords.append(pnt[0] * 1000.0);  // x: m -> mm
	coords.append((dim >= 2) ? pnt[1] * 1000.0 : 0.0);  // y: m -> mm
	coords.append((dim >= 3) ? pnt[2] * 1000.0 : 0.0);  // z: m -> mm

	// Call Radia with mm coordinates
	py::module_ rad = py::module_::import("radia");
	py::object B_result = rad.attr("Fld")(radia_obj, field_comp, coords);

	// Field values are in Tesla (no conversion needed)
	...
}
```

## Conversion Table

| NGSolve (m) | Radia (mm) | Example |
|-------------|------------|---------|
| 0.001 | 1 | 1mm |
| 0.010 | 10 | 1cm |
| 0.020 | 20 | 2cm |
| 1.000 | 1000 | 1m |

## Field Units

Field values (B, H, A) are in SI units and **do not require conversion**:
- **B-field**: Tesla (T)
- **H-field**: A/m
- **A-field**: T·m

## Verification

Run the unit conversion test:

```bash
python test_units_verification.py
```

Expected output:
```
[PASS] [0, 0, 0] m = [0, 0, 0] mm
  NGSolve CF: Bx=-0.000000, By=-0.000000, Bz=0.949718 T
  Radia:      Bx=-0.000000, By=-0.000000, Bz=0.949718 T
  Error:      3.777423e-17 T
...
[SUCCESS] Unit conversion working correctly!
```

## Important Notes

1. **Mesh Creation**: Create NGSolve meshes in meters:
   ```python
   # Correct: mesh from -50mm to +50mm
   geo.Add(OrthoBrick(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05)))
   ```

2. **Radia Objects**: Create Radia objects in millimeters (Radia's default):
   ```python
   # Correct: 20x20x30 mm magnet
   magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
   ```

3. **Point Evaluation**: Use meter coordinates with NGSolve:
   ```python
   # Correct: evaluate at 20mm = 0.02m
   mesh_pt = mesh(0, 0, 0.02)
   B = B_cf(mesh_pt)
   ```

4. **Automatic**: You don't need to do any manual conversion - `rad_ngsolve` handles it automatically.

## Example: Magnet at Different Scales

```python
# Small magnet (10mm = 0.01m)
magnet_small = rad.ObjRecMag([0, 0, 0], [10, 10, 15], [0, 0, 1.2])
geo_small = CSGeometry()
geo_small.Add(OrthoBrick(Pnt(-0.02, -0.02, -0.02), Pnt(0.02, 0.02, 0.02)))

# Large magnet (100mm = 0.1m)
magnet_large = rad.ObjRecMag([0, 0, 0], [100, 100, 150], [0, 0, 1.2])
geo_large = CSGeometry()
geo_large.Add(OrthoBrick(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2)))
```

Both will work correctly with automatic unit conversion.
