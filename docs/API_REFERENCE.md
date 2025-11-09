# Radia Python API Reference

Complete reference for Radia Python API including original ESRF functions and custom extensions.

**Original ESRF Documentation**: https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html

**Date**: 2025-11-08

---

## Table of Contents

- [Geometry Objects](#geometry-objects)
  - [Magnetized Objects](#magnetized-objects)
  - [Current-Carrying Objects](#current-carrying-objects)
  - [Field Source Objects](#field-source-objects)
  - [Object Containers](#object-containers)
- [Materials](#materials)
- [Transformations](#transformations)
- [Relaxation and Solving](#relaxation-and-solving)
- [Field Computation](#field-computation)
- [Utilities](#utilities)
- [Extensions](#extensions)

---

## Geometry Objects

### Magnetized Objects

#### ObjRecMag
```python
obj = rad.ObjRecMag(center, dimensions, magnetization)
```
Creates a rectangular parallelepiped with uniform magnetization.

**Parameters**:
- `center`: `[x, y, z]` - Center point (mm)
- `dimensions`: `[dx, dy, dz]` - Dimensions (mm)
- `magnetization`: `[Mx, My, Mz]` - Magnetization vector (T)

**Example**:
```python
# 10×10×10 mm cube centered at origin, magnetized in z-direction
magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1.2])
```

---

#### ObjThckPgn
```python
obj = rad.ObjThckPgn(x, dx, points, axis, magnetization)
```
Creates an extruded polygon (prism) with uniform magnetization.

**Parameters**:
- `x`: Coordinate along extrusion axis (mm)
- `dx`: Extrusion length (mm)
- `points`: `[[x1,y1], [x2,y2], ...]` - Polygon vertices in 2D
- `axis`: `'x'|'y'|'z'` - Extrusion axis direction
- `magnetization`: `[Mx, My, Mz]` - Magnetization vector (T)

**Example**:
```python
# Triangular prism
points = [[0,0], [10,0], [5,8.66]]
prism = rad.ObjThckPgn(0, 20, points, 'z', [0, 0, 1])
```

---

#### ObjMltExtRtg
```python
obj = rad.ObjMltExtRtg([[x1,dx1], [x2,dx2], ...])
```
Creates a convex polyhedron by recursive rectangular extrusion.

**Parameters**:
- List of `[position, thickness]` pairs for nested extrusions

**Example**:
```python
# Nested rectangular structure
obj = rad.ObjMltExtRtg([[0,10], [0,10], [0,10]])
```

---

#### ObjMltExtTri
```python
obj = rad.ObjMltExtTri([[x1,dx1], [x2,dx2], ...])
```
Creates a polyhedron by recursive triangular extrusion.

---

### Current-Carrying Objects

#### ObjRaceTrk
```python
coil = rad.ObjRaceTrk(center, radii, heights, current, n_segments)
```
Creates a racetrack coil (rectangular with rounded ends).

**Parameters**:
- `center`: `[x, y, z]` - Center point (mm)
- `radii`: `[r_straight, r_arc]` - Radii (mm)
- `heights`: `[h_inner, h_outer]` - Heights (mm)
- `current`: Current (A)
- `n_segments`: Number of segments for arc discretization

**Example**:
```python
# Racetrack coil: 30mm straight, 20mm arc radius, 1000A current
coil = rad.ObjRaceTrk([0,0,20], [30,20], [5,5], 1000, 3)
```

---

#### ObjFlmCur
```python
filament = rad.ObjFlmCur([[x1,y1,z1], [x2,y2,z2], ...], current)
```
Creates a filament polygonal line conductor.

**Parameters**:
- `points`: List of 3D points defining the path (mm)
- `current`: Current (A)

**Example**:
```python
# Rectangular loop
points = [[0,0,0], [10,0,0], [10,10,0], [0,10,0], [0,0,0]]
loop = rad.ObjFlmCur(points, 100)
```

---

### Field Source Objects

#### ObjBckg ⭐ EXTENDED
```python
field_src = rad.ObjBckg([Bx, By, Bz])
```
Creates a source of uniform background magnetic field.

**Parameters**:
- `[Bx, By, Bz]`: Uniform field vector (T)

**Example**:
```python
# 0.5 T uniform field in z-direction
bg = rad.ObjBckg([0, 0, 0.5])
```

**Documentation**: See [API_EXTENSIONS.md](API_EXTENSIONS.md#objbckg)

---

#### ObjBckgCF ⭐ EXTENDED
```python
field_src = rad.ObjBckgCF(callback_function)
```
Creates a source of arbitrary spatially-varying background field.

**Parameters**:
- `callback_function`: Python function `f([x,y,z]) -> [Bx,By,Bz]`
  - Input: Position in mm
  - Output: Field in T

**Example**:
```python
def gradient_field(pos):
	x, y, z = pos
	return [0.01*x, 0.01*y, 0.5]

bg = rad.ObjBckgCF(gradient_field)
```

**Documentation**: See [API_EXTENSIONS.md](API_EXTENSIONS.md#objbckgcf)

---

### Object Containers

#### ObjCnt
```python
group = rad.ObjCnt([obj1, obj2, obj3, ...])
```
Creates a container (group) of objects.

**Parameters**:
- List of object keys

**Returns**:
- Container object key

**Example**:
```python
magnet1 = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])
magnet2 = rad.ObjRecMag([20,0,0], [10,10,10], [0,0,1])
group = rad.ObjCnt([magnet1, magnet2])
```

---

#### ObjCutMag
```python
pieces = rad.ObjCutMag(obj, point, normal, option)
```
Cuts a magnetic object with a plane.

**Parameters**:
- `obj`: Object key to cut
- `point`: `[x, y, z]` - Point on cutting plane (mm)
- `normal`: `[nx, ny, nz]` - Plane normal vector
- `option`: `'Frame->Lab'|'Frame->Loc'` - Coordinate frame

**Returns**:
- List of resulting object keys

---

#### ObjDivMag
```python
pieces = rad.ObjDivMag(obj, subdivisions, option)
```
Subdivides a magnetic object.

**Parameters**:
- `obj`: Object key to subdivide
- `subdivisions`: `[nx, ny, nz]` or `[[kx1,q1],[ky1,q2],[kz1,q3]]`
  - Simple form: Number of divisions per axis
  - Advanced form: `[ki, qi]` = ki subdivisions with ratio qi
- `option`: Options string (e.g., `'Frame->Lab'`)

**Returns**:
- Container of subdivided elements

**Example**:
```python
# Uniform 3×3×3 subdivision
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])
pieces = rad.ObjDivMag(magnet, [3, 3, 3])

# Non-uniform subdivision with ratio
pieces = rad.ObjDivMag(magnet, [[3,1.5], [3,1.5], [3,1.5]])
```

---

## Materials

### MatLin
```python
material = rad.MatLin([ksi_parallel, ksi_perpendicular], remanence)
```
Creates a linear anisotropic magnetic material.

**Parameters**:
- `[ksi_parallel, ksi_perpendicular]`: Susceptibilities parallel/perpendicular to easy axis
- `remanence`: Scalar magnitude or `[Mx, My, Mz]` vector (T)
  - Scalar: Easy axis set by object's magnetization direction
  - Vector: Easy axis set explicitly

**Example**:
```python
# Isotropic linear, no remanence (soft iron)
mat = rad.MatLin([999, 999], [0, 0, 0])

# Anisotropic with remanence in z
mat = rad.MatLin([1000, 100], [0, 0, 1.2])
```

---

### MatSatIsoFrm
```python
material = rad.MatSatIsoFrm([ksi1,ms1], [ksi2,ms2], [ksi3,ms3])
```
Creates nonlinear isotropic material from formula.

**Formula**: M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3)

**Parameters**:
- `[ksi_i, ms_i]`: Susceptibility and saturation magnetization pairs (T)
- Terms can be omitted (2nd and 3rd parameters optional)

---

### MatSatIsoTab
```python
material = rad.MatSatIsoTab([[H1,M1], [H2,M2], ...])
```
Creates nonlinear isotropic material from M-H table.

**Parameters**:
- List of `[H, M]` pairs defining M(H) curve (T)

**Example**:
```python
# Soft iron M-H curve
MH_data = [
	[0, 0], [200, 0.7], [600, 1.2], [1200, 1.4],
	[2000, 1.5], [3500, 1.54], [6000, 1.56], [12000, 1.57]
]
mat = rad.MatSatIsoTab(MH_data)
```

---

### MatSatLamFrm / MatSatLamTab
Laminated materials (similar to `MatSatIsoFrm/Tab` with packing factor and lamination normal).

---

### MatSatAniso
```python
material = rad.MatSatAniso(data_parallel, data_perpendicular)
```
Creates nonlinear anisotropic material with separate parallel/perpendicular characteristics.

---

### Common Material Configuration Examples

Below are typical material parameters for common magnetic materials. These can be created using `MatLin` for permanent magnets or `MatSatIsoFrm` for soft magnetic materials.

#### Permanent Magnet Materials (Linear Anisotropic)

**NdFeB** (Neodymium-Iron-Boron):
```python
mat = rad.MatLin([0.06, 0.17], [0, 0, 1.2])  # ksi_par=0.06, ksi_perp=0.17, Mr=1.2 T
```

**SmCo5** (Samarium-Cobalt):
```python
mat = rad.MatLin([0.005, 0.04], [0, 0, 0.85])  # Mr=0.85 T
```

**Sm2Co17** (Samarium-Cobalt 2:17):
```python
mat = rad.MatLin([0.005, 0.04], [0, 0, 1.05])  # Mr=1.05 T
```

**Ferrite**:
```python
mat = rad.MatLin([0.07, 0.2], [0, 0, 0.35])  # Mr=0.35 T
```

#### Soft Magnetic Materials (Nonlinear Isotropic)

**Xc06** (Low Carbon Steel, C<0.06%):
```python
mat = rad.MatSatIsoFrm([2118., 1.362], [63.06, 0.2605], [17.138, 0.4917])
```

**Steel37** (C<0.13%):
```python
mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
```

**Steel42** (C<0.19%):
```python
mat = rad.MatSatIsoFrm([968.66, 1.441], [24.65, 0.2912], [8.3, 0.3316])
```

**AFK502** (Vanadium Permendur, Fe:49%, Co:49%, V:2%):
```python
mat = rad.MatSatIsoFrm([10485., 1.788], [241.5, 0.437], [7.43, 0.115])
```

**AFK1** (FeCo alloy, Fe:74.2%, Co:25%):
```python
mat = rad.MatSatIsoFrm([2001., 1.704], [38.56, 0.493], [1.24, 0.152])
```

**Note**: For soft magnetic materials, the magnetization is modeled as:
```
M(H) = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3)
```
where H is the magnitude of magnetic field strength in Tesla.

---

### MatApl
```python
rad.MatApl(obj, material)
```
Applies material to object.

**Parameters**:
- `obj`: Object key
- `material`: Material key

**Example**:
```python
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,0.1])
mat = rad.MatSatIsoTab(MH_data)
rad.MatApl(magnet, mat)
```

---

### MatMvsH
```python
M = rad.MatMvsH(obj, 'mx|my|mz'|'', [Hx, Hy, Hz])
```
Computes magnetization from field strength for object's material.

**Parameters**:
- `obj`: Object with material
- Component selector: `'mx'`, `'my'`, `'mz'`, or `''` (magnitude)
- `[Hx, Hy, Hz]`: Field strength vector (T)

**Returns**:
- Magnetization component or magnitude (T)

---

## Transformations

### TrfTrsl
```python
obj_copy = rad.TrfTrsl(obj, [dx, dy, dz])
```
Translates (moves) an object.

**Parameters**:
- `obj`: Object key
- `[dx, dy, dz]`: Translation vector (mm)

**Returns**:
- Transformed object key (same as input)

---

### TrfRot
```python
obj_copy = rad.TrfRot(obj, [x, y, z], [nx, ny, nz], angle)
```
Rotates an object around an axis.

**Parameters**:
- `obj`: Object key
- `[x, y, z]`: Point on rotation axis (mm)
- `[nx, ny, nz]`: Axis direction vector
- `angle`: Rotation angle (radians)

---

### TrfOrnt
```python
obj_copy = rad.TrfOrnt(obj, [x, y, z], [nz_old, nx_old], [nz_new, nx_new])
```
Reorients an object by specifying old and new coordinate frame directions.

---

### TrfMlt
```python
array = rad.TrfMlt(obj, transformation, n_copies)
```
Creates multiple copies with repeated transformation.

**Parameters**:
- `obj`: Object key to copy
- `transformation`: Transformation key (from TrfTrsl, TrfRot, etc.)
- `n_copies`: Number of copies

**Returns**:
- Container of copies

**Example**:
```python
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])
tr = rad.TrfTrsl(rad.ObjCnt([]), [15, 0, 0])  # Translation
array = rad.TrfMlt(magnet, tr, 5)  # 5 copies, 15mm apart
```

---

### TrfZerPara / TrfZerPerp
```python
mirror = rad.TrfZerPara(obj, [x,y,z], [nx,ny,nz])
mirror = rad.TrfZerPerp(obj, [x,y,z], [nx,ny,nz])
```
Creates mirror symmetry with field perpendicular/parallel to mirror plane.

---

## Relaxation and Solving

### RlxPre ⭐ EXTENDED
```python
intrc = rad.RlxPre(obj, srcobj=0)
```
Builds interaction matrix for relaxation.

**Parameters**:
- `obj`: Main object (magnetizable)
- `srcobj`: Optional external field source (default=0)

**Returns**:
- Interaction matrix key

**Time Complexity**: O(N²) where N = number of elements

**Example**:
```python
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,0.1])
rad.MatApl(magnet, material)

# Self-interaction only
intrc = rad.RlxPre(magnet, magnet)

# With external field source
coil = rad.ObjRaceTrk([0,0,20], [30,20], [5,5], 1000, 3)
intrc = rad.RlxPre(magnet, coil)
```

**Documentation**: See [API_EXTENSIONS.md](API_EXTENSIONS.md#rlxpre)

---

### SetRelaxSubInterval ⭐ NEW
```python
rad.SetRelaxSubInterval(intrc, start_idx, end_idx, relax_together=1)
```
Configures element grouping for LU decomposition solver.

**Parameters**:
- `intrc`: Interaction matrix key
- `start_idx`: Starting element index (0-based)
- `end_idx`: Ending element index (0-based, inclusive)
- `relax_together`: `1` = LU decomposition, `0` = Gauss-Seidel

**Example**:
```python
intrc = rad.RlxPre(grp, grp)
rad.SetRelaxSubInterval(intrc, 0, N-1, 1)  # All elements with LU
rad.RlxMan(intrc, 5, 100, 1.0)  # Method 5
```

**Documentation**: See [API_EXTENSIONS.md](API_EXTENSIONS.md#setrelaxsubinterval)

---

### RlxMan ⭐ EXTENDED
```python
rad.RlxMan(intrc, method, iter_num, relax_param)
```
Executes manual relaxation procedure.

**Parameters**:
- `intrc`: Interaction matrix key
- `method`: Solver method (0-5) ⭐ **Extended to support Method 5**
  - 0: Simple iteration
  - 1: Over-relaxation
  - 2: Under-relaxation
  - 3: Gauss-Seidel
  - 4: Gauss-Seidel
  - **5: LU decomposition** ⭐ NEW (requires `SetRelaxSubInterval`)
- `iter_num`: Number of iterations
- `relax_param`: Relaxation parameter (typically 1.0)

**Time Complexity**:
- Methods 0-4: O(N²) per iteration
- Method 5: O(N³) for LU decomposition

**Example**:
```python
# Gauss-Seidel
rad.RlxMan(intrc, 4, 100, 1.0)

# LU decomposition
rad.SetRelaxSubInterval(intrc, 0, N-1, 1)
rad.RlxMan(intrc, 5, 100, 1.0)
```

**Documentation**: See [API_EXTENSIONS.md](API_EXTENSIONS.md#rlxman---method-5-support)

---

### RlxAuto
```python
rad.RlxAuto(intrc, precision, max_iter, method=4, 'ZeroM->True|False')
```
Executes automatic relaxation with convergence criterion.

**Parameters**:
- `intrc`: Interaction matrix key
- `precision`: Convergence threshold (change in M)
- `max_iter`: Maximum iterations
- `method`: Solver method (0-4, default=4)
- `'ZeroM->True'`: Start with M=0 (default)
- `'ZeroM->False'`: Start with existing M values

**Example**:
```python
rad.RlxAuto(intrc, 0.001, 1000, 4)
```

---

### RlxUpdSrc
```python
rad.RlxUpdSrc(intrc)
```
Updates external field source data without rebuilding interaction matrix.

**Use case**: When coil currents change but geometry is fixed.

---

### Solve
```python
rad.Solve(obj, precision, max_iter, method=4)
```
High-level function combining matrix build and relaxation.

**Parameters**:
- `obj`: Object key
- `precision`: Convergence threshold
- `max_iter`: Maximum iterations
- `method`: Solver method (0-4, default=4)

**Equivalent to**:
```python
intrc = rad.RlxPre(obj, obj)
rad.RlxAuto(intrc, precision, max_iter, method)
```

---

## Field Computation

### Fld
```python
field = rad.Fld(obj, component, point_or_points)
```
Computes magnetic field at point(s).

**Parameters**:
- `obj`: Object key
- `component`: Field component selector
  - Magnetic flux density: `'bx'`, `'by'`, `'bz'`, `'b'` (magnitude)
  - Magnetic field: `'hx'`, `'hy'`, `'hz'`, `'h'` (magnitude)
  - Vector potential: `'ax'`, `'ay'`, `'az'`, `'a'` (magnitude)
  - Magnetization: `'mx'`, `'my'`, `'mz'`, `'m'` (magnitude)
- `point_or_points`: `[x,y,z]` or `[[x1,y1,z1], [x2,y2,z2], ...]` (mm)

**Returns**:
- Field value (T) or list of values

**Example**:
```python
# Single point
Bz = rad.Fld(magnet, 'bz', [0, 0, 10])

# Multiple points
points = [[0,0,0], [5,0,0], [10,0,0]]
B_values = rad.Fld(magnet, 'bz', points)
```

---

### FldLst
```python
field_list = rad.FldLst(obj, component, p1, p2, n_points, 'arg|noarg', start)
```
Computes field along a line segment.

**Parameters**:
- `obj`: Object key
- `component`: Field component selector
- `p1`, `p2`: `[x,y,z]` - Start and end points (mm)
- `n_points`: Number of points
- `'arg'|'noarg'`: Include longitudinal position in output
- `start`: Start value for longitudinal position (optional)

**Returns**:
- List of field values, or list of `[position, field]` pairs

---

### FldInt
```python
integral = rad.FldInt(obj, 'inf|fin', component, p1, p2)
```
Computes field integral along a line.

**Parameters**:
- `obj`: Object key
- `'inf'|'fin'`: Infinite or finite integral
- `component`: `'ibx'|'iby'|'ibz'` - Integral component
- `p1`, `p2`: `[x,y,z]` - Start and end points (mm)

**Returns**:
- Field integral (T·mm)

---

### FldPtcTrj
```python
trajectory = rad.FldPtcTrj(obj, energy, [x0,dxdy0,z0,dzdy0], [y0,y1], n_points)
```
Computes relativistic particle trajectory.

**Parameters**:
- `obj`: Object key (field source)
- `energy`: Particle energy (GeV)
- `[x0, dxdy0, z0, dzdy0]`: Initial transverse coordinates and angles
- `[y0, y1]`: Longitudinal coordinate range (mm)
- `n_points`: Number of integration steps

**Returns**:
- List of trajectory points

---

### FldEnr / FldEnrFrc / FldEnrTrq
```python
energy = rad.FldEnr(obj_dst, obj_src)
force = rad.FldEnrFrc(obj_dst, obj_src, 'fx|fy|fz')
torque = rad.FldEnrTrq(obj_dst, obj_src, 'tx|ty|tz', [x,y,z])
```
Computes energy, force, or torque on obj_dst in field of obj_src.

**Returns**:
- Energy (J), force (N), or torque (N·mm)

---

### FldCmpPrc
```python
rad.FldCmpPrc('PrcB->prb, PrcA->pra, ...')
```
Sets computation precision for various field quantities.

**Options**:
- `PrcB`: Magnetic flux density precision (T)
- `PrcA`: Vector potential precision
- `PrcBInt`: Field integral precision (T·mm)
- `PrcForce`: Force precision (N)
- `PrcTorque`: Torque precision (N·mm)
- `PrcEnergy`: Energy precision (J)

---

### FldUnits
```python
rad.FldUnits()
```
Returns current units in use.

**Default units**:
- Length: millimeters (mm)
- Magnetic flux density: Tesla (T)
- Current: Amperes (A)
- Force: Newtons (N)

---

## Utilities

### UtiDel / UtiDelAll
```python
rad.UtiDel(obj)      # Delete one object
rad.UtiDelAll()      # Delete all objects
```
Removes objects from memory.

---

### UtiVer
```python
version = rad.UtiVer()
```
Returns Radia library version string.

---

### UtiDmp / UtiDmpPrs
```python
# Serialize object to string
data = rad.UtiDmp(obj, 'asc'|'bin')

# Deserialize object from string
obj = rad.UtiDmpPrs(data)
```
Saves/loads objects to/from byte strings.

---

## Extensions

### SolverHMatrixDisable / SolverHMatrixEnable ⭐ NEW
```python
rad.SolverHMatrixDisable()  # Use dense matrix
rad.SolverHMatrixEnable()   # Use H-matrix (default)
```
Controls hierarchical matrix acceleration.

**When to use**:
- Disable: Small problems (N < 1000), benchmarking, debugging
- Enable: Large problems (N > 1000), production

**Documentation**: See [API_EXTENSIONS.md](API_EXTENSIONS.md#solverhmatrixdisableenable)

---

## Units Summary

| Quantity | Units |
|----------|-------|
| Length | millimeters (mm) |
| Magnetic flux density B | Tesla (T) |
| Magnetic field H | A/m |
| Magnetization M | T |
| Current | Amperes (A) |
| Force | Newtons (N) |
| Energy | Joules (J) |
| Torque | Newton·millimeters (N·mm) |

---

## References

1. **Original ESRF Radia**
   - Website: https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia.html
   - Reference Guide: https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html

2. **Extensions Documentation**
   - [API_EXTENSIONS.md](API_EXTENSIONS.md) - Detailed documentation of custom extensions

3. **Examples**
   - [`examples/solver_time_evaluation/README.md`](../examples/solver_time_evaluation/README.md) - Solver benchmarks
   - [`examples/simple_problems/`](../examples/simple_problems/) - Basic examples
   - [`examples/complex_coil_geometry/`](../examples/complex_coil_geometry/) - Complex geometries

4. **Integration**
   - [NGSOLVE_USAGE_GUIDE.md](NGSOLVE_USAGE_GUIDE.md) - NGSolve integration

---

**Last Updated**: 2025-11-08
**Maintained By**: Radia Development Team
**License**: LGPL-2.1 (modifications), BSD-style (original RADIA from ESRF)
