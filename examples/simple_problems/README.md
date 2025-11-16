# Radia Simple Examples

This folder contains basic Radia examples demonstrating fundamental features.

## Example Files

- [arc_current_with_magnet.py](arc_current_with_magnet.py) - Arc current with rectangular magnet
- [arc_current_dual_magnets.py](arc_current_dual_magnets.py) - Arc current with two magnets
- [chamfered_pole_piece.py](chamfered_pole_piece.py) - Multiple extrusion with chamfer
- [cubic_polyhedron_magnet.py](cubic_polyhedron_magnet.py) - Polyhedron (cube) magnet
- [compare_magpylib.py](compare_magpylib.py) - Comparison with magpylib library
- [hmatrix_update_magnetization.py](hmatrix_update_magnetization.py) - H-matrix magnetization update

## How to Run

Each Python file can be run independently:

```bash
cd examples/simple_problems
python arc_current_with_magnet.py
python arc_current_dual_magnets.py
python chamfered_pole_piece.py
python cubic_polyhedron_magnet.py
```

## Example Descriptions

### arc_current_with_magnet.py
- Creates arc current element and rectangular magnet
- Applies linear material properties
- Calculates magnetic field at origin

**Key Functions:**
- `rad.ObjArcCur()` - Create arc current element
- `rad.ObjRecMag()` - Create rectangular magnet
- `rad.MatLin([ksi_par, ksi_perp], [mx,my,mz])` - Define anisotropic linear material
- `rad.MatApl(obj, mat)` - Apply material to object
- `rad.Fld(obj, 'b', [x,y,z])` - Calculate magnetic field

### arc_current_dual_magnets.py
- Arc current with two magnets at different positions
- Manages multiple objects with container
- Applies linear material properties

**Key Functions:**
- `rad.ObjCnt([obj1, obj2, ...])` - Create object container

### chamfered_pole_piece.py
- Complex extrusion shape with chamfer
- Defines multiple cross-sections and extrudes
- Element subdivision

**Key Functions:**
- `rad.ObjMltExtRtg()` - Create multiple extrusion rectangle
- `rad.ObjDivMag()` - Subdivide magnet

### cubic_polyhedron_magnet.py
- Creates polyhedron from vertices and faces
- Cube example

**Key Functions:**
- `rad.ObjPolyhdr()` - Create polyhedron

## Radia Python API Basics

1. **Import module**
   ```python
   import radia as rad
   ```

2. **Create objects**
   ```python
   # Rectangular magnet: position, dimensions, magnetization
   mag = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1.0])

   # Arc current: center, [rmin,rmax], [phimin,phimax], height, segments, current
   arc = rad.ObjArcCur([0,0,0], [100,150], [0, 6.28], 20, 20, 10)
   ```

3. **Material properties**
   ```python
   # Anisotropic linear material: [ksi_parallel, ksi_perpendicular], magnetization vector
   mat = rad.MatLin([0.06, 0.17], [0,0,1])
   rad.MatApl(obj, mat)
   ```

4. **Field calculation**
   ```python
   # Magnetic flux density [T]
   B = rad.Fld(obj, 'b', [x, y, z])

   # Magnetic field strength [A/m]
   H = rad.Fld(obj, 'h', [x, y, z])
   ```

5. **3D Visualization**
   ```python
   # Export to VTK for ParaView visualization
   from radia_vtk_export import exportGeometryToVTK
   exportGeometryToVTK(obj, 'output.vtk')
   ```

## Requirements

- Python 3.12+
- Radia module (installed via pip or built from source)
- NumPy

## Troubleshooting

### ModuleNotFoundError: No module named 'radia'

Install Radia:
```bash
pip install radia
```

Or build from source and install locally:
```bash
cd ../..
python -m pip install .
```

### Field calculation returns zero

- Check that magnetization is set on the object
- Verify material properties are correctly applied
- Ensure object ID is valid

## See Also

- [Radia Official Documentation](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)
- [API Reference](../../docs/API_REFERENCE.md)
- [H-Matrix User Guide](../../docs/HMATRIX_USER_GUIDE.md)
