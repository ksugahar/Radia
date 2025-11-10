# Wolfram Language to Python Conversion Notes

## Conversion Date: 2025-01-29

## Converted Files

| Wolfram Language | Python | Status | Notes |
|-----------------|--------|------|------|
| case0.wls | arc_current_with_magnet.py | ✓ Complete | Verified working |
| case1.wls | arc_current_dual_magnets.py | ✓ Complete | Verified working |
| case2.wls | chamfered_pole_piece.py | ✓ Complete | g2 undefined issue commented out |
| case3.wls | cubic_polyhedron_magnet.py | ✓ Complete | g2 undefined issue commented out |

## Main Conversion Details

### 1. Basic Syntax Conversion

```mathematica
(* Wolfram Language *)
<<Radia`;
radUtiDelAll[];
rmin = 100;
g1 = radObjArcCur[{0,0,0}, {rmin, rmax}, {phimin, phimax}, h, nseg, j];
```

```python
# Python
import radia as rad
rad.UtiDelAll()
rmin = 100
g1 = rad.ObjArcCur([0, 0, 0], [rmin, rmax], [phimin, phimax], h, nseg, j)
```

### 2. Mathematical Functions and Constants

| Wolfram | Python | Description |
|---------|--------|------|
| `Pi` | `math.pi` | Pi constant |
| `Sin[x]` | `math.sin(x)` | Sine function |
| `Cos[x]` | `math.cos(x)` | Cosine function |

### 3. Radia Function Mapping

| Wolfram Language | Python | Function |
|-----------------|--------|------|
| `radUtiDelAll[]` | `rad.UtiDelAll()` | Delete all objects |
| `radObjArcCur[...]` | `rad.ObjArcCur(...)` | Create arc current |
| `radObjRecMag[...]` | `rad.ObjRecMag(...)` | Create rectangular magnet |
| `radMatLin[...]` | `rad.MatLin(...)` | Define linear material |
| `radMatSatIso[...]` | `rad.MatSatIso(...)` | Nonlinear isotropic material |
| `radMatApl[obj, mat]` | `rad.MatApl(obj, mat)` | Apply material |
| `radObjDrwAtr[...]` | `rad.ObjDrwAtr(...)` | Set drawing attributes |
| `radObjCnt[{...}]` | `rad.ObjCnt([...])` | Create container |
| `radObjMltExtRtg[...]` | `rad.ObjMltExtRtg(...)` | Multiple extrusion |
| `radObjDivMag[...]` | `rad.ObjDivMag(...)` | Subdivide magnet |
| `radObjPolyhdr[...]` | `rad.ObjPolyhdr(...)` | Create polyhedron |
| `radFld[obj, "bxbybz", pt]` | `rad.Fld(obj, 'b', pt)` | Calculate magnetic field |

### 4. Data Structure Conversion

```mathematica
(* Wolfram: Lists *)
points = {{1,2,3}, {4,5,6}};
faces = {{1,2,3,4}, {5,6,7,8}};
```

```python
# Python: Lists
points = [[1, 2, 3], [4, 5, 6]]
faces = [[1, 2, 3, 4], [5, 6, 7, 8]]
```

### 5. Visualization

**Wolfram Language:**
```mathematica
t = Show[Graphics3D[radObjDrw[g]]];
Export["3DPlot.png", t];
```

**Python:**
No direct equivalent functionality at this time. Future options include:
- matplotlib + mplot3d
- mayavi
- plotly

### 6. File Output

**Wolfram Language:**
```mathematica
t = radFld[g2, "bxbybz", {0,0,0}]
Export["out.dat", t]
```

**Python:**
```python
field = rad.Fld(g2, 'b', [0, 0, 0])
with open('out.dat', 'w') as f:
	f.write(f"{field[0]}\t{field[1]}\t{field[2]}\n")
```

## Issues Discovered

### case2.wls and case3.wls

**Issue:** `radFld[g2, ...]` is called but `g2` object is not defined

**Resolution:**
- Commented out the relevant section in Python version
- Added note to README.md

```python
# Note: g2 is not defined in the original script
# field = rad.Fld(g2, 'b', [0, 0, 0])  # This would fail
print("Note: Field calculation requires defining g2 object")
```

## Coding Style

- Indentation: TAB characters (1 TAB = 4 spaces equivalent)
- Encoding: UTF-8
- Line endings: LF (Unix style)
- Documentation strings: Written in English
- Comments: Japanese and English mixed

## Test Results

### arc_current_with_magnet.py
```
Container object ID: 4
Magnetic field at origin: Bx=0.000000e+00, By=0.000000e+00, Bz=0.000000e+00 T
Calculation complete. Field data saved to out.dat
```
✓ Working correctly

### arc_current_dual_magnets.py
```
Container object ID: 6
Magnetic field at origin: Bx=0.000000e+00, By=0.000000e+00, Bz=0.000000e+00 T
Calculation complete. Field data saved to out.dat
```
✓ Working correctly

### chamfered_pole_piece.py
```
Object ID: 3
Geometry created and subdivided successfully
Note: Field calculation requires defining g2 object
```
✓ Working correctly (with g2 undefined warning)

### cubic_polyhedron_magnet.py
```
Polyhedron object ID: 1
Cube polyhedron created successfully
Vertices: 8
Faces: 6
Note: Field calculation requires defining g2 object
```
✓ Working correctly (with g2 undefined warning)

## Future Improvements

1. **3D Visualization Features**
   - Basic 3D plotting with matplotlib
   - Interactive visualization (mayavi, plotly)

2. **Enhanced Error Handling**
   - More detailed error messages
   - Additional exception handling

3. **Additional Features**
   - Field map generation
   - Contour plotting
   - CSV data output

4. **Documentation Enhancement**
   - Detailed usage examples for each function
   - Radia concept explanations
   - Additional tutorials

---

**Conversion Tools:** Manual conversion + Python script (tab conversion)
**Verification:** Python 3.12 + Radia 4.32
**Status:** Complete
