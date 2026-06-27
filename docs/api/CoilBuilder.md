# CoilBuilder API Reference

CoilBuilder provides a fluent interface for constructing multi-segment coil geometries. It tracks position and orientation automatically, so you only describe the coil path.

**Module:** `src/radia/coil_builder.py`

## Quick Start

```python
from coil_builder import CoilBuilder

mm = 1e-3  # Radia always uses meters

# Simple racetrack coil
coil = (CoilBuilder(current=1000)
    .set_start([0, -25*mm, 0])
    .set_cross_section(width=10*mm, height=20*mm)
    .add_straight(50*mm)
    .add_arc(radius=15*mm, arc_angle=180)
    .add_straight(50*mm)
    .add_arc(radius=15*mm, arc_angle=180))

# Use it
radia_objects = coil.to_radia()       # Radia ObjRecCur/ObjArcCur list
coil.write_step("coil.step")          # STEP export for visualization
segments, I = coil.to_wire_segments() # Wire segments for Biot-Savart
```

## Constructor

```python
CoilBuilder(current)
```

| Arg | Type | Description |
|-----|------|-------------|
| `current` | float | Coil current in Amperes |

## Builder Methods

All builder methods return `self` for method chaining.

### set_start(position, orientation=None)

Set starting position and optional orientation.

| Arg | Type | Default | Description |
|-----|------|---------|-------------|
| `position` | list[3] | required | Starting point [x, y, z] in meters |
| `orientation` | ndarray(3,3) | identity | Row-vector orientation matrix |

### set_cross_section(width, height)

Set rectangular conductor cross-section.

| Arg | Type | Description |
|-----|------|-------------|
| `width` | float | Radial width in meters |
| `height` | float | Axial height in meters |

### add_straight(length, tilt=0)

Add a straight conductor segment.

| Arg | Type | Default | Description |
|-----|------|---------|-------------|
| `length` | float | required | Segment length in meters |
| `tilt` | float | 0 | Tilt angle in degrees (rotation around current direction) |

### add_arc(radius, arc_angle, tilt=0)

Add a toroidal arc segment.

| Arg | Type | Default | Description |
|-----|------|---------|-------------|
| `radius` | float | required | Centerline arc radius in meters |
| `arc_angle` | float | required | Arc angle in degrees |
| `tilt` | float | 0 | Tilt angle in degrees |

## Output Methods

### to_radia()

Convert all segments to Radia objects (`ObjRecCur`, `ObjArcCur`) with transformations.

Returns: `list[int]` — Radia object handles. Pass to `rad.ObjCnt(...)`.

```python
coil_objects = coil.to_radia()
container = rad.ObjCnt(coil_objects)
B = rad.Fld(container, 'b', [0, 0, 0.1])
```

### write_step(filename)

Export coil geometry to STEP file. Requires `netgen.occ`.

```python
coil.write_step("my_coil.step")
```

### to_occ()

Convert to `netgen.occ` shape for programmatic OCC operations.

Returns: OCC shape object.

### to_wire_segments(n_arc=20)

Extract centerline wire segments for Biot-Savart field computation.

| Arg | Type | Default | Description |
|-----|------|---------|-------------|
| `n_arc` | int | 20 | Number of straight segments per arc |

Returns: `(segments, current)` where segments is `list[((x1,y1,z1), (x2,y2,z2))]`.

```python
from biot_savart import h_filament
import numpy as np

segments, current = coil.to_wire_segments()
mu_0 = 4 * np.pi * 1e-7

H = sum(h_filament(p1, p2, obs, current) for p1, p2 in segments)
B = mu_0 * H
```

## Symmetry Methods

### mirror(plane)

Create a mirrored copy with reversed current.

| Arg | Values | Description |
|-----|--------|-------------|
| `plane` | `'xz'`, `'yz'`, `'xy'` | Mirror plane |

Returns: new `CoilBuilder` instance.

```python
upper = coil
lower = coil.mirror('xy')  # Mirrored coil with -I
```

### rotate_copies(axis, n_copies)

Create rotational copies around an axis.

| Arg | Type | Description |
|-----|------|-------------|
| `axis` | `'x'`, `'y'`, `'z'` | Rotation axis |
| `n_copies` | int | Number of copies (including original) |

Returns: list of `CoilBuilder` instances.

```python
quadrupole_coils = coil.rotate_copies('z', n_copies=4)
```

### close(tolerance=1e-6)

Automatically adjust arc angles to close the loop. Uses `scipy.optimize.minimize`.

```python
coil.close()
print(coil.gap)        # distance between end and start
print(coil.is_closed)  # True if gap < 1e-10
```

### combined_occ(others=None)

Combine multiple coils into a single OCC shape for STEP export.

```python
shape = upper.combined_occ([lower])
shape.WriteStep("both_coils.step")
```

## Properties

| Property | Type | Description |
|----------|------|-------------|
| `gap` | float | Distance between end position and start position (meters) |
| `is_closed` | bool | True if gap < 1e-10 |

## Cubit Integration

CoilBuilder works with Cubit via the `coil` APREPRO command:

```
coil "script.py" [output "path.step"] [noimport]
```

The script must define a `build_coil()` function returning a `CoilBuilder` instance:

```python
# my_coil.py
from coil_builder import CoilBuilder

def build_coil():
    mm = 1e-3
    return (CoilBuilder(current=2000)
        .set_start([0, 0.13, 0])
        .set_cross_section(35*mm, 105*mm)
        .add_straight(62.5*mm)
        .add_arc(22.5*mm, 180)
        .add_straight(62.5*mm)
        .add_arc(22.5*mm, 180))
```

In a Cubit journal file:

```
# ... mesh yoke geometry ...
coil "my_coil.py" output "coil.step"
export step "yoke_with_coil.step" overwrite
```

The `coil` command runs CoilBuilder in an external Python subprocess, generates a STEP file, and imports it into Cubit.

## Complex Example: 8-Segment Beam Steering Magnet

```python
mm = 1e-3

coil = (CoilBuilder(current=1265)
    .set_start([218*mm, -16.4*mm, -81*mm])
    .set_cross_section(122*mm, 122*mm)
    .add_straight(32.9*mm)
    .add_arc(121*mm, 64.6, tilt=90)
    .add_straight(1018.5*mm, tilt=90)
    .add_arc(121*mm, 115.4, tilt=-90)
    .add_straight(906.9*mm, tilt=90)
    .add_arc(121*mm, 115.4, tilt=-90)
    .add_straight(1018.5*mm, tilt=90)
    .add_arc(121*mm, 64.6, tilt=-90))

coil.write_step("beam_steering.step")
```

The `tilt` parameter rotates the cross-section around the current direction, enabling 3D coil paths that twist out of plane.

This 8-segment coil is built, assembled, and field-mapped end-to-end in the
runnable showcase notebook
[`docs/complex_coil_geometry/complex_coil.ipynb`](../complex_coil_geometry/complex_coil.ipynb)
(reproducible helper scripts kept beside the notebook in
`docs/complex_coil_geometry/`).

## Architecture

```
CoilSegment (ABC)
  StraightSegment  ->  ObjRecCur + ZXZ Euler transform
  ArcSegment       ->  ObjArcCur + ZXZ Euler transform

CoilBuilder
  State: position, orientation (3x3), cross-section
  add_straight/add_arc -> append segment, auto-update state
  to_radia()           -> convert all segments to Radia objects
  write_step()         -> OCC Loft + STEP export
  to_wire_segments()   -> centerline discretization for Biot-Savart
```

## Orientation Convention

Row-vector format: `orientation[i, :]` = i-th basis vector.

| Row | Direction | Description |
|-----|-----------|-------------|
| `orientation[0, :]` | X-axis | Radial (toward arc center) |
| `orientation[1, :]` | Y-axis | Current direction |
| `orientation[2, :]` | Z-axis | Cross-product (axial) |

Tilt rotates around the Y-axis (current direction) before segment construction.
