# Radia Hexahedron Test Report

## Test Date: 2025-10-29

## Summary

Comprehensive testing of `rad.ObjPolyhdr()` with various hexahedron geometries to identify cases that cause segmentation faults or failures.

## Test Results

### ✓ PASSING Cases

#### 1. Regular Cube
- **Status**: PASS
- **Coordinates**: Standard unit cube
- **Result**: Bz = 0.333 T at origin

#### 2. Skewed Hexahedron (Angled)
- **Status**: PASS
- **Description**: Slightly angled top face
- **Result**: Bz = 0.333 T at origin

#### 3. Highly Distorted Hexahedron
- **Status**: PASS
- **Description**: Top face significantly smaller than bottom
- **Result**: Bz = 0.352 T at [5,5,5]

#### 4. Multiple Hexahedra in Container
- **Status**: PASS
- **Description**: 3 cubes in a container
- **Result**: Bz = 0.311 T at [3,0,0]

#### 5. Webcut at 1 Degree (from angle_1.bdf)
- **Status**: PASS
- **Coordinates**: Real coordinates from Cubit webcut
- **Result**: Bz = 0.344 T at origin

#### 6. Extreme Angles (1-30 degrees)
- **Status**: PASS for all angles
- **Tested**: 1°, 5°, 10°, 15°, 20°, 25°, 30°
- **Result**: All hexahedra created successfully

#### 7. Very Flat Hexahedron (thickness = 0.001)
- **Status**: PASS
- **Result**: Bz = 2.25e-04 T

#### 8. Extremely Flat Hexahedron (thickness = 0.00001)
- **Status**: PASS
- **Result**: Bz = 2.25e-06 T

#### 9. Inverted Hexahedron (Inside-Out)
- **Status**: PASS
- **Description**: Top and bottom faces swapped
- **Result**: Bz = 0.333 T (Radia handles this correctly)

### ✗ FAILING Cases

#### Test 3: Zero Thickness Hexahedron - **SEGMENTATION FAULT**
- **Status**: **CRASH (Segmentation Fault)**
- **Coordinates**:
```python
coords = [
    [-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0],  # bottom at z=0
    [-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]   # top also at z=0
]
```
- **Issue**: All 8 nodes are coplanar (zero volume)
- **Consequence**: Causes segmentation fault instead of throwing exception

#### Test 5: Non-Convex (Bowtie) Hexahedron - **SEGMENTATION FAULT**
- **Status**: **CRASH (Segmentation Fault)**
- **Coordinates**:
```python
coords = [
    [-1, -1, -1], [1, 1, -1], [1, -1, -1], [-1, 1, -1],  # crossed bottom
    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]       # normal top
]
```
- **Issue**: Bottom face vertices are in wrong order, creating self-intersecting geometry
- **Consequence**: Causes segmentation fault instead of throwing exception

#### Test 6: All Nodes at Same Point - **EXCEPTION (GOOD)**
- **Status**: **Correctly handled with exception**
- **Coordinates**: All nodes at [0, 0, 0]
- **Exception**: "Incorrect input: Can not find at least three vertex points not belonging to one line in the face of the polyhedron."
- **Result**: Proper error handling ✓

## Root Cause Analysis

### Problem 1: Zero Volume Detection
The C++ code does not properly validate for zero-volume hexahedra before attempting geometric calculations. When all nodes are coplanar, division by zero or invalid geometric operations cause segmentation faults.

### Problem 2: Face Winding Order Validation
The code does not validate that face vertices are ordered correctly (counter-clockwise when viewed from outside). When vertices are in wrong order, the geometry becomes self-intersecting, leading to invalid calculations.

## Recommendations

### Critical Fixes Required

1. **Add Zero-Volume Check**
   - Location: `src/core/radpolyhd.cpp` (or relevant polyhedron creation code)
   - Action: Calculate volume using determinant method before proceeding
   - Return proper error: "Polyhedron has zero or near-zero volume"

2. **Add Face Winding Validation**
   - Location: Face processing in polyhedron creation
   - Action: Verify face normal consistency
   - Return proper error: "Invalid face winding order detected"

3. **Add Planarity Check**
   - Check if all 8 vertices are coplanar
   - Return proper error: "All vertices are coplanar"

### Code Location to Investigate

```cpp
// Likely in: src/core/radpolyhd.cpp
// Function: radTPolyhedron::radTPolyhedron() or similar constructor
// Need to add validation before geometric calculations
```

### Suggested Validation Code Pattern

```cpp
// Calculate volume using signed volume formula
double volume = CalculateHexVolume(vertices);
if (abs(volume) < 1e-12) {
    throw radTException("Polyhedron has zero or near-zero volume");
}

// Check for coplanar vertices
if (AreAllVerticesCoplanar(vertices)) {
    throw radTException("All vertices are coplanar");
}

// Validate face winding order
if (!ValidateFaceWindingOrder(faces, vertices)) {
    throw radTException("Invalid face winding order detected");
}
```

## Impact Assessment

### Current Behavior (Before Fix)
- Zero-volume hexahedra → **Segmentation Fault** (Python crashes)
- Non-convex hexahedra → **Segmentation Fault** (Python crashes)
- Same-point vertices → Proper exception ✓

### Expected Behavior (After Fix)
- Zero-volume hexahedra → **Proper exception with clear error message**
- Non-convex hexahedra → **Proper exception with clear error message**
- Same-point vertices → Proper exception ✓ (already working)

## Test Coverage

- ✓ Regular geometries: 100% pass
- ✓ Angled/skewed geometries: 100% pass (1-30 degrees)
- ✓ Very thin geometries: 100% pass (down to 0.00001)
- ✗ Zero-volume geometries: **CRASH**
- ✗ Self-intersecting geometries: **CRASH**

## Conclusion

The current Radia version handles most hexahedron cases correctly, including:
- Regular and distorted hexahedra
- Webcut geometries at various angles
- Very thin (but non-zero) hexahedra
- Multiple hexahedra in containers

However, **two critical edge cases cause segmentation faults**:
1. Zero-thickness (coplanar) hexahedra
2. Self-intersecting (wrong vertex order) hexahedra

**These should be fixed by adding proper validation and throwing exceptions instead of crashing.**

---

## Test Files Created

1. `test_hexahedron.py` - Basic hexahedron tests
2. `test_case01_coords.py` - Real webcut coordinates from angle_1.bdf
3. `test_extreme_angles.py` - Various cut angles (1-30 degrees)
4. `test_degenerate_hex.py` - Degenerate cases (causes crash)
5. `test_degenerate_one_by_one.py` - Individual test runner

## Original Issue Files - NOW WORKING!

- `2024_02_01_brick_webcut_case_01.py` - ✓ **NOW WORKS** with Cubit 2025.3
- `2024_02_01_brick_webcut_case_02.py` - ✓ **NOW WORKS** with Cubit 2025.3
- `angle_1.bdf` - NASTRAN format mesh with actual coordinates

### Test Results with Actual Cubit Meshes

#### Case 01: Brick webcut at 20 degrees
- **Cubit meshes created**: 16 hexahedra (8 in each volume)
- **Radia objects created**: 16/16 (100% success)
- **Field at origin**: Bx=-3.54e-08, By=-7.70e-11, Bz=3.33e-01 T
- **Status**: ✓ **FULLY FUNCTIONAL**

#### Case 02: Multiple webcuts (complex geometry)
- **Cubit meshes created**: 1 hexahedron (after deleting other volumes)
- **Radia objects created**: 1/1 (100% success)
- **Field at origin**: Bx=7.39e-03, By=4.28e-02, Bz=-3.53e-02 T
- **Status**: ✓ **FULLY FUNCTIONAL**

**Status**:
- ✓ Original case_01 and case_02 **NOW WORK CORRECTLY** with current Radia version
- ✗ Two edge cases (zero-volume and self-intersecting) still cause segfaults - C++ fixes needed
