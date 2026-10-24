# PEEC SIBC Validation: Circular Wire vs Bessel Analytical Solution

This validation workflow compares Radia PEEC SIBC implementation against exact Bessel function analytical solution for circular wires.

## Overview

**Goal**: Validate that PEEC with Surface Impedance Boundary Condition (SIBC) accurately captures skin effect in conductors.

**Method**:
1. Generate the declared 1D wire/coil mesh with native Netgen
2. Import mesh into Radia PEEC (filament model, like FastHenry/FastImp)
3. Apply SIBC for frequency-dependent resistance
4. Compare with exact Bessel function solution for circular wires

**Reference**:
- F.W. Grover, "Inductance Calculations", Dover, 1946
- Bessel function solution: `Z = (k*L)/(2πrσ) * J0(kr)/J1(kr)`

## File Structure

```
validation_test/peec_integration/verification/
├── README.md                           # This file
├── netgen_mesh_generation/
│   ├── generate_centerlines.py       # Netgen script for the analytic centerlines
│   ├── circular_wire_centerline.vol    # Native Netgen mesh
│   └── circular_wire_params.txt        # Wire parameters
├── validate_circular_wire_sibc.py      # Validation script (main)
└── validation_circular_wire_sibc.png   # Results plot
```

`centerline_vol.py` loads native Netgen 1D `.vol` files and checks their
versioned `.labels.json` contracts before returning points and edges. These
meshes contain actual centerline segments, no fabricated volume elements. The
existing 10-segment wire and 36-segment coil, physical parameters, references,
and error criteria are preserved. Their source is the analytic definitions in
`netgen_mesh_generation/generate_centerlines.py`; no Gmsh file is read.

## Step-by-Step Instructions

### Step 1: Generate native Netgen centerline meshes

Run the Netgen script to generate a 1m straight wire with 1mm radius, divided into 10 segments:

```bash
cd netgen_mesh_generation
python generate_centerlines.py
```

**Output**:
- `circular_wire_centerline.vol` - Native Netgen 1D .vol mesh with 1D edge elements
- `circular_wire_params.txt` - Wire parameters (radius, length, sigma, segments)

**Expected output**:
```
Wire parameters:
  Length: 1000 mm
  Radius: 1.00 mm (diameter: 2.00 mm)
  Segments: 10

Generated mesh:
  Nodes: 11
  Elements: 10
  Total length: 1.000000 m
```

### Step 2: Run Validation

Run the validation script to compare PEEC SIBC with Bessel analytical solution:

```bash
cd ..
python validate_circular_wire_sibc.py
```

**What it does**:
1. Reads the analytic Netgen centerline (`circular_wire_centerline.vol`)
2. Creates PEEC model with 10 filament segments
3. Sweeps frequency from 10 Hz to 1 MHz (30 points)
4. For each frequency:
   - Computes PEEC resistance with SIBC (rectangular approximation)
   - Computes exact Bessel function solution (circular wire)
   - Calculates error
5. Generates 4-panel validation plot

**Expected results**:
```
Resistance Error:
  Mean error: < 5%
  Max error: < 10%
  RMS error: < 6%

Reactance Error (f > 100 Hz):
  Mean error: < 5%
  Max error: < 10%
```

**Output**:
- `validation_circular_wire_sibc.png` - 4-panel plot showing:
  - Resistance vs Frequency
  - Reactance vs Frequency
  - Resistance Error vs Frequency
  - Impedance Magnitude vs Frequency

### Step 3: Interpret Results

**Acceptable Error**:
- Resistance error < 10% is acceptable
- Error comes from **rectangular approximation** of circular cross-section
- PEEC uses rectangular SIBC with equivalent square: `side = √π * radius`

**Key Observations**:
1. **Low frequency (< 1 kHz)**: Both methods converge to DC resistance
2. **Medium frequency (1-100 kHz)**: Good agreement (< 5% error)
3. **High frequency (> 100 kHz)**: Rectangular approximation introduces error

**For exact circular wire results**: Use `ComputeCircular()` method (implemented but not yet exposed to Python API).

## Theory: Rectangular vs Circular SIBC

### Circular Wire (Exact Bessel Solution)

```
Z = (k*L)/(2πrσ) * J0(kr)/J1(kr)
```

Where:
- `k = √(jωμσ)` - Propagation constant
- `J0, J1` - Bessel functions of first kind
- `r` - Wire radius
- `L` - Wire length
- `σ` - Conductivity

### Rectangular Approximation (PEEC)

```
Z_s = √(jωμ/σ) * (perimeter / area)
```

For equivalent square:
- `side = √π * radius`
- `perimeter = 4 * √π * radius`
- `area = π * radius²`

**Error source**: Rectangular SIBC assumes uniform current distribution on 4 flat sides, while circular wire has smooth current distribution.

## References

1. **SIBC Theory**:
   - R.F. Harrington, "Field Computation by Moment Methods", 1968
   - C.R. Paul, "Inductance: Loop and Partial", 2010

2. **Bessel Function Solution**:
   - F.W. Grover, "Inductance Calculations", Dover, 1946
   - E.B. Rosa, F.W. Grover, "Formulas and Tables for Calculation of Mutual and Self-Inductance", NBS Bulletin, 1916

3. **FastHenry/FastImp**:
   - M. Kamon et al., "FastHenry: A Multipole-Accelerated 3-D Inductance Extraction Program", IEEE TMTT, 1994
   - Z. Zhu et al., "FastImp: A Fast BEM Package for 3-D Impedance Extraction", IEEE TCAD, 2002

## Next Steps

After validating basic SIBC:

1. **Phase 1**: Analytical Bessel proximity effect for parallel wires
2. **Phase 2**: ESIM (Effective Surface Impedance Method) for complex geometries
3. **Phase 3**: Full 2D FEM cross-section solver

---

**Last Updated**: 2026-02-13
**Status**: ✅ Ready for validation
