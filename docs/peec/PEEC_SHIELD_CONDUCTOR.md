# PEEC Shield Conductor Modeling

**Status**: Implemented (2026-02-17)
**Module**: `src/radia/peec_shield.py`

---

## Overview

Shield conductor modeling for PEEC analysis. Aluminum (or other conductor) shield
plates are added as PEEC filament grids to PEECBuilder. Shield currents are induced
by magnetic coupling with coil segments -- no external port excitation.

### Physics

The shield effect on coil impedance appears as reflected impedance:

```
Delta_Z_shield(f) = Z_port(with_shield) - Z_port(without_shield)
Delta_R_shield = Re(Delta_Z_shield)    (eddy current loss)
Delta_L_shield = Im(Delta_Z_shield) / omega  (flux screening)
```

Physical behavior:
- Shield eddy currents oppose the coil field -> reduced flux linkage -> **lower L**
- Eddy current losses dissipate power -> reflected resistance -> **higher R**
- At low frequency, skin depth >> shield thickness -> weak shielding
- At high frequency, skin depth << shield thickness -> strong shielding

### Architecture

```
Coil (PEEC port)          Shield (no port, floating)
   N1 -- seg1 --> N2         Sn1 -- sh_seg --> Sn2
   N2 -- seg2 --> N3         Sn3 -- sh_seg --> Sn4
   ...                       ... (grid of filaments)

   Port: N1 - Nn             No port: floating conductor

MNA solver handles both as connected components:
  - Coil: one port, one connected component
  - Shield: no port, separate connected component
  - Coupling through mutual inductance L matrix
```

---

## API

### Functions

| Function | Purpose |
|----------|---------|
| `add_shield_mesh()` | Add rectangular conducting shield grid to PEECBuilder |
| `add_coil_from_segments()` | Add coil segments to PEECBuilder with connectivity |
| `compute_shield_effect()` | 1-port: compute shield effect on single coil Z(f) |
| `compute_shield_effect_2port()` | 2-port: compute shield effect on TX/RX coupling |

### Basic Usage (1-port)

```python
from peec_matrices import PEECBuilder
from peec_shield import add_shield_mesh, add_coil_from_segments, compute_shield_effect

# Square loop coil
p1, p2 = create_square_loop(side=0.05, z_m=0.0, n_seg_per_side=8)

# Compute shield effect over frequency range
result = compute_shield_effect(
    p1, p2, wire_w=1e-3, wire_h=1e-3,
    shield_center=[0, 0, -0.005],
    shield_size_x=0.08, shield_size_y=0.08,
    shield_thickness=1e-3, shield_sigma=3.5e7,  # aluminum
    shield_nx=6, shield_ny=6,
    freqs=np.logspace(2, 6, 20))

# Results
print(f"Delta_L = {result['Delta_L'][-1]*1e9:.1f} nH")  # Negative (screening)
print(f"Delta_R = {result['Delta_R'][-1]*1e3:.1f} mOhm")  # Positive (eddy loss)
```

### 2-Port Usage (TX/RX with shield)

```python
from peec_shield import compute_shield_effect_2port

result = compute_shield_effect_2port(
    tx_p1, tx_p2, tx_wire_w, tx_wire_h,
    rx_p1, rx_p2, rx_wire_w, rx_wire_h,
    shield_center=[0, 0, -0.005],
    shield_size_x=0.10, shield_size_y=0.10,
    shield_thickness=1e-3, shield_sigma=3.5e7,
    shield_nx=8, shield_ny=8,
    freq=100e3)

print(f"k without shield: {result['k_no']:.4f}")
print(f"k with shield:    {result['k_with']:.4f}")
print(f"Delta_R_tx: {result['Delta_R1']*1e3:.1f} mOhm")
print(f"Delta_R_rx: {result['Delta_R2']*1e3:.1f} mOhm")
```

---

## Shield Mesh Generation

`add_shield_mesh()` creates a grid of (nx+1)*(ny+1) nodes connected by
nx*(ny+1) horizontal + (nx+1)*ny vertical segments.

### Dual Mesh Width

Each filament's cross-section width is determined by the dual mesh:
- **Interior filaments**: width = full cell size (dx or dy)
- **Boundary filaments**: width = half cell size (dx/2 or dy/2)
- **Height**: shield thickness for all filaments

```
Dual mesh (6x4 example):

  dy/2 |---*---*---*---*---*---*---|
  dy   |---*---*---*---*---*---*---|
  dy   |---*---*---*---*---*---*---|
  dy   |---*---*---*---*---*---*---|
  dy/2 |---*---*---*---*---*---*---|
       dx/2 dx  dx  dx  dx  dx dx/2
```

---

## Closed-Loop Coil Handling

`add_coil_from_segments()` automatically detects closed-loop coils (all nodes have
degree 2) and handles them by:

1. Identifying the loop closure point
2. Duplicating the closing node (same position, new node ID)
3. Redirecting the last segment to the new node
4. Setting port between original and duplicated node

This "breaks" the loop at one point to allow port voltage measurement.

---

## Validation Results

### Test 1: Basic shield effect (3/3 PASS)

- L decreases with shield at high frequency
- R increases with shield (positive Delta_R)
- Shield effect increases with frequency

### Test 2: Shield effect vs distance (PASS)

- Closer shield -> larger Delta_L and Delta_R

### Test 3: No-shield baseline (PASS)

- Self-inductance matches analytical within 7.1%

---

## Bug Fix: Mutual Inductance Sign (2026-02-17)

### Problem

For anti-parallel filaments, the Rosa/Grover analytical path used
`std::abs(dot)` which discarded the sign, making all mutual inductances
positive regardless of current direction.

### Root Cause

In `rad_peec_matrices.cpp`, the parallel/anti-parallel analytical formula:
```cpp
// BEFORE (BUGGY):
return std::abs(dot) * M;  // Discards sign for anti-parallel

// AFTER (FIXED):
return dot * M;  // Preserves sign for anti-parallel filaments
```

### Impact

- Mutual inductance between anti-parallel segments was incorrectly positive
- For spiral coils (many anti-parallel segments), positive and negative terms
  canceled to near-zero instead of summing correctly
- Effect: M ~0 instead of correct value (e.g., -2.89 uH)

### Why F(x,d) Even-ness Matters

The Rosa/Grover formula uses F(x,d) = x*arsinh(x/d) - sqrt(x^2+d^2).
This is an even function: F(-x,d) = F(x,d). Therefore the
`if (dot < 0) p = -p;` correction has NO effect on the integral.
The sign MUST come from the dot product multiplication at the end.

---

## Source Files

| File | Description |
|------|-------------|
| `src/radia/peec_shield.py` | Shield mesh, coil segments, compute functions |
| `src/core/rad_peec_matrices.cpp` | C++ Neumann integral (sign fix line 776) |
| `validation_test/peec_integration/verification/validate_shield_delta_r.py` | Validation (3/3 PASS) |
