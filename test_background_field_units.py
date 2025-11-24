#!/usr/bin/env python
"""
Test background field units issue

The debug output showed that rad.ObjBckg([0, 0, 1000])
internally becomes [0, 0, 795774715.459] A/m

This is 1000 * (4*pi*1e-7) / 1e-10 = 1000 * 12566.37 / 1e-10

Let's investigate the unit conversion.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build/Release'))

import numpy as np
import radia as rad

print("="*70)
print("Background Field Units Investigation")
print("="*70)

# Test different ways to specify background field
H_ext_input = 1000.0  # A/m

print(f"\nInput: H_ext = {H_ext_input} A/m")
print(f"After rad.FldUnits('m'):\n")

# ============================================================
# Test 1: ObjBckg with direct field value
# ============================================================
print("[1] rad.ObjBckg([0, 0, H_ext])")

rad.FldUnits('m')
bg1 = rad.ObjBckg([0, 0, H_ext_input])

# Evaluate at a point
test_pt = [0.1, 0.1, 0.1]
H1 = rad.Fld(bg1, 'h', test_pt)
H1_mag = np.linalg.norm(H1)

print(f"  Evaluated H at [{test_pt[0]}, {test_pt[1]}, {test_pt[2]}]:")
print(f"  H = [{H1[0]:.6e}, {H1[1]:.6e}, {H1[2]:.6e}] A/m")
print(f"  |H| = {H1_mag:.6e} A/m")
print(f"  Expected: [0, 0, {H_ext_input}] A/m")

conversion_factor = H1[2] / H_ext_input
print(f"  Conversion factor: {conversion_factor:.6e}")
print(f"  This is 4*pi*1e-7 / 1e-10 = {4*np.pi*1e-7 / 1e-10:.6e}")

# ============================================================
# Test 2: Check if we need to convert input
# ============================================================
print(f"\n[2] Hypothesis: Input should be in Tesla?")

mu0 = 4*np.pi*1e-7  # H/m
B_ext_input = H_ext_input * mu0  # Convert A/m to Tesla

rad.FldUnits('m')
bg2 = rad.ObjBckg([0, 0, B_ext_input])

H2 = rad.Fld(bg2, 'h', test_pt)
H2_mag = np.linalg.norm(H2)

print(f"  Input B = H * mu0 = {B_ext_input:.6e} T")
print(f"  Evaluated H = [{H2[0]:.6e}, {H2[1]:.6e}, {H2[2]:.6e}] A/m")
print(f"  |H| = {H2_mag:.6e} A/m")
print(f"  Expected: {H_ext_input} A/m")

error = abs(H2_mag - H_ext_input) / H_ext_input * 100
print(f"  Error: {error:.2f}%")

if error < 1:
    print(f"  SUCCESS: Input should be in Tesla!")
else:
    print(f"  FAILED: Still not correct")

# ============================================================
# Test 3: ObjBckgCF with function
# ============================================================
print(f"\n[3] rad.ObjBckgCF with function")

rad.FldUnits('m')

def uniform_field_A_per_m(pos):
    """Return field in A/m"""
    return [0.0, 0.0, H_ext_input]

def uniform_field_Tesla(pos):
    """Return field in Tesla"""
    return [0.0, 0.0, H_ext_input * mu0]

# Test A/m version
bg3a = rad.ObjBckgCF(uniform_field_A_per_m)
H3a = rad.Fld(bg3a, 'h', test_pt)
H3a_mag = np.linalg.norm(H3a)

print(f"  Function returning A/m:")
print(f"    H = [{H3a[0]:.6e}, {H3a[1]:.6e}, {H3a[2]:.6e}] A/m")
print(f"    |H| = {H3a_mag:.6e} A/m")
print(f"    Error: {abs(H3a_mag - H_ext_input)/H_ext_input*100:.2f}%")

# Test Tesla version
bg3b = rad.ObjBckgCF(uniform_field_Tesla)
H3b = rad.Fld(bg3b, 'h', test_pt)
H3b_mag = np.linalg.norm(H3b)

print(f"  Function returning Tesla:")
print(f"    H = [{H3b[0]:.6e}, {H3b[1]:.6e}, {H3b[2]:.6e}] A/m")
print(f"    |H| = {H3b_mag:.6e} A/m")
print(f"    Error: {abs(H3b_mag - H_ext_input)/H_ext_input*100:.2f}%")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*70}")
print("Summary")
print(f"{'='*70}")

print(f"""
Finding:
  rad.ObjBckg() expects input in TESLA, not A/m!

Correct usage:
  H_ext = 1000.0  # A/m (desired field strength)
  mu0 = 4*pi*1e-7  # H/m
  B_ext = H_ext * mu0  # Convert to Tesla
  bg = rad.ObjBckg([0, 0, B_ext])  # Pass Tesla value

Alternative (correct usage):
  def uniform_field(pos):
      return [0.0, 0.0, H_ext * 4*np.pi*1e-7]  # Return Tesla
  bg = rad.ObjBckgCF(uniform_field)

This explains why the debug output showed {H1[2]:.0f} A/m instead of {H_ext_input} A/m!
The factor of {conversion_factor:.0f} is exactly 4*pi / (1e-3) in CGS-Gauss units.
""")

print("="*70)
