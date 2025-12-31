#!/usr/bin/env python
"""
Debug script: Compare A field from ObjRecMag vs ObjHexahedron
at various points including [0,0,0.05] on the symmetry axis.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'radia'))

import numpy as np
import radia as rad

print('=' * 70)
print('Debug: A field comparison ObjRecMag vs ObjHexahedron')
print('=' * 70)

rad.UtiDelAll()
rad.FldUnits('m')

# Magnet dimensions
dx, dy, dz = 0.02, 0.02, 0.03  # half dimensions
MU_0 = 4 * np.pi * 1e-7
Br = 1.2  # T
Mr = Br / MU_0  # A/m

# Create ObjRecMag
center = [0, 0, 0]
dimensions = [2*dx, 2*dy, 2*dz]  # [0.04, 0.04, 0.06]
rec_mag = rad.ObjRecMag(center, dimensions, [0, 0, Mr])

# Create ObjHexahedron
# Vertices 1-4 define the bottom face (counter-clockwise when viewed from below),
# vertices 5-8 are directly above vertices 1-4.
vertices = [
    [-dx, -dy, -dz],  # 1 (bottom, front-left)
    [ dx, -dy, -dz],  # 2 (bottom, front-right)
    [ dx,  dy, -dz],  # 3 (bottom, back-right)
    [-dx,  dy, -dz],  # 4 (bottom, back-left)
    [-dx, -dy,  dz],  # 5 (top, front-left)
    [ dx, -dy,  dz],  # 6 (top, front-right)
    [ dx,  dy,  dz],  # 7 (top, back-right)
    [-dx,  dy,  dz],  # 8 (top, back-left)
]
hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

print()
print('Magnet: %.1f x %.1f x %.1f mm, Mz = %.0f A/m' % (
    2*dx*1000, 2*dy*1000, 2*dz*1000, Mr))
print()

# Test points
test_points = [
    # Off-axis points (should work)
    [0.05, 0.0, 0.0],
    [0.0, 0.05, 0.0],
    [0.03, 0.03, 0.03],
    [0.04, 0.02, 0.05],
    # On-axis points (problem area)
    [0.0, 0.0, 0.05],   # z-axis
    [0.0, 0.0, 0.035],  # closer
    [0.0, 0.0, 0.10],   # farther
]

print('%-30s  %12s  %12s  %12s' % ('Point', '|A_rec|', '|A_hex|', 'Diff %'))
print('-' * 70)

for pt in test_points:
    A_rec = rad.Fld(rec_mag, 'a', pt)
    A_hex = rad.Fld(hex_mag, 'a', pt)

    A_rec_mag = np.sqrt(A_rec[0]**2 + A_rec[1]**2 + A_rec[2]**2)
    A_hex_mag = np.sqrt(A_hex[0]**2 + A_hex[1]**2 + A_hex[2]**2)

    if A_rec_mag > 1e-15:
        diff_pct = 100 * abs(A_rec_mag - A_hex_mag) / A_rec_mag
    else:
        diff_pct = 0 if A_hex_mag < 1e-15 else 100

    print('[%6.3f, %6.3f, %6.3f]  %12.4e  %12.4e  %10.2f%%' % (
        pt[0], pt[1], pt[2], A_rec_mag, A_hex_mag, diff_pct))

# Detailed analysis for [0,0,0.05]
print()
print('=' * 70)
print('Detailed analysis at [0, 0, 0.05]')
print('=' * 70)

pt = [0.0, 0.0, 0.05]
A_rec = rad.Fld(rec_mag, 'a', pt)
A_hex = rad.Fld(hex_mag, 'a', pt)
B_rec = rad.Fld(rec_mag, 'b', pt)
B_hex = rad.Fld(hex_mag, 'b', pt)

print()
print('ObjRecMag:')
print('  A = [%.6e, %.6e, %.6e]' % tuple(A_rec))
print('  |A| = %.6e' % np.sqrt(sum(x**2 for x in A_rec)))
print('  B = [%.6f, %.6f, %.6f] T' % tuple(B_rec))
print('  |B| = %.6f T' % np.sqrt(sum(x**2 for x in B_rec)))

print()
print('ObjHexahedron:')
print('  A = [%.6e, %.6e, %.6e]' % tuple(A_hex))
print('  |A| = %.6e' % np.sqrt(sum(x**2 for x in A_hex)))
print('  B = [%.6f, %.6f, %.6f] T' % tuple(B_hex))
print('  |B| = %.6f T' % np.sqrt(sum(x**2 for x in B_hex)))

# At z=0.05, the observation point is on the z-axis
# For a z-magnetized block, by symmetry, A should be perpendicular to z
# A = (1/4pi) * M x BufVect
# M = [0, 0, Mz], so A_x = (1/4pi) * Mz * BufVect_y, A_y = -(1/4pi) * Mz * BufVect_x
# On the z-axis, BufVect should also have only z-component by symmetry
# Therefore A should be zero on the z-axis!

print()
print('Expected behavior on z-axis:')
print('  M = [0, 0, Mz] -> A = (1/4pi) * M x BufVect')
print('  On z-axis, BufVect should have only z-component by symmetry')
print('  Therefore: A = (1/4pi) * [Mz*BufVect_y, -Mz*BufVect_x, 0]')
print('  If BufVect_x = BufVect_y = 0 (symmetry), then |A| = 0')
print()
print('However, ObjRecMag shows |A| = %.4e (NOT zero!)' % np.sqrt(sum(x**2 for x in A_rec)))
print('This suggests BufVect has x,y components even on z-axis.')

# Check if RecMag is using a perturbation to avoid singularities
print()
print('=' * 70)
print('Testing near-axis points to understand the behavior')
print('=' * 70)

epsilon_values = [0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3]
z_val = 0.05

print()
print('Testing pt = [eps, 0, 0.05] for various eps:')
print('%-12s  %12s  %12s' % ('eps', '|A_rec|', '|A_hex|'))
print('-' * 40)

for eps in epsilon_values:
    pt = [eps, 0, z_val]
    A_rec = rad.Fld(rec_mag, 'a', pt)
    A_hex = rad.Fld(hex_mag, 'a', pt)
    A_rec_mag = np.sqrt(sum(x**2 for x in A_rec))
    A_hex_mag = np.sqrt(sum(x**2 for x in A_hex))
    print('%12.2e  %12.4e  %12.4e' % (eps, A_rec_mag, A_hex_mag))

# Additional test: Check component-wise values
print()
print('=' * 70)
print('Component-wise comparison at [0, 0, 0.05]')
print('=' * 70)
pt = [0.0, 0.0, 0.05]
A_rec = rad.Fld(rec_mag, 'a', pt)
A_hex = rad.Fld(hex_mag, 'a', pt)
print(f'ObjRecMag: A = [{A_rec[0]:.10e}, {A_rec[1]:.10e}, {A_rec[2]:.10e}]')
print(f'ObjHexahedron: A = [{A_hex[0]:.10e}, {A_hex[1]:.10e}, {A_hex[2]:.10e}]')

# Check if RecMag uses perturbation
# The perturbation in RecMag is applied when BfSt.x[ii] == 0, etc.
# For point [0, 0, 0.05], P_min_CenPo = [0, 0, 0.05]
# BfSt.x[0] = 0 - (-0.02) = 0.02, BfSt.x[1] = 0 - 0.02 = -0.02
# Neither is zero, so no perturbation should be applied!

print()
print('Analysis of BfSt values for point [0, 0, 0.05]:')
P_min_CenPo = [0, 0, 0.05]
HalfDim = [dx, dy, dz]
for ii in [0, 1]:
    Eps = ii*2-1
    BfSt_x = -P_min_CenPo[0] + Eps*HalfDim[0]
    BfSt_y = -P_min_CenPo[1] + Eps*HalfDim[1]
    BfSt_z = -P_min_CenPo[2] + Eps*HalfDim[2]
    print(f'ii={ii}: BfSt = [{BfSt_x:.6f}, {BfSt_y:.6f}, {BfSt_z:.6f}]')
    print(f'       Any zero? x={BfSt_x==0}, y={BfSt_y==0}, z={BfSt_z==0}')

rad.UtiDelAll()
