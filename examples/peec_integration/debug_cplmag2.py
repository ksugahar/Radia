"""Debug script for CplMag solver - step by step."""

import sys
sys.path.insert(0, '../../src/radia')
import radia as rad
import numpy as np

rad.UtiDelAll()
rad.FldUnits('m')

# Create simple coil
loop_radius = 0.05  # 50 mm
coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
print(f'Coil handle: {coil}')

# Get coil geometry info - skip CndInfo since it doesn't exist

# Create a very simple magnetic core - small block inside the coil
# Use ObjRecMag to have a simpler case
core = rad.ObjRecMag([0, 0, 0], [0.02, 0.02, 0.02], [0, 0, 0])  # 20mm cube at origin
print(f'Core handle: {core}')

# Apply material
mat = rad.MatLin(1000)  # mu_r = 1000
rad.MatApl(core, mat)

# Test PEEC only
print("\n=== PEEC only (no core) ===")
rad.CndSetFrequency(coil, 1000)
rad.CndSetVoltage(coil, 1.0, 0.0)
rad.CndSolve(coil)
Z_peec = rad.CndGetImpedance(coil)
print(f'PEEC Z = {Z_peec.real*1000:.4f} + j{Z_peec.imag*1000:.4f} mOhm')
L_peec = Z_peec.imag / (2 * np.pi * 1000)
print(f'PEEC L = {L_peec * 1e9:.2f} nH')

# Before testing CplMag, let's verify the core works with standard Radia solve
print("\n=== Testing standard Radia Solve on core ===")
# Add a uniform background field
bkg = rad.ObjBckg(lambda p: [0, 0, 0.01])  # 0.01 T = 10 mT background
container = rad.ObjCnt([core, bkg])
result = rad.Solve(container, 0.0001, 1000, 0)  # LU solver
print(f'Radia Solve iterations: {result}')

# Get magnetization from the core
M = rad.Fld(core, 'm', [0, 0, 0])
print(f'Core M at center: [{M[0]:.1f}, {M[1]:.1f}, {M[2]:.1f}] A/m')

# Clean up for next test
rad.UtiDelAll()
rad.FldUnits('m')

# Recreate coil
coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

# Recreate core
core = rad.ObjRecMag([0, 0, 0], [0.02, 0.02, 0.02], [0, 0, 0])
mat = rad.MatLin(1000)
rad.MatApl(core, mat)

# Now test CplMag
print("\n=== CplMag (with core) ===")
solver = rad.CplMagCreate(coil, core)
print(f'Solver handle: {solver}')

freq = 1000  # 1 kHz
rad.CplMagSetFrequency(solver, freq)
rad.CplMagSetVoltage(solver, 1.0, 0.0)
rad.CplMagSetMu(solver, 1000, 50)

# Solve
result = rad.CplMagSolve(solver)
Z = result['Z']
omega = 2 * np.pi * freq
L = Z.imag / omega if omega > 0 else 0

print(f'Result:')
print(f'  Z = {Z.real*1000:.4f} + j{Z.imag*1000:.4f} mOhm')
print(f'  L = {L * 1e9:.2f} nH')
print(f'  P_cond = {result["P_conductor"]:.6f} W')
print(f'  P_mag = {result["P_magnet"]:.6f} W')

# Comparison
print(f'\n=== Comparison ===')
print(f'PEEC L (air core): {L_peec * 1e9:.2f} nH')
print(f'CplMag L (mu_r=1000 core): {L * 1e9:.2f} nH')
if L_peec > 0:
    print(f'Ratio: {L / L_peec:.3f} (should be > 1 for mu_r > 1)')

rad.CplMagDelete(solver)
