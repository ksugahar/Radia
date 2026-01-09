# ESIM + RWG-EFIE Examples

This folder contains examples demonstrating the integration of ESIM (Effective Surface Impedance Method) with RWG-EFIE (Rao-Wilton-Glisson Electric Field Integral Equation) solvers for induction heating analysis.

## Demo Scripts

### 1. demo_esim_rwg_nonlinear.py
Demonstrates the nonlinear coupled EFIE solver with ESIM for analyzing ferromagnetic workpieces:
- ESI table generation from B-H curve
- Nonlinear iteration with underrelaxation
- Power loss computation
- Field-dependent surface impedance

### 2. demo_pfft_acceleration.py
Demonstrates pFFT (Pre-corrected FFT) acceleration for large-scale problems:
- O(N log N) matrix-vector products instead of O(N^2)
- Memory reduction for large problems
- Integration with iterative solvers

### 3. demo_hybrid_peec_rwg.py
Demonstrates the hybrid PEEC + RWG solver:
- PEEC for wire/coil modeling (thin-wire approximation)
- RWG for workpiece surface currents
- Mutual coupling between coil and workpiece
- Optimal DOF allocation

### 4. demo_loop_star_decomposition.py
Demonstrates Loop-Star decomposition for low-frequency stability:
- Eliminates low-frequency breakdown in EFIE
- Separates inductive (Loop) and capacitive (Star) effects
- MQS mode for induction heating (Loop-only)
- Reduced DOF and faster computation

## Quick Start

```python
import sys
sys.path.insert(0, '../../src/radia')

from rwg_efie_solver import (
    RWGMesh,
    NonlinearCoupledEFIESolver,
    create_nonlinear_solver,
    create_induction_heating_model_with_esim,
)
from esim_cell_problem import ESITable

# Define B-H curve
STEEL_BH_CURVE = [
    [0, 0], [100, 0.2], [500, 0.9], [1000, 1.3],
    [5000, 1.8], [10000, 1.95], [50000, 2.1],
]

# Create meshes
coil_mesh = RWGMesh()
coil_mesh.create_circular_disk(center=[0,0,0.01], radius=0.03, normal=[0,0,1])

workpiece_mesh = RWGMesh()
workpiece_mesh.create_rectangular_plate(center=[0,0,0], Lx=0.1, Ly=0.1, normal=[0,0,1])

# Create solver using helper function
solver = create_nonlinear_solver(
    coil_mesh=coil_mesh,
    workpiece_mesh=workpiece_mesh,
    bh_curve=STEEL_BH_CURVE,
    sigma=2e6,       # Hot steel conductivity
    frequency=50000  # 50 kHz
)

# Solve
solver.set_voltage_excitation(10.0)  # 10V
result = solver.solve()

print(f"Heating power: {result['power_workpiece']:.2f} W")
```

## Key Features

| Feature | Description |
|---------|-------------|
| ESIM | Nonlinear surface impedance Z(H) from cell problem |
| RWG-EFIE | Accurate surface current modeling |
| pFFT | O(N log N) acceleration for large problems |
| Hybrid PEEC+RWG | Optimal coil+workpiece modeling |
| Nonlinear Iteration | Converges in 5-10 steps with underrelaxation |
| Loop-Star | Low-frequency stable MQS formulation |

## Requirements

- Python 3.8+
- NumPy
- SciPy (for linear algebra)
- Radia package

## References

1. K. Hollaus et al., "A Nonlinear Effective Surface Impedance," IEEE Trans. Magnetics, 2025
2. S.M. Rao et al., "Electromagnetic scattering by surfaces of arbitrary shape," IEEE TAP, 1982
3. J.R. Phillips, J.K. White, "A precorrected-FFT method," IEEE TCAD, 1997
4. A.E. Ruehli, "Equivalent Circuit Models for 3D Multiconductor Systems," IEEE Trans. MTT, 1974
5. G. Vecchi, "Loop-Star Decomposition of Basis Functions in EFIE," IEEE TAP, 1999
