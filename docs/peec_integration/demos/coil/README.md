# coil — Coil Impedance Analysis

PEEC-based coil impedance analysis demos.

## Files

| File | Description |
|------|-------------|
| `coil_impedance_peec.py` | Coil impedance frequency response via PEEC |
| `coil_on_magnetic_core_peec.py` | Coil on magnetic core analysis (CplMag solver) |
| `demo_circular_coil_gmsh.py` | Historical filename; now a Netgen/OCC `.vol` boundary-surface demo |

## Outputs

- `coil_impedance_peec.png` — Impedance frequency response
- `coil_impedance_peec_solver.png` — Solver details
- `coil_impedance_peec_coil.msh` — GMSH MSH v4.1 coil centerline
- `coil_on_magnetic_core_*.png` — Magnetic core analysis results
