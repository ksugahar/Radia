# Induction Heating

Examples demonstrating electromagnetic induction heating analysis using the ESIM (Effective Surface Impedance Method) and RWG-EFIE solvers in Radia, including nonlinear ferromagnetic materials, coil-workpiece coupling, and wireless power transfer.

## Scripts

| File | Description |
|------|-------------|
| `esim_demo.py` | Demonstrates ESIM fundamentals: BH-curve interpolation, 1D cell problem solver, ESI table generation, workpiece power loss computation, and frequency sweep analysis for a steel billet. |
| `esim_induction_heating_demo.py` | Complete coupled induction heating workflow with spiral/loop coil models and ESIM workpieces, including frequency sweep and current scaling with saturation effects. |
| `demo_esim_impedance.py` | Coil-workpiece coupled impedance calculation with complex permeability support, H-dependent nonlinear permeability, reflected impedance analysis, and resonance capacitor design. |
| `demo_rwg_efie_3d.py` | RWG-EFIE (Rao-Wilton-Glisson Electric Field Integral Equation) solver for 3D surface-element induction heating, demonstrating mesh creation (loop, spiral, plate, disk, cylinder), coupled solving, and frequency sweep. |
| `demo_wpt_coupling.py` | Wireless power transfer coil coupling analysis: mutual inductance via the Neumann integral, coupling coefficient vs. distance and misalignment, mutual resistance (proximity effect), WPT system efficiency, and S-S resonant topology design. |
| `test_esim_integration.py` | Integration tests validating the full ESIM workflow: BH-curve interpolation accuracy, cell problem solver, ESI table generation, coil field computation, coupled solver convergence, VTK export, and physical consistency checks. |

## Dependencies

- `numpy`
- `scipy` (used for elliptic integrals in the WPT coupling demo)
- `radia` (ESIM modules: `esim_cell_problem`, `esim_workpiece`, `esim_coupled_solver`, `esim_vtk_export`; RWG-EFIE C++ backend; WPT modules: `WPTCoupledSolver`, `compute_mutual_inductance`, `compute_coupling_coefficient`)

## References

- K. Hollaus, M. Kaltenbacher, J. Schoberl, "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation," IEEE Trans. Magnetics, 2025, DOI: 10.1109/TMAG.2025.3613932
- SAE J2954 - Wireless Power Transfer for Light-Duty Plug-In/Electric Vehicles and Alignment Methodology (referenced in WPT demos at 85 kHz)
