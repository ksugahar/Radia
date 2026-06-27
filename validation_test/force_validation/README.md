# Force and torque validation corpus

Small public-safe electromagnetic force validations that are not tied to a
machine topology.  These scripts are meant as readable checks before moving
the same force identities into FEM, BEM, or CAD/mesh workflows.

The human-facing, result-saved notebook layer is
`docs/force_validation/force_validation.ipynb`. This directory is the
executable validation surface; each script refreshes its
adjacent JSON summary with timestamp and runtime version metadata.

| Example | Shows |
|---|---|
| [`validation_parallel_wire_virtual_work_force.py`](validation_parallel_wire_virtual_work_force.py) | Ampere two-wire Lorentz force and fixed-current coenergy gradient give the same force per unit length |
| [`validation_virtual_work_force_sweep_audit.py`](validation_virtual_work_force_sweep_audit.py) | Coenergy-vs-displacement sweeps keep finite-difference stencils, force-gradient estimates, and reference-force errors |
| [`validation_coenergy_torque_table_consistency.py`](validation_coenergy_torque_table_consistency.py) | Torque-angle tables and coenergy-angle tables agree under fixed-current virtual work, including nonzero mean work terms |
| [`validation_torque_waveform_comparison.py`](validation_torque_waveform_comparison.py) | Periodic torque tables can be compared by mean drift, sample error, and harmonic ripple deltas |
| [`validation_maxwell_contour_segment_balance.py`](validation_maxwell_contour_segment_balance.py) | Closed 2D Maxwell stress contours expose large local segment forces while the symmetric net force cancels |
| [`validation_torque_waveform_health.py`](validation_torque_waveform_health.py) | A single periodic torque table is checked by mean torque, RMS ripple ratio, dominant harmonic, and harmonic variance budget |

```powershell
python validation_test/force_validation/validation_parallel_wire_virtual_work_force.py
python validation_test/force_validation/validation_virtual_work_force_sweep_audit.py
python validation_test/force_validation/validation_coenergy_torque_table_consistency.py
python validation_test/force_validation/validation_torque_waveform_comparison.py
python validation_test/force_validation/validation_maxwell_contour_segment_balance.py
python validation_test/force_validation/validation_torque_waveform_health.py
```
