# Force And Torque Validation Examples

Small public-safe electromagnetic force validations that are not tied to a
machine topology.  These examples are meant as readable checks before moving
the same force identities into FEM, BEM, or CAD/mesh workflows.

| Example | Shows |
|---|---|
| [`validation_parallel_wire_virtual_work_force.py`](validation_parallel_wire_virtual_work_force.py) | Ampere two-wire Lorentz force and fixed-current coenergy gradient give the same force per unit length |
| [`validation_virtual_work_force_sweep_audit.py`](validation_virtual_work_force_sweep_audit.py) | Coenergy-vs-displacement sweeps keep finite-difference stencils, force-gradient estimates, and reference-force errors |
| [`validation_coenergy_torque_table_consistency.py`](validation_coenergy_torque_table_consistency.py) | Torque-angle tables and coenergy-angle tables agree under fixed-current virtual work, including nonzero mean work terms |
| [`validation_torque_waveform_comparison.py`](validation_torque_waveform_comparison.py) | Periodic torque tables can be compared by mean drift, sample error, and harmonic ripple deltas |

```powershell
python validation_parallel_wire_virtual_work_force.py
python validation_virtual_work_force_sweep_audit.py
python validation_coenergy_torque_table_consistency.py
python validation_torque_waveform_comparison.py
```
