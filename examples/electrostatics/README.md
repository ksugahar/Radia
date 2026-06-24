# Electrostatics examples

Closed-form validation examples for capacitance, image charges, dielectric
response, and electrostatic force identities.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_parallel_plate_electrostatic_force.py`](validation_parallel_plate_electrostatic_force.py) | Validation-class parallel-plate force identity: capacitance energy, pressure, Maxwell traction, and fixed-voltage 1/gap^2 scaling | `parallel_plate_capacitor_energy_force`, `electrostatic_traction_summary`, `EPS0` |
| [`validation_capacitance_gradient_force.py`](validation_capacitance_gradient_force.py) | Validation-class capacitance-gradient electrostatic force with signed gap/closing coordinates and sphere-ground finite-difference attraction | `capacitance_gradient_force_summary`, `parallel_plate_capacitor_energy_force`, `sphere_above_plane_capacitance` |
| [`validation_sphere_ground_plane_capacitance.py`](validation_sphere_ground_plane_capacitance.py) | Validation-class conducting sphere above a grounded plane: image-series truncation, near-plane capacitance growth, far-field isolated-sphere limit, and finite-difference attraction check | `sphere_above_plane_capacitance`, `EPS0` |
| [`validation_electrostatic_maxwell_traction.py`](validation_electrostatic_maxwell_traction.py) | Validation-class electrostatic Maxwell stress / traction identities for normal, tangential, and oblique electric fields | `electrostatic_stress_tensor`, `electrostatic_traction_summary`, `EPS0` |

```powershell
python validation_parallel_plate_electrostatic_force.py
python validation_capacitance_gradient_force.py
python validation_sphere_ground_plane_capacitance.py
python validation_electrostatic_maxwell_traction.py
```

The examples use public textbook formulas and are suitable as solver regression
targets without depending on private benchmark data.
