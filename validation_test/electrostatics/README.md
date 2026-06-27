# Electrostatics validation corpus

Closed-form validation scripts for capacitance, image charges, dielectric
response, and electrostatic force identities.

The human-facing, result-saved notebook layer is
`docs/electrostatics/electrostatics.ipynb` plus the source archive
`docs/electrostatics/electrostatics_examples_archive.ipynb`.  This directory is
the executable validation surface; each script refreshes its adjacent JSON
summary with timestamp and runtime version metadata.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_parallel_plate_electrostatic_force.py`](validation_parallel_plate_electrostatic_force.py) | Validation-class parallel-plate force identity: capacitance energy, pressure, Maxwell traction, and fixed-voltage 1/gap^2 scaling | `parallel_plate_capacitor_energy_force`, `electrostatic_traction_summary`, `EPS0` |
| [`validation_coaxial_capacitor_force.py`](validation_coaxial_capacitor_force.py) | Validation-class coaxial capacitor force identity: cylindrical capacitance, Maxwell pressure, and inner/outer radius capacitance-gradient forces | `coaxial_capacitor_energy_force`, `capacitance_gradient_force_summary` |
| [`validation_capacitance_gradient_force.py`](validation_capacitance_gradient_force.py) | Validation-class capacitance-gradient electrostatic force with signed gap/closing coordinates and sphere-ground finite-difference attraction | `capacitance_gradient_force_summary`, `parallel_plate_capacitor_energy_force`, `sphere_above_plane_capacitance` |
| [`validation_sphere_ground_plane_capacitance.py`](validation_sphere_ground_plane_capacitance.py) | Validation-class conducting sphere above a grounded plane: image-series truncation, near-plane capacitance growth, far-field isolated-sphere limit, and finite-difference attraction check | `sphere_above_plane_capacitance`, `EPS0` |
| [`validation_electrostatic_maxwell_traction.py`](validation_electrostatic_maxwell_traction.py) | Validation-class electrostatic Maxwell stress / traction identities for normal, tangential, and oblique electric fields | `electrostatic_stress_tensor`, `electrostatic_traction_summary`, `EPS0` |

```powershell
python validation_test/electrostatics/validation_parallel_plate_electrostatic_force.py
python validation_test/electrostatics/validation_coaxial_capacitor_force.py
python validation_test/electrostatics/validation_capacitance_gradient_force.py
python validation_test/electrostatics/validation_sphere_ground_plane_capacitance.py
python validation_test/electrostatics/validation_electrostatic_maxwell_traction.py
```

The validations use public textbook formulas and are suitable as solver regression
targets without depending on private benchmark data.
