# Electrostatics examples

Closed-form validation examples for capacitance, image charges, dielectric
response, and electrostatic force identities.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_sphere_ground_plane_capacitance.py`](validation_sphere_ground_plane_capacitance.py) | Validation-class conducting sphere above a grounded plane: image-series truncation, near-plane capacitance growth, far-field isolated-sphere limit, and finite-difference attraction check | `sphere_above_plane_capacitance`, `EPS0` |
| [`validation_electrostatic_maxwell_traction.py`](validation_electrostatic_maxwell_traction.py) | Validation-class electrostatic Maxwell stress / traction identities for normal, tangential, and oblique electric fields | `electrostatic_stress_tensor`, `electrostatic_traction_summary`, `EPS0` |

```powershell
python validation_sphere_ground_plane_capacitance.py
python validation_electrostatic_maxwell_traction.py
```

The examples use public textbook formulas and are suitable as solver regression
targets without depending on private benchmark data.
