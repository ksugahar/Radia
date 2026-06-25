# Force Validation Examples

Small public-safe electromagnetic force validations that are not tied to a
machine topology.  These examples are meant as readable checks before moving
the same force identities into FEM, BEM, or CAD/mesh workflows.

| Example | Shows |
|---|---|
| [`validation_parallel_wire_virtual_work_force.py`](validation_parallel_wire_virtual_work_force.py) | Ampere two-wire Lorentz force and fixed-current coenergy gradient give the same force per unit length |

```powershell
python validation_parallel_wire_virtual_work_force.py
```
