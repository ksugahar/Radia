# TEAM Workshop Problem 36: radia-ih validation

This retained validation implements the public TEAM Workshop Problem 36
magneto-thermal induction-heating benchmark.  The benchmark definition is the
[original COMPUMAG problem statement](https://www.compumag.org/wp/wp-content/uploads/2021/07/problem-36.pdf).

The runner uses the radia axisymmetric Henrotte space for the time-harmonic
electromagnetic solve and a first-order axisymmetric thermal solve.  It includes
the original temperature tables for electrical resistivity, thermal
conductivity, and heat capacity; the field- and temperature-dependent
permeability law; convection; radiation; and the 250 s coupling history.

The electromagnetic and thermal meshes are independently generated and saved
as Netgen `.vol` files in `C:\temp`.  They intentionally use different radial
partitions.  Temperature is mapped from the thermal mesh to the electromagnetic
mesh, while volumetric Joule loss is mapped back and normalized by the full
axisymmetric `2*pi*r` integral.  The result artifact records both topology
hashes and the mapping power error.

Fast development smoke:

```powershell
python validation_test\induction_heating\team36_axisymmetric\run_radia_ih.py `
  --profile smoke --mesh-profile baseline
```

Full 250 s retained validation:

```powershell
python validation_test\induction_heating\team36_axisymmetric\run_radia_ih.py `
  --profile validation --mesh-profile baseline
```

Use `--mesh-profile refined` for the refinement companion run.  Long retained
runs belong on an MDX or Hibino worker; the 25 s smoke is suitable for a local
development machine.

The result has two distinct gates:

- `accepted_for_solver_execution`: the complete public benchmark contract,
  noncoincident mesh coupling, nonlinear convergence, and 250 s history pass.
- `accepted_for_cross_validation`: an independent result with the same
  geometry/material/excitation identity is also supplied and the selected
  named 250 s observables agree. The gate reads those values from the Radia
  history itself; an optional `radia_value` in the reference must match it.

A standalone radia-ih run never promotes itself to cross-validation or MCP
learning evidence.  This prevents a numerically plausible temperature history
from being mistaken for an independent comparison.
