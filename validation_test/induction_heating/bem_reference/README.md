# IH BEM Reference Validation

This directory owns runnable validation drivers and durable result data for
the induction-heating BEM coupling route. Numerical implementation remains in
`src/radia`; validation must not keep a private solver copy.

- `calc_inductance.py` is a compatibility launcher for
  `radia.panels.calc_inductance`.
- `sweep_results.json` is a historical workpiece mesh-convergence record.
  Its driver depended on retired panel research solvers and is no longer distributed;
  this record does not certify the current solver.
- Generated `.msh` and `.geo` visualization files belong to the run directory
  and are not repository evidence.

The human production interface is the masked IH block in the Radia Simulink
library. Cubit supplies CAD, labels, and `.vol` solver meshes; it is not the IH
application GUI.
