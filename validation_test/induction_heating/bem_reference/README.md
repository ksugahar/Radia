# IH BEM Reference Validation

This directory owns runnable validation drivers and durable result data for
the induction-heating BEM coupling route. Numerical implementation remains in
`src/radia`; validation must not keep a private solver copy.

- `calc_inductance.py` is a compatibility launcher for
  `radia.panels.calc_inductance`.
- `sweep_bem_coupled_maxh.py` generates `sweep_results.json` for the workpiece
  mesh-convergence study.
- Generated `.msh` and `.geo` visualization files belong to the run directory
  and are not repository evidence.

The human production interface is the masked IH block in the Radia Simulink
library. Cubit supplies CAD, labels, and `.vol` solver meshes; it is not the IH
application GUI.
