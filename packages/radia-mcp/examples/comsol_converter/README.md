# TEAM-20 cross-validation (P1)

Worked example: ONE magnetostatic problem (TEAM Workshop Problem 20, a DC
lifting electromagnet) solved with the radia-ngsolve stack and compared against
an independent reference solve, from the same geometry / material / excitation
specification.

## Files

| File            | Solver  | What it does |
|-----------------|---------|--------------|
| `t20_ngsolve.py`| NGSolve | Builds yoke + pole + circular coil + air, linear steel mur=1000, azimuthal coil current (NI=1000) via `radia_mcp.radia_ngsolve.solve_magnetostatic_Aform`; Maxwell-stress force on the pole surface. |

The reference solve (same spec, an independent solver) is kept lab-private; only
its result numbers are retained below as a stored regression reference.

## Result

| Quantity      | reference (83k tets) | NGSolve (57k tets) | agreement |
|---------------|----------------------|--------------------|-----------|
| `|B|` @ pole  | 0.681 T              | 0.690 T            | **1.3 %** |
| `Fz` on pole  | -7.76 N (mag 7.76)   | 9.24 N             | ~16 %     |

The **field** agrees to 1.3 % -- the two independent solvers reproduce the
same magnetostatic solution from one specification. The **force** differs ~16 %:
it is dominated by the 1.5 mm pole-yoke bottom gap, which neither mesh resolves
(steel maxh 4-8 mm >> 1.5 mm). Both bracket the published TEAM-20 measurement of
8.1 N. To tighten the force, refine the air gap on both sides (not the bulk iron).

## Notes

* This uses a LINEAR steel (mur=1000) and a CIRCULAR coil so the external
  current density stays solenoidal -- the cleanest first cross-check. The
  production TEAM-20 (nonlinear B-H, rectangular coil) is in
  `validation/force/validate_force_xval.py::validate_team20_static_force`
  (NGSolve, 7.9-8.23 N).
