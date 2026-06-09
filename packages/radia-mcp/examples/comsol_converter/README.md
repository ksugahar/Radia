# comsol-converter -- TEAM-20 COMSOL <-> NGSolve cross-validation (P1)

Worked example of the comsol-converter pipeline: ONE magnetostatic problem
(TEAM Workshop Problem 20, a DC lifting electromagnet) solved on BOTH
COMSOL (via LiveLink for MATLAB) and NGSolve, from the same geometry /
material / excitation specification, and compared.

Built 2026-06-05 driving COMSOL 6.4 from an MCP-attached MATLAB
(`mphstart(2036)` onto the running "COMSOL with MATLAB" server).

## Files

| File            | Solver  | What it does |
|-----------------|---------|--------------|
| `t20_comsol.m`  | COMSOL  | Builds yoke + pole + circular coil + air, linear steel mur=1000, azimuthal coil current (NI=1000), Ampere's Law on the steel, stationary solve, Maxwell force on the pole. Run from MATLAB with LiveLink connected. Saves `t20_comsol.mph`. |
| `t20_ngsolve.py`| NGSolve | Same geometry/material/current via `radia_mcp.radia_ngsolve.solve_magnetostatic_Aform`; Maxwell-stress force on the pole surface. |

## Result

| Quantity      | COMSOL (83k tets) | NGSolve (57k tets) | agreement |
|---------------|-------------------|--------------------|-----------|
| `|B|` @ pole  | 0.681 T           | 0.690 T            | **1.3 %** |
| `Fz` on pole  | -7.76 N (mag 7.76)| 9.24 N             | ~16 %     |

The **field** agrees to 1.3 % -- the two independent solvers reproduce the
same magnetostatic solution from one specification (the converter's goal).
The **force** differs ~16 %: it is dominated by the 1.5 mm pole-yoke bottom
gap, which neither mesh resolves (steel maxh 4-8 mm >> 1.5 mm). Both bracket
the published TEAM-20 measurement of 8.1 N. To tighten the force, refine the
air gap on both sides (not the bulk iron).

## Notes

* This uses a LINEAR steel (mur=1000) and a CIRCULAR coil so the External
  Current Density stays solenoidal -- the cleanest first cross-check. The
  production TEAM-20 (nonlinear B-H, rectangular coil) is in
  `tests/test_force_xval.py::test_team20_static_force` (NGSolve, 7.9-8.23 N).
* The non-obvious COMSOL-LiveLink trip-wires hit while building this
  (mphinterp coord unit, `FreeSpace` ignoring mur, solenoidal-Je
  requirement, ...) are recorded in
  `interop` lab tips, topic **`livelink_matlab`**.
