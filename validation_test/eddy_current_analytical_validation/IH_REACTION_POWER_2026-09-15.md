# IH SIBC reaction and heating: repository acceptance

Status: **implemented and accepted within the scope below**.
No Takahashi inputs or sixfold-ratio target are used. Existing installed
Radia environments were not repointed. Solver work ran on idle hibino with
staged Python sources and an identified native runtime.

## Accepted model and criteria

Prescribed circular current at 1 kHz, copper (sigma=5.8e7 S/m, mu_r=1),
uniform planar SIBC, cylinder radius 25 mm and height 25 mm, with/without a
7.5 mm axial bore. Coil centre radius is 30 mm, section radius 0.5 mm.
FEM integrates uniform current over the disk; BEM uses degree-five disk
cubature (seven rings of 720 straight segments each).

Independent **2D axisymmetric A-form FEM**, order 3, varies the finite
exterior from 150 to 300 mm. The BEM is a 3D **surface** discretization of
the same axisymmetric geometry, at mesh sizes 3, 2 and 1.5 mm. No 3D volume
FEM or VOL Save/reload/Curve workaround is involved. These are explicit
in-memory OCC/Netgen verification meshes, not production export fixtures.

The user's acceptance ceiling is **2%** for loss, reaction power and the
relative L2 error of the **complex tangential field**, not just its magnitude.
This is an IH-specific criterion, not a repository-wide tolerance. Accelerator
applications require their own observable-specific acceptance (0.5% or tighter).

| Finest check | Solid | Bored | Limit |
|---|---:|---:|---:|
| Surface loss versus FEM | 0.6032% | 0.6905% | 2% |
| Reaction power versus FEM | 0.5205% | 1.1781% | 2% |
| Local complex tangential field, relative L2 | 1.5573% | 1.5499% | 2% |
| BEM independent power imbalance | 0.0832% | 0.4909% | 1% |
| Last BEM loss refinement change | 0.5900% | 0.6712% | 1.5% |
| FEM exterior enlargement change | 0.3296% | 0.3295% | 0.5% |

Heat projection error is below 1.3e-14; FEM-SIBC power imbalance below
2.1e-11. All gates pass. The field limit was tightened from exploratory
5% to 2% using unchanged saved solver measurements. The --recheck route
checks numerical source hashes, refuses relaxed limits and records the
separate acceptance-gate hash.

## Production changes

- Complete SIBC reciprocity replaces the magnetic-only reaction. With peak
  exp(+i omega t), incident potential psi and total potential phi:

      delta_L = integral(phi n.B_inc)/I^2
                + Zs/(i omega I^2) * psi.T K phi

  No conjugation is used. Reaction is independent of dissipated power;
  resistance is not manufactured by reversing the loss.
- The genus-1 reaction and local heat use the **full loop field**. Its
  current.A contribution retains the cut term. Plain phi alone is not used
  as a substitute for the multivalued field.
- NGSolve assembles a positive lumped P1 heat projection from the solved
  field. The globally rescaled incident Biot-Savart heat pattern is removed.
- A runtime imbalance above 10%, nonfinite/negative power, inconsistent
  projection, incomplete vertex mapping, failed SOL export or failed
  requested GMSH export raises. The 10% safety ceiling is NOT the tighter
  acceptance accuracy. GMSH includes q_surf_W_per_m2.
- Higher-order weak BEM and per-panel ESIM postprocessing fail fast pending
  mapped weighted surface forms. This EM restriction does not restrict
  standard **thermal H1 order 2**. Unsupported handle modes remain rejected.
- MCP guidance describes the corrected route and withdraws the old calibrated
  heat-pattern and plain-phi loop-inductance recommendations.

The two-term reciprocity expression is also given in equation 21 of
[Luo and Di Rienzo, High-order surface impedance boundary conditions in
three-dimensional boundary element modeling of eddy-current problems](https://www.sciencedirect.com/science/article/pii/S0955799726002420).
The independent computations above, not the citation alone, establish acceptance.

## Regression and evidence

- validate_ih_sibc_closure.py --output result.json: FEM/BEM comparisons,
  mesh/exterior refinement, runtime and source hashes.
- test_ih_sibc_production.py: actual workpiece extraction, incident projection,
  BIE/loop solve, reaction, GMSH write and SOL write/read with integral check.
  Only the file loader is replaced by an in-memory fixture; solved fields
  are not substituted. Covers solid dense/HACApK and bored dense.
- validation_test/bem/test_loop_extension_ring.py: thin-ring analytic current,
  frozen-subsystem equivalence, passive full reaction and finite local heat.
- tests/test_sibc_reaction_contract.py: complex bilinear pairing, current and
  gauge invariance, invalid inputs and refusal of inconsistent power.
- Existing loop CLI, incident-potential, application/manifest and MCP heat
  contracts are included in focused regression checks.
- ih_sibc_closure_20260915.json contains accepted numerical evidence.
  ih_sibc_production_20260915.json records workpiece-stage output checks.
- Focused local contracts: 60 passed, one CAD-dependent case deselected;
  remote production/loop regressions: nine passed. GMSH 4.15.2 reopened all
  three saved meshes and verified node coverage, finite fields and positive heat.

The earlier full-conductor axisymmetric FEM baseline remains complementary:
solid/bored p=2,3 power imbalance is below 2.8e-9. It is not substituted for
the matching **FEM-SIBC** comparison.

## Boundaries

### Simulink runtime acceptance and remaining delivery gates

The existing `tests/matlab/test_native_ih_sfun_integration.m` ran through
Python MATLAB Engine on LAB R2026a Update 3: eight passed, zero failed or
incomplete, recorded in `ih_simulink_integration_20260915.json`. This covers
closed-loop heating, current scaling, angle-dependent heat, conservative
rotation, temperature feedback, distinct field spaces, configuration loading
and recovery after a singular Eddy solve. It uses scratch model harnesses,
not acceptance of the tracked `radia_ih.slx` or a complete geometry workflow.

Do not call the full IH application finished: the tracked model's official
Toolkit acceptance, an unmocked strict-label geometry/configuration run and
the native axisymmetric operator route remain open. The checked-in cylinder
volume fixture currently labels its boundary `outer`, not required `sibc`;
it cannot silently stand in for a production-contract fixture. Native geometry
assembly still requires P1 and a 3D volume. Main integration and distribution
have not occurred.

### Additional native thermal evidence

`validation_test/induction_heating/validate_ih_mex_radial.m` now independently
assembles 33-node radial axisymmetric P1 operators for solid and bored cylinders.
Both execute 1000 Eddy/Thermal MEX steps with prescribed uniform heating,
an insulated inner boundary (or symmetry axis) and outer convection. The
independent MATLAB backward-Euler solve differs by at most 2.7e-11 K;
the analytic steady temperature-rise relative L2 error is at most 0.002901%.
Power imbalance is below 4.3e-12, and native handles return from zero to zero.
The JSON `ih_mex_radial_20260915.json` records the actual binary build identity.
The binary SHA-256 was checked; native sources and CMake have no diff between
its build commit and the IH correction commit. Execution used the existing
LAB MATLAB R2026a Update 3 via the official Python Engine, without restarting it.

This radial test establishes spatial thermal MEX behavior. Its prescribed source does not validate EM mapping,
inner-wall convection or production VOL export. Reproduce after `radia.setup`
by calling `validate_ih_mex_radial` with an output JSON path from MATLAB.

### BEM-to-thermal-MEX connection evidence

`prepare_ih_mex_chain.py` executes the corrected weak BEM workpiece stage,
loads its actual saved heat field and calls the production
`_assemble_thermal_operators`. `validate_ih_mex_chain.m` consumes these operators
in native Eddy/Thermal handles. Solid/bored cylinders use 565/622 thermal DOFs
and 100 steps at 100 A. Heat-vector transfer error is zero, integrated power
is preserved to roundoff, and maximum temperature-rise relative error against
independently assembled NGSolve forms is 2.06e-7 (0.0000206%). Maximum absolute
temperature difference is 2.02e-9 K; handles return from zero to zero.

Reproduce with Python `prepare_ih_mex_chain.py --output C:/temp/ih-mex-chain`
then MATLAB `validate_ih_mex_chain('C:/temp/ih-mex-chain', outputJson)` after
configuring the selected MEX runtime. `ih_mex_chain_20260915.json` records results
and source/build identity. These are coarse coupling tests, not a new 2% EM
mesh-convergence certificate. The existing fine-mesh FEM-SIBC comparisons own
that accuracy claim. File loading is replaced by an in-memory fixture, and
native-config fields are packaged explicitly: full geometry/configuration
orchestration, CAD and strict VOL checks are not certified by this test.
The current production thermal assembler accepts only 3D volumes, so this
axisymmetric geometry uses a 3D discretization here; it does not imply that
3D physics is necessary or that the axisymmetric native assembly route exists.

This is not a certificate for arbitrary frequency, material or geometry.
Nonlinear BH accuracy, spatially varying impedance, higher-order/curved
mapped postprocessing, strong coil-current feedback, and native Simulink
MEX operators are outside this acceptance. The separate strong route's
post-convergence loop treatment must not inherit this claim. Coil CAD and
production VOL export/check gates are separate from the in-memory
workpiece-stage tests. The received case's exact sixfold ratio is neither
asserted fixed nor used as an acceptance input.

This change is a tested repository implementation, not a PyPI publication
or deployment into existing student environments.
