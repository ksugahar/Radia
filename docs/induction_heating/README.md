# Induction heating: electromagnetic loss and thermal response

Radia computes electromagnetic workpiece loss and transient temperature fields.
The **Radia IH MCP tool family** owns current input contracts and operating
instructions. Human–AI collaboration uses Cubit for geometry/mesh and the Radia
Simulink library for application models.
These executed notebooks explain numerical evidence, not production controls.
Numerical evidence and its checked JSON remain under `validation_test/`.

## Choose your next step

| Purpose | Entry point |
| --- | --- |
| Inspect temperature and conservative source transfer | [Thermal feature notebook](axisymmetric_p2_thermal.ipynb) |
| Operate the supported analysis with AI | [MCP setup and tools](../../packages/radia-mcp/README.md) |
| Prepare geometry and labeled mesh | [Cubit workflow](../cubit_mesh_export/README.md) |
| Collaborate on an application model | [Radia Simulink library](../../matlab/README.md) |
| Read executable field-transfer examples | [Thermal transfer tests](../../tests/test_ih_thermal_transfer.py) |
| Inspect electromagnetic/thermal operator checks | [IH physics validation](../../validation_test/induction_heating/README.md#simulink-operator-physics-golden) |

The manufactured thermal transient, field-transfer contracts and real
electromagnetic operator check answer different questions. Use the validation
case matching your quantity; their presence alone does not qualify a complete
coupled application on a new geometry.

## Executed demonstrations

- [Axisymmetric P2/Q2 thermal verification](axisymmetric_p2_thermal.ipynb):
  an axis-touching cylinder, quadratic manufactured transient, conservative
  P1-source/Q2-temperature transfer, and saved mesh/temperature WebGUI scenes.
  Separately recorded native MEX and tracked Simulink comparisons are displayed.
- [`induction_heating_demo_showcase.ipynb`](induction_heating_demo_showcase.ipynb)
  covers the ESIM/Bessel comparison and historical nonlinear surface-impedance
  study, not nonlinear BH support in the native runtime.
- [Eddy–thermal architecture and acceptance](../IH_THERMAL_WORKFLOW.md).

## Verification boundary

Thermal order 2 uses standard NGSolve H1 (P2 on triangles, Q2 on quads) and
the physical `2*pi*r` measure, not the magnetic axis-reduced Henrotte basis.
It does not request post-load geometry curving. Higher-order state contains
FE coefficients; monitor extrema are sampled physical temperatures.

The declared independent IH workpiece FEM/BEM comparison meets its 2% bound.
That physical accuracy criterion is application-specific and distinct from
the much tighter native/reference arithmetic check. See the
[numerical report](../../validation_test/eddy_current_analytical_validation/IH_REACTION_POWER_2026-09-15.md).

P2/Q2 and tracked-model changes are on main through
[PR #280](https://github.com/ksugahar/Radia/pull/280). This is not release
completion: unmocked production CAD/VOL-to-EM-to-heat acceptance, full-window
UI inspection and the four-host exact-package gate remain separate requirements.
The received-case sixfold discrepancy is not asserted resolved by these tests.

## Evidence ownership

`validation_test/induction_heating/` owns the native thermal records;
`validation_test/eddy_current_analytical_validation/` owns the independent
workpiece physical comparison; `validation_test/ih_esim_benchmark/` owns the
ESIM benchmark corpus. The notebooks show those claims and limits without
replacing the checked evidence or the IH MCP operating manual.
