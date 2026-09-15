# Induction heating: electromagnetic loss and thermal response

Radia computes electromagnetic workpiece loss and transient temperature fields.
The **Radia IH MCP tool family** owns current input contracts and operating
instructions; the Radia Simulink library is the formal human interface.
These executed notebooks explain numerical evidence, not production controls.
Numerical evidence and its checked JSON remain under `validation_test/`.

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
