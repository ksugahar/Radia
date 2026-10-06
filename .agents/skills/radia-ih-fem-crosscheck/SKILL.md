---
name: radia-ih-fem-crosscheck
description: Verify Radia-IH electromagnetic heating results by comparing PEEC-BEM with an independent FEM solution for the same case. Use for IH solver, material, mesh, Zs, excitation or student Simulink changes and reported IH powers; interface-only checks do not establish numerical accuracy.
---

# Radia-IH independent FEM comparison

For every new or changed IH case, run PEEC-BEM and an independent FEM reference
before calling its heating result verified. A passing unit test, power identity,
or a previously validated copper case does not certify another material or mesh.
If a comparable FEM route is unavailable, record `not-verified` and its reason;
do not substitute a different geometry or silently approve the result.

## Make the comparison physical

Record geometry dimensions and units, mesh/input/source SHA-256, frequency,
complex coil current and peak/RMS convention, conductivity, permeability or the
exact B-H curve, temperature, boundary conditions and external-domain treatment.
Use the same coil geometry and excitation. Resolve skin depth in a volumetric
FEM reference and test its mesh and exterior-domain convergence. A 2D reference
is admissible only for axisymmetric geometry, excitation and material data.

Distinguish two questions explicitly:
- FEM-SIBC with the same Zs checks the BEM discretization and coupling. It does
  not validate the SIBC or a nonlinear ESIM constitutive approximation.
- Volumetric eddy-current FEM with resolved skin and matching material data
  checks that approximation too. `ih_axisym_coupled` uses linear mu(T), not B-H
  saturation. Use `radia.panels.calc_axisym_volumetric` for its supported
  elementwise B-H iteration; do not relabel mu(T) as magnetic nonlinearity.

Arbitrary supplied element Zs has no automatic volumetric-material equivalent.
Require that equivalence or label the comparison as same-boundary-law only.
Frozen-field loss recalculation is sensitivity analysis, not a reference solve.

## Execute and decide

Develop on 100; run numerical tests in an idle LAB/compute-host job under
`C:\temp`. LAB uses fixed wheels, never editable installs. Hash transferred
inputs and recovered results. Preserve borrowed MATLAB/COMSOL sessions. Use
official MATLAB MCP for tracked Simulink model checks and edits.

Find the applicable current route in:
- `src/radia/panels/calc_inductance.py`: weak PEEC/BEM and local surface heat.
- `src/radia/panels/calc_fem_kelvin.py`: independent FEM-SIBC with coil excitation.
- `src/radia/panels/calc_axisym_volumetric.py`: supported volumetric axisym case.
- `validation_test/induction_heating/fem_sibc_through_hole.py` and
  `validation_test/eddy_current_analytical_validation/test_ih_sibc_production.py`:
  bounded examples and acceptance contracts, not substitutes for the user's case.
- `docs/induction_heating/PANEL_ESIM.md`: local Zs and unsupported topology.

Use at least two refinements for each numerical route. Set tolerances before
viewing results, justified by the case's error budget and refinement evidence.
Do not inherit the copper regression tolerance for magnetic materials.
Compare total workpiece power, reaction loss/port impedance, spatial q_surf and
complex tangential H on common physical regions or sample points. Preserve phase,
area/volume weights and peak-phasor factor 1/2. Do not compare raw DOF vectors
from different meshes. Check heat transfer preserves integrated electromagnetic
power and separately check each solver's residual and nonlinear convergence.

Write a case report with both commands, exact versions/SHAs, inputs, meshes,
convergence tables, predeclared tolerances, measured differences and per-gate
pass/fail. Mark unavailable comparisons `not-performed`, incompatible ones
`not-comparable`, and failures `failed`; none is verified. For discrepancies,
investigate scaling, material law, current distribution, topology and boundary
conditions before changing tolerances. Keep commercial reference provenance and
private numerical results outside public artifacts.

## Student Simulink acceptance

For element Zs verify mesh-bound JSON, finite passive values, frequency match,
P1 weak genus-0 support, weighted stiffness and locally integrated heat. Zs is
fixed during the current linear-response simulation; changing it requires a new
assembly and FEM crosscheck. Watch the Zs file content hash so edits cannot reuse
stale operators. Save and reopen the model, run the actual native Eddy/Thermal
path, and check integrated heat. Report interface acceptance separately from
the independent electromagnetic accuracy verdict.
