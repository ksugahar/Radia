# Electric-machine examples

Open, analytic-gated rotating-machine FEA with **radia-ngsolve** — the 2D→3D motor
post-processing distilled from the open ONELAB/GetDP electric-machine study, built on
NGSolve's high-order curved elements.

| Example | Shows | Capabilities used |
|---------|-------|-------------------|
| [`cogging_skew_demo.py`](cogging_skew_demo.py) | Reluctance/cogging torque `tau(theta)` of a PM rotor in a salient stator; conversion to physical N·m for a real stack; and **skew ripple cancellation** | `eggshell_torque_2d`, `MachineScaling` (axial length × symmetry), `skew_average` ↔ `skew_factor` |
| [`validation_machine_torque_scaling.py`](validation_machine_torque_scaling.py) | Validation-class 2D sector torque to whole-machine N·m scaling, including mesh unit, stack length, and symmetry factor | `MachineScaling`, `torque_scaling_summary` |
| [`validation_slot_winding_spectrum.py`](validation_slot_winding_spectrum.py) | Validation-class harmonic spectrum from explicit slot sign tables, including fractional-slot layouts | `slot_table_winding_factor`, `integral_slot_winding_factor`, `mmf_harmonic_direction` |
| [`validation_pm_drive_speed_map.py`](validation_pm_drive_speed_map.py) | Validation-class PM-machine speed map across MTPA/FW/MTPV/infeasible regions plus short-circuit demagnetising trend | `pm_drive_speed_sweep`, `pm_drive_operating_point`, `field_weakening_operating_point`, `short_circuit_operating_point` |
| [`validation_pm_emf_constant_table.py`](validation_pm_emf_constant_table.py) | Validation-class PM back-EMF / torque-constant table with phase/line and peak/RMS conversions | `pm_flux_linkage_constants`, `pm_no_load_back_emf`, `dq_voltages`, `dq_torque` |
| [`validation_lamination_mu_eff_sweep.py`](validation_lamination_mu_eff_sweep.py) | Validation-class laminated-steel complex-permeability sweep from static limit into deep skin effect | `laminated_mu_eff` |
| [`validation_cross_saturation_flux_map.py`](validation_cross_saturation_flux_map.py) | Validation-class d-q cross-saturation flux map with reciprocity and incremental inductance rolloff | `incremental_inductance_matrix`, `dq_flux_torque` |
| [`validation_carter_magnetizing_sweep.py`](validation_carter_magnetizing_sweep.py) | Validation-class Carter air-gap opening sweep linked to AC-machine magnetizing inductance | `carter_coefficient`, `effective_air_gap`, `slotted_air_gap_permeance_factor`, `magnetizing_inductance_per_phase` |
| [`validation_mtpa_saliency_sweep.py`](validation_mtpa_saliency_sweep.py) | Validation-class MTPA current-angle sweep across non-salient PM, IPM, SynRM, and salient-PM cases | `mtpa_operating_point`, `dq_torque`, `dq_torque_components` |
| [`validation_pm_loadline_demag_sweep.py`](validation_pm_loadline_demag_sweep.py) | Validation-class PM load-line and irreversible-demag margin sweep across gap, magnet length, knee, and shape factors | `pm_circuit_loadline_operating_point`, `demag_operating_field`, `demag_margin` |
| [`validation_cogging_skew_plan.py`](validation_cogging_skew_plan.py) | Validation-class cogging order and one-slot-pitch skew planning table across integer/fractional slot-pole layouts | `cogging_skew_plan`, `cogging_torque_order`, `machine_symmetry_sector`, `skew_factor` |
| [`validation_nonlinear_bh_circuit_sweep.py`](validation_nonlinear_bh_circuit_sweep.py) | Validation-class nonlinear B-H magnetic-circuit sweep with secant/incremental permeability and MMF split checks | `magnetic_circuit_bh_operating_summary` |
| [`validation_nonlinear_bh_inductance_sweep.py`](validation_nonlinear_bh_inductance_sweep.py) | Validation-class nonlinear B-H winding-inductance sweep with secant/incremental inductance and constant-mu limit checks | `magnetic_circuit_bh_inductance_summary` |
| [`validation_maxwell_stress_traction.py`](validation_maxwell_stress_traction.py) | Validation-class Maxwell stress tensor / traction identities for normal, tangential, and oblique fields | `maxwell_stress_tensor_air`, `maxwell_traction_air`, `maxwell_traction_summary`, `air_gap_maxwell_pressure` |
| [`validation_air_gap_force_sweep.py`](validation_air_gap_force_sweep.py) | Validation-class air-gap Maxwell pressure / holding-force sweep linked to nonlinear magnetic-circuit B | `air_gap_maxwell_pressure`, `air_gap_holding_force`, `air_gap_force_summary`, `magnetic_circuit_bh_operating_summary` |
| [`validation_air_gap_shear_torque.py`](validation_air_gap_shear_torque.py) | Validation-class air-gap Maxwell shear stress to motor torque identity | `air_gap_shear_stress`, `air_gap_shear_torque`, `air_gap_shear_torque_summary`, `maxwell_traction_summary` |
| [`validation_sampled_air_gap_shear_torque.py`](validation_sampled_air_gap_shear_torque.py) | Validation-class sampled air-gap `Br(θ), Bt(θ)` Maxwell shear torque integral, with uniform and harmonic gates | `air_gap_shear_torque_from_angle_samples`, `air_gap_shear_torque` |
| [`validation_force_resultant_torque.py`](validation_force_resultant_torque.py) | Validation-class force-row resultant and pivot torque identities for patch, element, pressure, or nodal loads | `force_moment_resultant_summary` |
| [`validation_planar_maxwell_contour_force.py`](validation_planar_maxwell_contour_force.py) | Validation-class 2D Maxwell-stress contour force, per unit depth, with pole-face pressure and closed-contour cancellation checks | `maxwell_line_segment_force_2d`, `maxwell_contour_force_2d`, `air_gap_maxwell_pressure` |
| [`validation_planar_lorentz_block_force.py`](validation_planar_lorentz_block_force.py) | Validation-class 2D Lorentz block force, per unit depth, checked against Ampere's two-wire force law | `planar_lorentz_force_summary`, `parallel_wire_lorentz_force_summary`, `two_wire_force_per_length` |
| [`validation_coenergy_torque_angle_sweep.py`](validation_coenergy_torque_angle_sweep.py) | Validation-class virtual-work torque from fixed-current coenergy angle samples | `coenergy_torque_from_angle_samples`, `coenergy_torque_summary` |
| [`validation_virtual_work_force_displacement_sweep.py`](validation_virtual_work_force_displacement_sweep.py) | Validation-class virtual-work force from displacement energy/coenergy samples, with the fixed-current vs fixed-flux sign gate | `virtual_work_force_from_displacement_samples`, `virtual_work_force_summary` |
| [`validation_virtual_work_symmetric_pair_force.py`](validation_virtual_work_symmetric_pair_force.py) | Validation-class virtual-work force from a matched `x0 +/- h` energy pair, including even-energy residual | `virtual_work_symmetric_pair_force_summary`, `virtual_work_force_summary` |
| [`validation_torque_ripple_harmonic_budget.py`](validation_torque_ripple_harmonic_budget.py) | Validation-class three-phase back-EMF harmonic pair budget for 6k torque ripple, including pitch/skew filtering | `three_phase_torque_ripple_pair_table`, `three_phase_torque_ripple_harmonics`, `skewed_winding_factor` |

```powershell
python cogging_skew_demo.py
python validation_machine_torque_scaling.py
python validation_slot_winding_spectrum.py
python validation_pm_drive_speed_map.py
python validation_pm_emf_constant_table.py
python validation_lamination_mu_eff_sweep.py
python validation_cross_saturation_flux_map.py
python validation_carter_magnetizing_sweep.py
python validation_mtpa_saliency_sweep.py
python validation_pm_loadline_demag_sweep.py
python validation_cogging_skew_plan.py
python validation_nonlinear_bh_circuit_sweep.py
python validation_nonlinear_bh_inductance_sweep.py
python validation_maxwell_stress_traction.py
python validation_air_gap_force_sweep.py
python validation_air_gap_shear_torque.py
python validation_sampled_air_gap_shear_torque.py
python validation_force_resultant_torque.py
python validation_planar_maxwell_contour_force.py
python validation_planar_lorentz_block_force.py
python validation_coenergy_torque_angle_sweep.py
python validation_virtual_work_force_displacement_sweep.py
python validation_virtual_work_symmetric_pair_force.py
python validation_torque_ripple_harmonic_budget.py
```

What it validates (no commercial tool needed to run or check):

1. **`tau(theta)`** via the weighted-Maxwell-stress (eggshell) torque on an air-gap
   band — the rotor is turned by rotating the magnetisation (no remesh).
2. **`MachineScaling`** — the 2D, per-unit-depth, per-sector torque becomes a physical
   whole-machine N·m via the ONELAB `AxialLength × SymmetryFactor` scaling.
3. **`skew_average`** (FE multi-slice skew) reduces a harmonic of order *n* by exactly
   `sinc(n·skew/2) = skew_factor(n, skew)` — checked here against the analytic factor on
   the **real FE curve**; a half-period skew nulls the order-2 ripple.

The public example is analytic / self-consistency gated; the same geometry is
additionally cross-checked against independent lab tooling internally (not shipped).

## Related

- Iron / eddy / excess loss from a rotor sweep: `radia_mcp.radia_ngsolve.coreloss`
  (classical-eddy term = the ONELAB lamination homogenisation; full Bertotti + harmonic).
- Laminated-steel AC homogenisation: `radia_mcp.radia_ngsolve.solve.laminated_mu_eff`
  (complex in-plane permeability with fill factor and lamination skin effect).
- Saturated-machine small-signal maps:
  `radia_mcp.radia_ngsolve.solve.incremental_inductance_matrix`
  (d-q tangent matrix, reciprocity, cross-saturation).
- Winding / leakage / magnetising analytics: `radia_mcp.radia_ngsolve.solve`
  (`winding_factor`, `skew_factor`, `slot_leakage_inductance`,
  `magnetizing_inductance_per_phase`, `effective_air_gap` = Carter).
- The open GetDP reference these mirror: `radia_mcp.motor.onelab_knowledge`
  (`magstadyn_source`, `twod_corrections`).
