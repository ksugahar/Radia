# Electric-machine examples

Open, analytic-gated rotating-machine FEA with **radia-ngsolve** — the 2D→3D motor
post-processing distilled from the open ONELAB/GetDP electric-machine study, built on
NGSolve's high-order curved elements.

| Example | Shows | Capabilities used |
|---------|-------|-------------------|
| [`cogging_skew_demo.py`](cogging_skew_demo.py) | Reluctance/cogging torque `tau(theta)` of a PM rotor in a salient stator; conversion to physical N·m for a real stack; and **skew ripple cancellation** | `eggshell_torque_2d`, `MachineScaling` (axial length × symmetry), `skew_average` ↔ `skew_factor` |
| [`validation_slot_winding_spectrum.py`](validation_slot_winding_spectrum.py) | Validation-class harmonic spectrum from explicit slot sign tables, including fractional-slot layouts | `slot_table_winding_factor`, `integral_slot_winding_factor`, `mmf_harmonic_direction` |
| [`validation_pm_drive_speed_map.py`](validation_pm_drive_speed_map.py) | Validation-class PM-machine speed map across MTPA/FW/MTPV/infeasible regions plus short-circuit demagnetising trend | `field_weakening_operating_point`, `dq_operating_point`, `short_circuit_operating_point` |
| [`validation_lamination_mu_eff_sweep.py`](validation_lamination_mu_eff_sweep.py) | Validation-class laminated-steel complex-permeability sweep from static limit into deep skin effect | `laminated_mu_eff` |
| [`validation_cross_saturation_flux_map.py`](validation_cross_saturation_flux_map.py) | Validation-class d-q cross-saturation flux map with reciprocity and incremental inductance rolloff | `incremental_inductance_matrix`, `dq_flux_torque` |
| [`validation_carter_magnetizing_sweep.py`](validation_carter_magnetizing_sweep.py) | Validation-class Carter air-gap opening sweep linked to AC-machine magnetizing inductance | `carter_coefficient`, `effective_air_gap`, `slotted_air_gap_permeance_factor`, `magnetizing_inductance_per_phase` |
| [`validation_mtpa_saliency_sweep.py`](validation_mtpa_saliency_sweep.py) | Validation-class MTPA current-angle sweep across non-salient PM, IPM, SynRM, and salient-PM cases | `mtpa_operating_point`, `dq_torque`, `dq_torque_components` |
| [`validation_pm_loadline_demag_sweep.py`](validation_pm_loadline_demag_sweep.py) | Validation-class PM load-line and irreversible-demag margin sweep across gap, magnet length, knee, and shape factors | `pm_circuit_loadline_operating_point`, `demag_operating_field`, `demag_margin` |
| [`validation_cogging_skew_plan.py`](validation_cogging_skew_plan.py) | Validation-class cogging order and one-slot-pitch skew planning table across integer/fractional slot-pole layouts | `cogging_skew_plan`, `cogging_torque_order`, `machine_symmetry_sector`, `skew_factor` |
| [`validation_nonlinear_bh_circuit_sweep.py`](validation_nonlinear_bh_circuit_sweep.py) | Validation-class nonlinear B-H magnetic-circuit sweep with secant/incremental permeability and MMF split checks | `magnetic_circuit_bh_operating_summary` |

```powershell
python cogging_skew_demo.py
python validation_slot_winding_spectrum.py
python validation_pm_drive_speed_map.py
python validation_lamination_mu_eff_sweep.py
python validation_cross_saturation_flux_map.py
python validation_carter_magnetizing_sweep.py
python validation_mtpa_saliency_sweep.py
python validation_pm_loadline_demag_sweep.py
python validation_cogging_skew_plan.py
python validation_nonlinear_bh_circuit_sweep.py
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
