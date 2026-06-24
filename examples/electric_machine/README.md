# Electric-machine examples

Open, analytic-gated rotating-machine FEA with **radia-ngsolve** — the 2D→3D motor
post-processing distilled from the open ONELAB/GetDP electric-machine study, built on
NGSolve's high-order curved elements.

| Example | Shows | Capabilities used |
|---------|-------|-------------------|
| [`cogging_skew_demo.py`](cogging_skew_demo.py) | Reluctance/cogging torque `tau(theta)` of a PM rotor in a salient stator; conversion to physical N·m for a real stack; and **skew ripple cancellation** | `eggshell_torque_2d`, `MachineScaling` (axial length × symmetry), `skew_average` ↔ `skew_factor` |
| [`validation_slot_winding_spectrum.py`](validation_slot_winding_spectrum.py) | Validation-class harmonic spectrum from explicit slot sign tables, including fractional-slot layouts | `slot_table_winding_factor`, `integral_slot_winding_factor`, `mmf_harmonic_direction` |

```bash
python cogging_skew_demo.py
python validation_slot_winding_spectrum.py
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
- Winding / leakage / magnetising analytics: `radia_mcp.radia_ngsolve.solve`
  (`winding_factor`, `skew_factor`, `slot_leakage_inductance`,
  `magnetizing_inductance_per_phase`, `effective_air_gap` = Carter).
- The open GetDP reference these mirror: `radia_mcp.motor.onelab_knowledge`
  (`magstadyn_source`, `twod_corrections`).
