# ESIM spatial example

Open [the executed notebook](esim_spatial_demo.ipynb) to inspect the coil,
workpiece mesh, tangential field, surface resistance and two different loss
distributions at 10, 50 and 100 kHz.

These are retained results from 2026-09-15. Recovery preserves the calculated
values and saved WebGUI states; it does not claim a new calculation or rendered
visual acceptance. The input meshes and result files are in
[`validation_test/esim_spatial`](../../validation_test/esim_spatial).

The local diagnostic `0.5 Re(Zs) |Ht|²` and the power-normalized transfer
artifact are different quantities. Their integral agreement with a prescribed
power is a conservation check, not independent validation of the spatial loss.
The example does not certify a temperature solve, strong coupling or a
converged volume eddy-current solution.

Run the notebook from this directory. Its default inputs are the checked,
bundled validation fixtures; `RADIA_ESIM_INPUTS` can select an alternative
directory only when all recorded input hashes match.
