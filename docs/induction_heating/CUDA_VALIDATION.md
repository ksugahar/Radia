# CUDA development and validation

Develop the reusable CPU and shared C++/CUDA kernels in `radia`, first exposed
through Python/pybind11. `radia-simulink` exposes them through thin MATLAB
mexcuda gateways and Simulink blocks. Follow the repository's current compiler,
runtime, explicit-backend, TaskManager and package-ownership policies.

## Physical accuracy and backend equivalence

Physical acceptance compares declared engineering quantities with an independent
converged reference. A tolerance on terminal resistance or inductance does not
imply the same tolerance at every field/current-density point. Check field
phase, conservation, boundary conditions and mesh convergence separately:
agreement of integrated quantities alone does not certify the entire field.

Backend equivalence compares identical declared inputs against the validated
Radia Python CPU reference, using separately user-approved numerical tolerances.
A physical 1% target is not automatically a CPU/GPU parity tolerance. On a shared
kernel, pybind11/MEX parity remains bitwise. Declare Simulink-only behavior
(ports, sample time, S-Function lifecycle or LTspice) and GPU-only kernels without
a CPU counterpart in the parity manifest, with reasons and separate validation.

GPU development can proceed alongside CPU validation; Python CPU accuracy
validation is not a prerequisite imposed on GPU users. Adoption must identify
the validated reference and its qualified scope. A backend name is not speedup
evidence: measure setup, transfers, solve, total time and memory, preserving
source/artifact identities. Do not claim an unexecuted backend is qualified.

## Fin observables

Retain the user-approved 1% targets on declared R/L observables and region-
integrated root/tip Joule heat. There is no blanket all-point 1% current-density
requirement. Freeze regions, frequency, excitation, current normalization,
ports and the inductive return/partial-inductance convention before judging.

For mean regional Joule power and peak drive-current phasors,
`R_region,eq = 2 P_region / |I_drive,peak|^2`; with RMS phasors the factor is one.
This is a loss contribution under the full-model excitation, not the intrinsic
two-terminal resistance of an isolated fin. A window-only sum is not a whole-
device resistance, and the full drive current need not pass through the fin.

The selected fin-L direction is an in-situ partial inductance including coupling
to the body, not an isolated-fin inductance. Freeze the self/mutual allocation
and matching reference convention before its 1% verdict; report both parts and
verify allocated contributions sum to the declared whole-model L, including
any region outside the measurement window. Regional magnetic energy alone
does not define a self-contained fin inductance. Do not double-count internal
reactive terms when comparing volume-current and surface-impedance models.

These criteria do not change the linear true-residual gate, reference
convergence requirements, baseline-regeneration or mutation procedures.
