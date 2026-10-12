# IH resident CUDA development validation

`run_resident_smoke.py` builds a validation-only mexcuda gateway and checks the
resident complex CSR action, scalar-drive RHS and coherent quadratic loss
against independent NumPy complex algebra. Ten iterations retain operators and
state on device, with no intermediate vector transfers or explicit stage
synchronization. The loss includes basis cross terms. This is not a physical
PEEC-BEM(ESIM) solve, pybind11 physics parity, thermal time integration or
Simulink acceptance. Those remain mandatory before the radia-gpu release.

Small builds and correctness checks may run on 100 under low load with a
below-normal-priority owned MATLAB session and one compiler invocation, without
replacing installed runtimes. Do not compete with active heavy CI computation.
Formal transfer-inclusive time and memory comparisons require an otherwise idle
100 and are diagnostic provenance, not formal idle-mdx performance acceptance.
Results are always written as JSON. A100 is never started or connected to.

Next connect checked persistent context ownership, complete complex PEEC-BEM
solve/residual checks, thermal state updates and shared Simulink block state.
Keep electromagnetic/thermal vectors resident between updates; transfer only
required controls and observed outputs. Preserve block output/update ordering,
SI units, ESIM physical contracts and independently validated CPU results.

Device-scalar BiCGStab stages now compile and pass the identity-system first-
correction convergence check, including bypass of the zero omega denominator.
Full cuBLAS device reductions, the persistent solve loop, true-residual gate,
preconditioning, thermal evolution and Simulink context are still pending.
The stage check is not iterative-solver or physical IH acceptance.

The next compiled validation connects complete complex BiCGStab (cuBLAS device
reductions, eight-iteration CUDA graph batches, final true residual) to coherent
loss, device heat-load projection and a backward-Euler temperature increment
solve. The 3-unknown synthetic system agrees with NumPy direct algebra and
passes true-residual gates for zero, positive and negative excitation. This
checks numerical plumbing, not a PEEC-BEM physical model. The context is currently
retained only during the validation call; checked persistent MATLAB ownership,
Simulink output/update lifecycle and actual CPU/pybind11 physics comparison are
still required. No timing/performance claim follows from this result.

`persistent-state-20261012.json` checks ten coupled updates with one uploaded
operator set, device-only electromagnetic/loss/thermal stage exchange, accepted
state preservation on invalid controls, reset and rejection of a missing dynamic
material implementation. The device state is currently scoped to a single
validation MEX invocation. Checked persistent MATLAB handles, actual ESIM cell
and reaction-field updates, independent physical CPU parity and Simulink
acceptance remain pending. No full IH or release acceptance is claimed.
