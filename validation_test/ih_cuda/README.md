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
