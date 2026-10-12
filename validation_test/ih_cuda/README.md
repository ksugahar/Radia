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

## Native nonlinear ESIM cells

`run_esim_cell_smoke.py` differentially checks the unmodified Python
`ESIMFiniteSlabSolver` against device finite differences for slab and cylinder
cells, constant real/complex permeability and monotone B(H) piecewise cubics.
Each CUDA thread solves one cell with the original boundary stencils, Picard
relaxation, convergence criterion and midpoint active/reactive power integration.
This is finite-difference cell arithmetic, not a replacement for NGSolve FE
semantics. Nonconvergence and linear true residual above 1e-6 fail loudly.

The validation supplies interpolation coefficients from the independent CPU
oracle. A MATLAB/native B(H) coefficient constructor, full material-state and
reaction-field coupling, checked persistent MATLAB handles, pybind11 model
parity and Simulink acceptance are still required. This validation is not a
Python prerequisite for the eventual GPU product, and its gateway is not shipped.
The official MATLAB MCP currently reports `failed to attach to MATLAB session`;
no production SLX acceptance has been performed. Engine-owned sessions were used
only for these standalone CUDA development checks.

## Current checkpoint: FP32 preconditioning and persistent Simulink state

This checkpoint supersedes the earlier pending implementation notes above.
The checked persistent mexcuda context now owns electromagnetic, loss and
thermal device state across MATLAB calls. Fixed right preconditioning can store
and apply its inverse-action matrix in FP32 while the original operator, outer
BiCGStab vectors/reductions and final true residual remain FP64. The supplied
inverse-action matrix is not yet a native IC/AMG builder.

`results/fp32-preconditioner-20261012.json` records nine algebra cases (no,
FP64 and FP32 preconditioning, three drives, ten thermal updates), independent
NumPy comparisons, true residuals, iteration counts and owned buffer sizes.
It does not contain timing acceptance or sampled peak memory.

Official MATLAB MCP attachment was restored in an operation-owned session.
The official Simulink Toolkit checked the disposable resident-state model;
save/close/reopen and repeated simulation passed with zero live contexts after
simulation. `results/fp32-simulink-algebra-20261012.json` records the FP32
preconditioned temperature and heat trajectories against independent NumPy
backward Euler algebra. This is synthetic orchestration validation, not physical
PEEC-BEM/ESIM parity. Actual native geometry assembly, nonlinear ESIM/reaction
coupling and independent pybind11 physical acceptance still block release.
A100 was not used. No speedup is claimed from these checks.
