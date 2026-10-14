# PEEC-BEM precision qualification

RTX 5000 development will evaluate FP32 and mixed precision. The A100 path
retains FP64; A100 is not started, connected to or operated for this work.
Hardware selects a candidate, not numerical acceptance. FP64 remains the
reference and the current released precision until the complete physical case
passes the gates below.

## First candidate: FP32 preconditioning, FP64 outer iteration

Qualify FP32 preconditioner storage and application first, on both intended GPU
families. Keep the physical operator, solution, outer reductions and original
true-residual evaluation in FP64. Compare iteration count and whole-workflow
cost with FP64 preconditioning; half-size preconditioner values alone do not
establish a speedup or half the total memory use.

Use IC/PCG only for blocks satisfying the required Hermitian positive-definite
contract. The full complex PEEC-BEM system must use a suitable nonsymmetric
outer solver. If preconditioner application varies between iterations (adaptive
inner tolerance, precision or cycles), use a flexible outer method such as
FGMRES. Preserve sparsecholesky ownership for required direct FE subsolves.
Detect FP32 factorization breakdown, loss of definiteness, nonfinite values and
outer stagnation explicitly. Do not relax residual limits or silently change
backend/precision. Record preconditioner setup, apply time, outer iterations,
transfer-inclusive total time and actual memory separately in the same JSON.

## Candidate boundaries

- Evaluate FP32 storage/multiplication for repeated PEEC/BEM operator actions
  and well-scaled interaction blocks. Keep original FP64 operators for checks.
- Keep dot products, convergence decisions, true residuals, coherent loss
  integration and thermal energy balance in FP64 initially.
- Keep singular/near-field integration, poorly conditioned blocks and native
  sparsecholesky subsolves in FP64 until independently qualified.
- If needed, evaluate explicit mixed-precision correction with FP64 residuals
  and FP32 correction solves. Report correction counts and the total cost.
- Precision is an explicit configuration and evidence field. Do not silently
  switch solvers, freeze ESIM, relax tolerances or claim FP32 acceleration from
  compilation alone. Failure of a declared precision mode raises an error.

## Required full-model evidence

Use identical public PEEC-BEM geometry, materials, discretization, excitation
and thermal timeline for Python/pybind11 CPU, GPU FP64 and candidate precision.
Compare complex filament currents, terminal impedance, total tangential vector
fields, local/integrated heat and temperature, including nonlinear ESIM and
reaction-field coupling. Apply each physical case's existing acceptance limits;
do not replace them with the small algebra-test tolerance.

Compute the true residual from the original FP64 operator and RHS, never only
from the rounded FP32 system. It must satisfy the requested tolerance and the
project's 1e-6 ceiling. Check nonlinear convergence and reaction/energy balance
separately. Include cancellation, small drive, skin-depth/material variation and
mesh/size scaling cases; zero-drive cases use the declared absolute convention.

JSON records precision of storage, products, accumulation, residual and state;
FP64-relative differences and physical CPU errors; requested and measured
residuals; convergence/correction counts; input/source/binary hashes; versions;
and failures. Only after full-model correctness, compare transfer-inclusive
cold and warm whole-workflow time and memory on otherwise idle 100. Label 100
measurements diagnostic, separately from formal idle-development/build host acceptance. Include
FP64 validation storage in memory budgeting and retain max(2 GiB,20%) free.
For transfers above 2 GB use HDF5.

The existing synthetic coupled algebra and standalone ESIM cell results do not
qualify FP32 PEEC-BEM, Simulink physics, performance or release acceptance.
