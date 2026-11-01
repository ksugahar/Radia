# radia-optuna execution policy

Develop and verify the shared CUDA kernels through the Python/pybind11 binding
first. Compile the common .cu source once with MATLAB-bundled nvcc and place
its object in a CMake static library. Link that same verified artifact into
pybind11 and the thin plain MSVC MEX adapter using matlab_add_mex. The mexcuda
build route is retired; do not retain it as a fallback.

TPE supports explicit CPU (OpenMP host kernels) and GPU (CUDA kernels)
selection. The MEX adapter must not include or link NGSolve TaskManager.
Implicit GPU-to-CPU fallback is prohibited. Verify CPU execution creates no
CUDA context and verify same-kernel Python/MEX parity before qualification.

Preserve Optuna 5.0.0 equations, ordered RNG consumption, first-maximum ties
and failure-state restoration on both routes. Validate both modes. Treat
execution selection as process-local configuration, separate from Study
storage; parallel workers require their own explicit selection.
