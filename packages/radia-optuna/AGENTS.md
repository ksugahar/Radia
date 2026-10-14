# radia-optuna execution policy

Use MATLAB mexcuda for the single native module. TPE supports explicit CPU
(OpenMP host kernels) and GPU (CUDA kernels) selection, following the user
policy of 2026-10-13. This Optuna-specific CPU route does not permit ordinary
mex builds, NGSolve TaskManager linkage or implicit GPU-to-CPU fallback.

Preserve Optuna 5.0.0 equations, ordered RNG consumption, first-maximum ties
and failure-state restoration on both routes. Validate both modes. Treat
execution selection as process-local configuration, separate from Study
storage; parallel workers require their own explicit selection.
