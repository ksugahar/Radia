# MATLAB-owned AMS runtime migration

The product constructs AMS within MATLAB/mexcuda. Python generates independent
CPU validation oracles only; GPU users do not supply a Python-built solver
snapshot or run a Python adapter. No standalone preprocessing executable is used.

## Existing implementation reused

The retired MEX gateway at commit `cf45378f0`, `src/matlab/radia_mex.cpp`,
`SparseSolvAMS`, directly constructs the shared `ngla::HypreBasedAMS` from the
real HCurl matrix, discrete gradient, vertex coordinates and free-DOF mask.
It deliberately leaves TaskManager inactive during hierarchy construction.
The new `ams_matlab_setup.cpp` uses the same constructor and exports the same
hierarchy to the existing FP64 CUDA cycle/PCG implementation. NGSolve continues
to own FE assembly/orientation and sparsecholesky coarse factorization.

## Clarified dependency boundary

The user permits `python312.dll` dependency and indirect TaskManager dependencies
through standard NGSolve headers/DLLs. TaskManager activation and parallel
execution remain forbidden. `SPARSESOLV_SERIAL_AMS_SETUP` removes internally
owned regions; host row/chunk work uses explicitly selected OpenMP threads.
Construction rejects an active external TaskManager. MATLAB controls invocation;
CUDA executes vector, sparse matrix and coarse-factor application operations.

The pinned upstream DLL import observations in
[dependency evidence](python_free_dependency_audit.json) remain valid. They no
longer imply a requirement to maintain a Python-disabled NGSolve build. Neither
Python interpreter initialization nor Python script calls are part of setup.

## Validation

`native_setup_validation.py` creates publishable manufactured HCurl problems,
computes independent Python CPU references, and compares every native hierarchy
array bitwise with the reference before testing GPU probe and solve accuracy.
The native command temporarily sets NGSolve's thread limit to one (without
activating TaskManager), then restores it. This is required both for matching
factor arithmetic and avoiding severe setup overhead from the default limit.
It saves JSON on success/failure, including serial/OpenMP native setup,
transfer/solve/download, cold and warmed solve times, memory diagnostics and
versions. It uses HDF5 input exchange, including for inputs exceeding 2 GB.

Build through `Build.ps1 -AmsCudaOnly`. Staged artifacts remain isolated from the
published runtime. Measurements on the GPU build host are diagnostics; formal
performance acceptance requires an idle GPU-equipped timing host. A100 stays off.
