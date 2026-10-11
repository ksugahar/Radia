# Python-free AMS runtime investigation

The requested product boundary is MATLAB/mexcuda execution without a Python
snapshot adapter. Python CPU reference calculations remain useful development
and differential-validation tools; they do not establish a Python-free runtime.

## Observed dependency

The current snapshot module creates the AMS hierarchy with NGSolve and exports
its sparsecholesky factors. The staged CUDA gateway consumes that hierarchy;
it does not yet construct it directly from the FEM matrix and discrete gradient.

Inspection of the pinned Windows NGSolve 6.2.2607 wheel's PE import tables found
that both `libngsolve.dll` and `ngcore.dll` directly import `python312.dll`.
[Dependency evidence](python_free_dependency_audit.json) records library hashes.
Replacing the Python wrapper with a native executable linked against these
same libraries would therefore still require the Python runtime. No execution
or performance certification of a Python-free implementation has been made.

## Constraints and candidate boundary

Preserve the NGSolve sparsecholesky coarse factorization. MATLAB `chol`, PARDISO,
another direct backend, and implicit CPU solver fallback are not substitutes.
The mexcuda binary must neither include nor link NGSolve TaskManager.

A candidate is an independent native preprocessing executable built against
a Python-disabled NGSolve/Netgen configuration. It would construct the shared
AMS hierarchy and export the same checked numeric factors. MATLAB would invoke
it at initialization/update boundaries and pass the result to mexcuda. This
keeps NGSolve dependencies out of the mexcuda process; it introduces a native
helper rather than eliminating every adapter. NGSolve remains responsible for
FE spaces, orientation and assembly; inputs at the solver boundary are checked
matrices, discrete gradient, coordinates and free-DOF mask.

Upstream [Netgen CMake](https://github.com/NGSolve/netgen/blob/master/CMakeLists.txt)
provides `USE_PYTHON`; [NGSolve CMake](https://github.com/NGSolve/ngsolve/blob/master/CMakeLists.txt)
conditions Python linking on `NETGEN_USE_PYTHON`. These are source-level
configuration observations, not proof that the required standalone build and
distribution already work for the pinned release.

## Required evidence before adoption

1. Build the pinned upstream revision without Python and GUI support, using
   the project build entry point on the build host; no PARDISO route.
2. Inspect the complete native dependency closure and demonstrate operation
   where Python is unavailable. Record exact sources and artifact hashes.
3. Preserve sparsecholesky numerical contracts, singular coarse-space handling,
   factor ordering and input validation; verify CPU reference parity.
4. Define checked/versioned initialization and update artifacts. Use HDF5 for
   exchanges exceeding 2 GB. Measure helper startup and file exchange costs.
5. Validate MATLAB-only initialization, repeated GPU solves, failures, memory
   admission and cleanup. Save before/after transfer-inclusive timing JSON.
6. Resolve native-library packaging/license ownership and the existing Python
   CPU execution boundary before treating the helper as a supported product.

The Python/OpenMP snapshot experiment remains development evidence, not the
requested product architecture. Its larger-size timing run was not started
after the request to remove the Python adapter. A100 remains off.
