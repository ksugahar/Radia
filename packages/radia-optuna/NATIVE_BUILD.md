Shared native build
===================

`native/optuna_cuda.cu` contains the existing CPU/OpenMP and CUDA equations.
The package entry compiles it once with explicit MATLAB-bundled nvcc12.8 /
MSVC14.44 paths; CMake(CXX) archives that external object in `optuna_kernels.lib`.
The private Python
binding and ordinary MSVC MEX link that same archive. There is no executable
mexcuda branch or alternate compiler fallback.

On reviewed source `0acab7dec`, the Windows RTX 5000 build consumed identical
CUDA object and archive bytes in Python and MEX. Stored CPU/GPU oracle cases,
ten density arrays per backend in the bitwise binding comparison, and 37
existing MATLAB tests passed. Loaded module observations confirmed the Python
pip runtime and MATLAB-bundled runtime separately; CPU execution created no
CUDA context in the observed fresh processes. These results belong to that
source. Final 0.5.0 artifacts and installed-wheel acceptance require their own
source-bound evidence; A100 has not qualified the new route. Existing oracle
bytes and numerical limits remain unchanged.

The package build entry is `tools/build_native.ps1`, with explicit `-Python`,
`-MatlabRoot` and `-BuildDirectory` paths. Stage `python` builds the private
kernel module. Run `validation_test/validate_kernels.py --backend cpu` and then
`--backend gpu` with that `python-build.json` receipt; GPU execution requires
its own resource authorization. Stage `matlab` requires the resulting
`-PythonVerification` receipt and rejects changed archive/objects/source.
`buildOptunaMex(python,buildDirectory,pythonVerification)` is only a thin MATLAB
launcher for this plain MEX target. The root Optuna switch is a delegate.

`validation_test/validate_binding_parity.m` compares same-input density bits,
not a tolerance, with the Python output. Full existing MATLAB upstream tests
still guard model construction, host RNG state/consumption, candidate ordering
and first-maximum ties. A primitive-binding comparison alone does not certify
the whole optimizer. Mathematical rules and thresholds must not be relaxed to
accommodate a changed compiler, layout or lifecycle behavior.

No MATLAB/NGSolve/TaskManager type enters the common library. MATLAB GPU-array
allocation validity stays in its adapter. Python allocation owners retain a
driver context lifetime ID; pointer address equality is not reset detection.
The system driver is loaded only on explicit GPU initialization, and missing
context identity API exports fail loudly rather than using pointer attributes.
Prepared buffers keep shared workspace ownership, and shutdown refuses live
buffers before discarding callbacks. Owned-process reset checks rejected stale
allocations explicitly. MATLAB host RNG state was unchanged, and an explicit
clear/rebuild preserved density bits. This is not silent recovery: callers must
discard stale models and rebuild explicitly after reset. Physical address
reuse, every mexAtExit exception path and CUDA destructor ordering were not
certified by normal exit or these reset checks.

`radia_optuna.kernels.load(backend)` is private and requires `cpu` or `gpu`.
It never guesses another backend. Ordinary package import does not import the
native module. The optional `native` extra supplies the already approved
CUDA runtime12.8.57; nvcc fatbinary registration may map cudart even on CPU
native-module load. DLL mapping is distinct from context creation. Explicit
CPU must still prove zero CUDA contexts/launches/transfers in a fresh process.
Python does not search MATLAB DLL directories. MATLAB-only consumers do not
need the optional Python runtime extra or a Python solver process.

The production extension uses its exact CPython ABI, so wheels have interpreter
and ABI tags rather than `py3-none`. The build receipt binds the Python module
and MEX hashes to the archive/object hashes, without private host paths in the
wheel. Detailed compiler paths/commands remain private build evidence. Same
compiled inputs do not imply cross-process runtime/context sharing.

Pin provenance: pybind11==3.0.2 is the existing repository CI baseline; the
CUDA runtime pin matches reviewed root/radia-simulink metadata. CMake>=3.20 and
ninja follow existing build metadata/CI. No new pin was selected independently.
CuPy is neither a production dependency nor a fallback. CUB/pybind notices are
retained; complete SDK/dependency redistribution closure remains a publication
audit, separate from a successful build or runtime check.

Private primitive inputs are host vectors in row-major sample order; categorical
samples use the existing kernel's explicit one-based convention.
`best` returns a zero-based candidate index; the public MATLAB API is unchanged.
Invalid backend/sample/shape inputs are rejected. The Python-first stored-oracle runner
currently covers density and selection; existing MATLAB gates cover model
construction and sequential RNG. The executed lifecycle checks cover pinned
model retention, live-owner shutdown rejection, closed-model rejection and
owned reset behavior. Primitive density parity alone does not certify every
wrapper or complete upstream API coverage.
