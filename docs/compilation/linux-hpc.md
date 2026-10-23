# Linux HPC source builds

The Linux entry point is `tools/build_linux.py`. It builds CPU Python extensions
with Intel LLVM (`icx`/`icpx` by default), CMake and Ninja. It never installs
dependencies, loads environment modules, selects GCC as a fallback, submits a
scheduler job or modifies a shared Python installation. Windows continues to use
`Build.ps1` and Visual Studio Build Tools 2022 / MSVC 14.44.

## Status and dependency contract

The current NGSolve and Netgen pins are **6.2.2608**, owned by
`pyproject.toml`. The Linux builder reads those pins and rejects a mismatched
environment before configuring CMake. Rebuild every extension when the pins
change; an import of a previously built extension is not build validation.

An earlier external Linux trial reported successful IH/HACApK checks using
NGSolve 6.2.2607 and both GCC and Intel LLVM. That establishes a starting point,
not certification of 6.2.2608 or of the current source. The 6.2.2608 Intel LLVM
HPC path must pass an actual build, fresh imports and affected numerical tests
on the target platform before it is described as verified. Linux wheel
publication and full-package acceptance are separate from this source build.

Avoiding the GCC compiler does not remove GNU C++ runtime dependencies.
Intel LLVM on Linux can use GNU standard-library headers and `libstdc++` to
match the NGSolve wheel's C++ ABI. A sufficiently recent C++20 standard library
may be needed even with `icpx`; do not equate that with running `g++`.
Do not switch only Radia to LLVM `libc++` while linking a `libstdc++` NGSolve
wheel. Inspect the actual dynamic dependencies rather than assuming a compiler
name proves runtime independence.

## Environment setup

Use a fresh isolated job venv and an isolated checkout. Do not reuse a Windows
build directory, an existing CMake cache from another compiler, or another
job's environment. Load the site's Intel LLVM and oneMKL modules before running
these commands. Module names and filesystem paths are site-owned settings;
keep those in private operations notes.

```bash
python3.12 -m venv "$JOB_ROOT/venv"
. "$JOB_ROOT/venv/bin/activate"
python -m pip install 'ngsolve==6.2.2608' 'netgen-mesher==6.2.2608' \
  'pybind11==3.0.2' cmake ninja numpy scipy h5py threadpoolctl pytest
# MKLROOT must point to the oneMKL installation selected for this job.
test -f "$MKLROOT/include/mkl.h"
test -f "$MKLROOT/lib/libmkl_rt.so"
export LD_LIBRARY_PATH="$MKLROOT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
python tools/build_linux.py --build-dir "$JOB_ROOT/build-icpx" --jobs 1 --require-clean-source
```

Increase `--jobs` only within the admitted allocation. The default is one job;
the builder does not infer permission to consume all visible CPUs. Run long
builds in a scheduler allocation, not on a shared login node. The toolchain and
standard-library paths selected for the job must also be visible there.

`CC`/`CXX` or `--cc`/`--cxx` select an explicit LLVM toolchain. Clang is an
explicit alternative requiring its own platform verification; there is no
automatic fallback. All four CPU modules (`_radia_pybind`, `peec_matrices`,
`sparsesolv_ngsolve`, `axifem`) build by default. `--targets` limits an iteration
build; it does not establish acceptance of the omitted modules.

Use one source checkout per active build: CMake copies the resulting modules
to `src/radia`. CMake discovers the Linux wheel's `lib/cmake` and bundled
Netgen/pybind11 headers, links MKL and uses `.so` output names. On Windows those
outputs remain `.pyd`. Always use the same Python for configuring and importing.

## Execution and evidence

The build command imports the selected modules in a fresh process with this
checkout's `src` on `PYTHONPATH`. It writes `linux-build.json` in the build
directory on success or failure, including source SHA, source-state fingerprint, interpreter,
compiler, dependency versions, CPU affinity, job count, artifact SHA-256 hashes
and `ldd` output. This report records a source build, not release approval.
Retain the console log and the report with the job's private evidence.

For a clean checkout whose identity is unchanged during the build, each selected
module also receives a `.so.build.json` sidecar using the same
`radia.native-build-provenance.v1` contract as Windows. Dirty development builds
can run without `--require-clean-source`, but receive no acceptance sidecars.
Stale sidecars for selected targets are cleared before configuration.

For a separate computation, restore the same checked runtime environment:

```bash
export PYTHONPATH="$RADIA_SRC/src"
export NGS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -c 'import radia, ngsolve; print(radia.__file__); print(ngsolve.__version__)'
python -m pytest tests/test_hacapk_architecture_policy.py tests/test_peec_hacapk_smoke.py
```

Choose and record the MKL threading layer explicitly for the workload. The
earlier Linux trial reported that TBB did not obey `MKL_NUM_THREADS` or
`mkl_set_num_threads` as expected. Use scheduler CPU binding/affinity to bound
the process; environment variables alone are not proof of thread limits.
Record the actual backend, affinity and threads in numerical result JSON.

Existing IH/HACApK validation tests retain their current tolerances and reference
values. Run the relevant independent validation after the native build; a
successful import or regression comparison alone does not certify physics.

## CAD and loader issues

An earlier Linux trial found symbol collisions between different OpenCascade
versions loaded by Netgen and OCP in one process. A 6.2.2607-specific loader
wrapper is not certified for 6.2.2608 and is not installed by this builder.
Do not bypass its version guards or assume a compiler change fixes the issue.
Use checked exported `.vol` meshes for solver acceptance. If STEP processing
requires OCP, verify the 6.2.2608 loader behavior separately or isolate that CAD
stage in another process. Keep canonical CAD, mesh and label contracts.

The Linux build covers CPU Radia only. It does not build or deploy Cubit's
Windows plugin, MATLAB MEX, Simulink or GPU artifacts. Those distributions
retain their own toolchain and acceptance contracts.

## Troubleshooting

| Symptom | Check |
|---|---|
| Missing `icpx` | Load the intended Intel compiler module; inspect `icpx --version`. |
| C++20 header error | Inspect the standard-library headers selected by `icpx -v`; compiler and standard library have separate versions. |
| `GLIBCXX_*` missing | Compare `ldd` and the runtime library paths in the allocation; avoid mixing incompatible environments. |
| Unknown pybind11 base type | Confirm bundled Netgen headers precede the standalone pybind11 include path; rebuild with a fresh cache. |
| Missing MKL/OpenBLAS shared object | Check the selected MKL runtime paths and NGSolve wheel dependencies in the actual compute allocation. |
| STEP import crash | Investigate OpenCascade symbol collisions separately from solver or compiler behavior. |
| Slow first Python import | Check shared-filesystem metadata latency; retain source identity if staging dependencies onto job-local storage. |

For toolchain details see [Clang's toolchain documentation](https://clang.llvm.org/docs/Toolchain.html).
