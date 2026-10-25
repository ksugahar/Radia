---
name: build
description: Build Radia Python or standalone MATLAB MEX artifacts from the selected checkout using current Build.ps1, without changing running editable installations.
---

# Build Radia

## Select the source and environment

- The development/build host owns candidate/release builds and isolated heavy
  CPU numerical validation. Transfer reviewed source to an isolated checkout;
  recover artifacts/provenance to durable storage. Optional GPU evaluation uses
  hash-verified artifacts on the on-demand GPU worker; never replace user runtimes.
- Read the selected checkout's AGENTS.md, Build.ps1 parameters, and relevant
  CMake targets. Do not redirect a clean worktree build into the shared WIP tree.
- Use PowerShell 7 (`pwsh`), not Windows PowerShell 5. Scratch belongs in
  `C:\temp`; build outputs remain in the chosen checkout's standard build dirs.
- For release verification, use a dedicated Python venv. Put its Scripts
  directory first on PATH for this invocation, and verify `python` resolves
  there. Install exact NGSolve/Netgen pins from that checkout's pyproject.toml,
  plus its current build dependencies. Obtain MKL headers/libraries from
  `mkl-devel`; do not bundle or replace the user's NumPy.
- Inspect CMakeCache.txt before reusing a build directory: source path,
  Python, NGSolve/Netgen and MKL must match this checkout and environment.
  Do not copy a cache or native binary from another branch to make it pass.
- Coordinate heavy compilation with active performance measurements. mdx1/mdx2
  own CI; do not create a LAB runner or stop another task to free resources.

## Choose the narrow target

Run from the selected repository root:

```powershell
pwsh -NoProfile -File ./Build.ps1 -RadiaOnly
pwsh -NoProfile -File ./Build.ps1 -CudaMexOnly
pwsh -NoProfile -File ./Build.ps1 -OptunaMexOnly
pwsh -NoProfile -File ./Build.ps1 -AxiFemOnly
```

Radia CPU execution uses Python/pybind11; GPU execution uses MATLAB mexcuda.
The independent Optuna MEX target is not a Radia CPU solver build.
These switches are mutually exclusive. Use the one needed by the change.
Omitting them invokes the full build. Read current implementation before using
`-Rebuild` or `-Test`; neither is the default for a focused maintenance fix.
Use Build.ps1's current propagation paths, not historical .ccl/plugin recipes.

## Protect development installs

- Do not pass `-InstallToSitePackages` during isolated verification. Build
  success is not authorization to repoint 100号機 editable sources or replace
  LAB verified wheels.
- Never bulk-kill Python, MCP, MATLAB or Cubit to unlock artifacts. Build into
  an isolated checkout and defer replacement of an in-use file.
- Check import origin, version and native ABI in a fresh process from the
  selected source. A plain global `import radia` may load another worktree.
- Run focused regression tests for the changed boundary. Native numerical and
  performance evidence belongs in validation_test with result JSON; do not
  launch the entire suite merely because compilation finished.
- Record source SHA, interpreter/runtime versions, target and artifact hashes.
  Follow release-quad for deployment; LAB uses verified wheels and only
  100号機 uses an editable development install.
  MCP client reconnection is separate from file installation.

## Other distributions

- Standalone MATLAB Optuna uses `-OptunaMexOnly` and its independent mexcuda workflow;
  it must not acquire Radia/NGSolve/MKL dependencies. Use the package mexcuda recipe for GPU numerical kernels and CPU host control.
  Use bounded OpenMP for suitable parallel host work and retain serial
  state management and ordered random consumption.
- Eqnedit64 is maintained and released in its independent repository.
- For sparsesolv inspect the current in-repository CMake integration. Do not
  use the retired standalone repository or copy DLLs into global site-packages.
