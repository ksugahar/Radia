# MATLAB Gateway Setup Overhead

`callMex` checks the native runtime on every call, but does not configure
Simulink file-generation folders. Run `radia.setup()` explicitly at application
or model initialization to configure those folders. Existing production MEX
S-Functions already do this during setup.

The recovered WIP's `persistent configured` flag is intentionally not adopted:
it would bypass setup after a MATLAB path change or loss of the MEX gateway.
`radia.setup` already caches expensive Python runtime discovery. Its existing
PythonExecutable, Force and MEX-availability checks remain in effect. MATLAB
`pyenv` is not the selector for this standalone MEX runtime discovery; do not
claim Python-bridge lifecycle coverage from these gateway tests.

## Reproduce

Run from a checkout with a built MEX, using Python with MATLAB Engine installed:

```powershell
python validation_test/matlab_runtime/benchmark_setup.py --tests --output C:/temp/matlab-setup.json
```

An explicit `--mex-dir` can select an existing artifact for wrapper-only
comparison. The runner checks the resolved setup/MEX paths and records MEX and
source hashes. It owns and closes a private Engine; it never reconnects to a
user session or changes an editable install. This Windows harness is light
enough for LAB; it is not a solver-scale benchmark.

## Recorded Comparison

`results/lab_before.json` and `results/lab_after.json` use the same native
artifact and seven batches of 100 calls, in separate private MATLAB Engines.
The baseline temporarily restores only the original `callMex.m` in the isolated
worktree. Source hashes identify each measured snapshot; the base commit alone
does not describe the uncommitted changes. Only the after record includes the
four runtime regression tests. No numerical kernel was changed or rebuilt for
this wrapper-overhead measurement.

Median `api.commands` gateway time was 2.343 ms before and 0.819 ms after
(about 65% lower). Direct MEX time was 0.046/0.047 ms. Cold Engine/setup times
and all raw samples are retained separately. This does not demonstrate a 65%
reduction in solver runtime, and these timings are not a universal CI limit.

The runtime tests live in `tests/matlab/test_mex_runtime_setup.m`. The existing
mdx SparseSolv MEX lane runs them alongside numerical AMS parity, including on
changes to the shared setup/gateway files. Benchmark timing is not run in
normal CI. Original shared WIP remains preserved, not overwritten or deleted.
