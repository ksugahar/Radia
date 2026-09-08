# HEX BDM2 far-product candidate

Baseline: `703cb6f7400dada34741a594f623276d7e54fb6d` on the ongoing
`claude/hdiv-hex-gram-psd` development line. This change is intentionally
separate from main integration and the ongoing magnet certification jobs.

The baseline records that after shared near-block caching, 7,457 of 10,216
quadrature thread-seconds in the BDM2 quadrupole build belong to the far
family (commit `453d9efc8`). The candidate replaces its scalar contractions
with two tiled, single-thread MKL DGEMMs. Points, weights, image mapping,
near/far dispatch and final Coulomb scaling are unchanged. NGSolve retains
worker scheduling; the caller's MKL local thread count is restored.
Temporary arrays scale with 32 target points, not the full point-pair matrix.
BDM1 and builds without HAVE_LAPACK retain their existing numerical path.

## Evidence

`results/hex_far_product_blas_lab_20260908.json` is a synthetic-kernel
benchmark, not a complete Gram-build or nonlinear-solve performance claim.
Six rectangular cell/face-like dimensions include a partial tile and a
coincident-point guard. Each route uses 7 batches of 200 calls, reporting the
median warmed per-call time including scratch allocation. At production-like
MSVC `/O2 /Ob2 /fp:fast /arch:AVX2 /GL /MD` settings, speedups are 3.11-8.78x;
maximum block-relative difference is 6.29e-16. Inputs include signed mode
values. The native harness also checks MKL thread-setting restoration.

The isolated Radia build passed with NGSolve 6.2.2606. Existing HEX/WEDGE
tests passed (11 cases); the extended far-block test explicitly exercises a
quadratic charge against an independent Python product-quadrature reference.
The quadratic oracle uses eight rather than five Gauss points to attain the
same 2e-12 tolerance; this does NOT alter production quadrature defaults.

## Reproduction

In an x64 MSVC developer shell, compile `benchmark_hex_far_product.cpp` as
a DLL with `/O2 /Ob2 /DNDEBUG /MD /EHsc /fp:fast /arch:AVX2 /GL /std:c++20 /LD`.
Use the selected Python installation's `Library/include` and `Library/lib`
from pip `mkl-devel`, link `mkl_rt.lib`, and put outputs under `C:/temp`.
Then run:

```powershell
python validation_test/feec/run_hex_far_product_benchmark.py --library C:/temp/hex_far_product.dll --output C:/temp/hex_far_product.json --compiler-flags "<actual compiler and flags>"
python -m pytest tests/test_hdiv_vim_hex_wedge_rt2.py -q
```

The runner loads the DLL in a separate process using Python's MKL DLL
directory, with a 120-second timeout, and records machine/runtime/source hashes.
It does not change installed packages or shared native binaries.

## Remaining Acceptance

Run the existing quadrupole BDM2 timing case before/after on one idle host
with identical inputs, threads, accuracy and native-build provenance. Retain
Gram phase times, solve iterations, fields, total wall time and peak memory
in JSON. Repeat the relevant image/symmetry and physical-spectrum gates.
Do not multiply total-solve performance by this kernel speedup, declare the
HEX/TET gap closed, or deploy the candidate before that acceptance.
