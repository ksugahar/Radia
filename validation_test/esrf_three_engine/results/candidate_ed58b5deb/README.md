# Official wheel scalar PCG acceptance

The official CI wheel passed the strong-field native smoke and focused scalar
PCG regressions on LAB and hibino. This closes the reproduced straight WEDGE
order-2 inner-solve failure without changing tolerances, iteration limits,
default preconditioners, or global MKL dispatch settings.

## Identity

- Source: `ed58b5debb96efd3d4dc5df95ac6a5cd4056a987`
- CI: https://github.com/ksugahar/Radia/actions/runs/34942185100
- Distribution: Radia 5.0.0, CPython 3.12 Windows AMD64
- Wheel SHA-256: `819389941b4e48f11d91a23ccc5b4d31378c027290b7ff484a9562d4b27562af`
- Native SHA-256: `9174b12f02e0f156de72eb6bb659f85d53c604bd58c0f1a121c1a0c6341b76d3`

Each host used an isolated installed-wheel environment. The runner checked
the native and all 303 Python files against the wheel; `pip check` passed.
Commands, import paths, dependency inventories and exit codes are retained.

## Results

| Check | LAB | hibino |
| --- | --- | --- |
| Strong-field TET/HEX/WEDGE, orders 1 and 2 | 6 passed | 6 passed |
| Explicit frozen-PCG and diagonal/nonfinite tests | 7 passed | 7 passed |
| WEDGE order-2 nonlinear relative residual (target 2e-5) | 1.288040e-8 | 1.700737e-5 |
| Frozen-PCG iterations | 1881 | 1862 |
| Independently evaluated frozen-PCG residual (target 2.744577167e-4) | 1.563726960e-4 | 1.845639614e-4 |

The seven focused tests include four NaN/Inf cases and an exhausted-iteration
control. The frozen test rebuilds the demagnetization operator and checks it
against the saved reference within a stated tolerance; it is not a claim of
bit-identical reconstructed matrices. See the preceding candidate_d0d0bc4b5
evidence for the exact saved dense-system restart experiment.

The fix removes unconditional periodic restarts from scalar PCG only. True
residual verification, necessary corrective restarts, the iteration cap and
nonpositive/nonfinite curvature rejection remain. Batched PCG is unchanged.
Old-wheel host differences and tiny old/new operator differences are not fully
attributed; these results do not establish a single cause for every difference.

## Scope and recovery

### Cubit package combination

The LAB isolated environment also passed a combination check with official
CME 1.0.1 wheel SHA-256
`859211e96343cd1e97b202f22ea5d008997406de35459bd200130922e286bf73`.
All 344 installed Radia package files and 11 CME package files matched their
respective wheels. The standalone Cubit sphere export (curved order 2) passed
the structural/label checker and was consumed by HDiv `Solve(order=1,
mu_r=100, H_ext=(0,0,1000), gram_eps=1e-7, tol=1e-8)` with finite magnetization.
See `lab/cme-radia-combination.json`. This is interoperability evidence, not
a separate accuracy or nonlinear qualification.

CME 1.0.1 deliberately reuses the existing manifest-pinned native payloads;
no newly rebuilt payload is being published. Source-and-payload provenance
verification passed against the installed official CME wheel and the current
checkout. Private local rebuild outputs were retained for audit under
`S:/Radia/validation_artifacts/hdiv_cleanup_20260915/cme-local-rebuild`.
The release checkout was restored to its tracked CCM and the same hash-checked
curver release asset used by CI. QUAD preflight passed the version and native
provenance gates; its unsynchronized-main warning remains a final release gate.

These are strong-field straight-element checks, not a blanket qualification of
curved elements, IMA, pyramid elements or every application. Release publication
and deployment remain separate gates; no tag is authorized by this record alone.

The completed hibino job is recovered in
`S:/Radia/validation_artifacts/hdiv_cleanup_20260915/hdiv-pcg-ed58-20260915-recovery.tar.gz`.
Its SHA-256 is
`2be9b3518c3c9bc6a9d566ebf90c5310000ac885bef06cecf5752dbd54ea9a65`.
All 22 archived files were checked against `ed58-recovery-manifest.json` in the
same directory. The archive preserves inputs, scripts, wheel and results; the
disposable virtual environment is reproducible and is not archived.

After evidence commit `875ab988d`, the completed hibino job directory
`C:/temp/hdiv-pcg-ed58-20260915` and its sibling recovery archive were deleted
at 2026-09-15T07:53:43Z. Path, process, link and archive-hash checks passed;
both targets were verified absent. See `hibino/cleanup.json`. Unrelated jobs
and the installed global Radia environment were not modified.
