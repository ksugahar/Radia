# `_FieldEvalHMatrix` crossover benchmarks (step-2 de-risk)

Timing measured on **mdx** (idle, 38 logical cores) per the Benchmark Policy -- LAB timing is
codex-contaminated and must not back a decision.  Isolated env: the HEAD `src/radia` tree copied to
`C:\temp\radiabench\src` and put on `PYTHONPATH` over mdx's system Python (ngsolve 6.2.2604, MKL); the
shared machine-wide site-packages was NOT touched.  `MKL_NUM_THREADS=OMP_NUM_THREADS=1` so NGSolve
TaskManager (`SetNumThreads`) alone owns the parallelism swept.

`_FieldEvalHMatrix` embeds the rectangular obs x src field operator in a SQUARE SYMMETRIC HACApK matrix
`A = [[0,K],[K^T,0]]` over the combined `[obs; src]` points (the base `RadHACApKBase` is square-only).

## Verdict

**The H-matrix field eval wins ONLY when obs and src are spatially SEPARABLE into distinct clusters
(far-field-dominated maps).  For obs interleaved with the source (near / co-located) it is BOTH slower
than direct AND ~4-9% inaccurate.**

### `bench_fieldeval_hmatrix_far.py` -- far-separated (obs box shifted +5S off the src box): H-matrix WINS

| N (=N_obs=N_src) | direct rad.Fld | H-matrix (build+matvec) | speedup | err @ eps=1e-6 |
|---|---|---|---|---|
| 512  | 4.6 s   | 11.5 s | 0.40 | 1.1e-8 |
| 1000 | 12.7 s  | 16.7 s | 0.76 | 6.4e-7 |
| 1728 | 38.6 s  | 32.8 s | **1.18** | 1.5e-7 |
| 2744 | 109.7 s | 52.6 s | **2.08** | 2.1e-8 |

- **Crossover at N ~= 1500**; the win keeps growing (direct is O(N^2); the H-build is ~O(N^1.3)).  The
  matvec is ~free (6 ms) -- the whole cost is the ACA build (expensive per-entry polyhedron `B_genComp`).
- **Accuracy tracks eps** when obs/src separate: eps 1e-3 -> 1.9e-4, 1e-6 -> 1.1e-8, 1e-9 -> 2.5e-10.

### `bench_fieldeval_hmatrix_near.py` -- co-located (obs just outside each source cube): H-matrix LOSES

| N | direct | H-matrix | speedup | err @ eps=1e-6 |
|---|---|---|---|---|
| 64  | 0.06 s | 1.1 s  | 0.05 | 3e-16 (single dense block) |
| 216 | 0.66 s | 13.1 s | 0.05 | 6e-16 (single dense block) |
| 512 | 4.1 s  | 60.9 s | 0.07 | **9.3e-2** (multi-cluster) |

- **Root cause of the ~4-9% floor: the symmetric embed's checkerboard.**  When a cluster mixes obs and
  src points (co-located geometry), an admissible off-diagonal block `A[C1,C2]` is structurally
  `obs-obs = 0`, `src-src = 0`, `obs-src = kernel` -- a checkerboard that ACA cannot cross-approximate.
  Tightening eps barely helps (diag probe: eps 1e-3 -> 6.7%, 1e-6 -> 5.8%, 1e-9 -> 4.4%, 1e-12 -> 3.8%;
  it plateaus, it does NOT -> 0).  A single dense block (`leaf > 2N`) is exact (5e-16) -> the KERNEL is
  correct; only the multi-cluster ACA on co-located points fails.
- The build is also ~fully dense here (no admissible compression) AND pays clustering overhead -> 12-17x
  SLOWER than direct.

## Consequence for the two use cases

- **Path B (`reconstruct_field(backend='hmatrix')`, `prepare_cache_hmatrix`)** is SAFE + beneficial for
  FAR field maps (obs well separated from the magnetized body: stray-field maps at distance, particle
  trajectories, a coupling air region that stands off from the magnet).  For obs points hugging / inside
  the body it is slower AND silently ~4-9% wrong -- use `backend='direct'` there.  The default is
  `'direct'` (safe); `'hmatrix'` is an explicit opt-in for the far regime.
- **Step 3 (H-matrix the demag field at the body's OWN quad points, `assemble_demag_field`) is NOT
  worthwhile.**  Those quad points are co-located with the charges -> the embed hits the same checkerboard
  wall.  And `bench_hdiv_demag_field_batch.py` shows the dense C++ `_hdiv_demag_field_batch` (cheap charge
  closed-forms, ~170 ns/pair, 18x thread-scaling on 38 cores) is already sub-second to a few seconds and
  only becomes a bottleneck beyond ~45000 elements (far above current HDiv-VIM sizes); the Python packing
  (`t_pack`) is comparable to the C++ batch anyway.  A genuine rectangular H-matrix (separate obs/charge
  cluster trees) or FMM would be needed for the co-located case -- both out of current scope/policy.
