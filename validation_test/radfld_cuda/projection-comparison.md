# RadFld acceleration at the FEM projection boundary

Status: investigation contract, no backend selected or measured speedup.

The user-identified bottleneck is CoefficientFunction evaluation during
GridFunction.Set in FEM-MMM coupling, not an isolated arbitrary point cloud.
NGSolve must retain FE space, mapped quadrature, orientation, Piola maps and
projection semantics. Do not replace Set with nodal sampling.

## Source observations

- src/core/rad_ngsolve_field_coefficients.h: HDivFieldCoefficient evaluates
  individual mapped points or one mapped integration rule via EvaluateSerial.
- src/core/rad_ngsolve_radia_field.h: RadiaFieldCoefficient batches each mapped
  integration rule, offers explicit precomputation/cache machinery and calls
  RadFldBatchSerial for vector fields. Rule batching is not mesh-wide batching.
- These are distinct evaluators; identify the actual production CF and field
  (H, B, A or scalar potential) before choosing a kernel or reuse strategy.

## Comparison contract

First profile the unchanged projection: total Set time, CF evaluation time,
scalar/rule calls, target counts and batch histogram, repeated coordinates,
source updates and cache lifetime. Record mesh/FES/order/region/quadrature,
thread count, units, source family and input/source/binary hashes.

Compare direct, existing tree, candidate KIFMM and HACApK using the same FE
projection and source history. Unsupported backends remain explicitly pending.
Never label generic ACA or a treecode as HACApK or KIFMM. KIFMM is an evaluation
candidate here; core compact interactions retain the existing HACApK policy.

Gate both vector field values at held-out physical points and the projected
GridFunction against the unchanged FP64 direct evaluator. Relative vector and
integrated FE L2 error must each be <=1e-4 on the declared domain. Record
maximum pointwise errors with an explicit dimensional near-zero floor,
absolute errors and relevant integrated-field quantities. Keep original
same-algorithm CPU/CUDA parity gates; a 1e-4 approximation allowance does not
relax implementation parity. Projection/discretization error against physical
truth is separate from backend error against the same direct projection.

Measure cold construction and complete first projection, then repeated
projection with fixed geometry and changed source strengths. Include batching,
compression/equivalent-source setup, upload, evaluation, download and NGSolve
projection work in continuous end-to-end time. Report component times only as
diagnostics, together with reuse count and rebuild conditions. Do not launch
one transfer/kernel/synchronization for every integration point by design.

JSON must contain acceptance/pending reasons, accuracy metrics, all timing
samples, GPU owned bytes, sampled device/host memory and peak limitations,
precision, algorithm parameters, hashes and runtime versions. Use warmup and
ABBA repetitions only under idle conditions. 100 timings are diagnostic;
formal timing acceptance remains idle-mdx-only. Reserve max(2 GiB,20%) of RTX
5000 device capacity, HDF5 for >2 GB interchange. A100 is not operated.

## Primary method references

- Ying, Biros, Zorin (2004), kernel-independent adaptive FMM:
  https://cs.nyu.edu/~dzorin/papers/ying2004kia.pdf
- GPU H-matrix/ACA implementation study (2018):
  https://arxiv.org/abs/1806.11558

The papers motivate candidates, not evidence of speed or accuracy for this
RadFld projection. Backend integration and measurements are still pending.
