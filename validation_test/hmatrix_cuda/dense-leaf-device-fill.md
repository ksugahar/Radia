# Deferred dense leaf boundary: staged, GPU not connected

The current CUDA assembly evaluates dense leaf entries on the device, copies
them to native C `a1`, packs them into a CPU `Plan`, then uploads them again.
This note defines the ownership boundary for removing those round trips.
It is an implementation investigation, not a completed GPU backend.

## C boundary and current acceptance

`RAD_HACAPK_DENSE_LEAF_SINK` is **undefined** in ordinary CPU and current
mexcuda builds. When explicitly defined, `fill_dense_leaf_exact` offers the
already selected dense leaf to a native consumer before allocating `a1`.
Zero deliberately retains the original fill. Nonzero acknowledges ownership
of materialization; the C leaf retains its dimensions, permutation ranges,
kind 2 and rank zero, with a null `a1`.

This occurs after the existing ACA decision, including rank-cap dense
promotion. It changes neither pivot selection nor the update arithmetic.
The wrapper's post-fill work inspects leaf kinds/ranks, not dense values.
Existing `fromHacapk` rejects a null dense payload. No supported command
publishes or applies such an incomplete matrix.

The optional `check_full_host.py --deferred-dense` C test covers BBOX/PCA,
rank caps 32/1, full/symmetric fills, nonidentity permutation, exact metadata
and low-rank factor preservation, query-count conservation and destruction
before materialization. Its independent synthetic direct matrix checks the
fully materialized result. The test materializes on CPU solely to exercise
the C contract; it does **not** certify GPU filling or performance.

## Required device consumer

1. Use an owner-scoped, module-local thread binding, without a process-global
   setter. Copy or own all borrowed metadata; no pointers cross MATLAB.
   Capture every callback error without throwing through C, and rethrow it
   before packing or registry publication.
2. Add an explicit deferred segment description to the checked native plan.
   Record logical coefficient offsets and lengths separately from uploaded
   low-rank data. Merely replacing dense values with zeros and uploading the
   full array would retain the transfer being removed.
3. Fill row-major dense segments using the actual cluster permutation and
   physical scaling. Keep external low-rank factors column-major, the ordered
   symmetric entry average, original normalization multiplication order and
   the existing touching/non-touching integral rules.
4. Ensure oracle/cache lifetime extends through filling. Current assembly
   destroys the oracle before constructing `Resident`; the new path must
   account for **simultaneous** cache, resident coefficients, maps and
   workspace. Existing maximum-of-stage counters would undercount that peak.
   Precheck all dimensions, products, total bytes and the GPU reserve.
5. Keep the context unpublishable until every deferred segment is filled,
   device error flags checked and completion synchronized. Exceptions must
   release the C tree, plan, device arrays and cache. No CPU or full-Gram
   fallback is allowed.

## Remaining gates

- Whole native CPU build/default scalar and batched routes after shared edits.
- Checked deferred layout, scale/permutation, failure transaction and memory
  accounting tests; the current C boundary test is only the first layer.
- RTX operator, eviction, invalid input and handle lifecycle gates.
- Independent CPU actions, whole solution and original true residual for
  affine/curved cases and rank-cap promotion.
- Same-artifact old/new transfer-inclusive construction/solve timings and
  exact owned versus sampled memory. A lower query count alone is not a
  speedup, and earlier diagnostic samples cannot be relabelled.

The dense mass factor and full-hierarchy capacity still block a general
100k-DoF claim. This boundary does not change either limit.
