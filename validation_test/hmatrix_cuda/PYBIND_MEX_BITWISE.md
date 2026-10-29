# Pending HDiv same-kernel frontend parity

Status: **PENDING**, including both record producers and the admitted Windows GPU
run. The comparison reader below exists; it is not a completed execution harness.
Linux Flow2 smoke, CPU oracle agreement and source-only tests do not satisfy this
gate. No change to CPU/GPU numerical tolerances is proposed.

Build pybind11 and mexcuda from the same exact commit for the same Windows GPU,
CUDA compiler version and architecture. Record the full commands, artifact hashes
and runtime identity. Both builds must use `--fmad=false --ftz=false
--prec-div=true --prec-sqrt=true`, without fast math. Matching flags alone do not
establish parity; any byte difference fails and must be investigated.

Feed both frontends identical FP64 factors, permutation, scale, charge-map data
and vectors, with explicit index-base conversion at the thin MATLAB boundary.
Hash the canonical zero-based SI input dataset before frontend conversion.

- `normalized-leaf-resident-boundary`: use the nontrivial scale/permutation and
  low-rank leaf fixture of `tests/matlab/test_hmatrix_normalized_leaves.m` on both
  shared resident routes. Capture the `apply` vector.
- `explicit-ng-host-hdiv-bh`: use the same two connected affine hex/order1 host
  payload, linear and analytic cubic material laws, RHS, initial state and solver
  controls through `radia.hdiv_gpu.PreparedHDivNonlinear` and
  `radia.gpu.PreparedHDivNonlinear`. Capture `linear_solution`, `cubic_solution`
  and `cubic_residual_history`. CPU independent-oracle certification is separate.

Each frontend must independently emit HDF5 with a UTF-8 `receipt_json` attribute:
`frontend` (`pybind11` or `mexcuda`), `source_sha`, `input_sha256`, `cuda_compiler`,
`architecture`, `device_uuid`, `deterministic_flags` (array), `case`, `output_names`
(array), plus command/artifact provenance. Outputs are nonempty FP64 datasets in
`/outputs`, same SI units and canonical ordering. The comparator checks receipt
agreement and exact shape/dtype/bytes, including signed zero; no tolerance.

After separately admitted producers finish, run:

```sh
python validation_test/hmatrix_cuda/check_pybind_mex_bitwise.py \
  --python-record python.h5 --mex-record mex.h5 --report parity.json
```

Record both input files, artifact/runtime identities and report hashes in private
evidence. Change the pending status only after reviewed evidence from that run.
No performance measurement or full-package acceptance is included.
