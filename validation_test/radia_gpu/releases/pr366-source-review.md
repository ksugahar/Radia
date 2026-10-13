# PR 366 independent source review

An independent read-only Codex reviewer examined the AMS/FEM changes against
main source `88049e80944ab3c9de38d4a320d74a1b14f051bf`. Scope included GPU
curvature control and small coarse actions, factor/index validation, CUDA
stream and graph use, fixed affine HCurlSystem ownership and lifecycle,
independent numerical validation and exact-wheel release requirements.

One P2 finding was confirmed: NaN values in the independent numerical
validation could bypass tolerance comparisons. Commit `a102a8d00` adds
finite checks for solution, reference, applied operator, RHS, norms and
computed metrics. Fifteen direct regression cases cover NaN/Inf in each
input, overflowing reference norms, zero RHS and original-equation residual
failure. All passed. The independent reviewer re-read the fix and reported
no unresolved actionable findings. No additional concrete CUDA cache,
synchronization or HCurlSystem lifecycle defect was identified.

The reviewed fixed validation source SHA256 is
`3cffba3279c54098dd211fcdcaea628ed60eed5e19ee8ad3e17a97265321a9c5`;
the regression source SHA256 is
`fe94a2ad98c3d0085dc360133a33c708b8c7be5d06567414ddb7b1449a27e8dc`.
The reviewer performed source review only; GPU execution and independent
CPU numerical acceptance are recorded separately. This report does not
qualify a stable release by itself. A stable source/version and its exact
main CI wheel need fresh installed acceptance and final-source review.
