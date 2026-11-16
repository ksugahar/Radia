# Fin PEEC against the adaptive P2 reference

This record covers the axisymmetric volume PEEC
(`radia.peec_axisym.solve_axisym_ring_peec`) on the beak-fin section, driven
axisymmetrically. It does not cover the 3D tetrahedral PEEC, a
non-axisymmetric port, temperature dependence or any other geometry.

## Reference

The reference comes from the verdict script
`validation_test/induction_heating/axisym_volume_peec_fin.py`. That script runs
four adaptive sequences of axisymmetric A_phi finite elements:

- orders 1 and 2, at exterior factors 1.5 and 3;
- ZZ indicator, Doerfler theta 0.5.

On 2026-11-16 the user chose the final order-2, exterior-3 level (194275 DOFs)
as the reference. All four sequences converged inside their budget. They agree
within one tenth of every 5 % limit. The largest spread is 0.31 % in R, 0.03 % in
L, 0.31 % in root heat and 0.20 % in tip heat.

## Series and criterion

Each quantity (R, partial L, root heat, tip heat) forms one series over the
boundary edge sizes 0.06, 0.04 and 0.028 mm. The interior size (0.5 mm) and
grading (0.25) stay fixed. The criterion has four parts:

- errors decrease strictly;
- the log-ratio rate is at least 1, a lower bound the user approved on 2026-11-16;
- the finest error is at most 0.05;
- every true residual is at most 1e-6.

The observed rates are about 1.4 to 2.9. The finest errors are 0.073 % in R,
0.0025 % in L, 0.083 % in root heat and 0.21 % in tip heat.

## Files

- `entry.py` derives `adoption-result.json` from the saved verdict JSON and the
  hashes of the three meridian meshes. It solves nothing, and the checker never
  executes it.
- The original verdict and the four reference-sequence JSON files are kept in
  private storage. `adoption-result.json` records the verdict's hash.
- `mesh_sha256` is the sha256 of each mesh in the encoding of the tests fixture.
  The 0.06 mm mesh is `tests/fixtures/ih_fin/beak_section_r200_h060_mesh.json`.

## Surrogate test

`tests/test_ih_fin_peec_axisym_surrogate.py` locks the output of the 0.06 mm
level at the recorded source, with rtol 1e-6 and atol 0. The user approved that
tolerance on 2026-11-16 (APPROVAL-20261116-FIN-PEEC-SURROGATE-TOL-01). Two
temporary production mutations were tried at two sizes and then restored:

- scaling the ring inductance by (1 + eps);
- scaling the conductivity by (1 + eps).

With eps = 1e-5 the test fails on both (R moves 4.8e-6 and 5.2e-6, against
an allowed 1e-6). With eps = 1e-6 both survive (R moves about 5e-7); this gap
follows from the tolerance and stays disclosed. The restored source passes.

## Review

The threshold wiring is verified by a reviewer
(`APPROVAL-20261116-FIN-PEEC-P2-01-INDEPENDENT-REVIEW`). The hashes bind the
serialized criterion, rationale and scope. They do not authenticate consent.
