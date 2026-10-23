# Recorded ESIM slab adoption

This record excludes nonlinear BH, full-coil models, GPU execution and independent FEM validation.

This is source-only numerical evidence using the Python ESIM source, not an
exact-wheel or native-artifact qualification. The runtime package version was
5.3.3 with the identified source overlaid. The exact original measurement entry
is preserved as entry.py; its owned execution layout contained source/ (the
recorded checkout) and venv/. The impedance entry also reads the then-existing
historical B24 fixture; the phase entry reads historical-A25.json for the
candidate historical observation. These optional observations do not define
the analytic convergence criterion. The checker never executes this entry.

result.json is a public scientific summary derived from immutable originals,
not a rewritten raw experiment. It records the original hashes. Normal and
restored runs agreed; temporary production mutations were restored and never
committed. Mutation survivors in the impedance-only record remain evidence of
a gap; the companion phase record closes that detection gap.

The regime a/delta about 105 is effectively deep-skin. It does not establish
finite-thickness or midplane behaviour. h_max is in negligible field: the rate
is interpreted with uniform refinement of the entire fixed grading, and
h_min/delta is recorded alongside it. The tests surrogate uses the same N25
production path; its coarse analytic error is deliberately not an assertion.

The v3 wrapper points to adoption-result.json, a separately identified metadata
export of the approved A ladder. Historical result.json and entry.py bytes
are unchanged. reference-verification.md supplies the algebraic specialization
of the published diffusion equation; independent scientific review remains
required. The A25 surrogate tolerance was approved in
APPROVAL-20261023-ESIM-SURROGATE-TOL-01 (rtol1e-9,atol0,negative phase check).
Schema hashes bind the serialized criterion, rationale and feature/record scope;
they do not authenticate consent or replace review of its original scope.

The canonical reference key is haus1989electromagnetic. The hash-bound
symbolic-verification.json receipt records eight exact SymPy identities and
a conjugate-profile negative control, with local derivation/script hashes.
This algebra check does not replay the saved convergence experiment.
