# CoilBuilder.mirror() straight-heading fix (true geometric mirror)

2026-07-10 LAB report, fixed 2026-07-11.  `CoilBuilder.mirror()` emitted
mirrored StraightSegment objects that extended in the WRONG direction,
outside the loop outline.

Root cause (src/radia/coil_builder.py, mirror()):

- Orientation rows are the local axes in world components; row 1 is the
  straight heading and rows 0/1 are the arc in-plane axes.
- The old code used `new_orient = M @ O` (which SCALES ROWS by diag(M),
  not a per-axis component mirror) and restored handedness by negating
  the WHOLE matrix when det < 0.  The global negation reversed row 1, so
  every mirrored straight walked from M @ start AWAY from M @ end
  (verified dump: mirror('xy') of (0.08,-0.10,+z)->(0.08,+0.10,+z) ran
  to (0.08,-0.30,-z)).  Arc endpoint dumps LOOKED right only because
  90-deg arcs cancel the broken cos-term; arc interiors swept the wrong
  quadrant (wire-path error 0.4 m on the report racetrack, 8.7e-2 on a
  general tilted coil, for ALL planes).
- Field consequence: the mirror('xy') racetrack partner delivered only
  ~24% of the correct central Bz (measured 1.020824e-03 T vs 4.324280e-03 T).

The fix (verified 2026-07-11, all values locked in
tests/test_coil_builder_mirror.py):

- Mirror each axis row component-wise (`new_orient = O @ M`), then negate
  ONLY row 2 (local Z) to restore det = +1.  Rows 0/1 keep their true
  mirrored directions, so straight AND arc paths map point-for-point to
  M @ p with length / radius / arc_angle / current UNCHANGED (the
  transform commutes with the arc frame rotation, so chaining stays
  consistent).  Row 2 only spans the symmetric cross-section height /
  arc axis, where the sign is immaterial.
- The old `current = -I` bookkeeping is GONE: with traversal order
  preserved, the true mirror needs the SAME current.  mirror() is now
  the true geometric mirror of the current distribution (pseudovector
  m -> -M @ m): an 'xy' pair of a z-axis loop ADDS its main field;
  'yz' / 'xz' mirrors reverse the circulation seen from +z (anti pair).
- VERIFIED: per-segment end_pos == M @ end_pos and wire-path == M @ path
  exact (0.0 error) for xy/yz/xz on the report racetrack and a tilted
  mixed coil (60/-75/120-deg arcs, tilt 15); mirror-lower central Bz ==
  upper == natively rebuilt lower to 7.4e-9 relative (ZXZ Euler beta=pi
  vs beta=0 trig round-trip noise); pair midplane Bz odd-in-x 9.8e-11,
  odd-in-y 1.4e-10 relative to B0.

Traps recorded for the next session:

- The buggy racetrack mirror was 180-deg rotation symmetric about z, so
  a PAIR-level odd-in-x field test passes even with the bug (8e-15).
  The sharp detectors are the per-segment end_pos invariant and the
  lower-only == upper-only central-Bz equality; the regression test
  asserts all of them.
- `rotate_copies()` had the same frame-transform disease (`R @ O`
  instead of `O @ R.T`): a 90-deg z-copy showed 6e-2 wire-path error vs
  `R @ p`.  The integration review fixed it and added pointwise x/y/z
  rotation coverage in `tests/test_coil_builder_rotate_copies.py`.
- scipy `Gimbal lock detected` warnings from `CoilSegment.__init__`
  (ZXZ `as_euler` on axis-aligned frames) are benign because the returned
  representative is valid.  The implementation now suppresses only that
  known warning locally instead of requiring pytest-wide filters.
- docs/clebsch_hodograph demo result JSONs (accel_pole_ends_fem,
  leaf_coupling_perturbation_3d) were generated with the BUGGY mirror
  lower coil and are stale until the demos are rerun.  Several docs
  demos also still use the bare `from coil_builder import ...` import
  that tests/conftest.py retired (pre-existing breakage of their
  validation node ids; the two mirror-consuming demos were switched to
  `from radia.coil_builder import ...` with this fix).
- validation_test/feec/test_clebsch_hodograph_research.py bands were
  authored against the buggy-mirror physics and were re-centred to the
  corrected physics in the same change (CI-mesh measurements with the
  fixed mirror): bz_body_T 0.313 (band 0.10..0.45, was ..0.30, also in
  test_clebsch_dipole_workflow_fem), end_overshoot 0.0168 (band > 0.01,
  was > 0.02, forward + curved_chamfer), leaf fringe fr = [1.084, 0.385]
  (fr[0] > 1.0, was > 1.2), open-boundary integrated dipole bbar1
  ~0.050 T*m (band 0.03..0.08, was 0.015..0.035).  The corrected lower
  coil roughly doubles pair-integrated quantities; end-overshoot
  RATIOS drop slightly because the body field rose.
