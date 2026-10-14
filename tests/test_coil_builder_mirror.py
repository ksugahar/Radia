#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regression test: CoilBuilder.mirror() geometry and field symmetry.

Locks the 2026-07-10 fix: mirror()'s handedness-restoration step used to
negate the whole mirrored frame (new_orient = -(M @ O)), which reversed
the local heading (orientation row 1) of StraightSegment objects while
keeping new_start = M @ start_pos.  Every mirrored straight then extended
backwards out of the loop outline (verified segment dump: the mirror('xy')
of a straight (0.08, -0.10, +z) -> (0.08, +0.10, +z) ran to
(0.08, -0.30, -z)), and mirrored arc interiors swept the wrong quadrant.
Field consequence: the mirror('xy') partner of a racetrack loop delivered
only ~24% of the correct central Bz.

Invariants locked here, for each mirror plane M in (xy, yz, xz):
  1. mirrored.segments[k].end_pos == M @ original.segments[k].end_pos
     (and the same for start_pos), segment order preserved;
  2. the mirrored segment chain is continuous (start[k+1] == end[k]);
  3. every mirrored orientation is a proper rotation (det = +1);
  4. the full wire path (to_wire_segments, arc interiors included) maps
     point-for-point to M @ p;
  5. mirror() preserves the current value (true geometric mirror -- the
     pseudovector behavior m -> -M @ m comes from the geometry alone);
  6. for an x-symmetric racetrack pair upper + upper.mirror('xy'), the
     mirrored lower coil alone reproduces the upper coil's central Bz,
     and the pair's midplane Bz has odd-in-x (and odd-in-y) parts
     < 1e-8 relative to the central field B0.
"""
import os
import sys

import numpy as np
import pytest

# Canonical package import (tests/conftest.py policy: from radia.X import Y);
# the src insert only serves standalone `python tests/test_...py` runs.
_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
	sys.path.insert(0, _SRC)

from radia import coil_builder as cb

MIRROR_MATRICES = {
	'xz': np.diag([1.0, -1.0, 1.0]),
	'yz': np.diag([-1.0, 1.0, 1.0]),
	'xy': np.diag([1.0, 1.0, -1.0]),
}

GEOM_TOL = 1e-12       # exact +-1 sign flips; allow only rounding noise
# Central-Bz equality between the mirrored and the natively built lower
# coil: the ZXZ Euler placement of mirrored frames takes a different trig
# round-trip (beta=pi instead of beta=0), giving ~7e-9 relative rounding
# noise; the locked bug produced a 7.6e-1 relative error.
FIELD_REL_TOL = 1e-7
ODD_PART_REL_TOL = 1e-8


def build_racetrack(z0, current=1000.0):
	"""Rounded-rectangle loop at z=z0 (the 2026-07-10 bug-report geometry).

	Side straights at x=+-0.08 (y in [-0.10, 0.10]), top/bottom straights
	at y=+-0.12 (x in [-0.06, 0.06]), 90-deg corner arcs r=0.02.
	"""
	r = 0.02
	return (cb.CoilBuilder(current=current)
		.set_start([0.08, -0.10, z0])
		.set_cross_section(0.004, 0.004)
		.add_straight(0.20)
		.add_arc(r, 90)
		.add_straight(0.12)
		.add_arc(r, 90)
		.add_straight(0.20)
		.add_arc(r, 90)
		.add_straight(0.12)
		.add_arc(r, 90))


def build_general_open_coil(current=500.0):
	"""Non-axis-aligned mixed open coil for the general-orientation check.

	Tilted start frame (30 deg about x), non-90 arc angles including a
	negative one, and a tilted straight.  90-deg arcs alone would hide
	arc endpoint errors (cos(90) = 0 cancels the broken term), so the
	angles here are deliberately not 90.
	"""
	th = np.deg2rad(30.0)
	rot_x = np.array([
		[1.0, 0.0, 0.0],
		[0.0, np.cos(th), np.sin(th)],
		[0.0, -np.sin(th), np.cos(th)],
	])
	return (cb.CoilBuilder(current=current)
		.set_start([0.03, -0.02, 0.01], orientation=rot_x)
		.set_cross_section(0.003, 0.002)
		.add_straight(0.05)
		.add_arc(0.02, 60)
		.add_straight(0.04, tilt=15)
		.add_arc(0.015, -75)
		.add_straight(0.03)
		.add_arc(0.025, 120))


def _check_mirror_geometry(coil, plane):
	"""Assert invariants 1-5 for coil.mirror(plane)."""
	M = MIRROR_MATRICES[plane]
	mirrored = coil.mirror(plane)

	assert mirrored.current == coil.current, (
		f"[{plane}] mirror() must preserve the builder current "
		f"(got {mirrored.current}, want {coil.current})")
	assert len(mirrored.segments) == len(coil.segments)

	for k, (seg, mseg) in enumerate(zip(coil.segments, mirrored.segments)):
		assert type(mseg) is type(seg), f"[{plane}] seg{k} type changed"
		assert mseg.current == seg.current, (
			f"[{plane}] seg{k} current changed")

		err_start = np.max(np.abs(
			np.asarray(mseg.start_pos) - M @ np.asarray(seg.start_pos)))
		assert err_start < GEOM_TOL, (
			f"[{plane}] seg{k} start_pos != M @ start_pos (err={err_start:.3e})")

		err_end = np.max(np.abs(
			np.asarray(mseg.end_pos) - M @ np.asarray(seg.end_pos)))
		assert err_end < GEOM_TOL, (
			f"[{plane}] seg{k} end_pos != M @ end_pos (err={err_end:.3e})")

		det = np.linalg.det(np.asarray(mseg.orientation))
		assert abs(det - 1.0) < GEOM_TOL, (
			f"[{plane}] seg{k} orientation is not a proper rotation "
			f"(det={det:.6f})")

		if k + 1 < len(mirrored.segments):
			gap = np.max(np.abs(
				np.asarray(mirrored.segments[k + 1].start_pos)
				- np.asarray(mseg.end_pos)))
			assert gap < GEOM_TOL, (
				f"[{plane}] chain gap after seg{k} (gap={gap:.3e})")

	# Invariant 4: the whole wire path (arc interiors included) maps to M @ p.
	wires_orig, i_orig = coil.to_wire_segments(n_arc=16)
	wires_mir, i_mir = mirrored.to_wire_segments(n_arc=16)
	assert i_mir == i_orig
	assert len(wires_mir) == len(wires_orig)
	worst = 0.0
	for (p1, p2), (q1, q2) in zip(wires_orig, wires_mir):
		worst = max(worst,
			np.max(np.abs(np.asarray(q1) - M @ np.asarray(p1))),
			np.max(np.abs(np.asarray(q2) - M @ np.asarray(p2))))
	assert worst < GEOM_TOL, (
		f"[{plane}] wire path != M @ path (worst={worst:.3e}); "
		f"arc interiors or straight headings are misplaced")
	print(f"  [OK] plane={plane}: {len(coil.segments)} segments, "
		f"wire-path error {worst:.1e}")


def test_mirror_segment_geometry_racetrack():
	"""Bug-report racetrack: per-segment end_pos == M @ end_pos, all planes."""
	for plane in ('xy', 'yz', 'xz'):
		_check_mirror_geometry(build_racetrack(z0=0.05), plane)


def test_mirror_segment_geometry_general():
	"""General tilted mixed coil (non-90 arcs, negative arc, tilt), all planes."""
	for plane in ('xy', 'yz', 'xz'):
		_check_mirror_geometry(build_general_open_coil(), plane)


def test_mirror_xy_pair_midplane_field():
	"""x-symmetric racetrack pair upper + upper.mirror('xy'):

	the mirrored lower coil alone must reproduce the upper coil's central
	Bz (the buggy mirror delivered only ~24% of it), and the pair midplane
	Bz must be even in x and y to < 1e-8 relative to B0.
	"""
	import radia as rad

	rad.UtiDelAll()
	try:
		z0 = 0.05
		upper = build_racetrack(+z0)
		lower = upper.mirror('xy')
		explicit_lower = build_racetrack(-z0)

		bz_up = rad.Fld(rad.ObjCnt(upper.to_radia()), 'bz', [0, 0, 0])
		bz_lo = rad.Fld(rad.ObjCnt(lower.to_radia()), 'bz', [0, 0, 0])
		bz_ex = rad.Fld(rad.ObjCnt(explicit_lower.to_radia()), 'bz', [0, 0, 0])
		print(f"  Bz(0,0,0): upper={bz_up:.6e}  mirror-lower={bz_lo:.6e}  "
			f"explicit-lower={bz_ex:.6e}")
		assert abs(bz_up) > 1e-6, "upper coil field unexpectedly small"
		assert abs(bz_lo - bz_up) < FIELD_REL_TOL * abs(bz_up), (
			f"mirror('xy') lower coil central Bz {bz_lo:.6e} != upper "
			f"{bz_up:.6e} (the 2026-07-10 straight-heading bug gave ~24%)")
		assert abs(bz_lo - bz_ex) < FIELD_REL_TOL * abs(bz_ex), (
			"mirror('xy') lower != explicitly rebuilt lower loop at -z")

		pair = rad.ObjCnt(upper.to_radia() + lower.to_radia())
		b0 = rad.Fld(pair, 'bz', [0, 0, 0])
		assert abs(b0 - 2.0 * bz_up) < FIELD_REL_TOL * abs(b0), (
			"pair central Bz must be twice the single-coil value")

		worst_odd_x = 0.0
		worst_odd_y = 0.0
		for d in (0.01, 0.02, 0.03):
			bp = rad.Fld(pair, 'bz', [+d, 0.0, 0.0])
			bm = rad.Fld(pair, 'bz', [-d, 0.0, 0.0])
			worst_odd_x = max(worst_odd_x, abs(bp - bm) / 2.0)
			bp = rad.Fld(pair, 'bz', [0.0, +d, 0.0])
			bm = rad.Fld(pair, 'bz', [0.0, -d, 0.0])
			worst_odd_y = max(worst_odd_y, abs(bp - bm) / 2.0)
		rel_x = worst_odd_x / abs(b0)
		rel_y = worst_odd_y / abs(b0)
		print(f"  pair B0={b0:.6e}  odd-x/B0={rel_x:.3e}  odd-y/B0={rel_y:.3e}")
		assert rel_x < ODD_PART_REL_TOL, (
			f"midplane Bz odd-in-x part {rel_x:.3e} exceeds {ODD_PART_REL_TOL}")
		assert rel_y < ODD_PART_REL_TOL, (
			f"midplane Bz odd-in-y part {rel_y:.3e} exceeds {ODD_PART_REL_TOL}")
	finally:
		rad.UtiDelAll()


def main():
	print("=" * 70)
	print("CoilBuilder.mirror() regression test")
	print("=" * 70)
	print("[1] racetrack geometry invariants (xy / yz / xz)")
	test_mirror_segment_geometry_racetrack()
	print("[2] general tilted mixed-coil geometry invariants (xy / yz / xz)")
	test_mirror_segment_geometry_general()
	print("[3] mirror('xy') pair midplane field symmetry")
	test_mirror_xy_pair_midplane_field()
	print("\n[OK] mirror() traces M @ p for every segment and the racetrack "
		"pair midplane field is even in x and y.")


if __name__ == "__main__":
	main()
