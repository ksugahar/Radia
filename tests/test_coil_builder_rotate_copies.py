"""Pointwise geometry contract for CoilBuilder.rotate_copies()."""

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from radia.coil_builder import CoilBuilder


GEOM_TOL = 2e-12


def _build_general_coil():
	angle = np.deg2rad(30.0)
	orientation = np.array([
		[1.0, 0.0, 0.0],
		[0.0, np.cos(angle), np.sin(angle)],
		[0.0, -np.sin(angle), np.cos(angle)],
	])
	return (CoilBuilder(current=500.0)
		.set_start([0.03, -0.02, 0.01], orientation=orientation)
		.set_cross_section(0.003, 0.002)
		.add_straight(0.05)
		.add_arc(0.02, 60)
		.add_straight(0.04, tilt=15)
		.add_arc(0.015, -75)
		.add_straight(0.03)
		.add_arc(0.025, 120))


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_rotate_copies_maps_complete_wire_path(axis):
	coil = _build_general_coil()
	copies = coil.rotate_copies(axis=axis, n_copies=4)
	axis_vec = {
		"x": np.array([1.0, 0.0, 0.0]),
		"y": np.array([0.0, 1.0, 0.0]),
		"z": np.array([0.0, 0.0, 1.0]),
	}[axis]
	original_wires, original_current = coil.to_wire_segments(n_arc=17)

	assert copies[0] is coil
	for index, rotated in enumerate(copies[1:], start=1):
		R = Rotation.from_rotvec(index * np.pi / 2 * axis_vec).as_matrix()
		assert rotated.current == coil.current
		assert len(rotated.segments) == len(coil.segments)

		for segment, rotated_segment in zip(coil.segments, rotated.segments):
			assert type(rotated_segment) is type(segment)
			assert rotated_segment.current == segment.current
			np.testing.assert_allclose(
				rotated_segment.start_pos, R @ segment.start_pos,
				atol=GEOM_TOL, rtol=0.0)
			np.testing.assert_allclose(
				rotated_segment.end_pos, R @ segment.end_pos,
				atol=GEOM_TOL, rtol=0.0)
			assert abs(np.linalg.det(rotated_segment.orientation) - 1.0) < GEOM_TOL

		rotated_wires, rotated_current = rotated.to_wire_segments(n_arc=17)
		assert rotated_current == original_current
		assert len(rotated_wires) == len(original_wires)
		for (p1, p2), (q1, q2) in zip(original_wires, rotated_wires):
			np.testing.assert_allclose(q1, R @ p1, atol=GEOM_TOL, rtol=0.0)
			np.testing.assert_allclose(q2, R @ p2, atol=GEOM_TOL, rtol=0.0)


@pytest.mark.parametrize(
	("kwargs", "message"),
	[
		({"axis": "q"}, "Unknown axis"),
		({"n_copies": 0}, "positive integer"),
		({"n_copies": 2.5}, "positive integer"),
	],
)
def test_rotate_copies_rejects_invalid_arguments(kwargs, message):
	with pytest.raises(ValueError, match=message):
		_build_general_coil().rotate_copies(**kwargs)
