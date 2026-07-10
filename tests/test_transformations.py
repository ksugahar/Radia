"""
Unit tests for radtransform.cpp - Transformation operations

Tests transformation creation and application:
- Translation (TrfTrsl)
- Rotation (TrfRot)
- Combined transformations (TrfCmbL)
- Transformation application (TrfOrnt)

TrfOrnt(obj, tr) MOVES the field source: subsequent field evaluations see the
object at the transformed position (classic Radia semantics; the field of the
moved object, not an extra copy).  Until 2026-07-10 polyhedron elements
(ObjHexahedron etc., incl. rad.magnet_box) silently IGNORED their transform
list, which is why older revisions of this file claimed TrfOrnt "adds a
symmetry copy" -- that description matched the bug, not the intended
behavior.  See tests/test_radiafield_transformed_container.py for the full
regression suite of that fix.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

import pytest
import radia as rad
import numpy as np


class TestTranslation:
	"""Test translation transformations"""

	def test_create_translation(self):
		"""Test creating translation transformation"""
		rad.UtiDelAll()

		# Create translation vector (in meters)
		tr = rad.TrfTrsl([0.010, 0.020, 0.030])
		assert tr > 0, "Translation should have valid index"

	def test_apply_translation(self):
		"""TrfOrnt translation moves the source: the field at the translated
		observation point equals the field of a magnet built at the
		translated position."""
		rad.UtiDelAll()

		# Magnet built directly at the translated position (reference)
		mag_ref = rad.magnet_box([0, 0, 0.05], [0.01, 0.01, 0.01], [0, 0, 954930])
		H_ref = rad.Fld(mag_ref, 'h', [0.05, 0, 0.05])

		rad.UtiDelAll()

		# Magnet at origin, then TrfOrnt-translated by the same offset
		mag = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])
		rad.TrfOrnt(mag, rad.TrfTrsl([0, 0, 0.05]))
		H_moved = rad.Fld(mag, 'h', [0.05, 0, 0.05])

		assert np.linalg.norm(np.subtract(H_moved, H_ref)) <= \
			1e-10 * np.linalg.norm(H_ref)

	def test_multiple_translations(self):
		"""Two successive TrfOrnt translations compose (total offset)."""
		rad.UtiDelAll()

		mag_ref = rad.magnet_box([0.010, 0.020, 0], [0.01, 0.01, 0.01],
		                         [0, 0, 954930])
		H_ref = rad.Fld(mag_ref, 'h', [0.010, 0.020, 0.030])

		rad.UtiDelAll()

		mag = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])
		rad.TrfOrnt(mag, rad.TrfTrsl([0.010, 0, 0]))
		rad.TrfOrnt(mag, rad.TrfTrsl([0, 0.020, 0]))

		H = rad.Fld(mag, 'h', [0.010, 0.020, 0.030])
		assert len(H) == 3
		assert np.linalg.norm(np.subtract(H, H_ref)) <= \
			1e-10 * np.linalg.norm(H_ref)


class TestRotation:
	"""Test rotation transformations"""

	def test_create_rotation(self):
		"""Test creating rotation transformation"""
		rad.UtiDelAll()

		# Rotation around z-axis by 90 degrees
		tr = rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi/2)
		assert tr > 0

	def test_apply_rotation_90deg(self):
		"""TrfOrnt 90-degree rotation moves the source: the rotated model
		reproduces the field of a magnet built directly at the rotated
		position/orientation."""
		rad.UtiDelAll()

		# Reference: magnet along y-axis with y-magnetization
		mag_ref = rad.magnet_box([0, 0.01, 0], [0.005, 0.005, 0.005], [0, 1, 0])
		H_ref = rad.Fld(mag_ref, 'h', [0, 0.015, 0])

		rad.UtiDelAll()

		# Magnet along x-axis with x-magnetization, rotated +90deg about z
		mag = rad.magnet_box([0.01, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])
		rad.TrfOrnt(mag, rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi / 2))
		H_rot = rad.Fld(mag, 'h', [0, 0.015, 0])

		assert np.linalg.norm(np.subtract(H_rot, H_ref)) <= \
			1e-9 * np.linalg.norm(H_ref)

	def test_rotation_180deg(self):
		"""TrfOrnt 180-degree rotation about z moves the magnet to -x."""
		rad.UtiDelAll()

		mag_ref = rad.magnet_box([-0.010, 0, 0], [0.005, 0.005, 0.005],
		                         [-1, 0, 0])
		H_ref = rad.Fld(mag_ref, 'h', [-0.015, 0, 0])

		rad.UtiDelAll()

		mag = rad.magnet_box([0.010, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])
		tr = rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi)
		rad.TrfOrnt(mag, tr)

		H = rad.Fld(mag, 'h', [-0.015, 0, 0])
		assert len(H) == 3
		assert np.linalg.norm(np.subtract(H, H_ref)) <= \
			1e-9 * np.linalg.norm(H_ref)

	def test_rotation_around_arbitrary_point(self):
		"""Test rotation around non-origin point"""
		rad.UtiDelAll()

		mag = rad.magnet_box([0.010, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])

		# Rotate around point [0.010, 0, 0] (magnet center)
		tr = rad.TrfRot([0.010, 0, 0], [0, 0, 1], np.pi/2)
		rad.TrfOrnt(mag, tr)

		# Verify field computation succeeds
		H = rad.Fld(mag, 'h', [0.010, 0, 0.020])
		assert len(H) == 3


class TestCombinedTransformations:
	"""Test combining multiple transformations"""

	def test_combine_two_translations(self):
		"""Test TrfCmbL with two translations: creation succeeds and the
		combined transform can be applied and evaluated."""
		rad.UtiDelAll()

		# Create two translations (in meters)
		tr1 = rad.TrfTrsl([0.010, 0, 0])
		tr2 = rad.TrfTrsl([0, 0.020, 0])

		# Combine them
		tr_combined = rad.TrfCmbL(tr1, tr2)
		assert tr_combined > 0

		# Apply to magnet (moves the source by the combined offset)
		mag = rad.magnet_box([0, 0, 0], [0.005, 0.005, 0.005], [0, 0, 954930])
		rad.TrfOrnt(mag, tr_combined)

		# Verify field computation succeeds
		H = rad.Fld(mag, 'h', [0.010, 0.020, 0.020])
		assert len(H) == 3

	def test_combine_rotation_and_translation(self):
		"""Test combining rotation and translation: the combined transform
		can be created, applied, and evaluated."""
		rad.UtiDelAll()

		# First rotate, then translate (in meters)
		tr_rot = rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi/2)
		tr_trsl = rad.TrfTrsl([0.050, 0, 0])

		# Combine
		tr_combined = rad.TrfCmbL(tr_rot, tr_trsl)

		# Apply to magnet (moves the source)
		mag = rad.magnet_box([0.010, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])
		rad.TrfOrnt(mag, tr_combined)

		# Verify field computation succeeds
		H = rad.Fld(mag, 'h', [0.050, 0.010, 0])
		assert len(H) == 3

	# test_multiply_transformation REMOVED (2026-01-31)
	# TrfMlt has been removed - use IMA symmetry instead
	# See docs/solver/IMA_SYMMETRY_DESIGN.md for the correct approach


# TestInversion REMOVED (2026-03-06)
# TrfInv() takes no arguments in current API - it creates an identity/inversion
# transform, not an inverse of a specific transform handle. The previous test
# incorrectly passed a transform handle to TrfInv().


class TestTransformationOnGroups:
	"""Test transformations applied to groups (containers).

	TrfOrnt on a group moves ALL members together; field evaluation applies
	the group transform to the observation point / field without mutating
	any member (thread-safe for batch evaluation -- see
	tests/test_radiafield_transformed_container.py).
	"""

	def test_transform_container(self):
		"""Translating a container moves the fields of all its members."""
		rad.UtiDelAll()

		mag1r = rad.magnet_box([0, 0.050, 0], [0.005, 0.005, 0.005],
		                       [0, 0, 954930])
		mag2r = rad.magnet_box([0.010, 0.050, 0], [0.005, 0.005, 0.005],
		                       [0, 0, 954930])
		group_ref = rad.ObjCnt([mag1r, mag2r])
		H_ref = rad.Fld(group_ref, 'h', [0.005, 0.050, 0.020])

		rad.UtiDelAll()

		mag1 = rad.magnet_box([0, 0, 0], [0.005, 0.005, 0.005], [0, 0, 954930])
		mag2 = rad.magnet_box([0.010, 0, 0], [0.005, 0.005, 0.005], [0, 0, 954930])
		group = rad.ObjCnt([mag1, mag2])
		rad.TrfOrnt(group, rad.TrfTrsl([0, 0.050, 0]))

		H = rad.Fld(group, 'h', [0.005, 0.050, 0.020])
		assert len(H) == 3
		assert np.linalg.norm(np.subtract(H, H_ref)) <= \
			1e-10 * np.linalg.norm(H_ref)

	def test_rotate_container(self):
		"""Rotating a container moves the fields of all its members."""
		rad.UtiDelAll()

		mag1r = rad.magnet_box([0, 0.010, 0], [0.005, 0.005, 0.005], [0, 1, 0])
		mag2r = rad.magnet_box([0, 0.020, 0], [0.005, 0.005, 0.005], [0, 1, 0])
		group_ref = rad.ObjCnt([mag1r, mag2r])
		H_ref = rad.Fld(group_ref, 'h', [0, 0.015, 0])

		rad.UtiDelAll()

		mag1 = rad.magnet_box([0.010, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])
		mag2 = rad.magnet_box([0.020, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])
		group = rad.ObjCnt([mag1, mag2])
		rad.TrfOrnt(group, rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi/2))

		H = rad.Fld(group, 'h', [0, 0.015, 0])
		assert len(H) == 3
		assert np.linalg.norm(np.subtract(H, H_ref)) <= \
			1e-9 * np.linalg.norm(H_ref)


# TestTransformationSymmetry REMOVED (2026-01-31)
# TrfMlt has been removed - use IMA symmetry instead
# See docs/solver/IMA_SYMMETRY_DESIGN.md for the correct approach
#
# For symmetric structures (quadrupole, arrays), use:
# 1. Explicit element creation (recommended)
# 2. IMA symmetry for plane symmetry with SetIMASymmetry() and BuildIMAMatrix()


if __name__ == "__main__":
	pytest.main([__file__, "-v"])
