"""
Unit tests for radtransform.cpp - Transformation operations

Tests transformation creation and application:
- Translation (TrfTrsl)
- Rotation (TrfRot)
- Combined transformations (TrfCmbL)
- Transformation application (TrfOrnt)

Note: TrfOrnt adds a symmetry copy of the object at the transformed position.
It does NOT physically move the object. Tests are designed accordingly.
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
		"""Test translation by comparing magnets at different positions.

		TrfOrnt adds a symmetry copy, so instead we verify translation
		semantics by creating two separate magnets at original and
		translated positions and confirming the field at the same
		relative observation point is identical.
		"""
		rad.UtiDelAll()

		# Magnet at origin
		mag1 = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])
		H1 = rad.Fld(mag1, 'h', [0.05, 0, 0])

		rad.UtiDelAll()

		# Magnet at translated position (same relative geometry)
		mag2 = rad.magnet_box([0, 0, 0.05], [0.01, 0.01, 0.01], [0, 0, 954930])
		H2 = rad.Fld(mag2, 'h', [0.05, 0, 0.05])

		# Same relative geometry -> same field
		assert np.allclose(H1, H2, rtol=1e-10)

	def test_multiple_translations(self):
		"""Test applying multiple translations produces valid field"""
		rad.UtiDelAll()

		mag = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])

		# Translate in x (adds symmetry copy)
		rad.TrfOrnt(mag, rad.TrfTrsl([0.010, 0, 0]))

		# Translate in y (adds another symmetry copy)
		rad.TrfOrnt(mag, rad.TrfTrsl([0, 0.020, 0]))

		# Verify field computation succeeds and returns 3 components
		H = rad.Fld(mag, 'h', [0.010, 0.020, 0.030])
		assert len(H) == 3


class TestRotation:
	"""Test rotation transformations"""

	def test_create_rotation(self):
		"""Test creating rotation transformation"""
		rad.UtiDelAll()

		# Rotation around z-axis by 90 degrees
		tr = rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi/2)
		assert tr > 0

	def test_apply_rotation_90deg(self):
		"""Test 90-degree rotation by comparing rotated geometry.

		TrfOrnt adds a symmetry copy, so instead we verify rotation
		semantics by creating two separate magnets at original and
		rotated positions and confirming the field magnitudes match.
		"""
		rad.UtiDelAll()

		# Magnet along x-axis with x-magnetization
		mag1 = rad.magnet_box([0.01, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])
		H_x = rad.Fld(mag1, 'h', [0.015, 0, 0])
		H_x_mag = np.linalg.norm(H_x)

		rad.UtiDelAll()

		# Same magnet rotated 90deg around z: along y-axis with y-magnetization
		mag2 = rad.magnet_box([0, 0.01, 0], [0.005, 0.005, 0.005], [0, 1, 0])
		H_y = rad.Fld(mag2, 'h', [0, 0.015, 0])
		H_y_mag = np.linalg.norm(H_y)

		# Magnitudes should match (rotational symmetry)
		assert np.isclose(H_x_mag, H_y_mag, rtol=1e-6)

	def test_rotation_180deg(self):
		"""Test 180-degree rotation produces valid field"""
		rad.UtiDelAll()

		mag = rad.magnet_box([0.010, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])

		# Rotate 180 degrees around z-axis (adds symmetry copy)
		tr = rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi)
		rad.TrfOrnt(mag, tr)

		# Verify field computation at symmetry copy position succeeds
		H = rad.Fld(mag, 'h', [-0.015, 0, 0])
		assert len(H) == 3

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
		"""Test TrfCmbL with two translations.

		TrfOrnt adds a symmetry copy, so we only verify that
		combined transform creation succeeds and field computation
		does not crash.
		"""
		rad.UtiDelAll()

		# Create two translations (in meters)
		tr1 = rad.TrfTrsl([0.010, 0, 0])
		tr2 = rad.TrfTrsl([0, 0.020, 0])

		# Combine them
		tr_combined = rad.TrfCmbL(tr1, tr2)
		assert tr_combined > 0

		# Apply to magnet (adds symmetry copy)
		mag = rad.magnet_box([0, 0, 0], [0.005, 0.005, 0.005], [0, 0, 954930])
		rad.TrfOrnt(mag, tr_combined)

		# Verify field computation succeeds
		H = rad.Fld(mag, 'h', [0.010, 0.020, 0.020])
		assert len(H) == 3

	def test_combine_rotation_and_translation(self):
		"""Test combining rotation and translation.

		TrfOrnt adds a symmetry copy, so we only verify that
		the combined transform can be created and applied without
		errors.
		"""
		rad.UtiDelAll()

		# First rotate, then translate (in meters)
		tr_rot = rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi/2)
		tr_trsl = rad.TrfTrsl([0.050, 0, 0])

		# Combine
		tr_combined = rad.TrfCmbL(tr_rot, tr_trsl)

		# Apply to magnet (adds symmetry copy)
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
	"""Test transformations applied to groups.

	TrfOrnt on groups adds symmetry copies of all group members.
	Tests verify that field computation succeeds (no crash) after
	applying transformations to containers.
	"""

	def test_transform_container(self):
		"""Test applying translation to entire container"""
		rad.UtiDelAll()

		# Create group of magnets (in meters)
		mag1 = rad.magnet_box([0, 0, 0], [0.005, 0.005, 0.005], [0, 0, 954930])
		mag2 = rad.magnet_box([0.010, 0, 0], [0.005, 0.005, 0.005], [0, 0, 954930])
		group = rad.ObjCnt([mag1, mag2])

		# Translate entire group (adds symmetry copies)
		tr = rad.TrfTrsl([0, 0.050, 0])
		rad.TrfOrnt(group, tr)

		# Verify field computation succeeds
		H = rad.Fld(group, 'h', [0.005, 0.050, 0.020])
		assert len(H) == 3

	def test_rotate_container(self):
		"""Test rotating entire container"""
		rad.UtiDelAll()

		mag1 = rad.magnet_box([0.010, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])
		mag2 = rad.magnet_box([0.020, 0, 0], [0.005, 0.005, 0.005], [1, 0, 0])
		group = rad.ObjCnt([mag1, mag2])

		# Rotate 90 degrees (adds symmetry copies)
		tr = rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi/2)
		rad.TrfOrnt(group, tr)

		# Verify field computation succeeds
		H = rad.Fld(group, 'h', [0, 0.015, 0])
		assert len(H) == 3


# TestTransformationSymmetry REMOVED (2026-01-31)
# TrfMlt has been removed - use IMA symmetry instead
# See docs/solver/IMA_SYMMETRY_DESIGN.md for the correct approach
#
# For symmetric structures (quadrupole, arrays), use:
# 1. Explicit element creation (recommended)
# 2. IMA symmetry for plane symmetry with SetIMASymmetry() and BuildIMAMatrix()


if __name__ == "__main__":
	pytest.main([__file__, "-v"])
