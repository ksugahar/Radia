"""
Unit tests for radgroup.cpp - Object grouping operations

Tests ObjCnt (Container/Group) functionality:
- Group creation
- Adding/removing objects
- Group transformations
- Nested groups
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../build/Release"))

import pytest
import radia as rad
import numpy as np


class TestGroupCreation:
	"""Test basic group creation and structure"""

	def test_empty_group(self):
		"""Test creating an empty group"""
		rad.UtiDelAll()
		group = rad.ObjCnt([])
		assert group > 0, "Empty group should have valid index"

	def test_single_object_group(self):
		"""Test group with single object"""
		rad.UtiDelAll()
		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		group = rad.ObjCnt([mag])
		assert group > 0

	def test_multiple_objects_group(self):
		"""Test group with multiple objects"""
		rad.UtiDelAll()

		mags = []
		for i in range(5):
			mag = rad.ObjRecMag([i*20, 0, 0], [10, 10, 10], [0, 0, 1])
			mags.append(mag)

		group = rad.ObjCnt(mags)
		assert group > 0

		# Test field evaluation on group
		H = rad.Fld(group, 'h', [50, 0, 0])
		assert len(H) == 3, "Should return 3D field vector"

	def test_nested_groups(self):
		"""Test creating nested groups"""
		rad.UtiDelAll()

		# Create two sub-groups
		group1_mags = [
			rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1]),
			rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])
		]
		group1 = rad.ObjCnt(group1_mags)

		group2_mags = [
			rad.ObjRecMag([0, 20, 0], [10, 10, 10], [0, 0, 1]),
			rad.ObjRecMag([20, 20, 0], [10, 10, 10], [0, 0, 1])
		]
		group2 = rad.ObjCnt(group2_mags)

		# Create parent group
		parent_group = rad.ObjCnt([group1, group2])
		assert parent_group > 0

		# Test field from nested group
		H = rad.Fld(parent_group, 'h', [10, 10, 0])
		assert len(H) == 3


class TestGroupTransformations:
	"""Test transformations applied to groups"""

	def test_translate_group(self):
		"""Test translating entire group"""
		rad.UtiDelAll()

		# Create group at origin
		mags = []
		for i in range(3):
			mag = rad.ObjRecMag([i*10, 0, 0], [5, 5, 5], [0, 0, 1])
			mags.append(mag)
		group = rad.ObjCnt(mags)

		# Field at origin before translation
		H_before = rad.Fld(group, 'h', [0, 0, 20])

		# Translate group
		rad.TrfOrnt(group, rad.TrfTrsl([0, 100, 0]))

		# Field at translated position
		H_after = rad.Fld(group, 'h', [0, 100, 20])

		# Fields should be similar (geometry moved)
		assert np.allclose(H_before, H_after, rtol=1e-6)

	def test_rotate_group(self):
		"""Test rotating entire group"""
		rad.UtiDelAll()

		# Create asymmetric group
		mag1 = rad.ObjRecMag([10, 0, 0], [5, 5, 5], [1, 0, 0])
		mag2 = rad.ObjRecMag([20, 0, 0], [5, 5, 5], [1, 0, 0])
		group = rad.ObjCnt([mag1, mag2])

		# Rotate 90 degrees around z-axis
		rad.TrfOrnt(group, rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi/2))

		# After rotation, magnets should be along y-axis
		# Field at rotated position should make sense
		H = rad.Fld(group, 'h', [0, 15, 0])
		assert len(H) == 3

	def test_multiple_transformations(self):
		"""Test combining multiple transformations"""
		rad.UtiDelAll()

		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		group = rad.ObjCnt([mag])

		# Translate
		rad.TrfOrnt(group, rad.TrfTrsl([10, 0, 0]))

		# Rotate
		rad.TrfOrnt(group, rad.TrfRot([0, 0, 0], [0, 0, 1], np.pi/4))

		# Should still compute field
		H = rad.Fld(group, 'h', [20, 0, 0])
		assert len(H) == 3


class TestGroupMaterialApplication:
	"""Test applying materials to groups"""

	def test_material_to_group(self):
		"""Test applying material to entire group"""
		rad.UtiDelAll()

		# Create group without material
		mags = []
		for i in range(3):
			mag = rad.ObjRecMag([i*15, 0, 0], [10, 10, 10], [0, 0, 1])
			mags.append(mag)
		group = rad.ObjCnt(mags)

		# Apply material to group (anisotropic material with easy axis in [1,1,1] direction)
		mat = rad.MatLin([1000, 0], [1, 1, 1])  # Easy axis direction
		rad.MatApl(group, mat)

		# Field should be computed
		H = rad.Fld(group, 'h', [20, 0, 0])
		assert len(H) == 3
		assert not np.allclose(H, [0, 0, 0])  # Should have non-zero field


class TestGroupFieldEvaluation:
	"""Test field evaluation from groups"""

	def test_field_from_group_vs_individual(self):
		"""Compare field from group vs sum of individual objects"""
		rad.UtiDelAll()

		# Create individual magnets
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		rad.MatApl(mag1, rad.MatLin([1000, 0], [1, 1, 1]))  # Anisotropic, easy axis [1,1,1]

		mag2 = rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])
		rad.MatApl(mag2, rad.MatLin([1000, 0], [1, 1, 1]))  # Anisotropic, easy axis [1,1,1]

		# Field from individual objects
		H1 = np.array(rad.Fld(mag1, 'h', [30, 0, 0]))
		H2 = np.array(rad.Fld(mag2, 'h', [30, 0, 0]))
		H_sum = H1 + H2

		# Field from group
		group = rad.ObjCnt([mag1, mag2])
		H_group = np.array(rad.Fld(group, 'h', [30, 0, 0]))

		# Should be equal (superposition principle)
		assert np.allclose(H_sum, H_group, rtol=1e-10)

	def test_batch_field_from_group(self):
		"""Test FldBatch with groups"""
		rad.UtiDelAll()

		mags = []
		for i in range(5):
			mag = rad.ObjRecMag([i*15, 0, 0], [10, 10, 10], [0, 0, 1])
			rad.MatApl(mag, rad.MatLin([1000, 0], [1, 1, 1]))  # Anisotropic, easy axis [1,1,1]
			mags.append(mag)

		group = rad.ObjCnt(mags)

		# Evaluate at multiple points
		points = [[x, 0, 0] for x in [10, 20, 30, 40]]
		H_batch = rad.FldBatch(group, 'h', points, 0)  # Direct calculation

		assert len(H_batch) == 4, "Should return 4 field vectors"
		for H in H_batch:
			assert len(H) == 3, "Each field should be 3D"


class TestGroupEdgeCases:
	"""Test edge cases and error handling"""

	def test_group_with_duplicate_objects(self):
		"""Test group containing same object multiple times"""
		rad.UtiDelAll()

		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])

		# Create group with duplicates (should handle gracefully)
		group = rad.ObjCnt([mag, mag, mag])
		assert group > 0

	def test_deeply_nested_groups(self):
		"""Test deeply nested group hierarchy"""
		rad.UtiDelAll()

		# Create deep nesting: group -> group -> group -> magnet
		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		level1 = rad.ObjCnt([mag])
		level2 = rad.ObjCnt([level1])
		level3 = rad.ObjCnt([level2])

		# Should still compute field
		H = rad.Fld(level3, 'h', [20, 0, 0])
		assert len(H) == 3


if __name__ == "__main__":
	pytest.main([__file__, "-v"])
