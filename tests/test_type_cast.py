"""
Unit tests for radcast.cpp - Type casting operations

Tests radTCast functionality:
- InteractCast: Cast to radTInteraction
- g3dCast: Cast to radTg3d
- GroupCast: Cast to radTGroup
"""

import sys
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import pytest
import radia as rad
import numpy as np


class TestTypeCasting:
	"""Test type casting operations"""

	def test_cast_magnet_to_g3d(self):
		"""Test that rectangular magnet can be used as g3d object"""
		rad.UtiDelAll()

		# Create magnet - should be radTg3d type
		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		assert mag > 0

		# Should work with g3d operations (field computation)
		H = rad.Fld(mag, 'h', [20, 0, 0])
		assert len(H) == 3
		assert not np.allclose(H, [0, 0, 0])

	def test_cast_group_operations(self):
		"""Test that groups can be operated on"""
		rad.UtiDelAll()

		# Create multiple magnets and group them
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mag2 = rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])
		group = rad.ObjCnt([mag1, mag2])

		assert group > 0

		# Group should support field computation
		H = rad.Fld(group, 'h', [10, 0, 0])
		assert len(H) == 3

	def test_different_object_types(self):
		"""Test different geometry types all work as g3d objects"""
		rad.UtiDelAll()

		# Rectangular block
		rec = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		H_rec = rad.Fld(rec, 'h', [20, 0, 0])
		assert len(H_rec) == 3

		# Polyhedron (via extrusion)
		rad.UtiDelAll()
		points = [[0, 0], [10, 0], [10, 10], [0, 10]]
		poly = rad.ObjThckPgn(0, 10, points, "z", [0, 0, 1])
		H_poly = rad.Fld(poly, 'h', [20, 0, 0])
		assert len(H_poly) == 3

		# Container/Group
		rad.UtiDelAll()
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mag2 = rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])
		cnt = rad.ObjCnt([mag1, mag2])
		H_cnt = rad.Fld(cnt, 'h', [30, 0, 0])
		assert len(H_cnt) == 3

	def test_nested_container_casting(self):
		"""Test that nested containers work properly"""
		rad.UtiDelAll()

		# Create nested structure
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mag2 = rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])
		group1 = rad.ObjCnt([mag1, mag2])

		mag3 = rad.ObjRecMag([0, 20, 0], [10, 10, 10], [0, 0, 1])
		mag4 = rad.ObjRecMag([20, 20, 0], [10, 10, 10], [0, 0, 1])
		group2 = rad.ObjCnt([mag3, mag4])

		# Container of containers
		parent = rad.ObjCnt([group1, group2])
		assert parent > 0

		# Should compute field from nested structure
		H = rad.Fld(parent, 'h', [10, 10, 0])
		assert len(H) == 3

	def test_transformation_on_casted_object(self):
		"""Test that transformations work on casted objects"""
		rad.UtiDelAll()

		# Create object and transform it
		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])

		# Apply transformation (requires g3d cast internally)
		rad.TrfOrnt(mag, rad.TrfTrsl([10, 0, 0]))

		# Field should reflect translated position
		H = rad.Fld(mag, 'h', [15, 0, 0])
		assert len(H) == 3

	def test_material_application_casting(self):
		"""Test that material application works (requires casting)"""
		rad.UtiDelAll()

		# Create magnet
		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])

		# Apply material (internally casts to appropriate type)
		mat = rad.MatLin([1000, 0], [1, 1, 1])  # Anisotropic, easy axis [1,1,1]
		rad.MatApl(mag, mat)

		# Field should be computed with material
		H = rad.Fld(mag, 'h', [20, 0, 0])
		assert len(H) == 3
		assert not np.allclose(H, [0, 0, 0])


class TestInvalidCasting:
	"""Test error handling for invalid casting"""

	def test_invalid_object_index(self):
		"""Test field computation with invalid object index"""
		rad.UtiDelAll()

		# Try to compute field for non-existent object
		# Should fail gracefully (not crash)
		try:
			H = rad.Fld(99999, 'h', [0, 0, 0])
			# If it doesn't raise, it might return empty or zero
		except:
			# Expected to fail
			pass


if __name__ == "__main__":
	pytest.main([__file__, "-v"])
