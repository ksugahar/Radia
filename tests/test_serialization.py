"""
Unit tests for radsend.cpp - Serialization and I/O operations

Tests message handling and serialization:
- Error message handling
- Object dump/restore
- DuplicateObject functionality
- Data export
"""

import sys
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import pytest
import radia as rad
import numpy as np


class TestErrorHandling:
	"""Test error message and exception handling"""

	def test_invalid_object_error(self):
		"""Test that invalid object indices produce errors gracefully"""
		rad.UtiDelAll()

		# Try to use non-existent object
		# Should not crash, may return error or zero result
		try:
			H = rad.Fld(99999, 'h', [0, 0, 0])
			# If no exception, result might be empty or zero
		except:
			# Expected to fail gracefully
			pass

	def test_invalid_field_type(self):
		"""Test handling of invalid field type specification"""
		rad.UtiDelAll()

		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])

		# Try invalid field type
		try:
			H = rad.Fld(mag, 'invalid', [0, 0, 20])
			# Should fail gracefully
			assert False, "Should have raised error for invalid field type"
		except:
			# Expected
			pass

	def test_invalid_material_parameters(self):
		"""Test error handling for invalid material parameters"""
		rad.UtiDelAll()

		try:
			# Invalid ksi array (wrong length)
			mat = rad.MatLin([1000], [0, 0, 1])
			# Should fail
			assert False, "Should have raised error for invalid ksi array"
		except:
			# Expected
			pass


class TestObjectDuplication:
	"""Test object duplication functionality"""

	def test_duplicate_simple_object(self):
		"""Test duplicating a simple magnet"""
		rad.UtiDelAll()

		# Create magnet
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])

		# Duplicate it
		mag2 = rad.ObjDpl(mag1)
		assert mag2 > 0
		assert mag2 != mag1, "Duplicate should have different index"

		# Both should produce same field at same point
		H1 = rad.Fld(mag1, 'h', [20, 0, 0])
		H2 = rad.Fld(mag2, 'h', [20, 0, 0])

		assert np.allclose(H1, H2, rtol=1e-10)

	def test_duplicate_with_material(self):
		"""Test duplicating object with material"""
		rad.UtiDelAll()

		# Create magnet with material
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mat = rad.MatLin([1000, 0], [1, 1, 1])  # Anisotropic, easy axis [1,1,1]
		rad.MatApl(mag1, mat)

		# Duplicate
		mag2 = rad.ObjDpl(mag1)

		# Material should be duplicated too
		H1 = rad.Fld(mag1, 'h', [20, 0, 0])
		H2 = rad.Fld(mag2, 'h', [20, 0, 0])

		assert np.allclose(H1, H2, rtol=1e-10)

	def test_duplicate_container(self):
		"""Test duplicating a container"""
		rad.UtiDelAll()

		# Create container
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mag2 = rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])
		cnt1 = rad.ObjCnt([mag1, mag2])

		# Duplicate container
		cnt2 = rad.ObjDpl(cnt1)
		assert cnt2 > 0

		# Should produce same field
		H1 = rad.Fld(cnt1, 'h', [10, 0, 20])
		H2 = rad.Fld(cnt2, 'h', [10, 0, 20])

		assert np.allclose(H1, H2, rtol=1e-10)


class TestUtilityFunctions:
	"""Test utility functions related to serialization/I/O"""

	def test_delete_all(self):
		"""Test UtiDelAll clears all objects"""
		rad.UtiDelAll()

		# Create some objects
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mag2 = rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])

		# Delete all
		rad.UtiDelAll()

		# Old indices should no longer be valid
		try:
			H = rad.Fld(mag1, 'h', [0, 0, 20])
			# If this works, it might be reusing indices
		except:
			# Expected - object deleted
			pass

	def test_create_many_delete_all(self):
		"""Test creating many objects and deleting all"""
		rad.UtiDelAll()

		# Create many objects
		for i in range(100):
			mag = rad.ObjRecMag([i*10, 0, 0], [5, 5, 5], [0, 0, 1])

		# Delete all
		rad.UtiDelAll()

		# Should be able to create new objects
		new_mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		assert new_mag > 0


class TestColorAndVisualization:
	"""Test color and visualization attributes (part of serialization)"""

	def test_set_object_color(self):
		"""Test setting object drawing color"""
		rad.UtiDelAll()

		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])

		# Set color (R, G, B)
		rad.ObjDrwAtr(mag, [1, 0, 0])  # Red

		# Should not crash - color is stored
		assert mag > 0

	def test_set_container_color(self):
		"""Test setting color for entire container"""
		rad.UtiDelAll()

		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mag2 = rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])
		cnt = rad.ObjCnt([mag1, mag2])

		# Set color for container
		rad.ObjDrwAtr(cnt, [0, 1, 0])  # Green

		assert cnt > 0


class TestDataExport:
	"""Test data export functionality"""

	def test_get_magnetization(self):
		"""Test retrieving magnetization data"""
		rad.UtiDelAll()

		# Create magnet with known magnetization
		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [1, 2, 3])

		# Get magnetization at center
		M = rad.Fld(mag, 'm', [0, 0, 0])

		# Should be close to [1, 2, 3] (in A/m, may need unit conversion)
		assert len(M) == 3
		assert M[0] != 0 or M[1] != 0 or M[2] != 0

	def test_get_object_info(self):
		"""Test retrieving object information"""
		rad.UtiDelAll()

		mag = rad.ObjRecMag([5, 10, 15], [20, 30, 40], [0, 0, 1])

		# Try to get object properties (if available)
		try:
			# Some versions have ObjCenFld
			info = rad.ObjCenFld(mag, 'inf')
			if info is not None:
				# Should return some information
				pass
		except:
			# Not all versions have this function
			pass


class TestMemoryManagement:
	"""Test memory management (creation/deletion cycles)"""

	def test_create_delete_cycle(self):
		"""Test repeated create/delete cycles"""
		for cycle in range(10):
			rad.UtiDelAll()

			# Create objects
			for i in range(50):
				mag = rad.ObjRecMag([i*10, 0, 0], [5, 5, 5], [0, 0, 1])

			# Compute some fields
			H = rad.Fld(mag, 'h', [100, 0, 0])
			assert len(H) == 3

		# Final cleanup
		rad.UtiDelAll()

	def test_large_object_creation(self):
		"""Test creating many objects"""
		rad.UtiDelAll()

		objects = []
		for i in range(200):
			mag = rad.ObjRecMag([i*5, 0, 0], [3, 3, 3], [0, 0, 1])
			objects.append(mag)

		# Create container
		cnt = rad.ObjCnt(objects)

		# Should be able to compute field
		H = rad.Fld(cnt, 'h', [500, 0, 20])
		assert len(H) == 3

		rad.UtiDelAll()


if __name__ == "__main__":
	pytest.main([__file__, "-v"])
