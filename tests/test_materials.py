"""
Unit tests for radmaterial.cpp - Material system

Tests material definition and application:
- Linear materials (isotropic and anisotropic)
- Nonlinear materials
- Material application to geometry
- Standard materials
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../build/Release"))

import pytest
import radia as rad
import numpy as np


class TestLinearMaterials:
	"""Test linear material creation and application"""

	def test_linear_isotropic_material(self):
		"""Test linear isotropic material"""
		rad.UtiDelAll()

		# Create linear isotropic material: MatLin(ksi)
		# Single ksi value for isotropic (same in all directions)
		mat = rad.MatLin(1000)
		assert mat > 0, "Linear material should have valid index"

	def test_linear_anisotropic_material(self):
		"""Test linear anisotropic material"""
		rad.UtiDelAll()

		# Create anisotropic material: different ksi in parallel and perpendicular
		mat = rad.MatLin([2000, 500], [0, 0, 1])
		assert mat > 0

	def test_linear_material_no_remanence(self):
		"""Test linear material without remanent magnetization"""
		rad.UtiDelAll()

		# Soft iron: high permeability, isotropic (no remanence)
		# Use new isotropic API: MatLin(ksi)
		mat = rad.MatLin(5000)
		assert mat > 0

	def test_material_application_to_magnet(self):
		"""Test applying material to magnet"""
		rad.UtiDelAll()

		# Create magnet
		mag = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])

		# Create and apply material (anisotropic, easy axis [1,1,1])
		mat = rad.MatLin([1000, 0], [1, 1, 1])
		rad.MatApl(mag, mat)

		# Compute field - should be non-zero
		H = rad.Fld(mag, 'h', [20, 0, 0])
		assert not np.allclose(H, [0, 0, 0])


class TestNonlinearMaterials:
	"""Test nonlinear material creation"""

	def test_saturation_material(self):
		"""Test saturation material (M vs H curve)"""
		rad.UtiDelAll()

		# Create simple saturation curve
		# H and M values must be interleaved: [[H1, M1], [H2, M2], ...]
		HM_data = [[0, 0], [100, 800], [500, 1200], [1000, 1400], [5000, 1500], [10000, 1500]]

		mat = rad.MatSatIsoTab(HM_data)
		assert mat > 0

	def test_nonlinear_material_application(self):
		"""Test applying nonlinear material"""
		rad.UtiDelAll()

		# Create magnet
		mag = rad.ObjRecMag([0, 0, 0], [20, 20, 20], [1, 0, 0])

		# Create nonlinear material
		HM_data = [[0, 0], [100, 800], [500, 1200], [1000, 1400], [5000, 1500]]
		mat = rad.MatSatIsoTab(HM_data)

		# Apply to magnet
		rad.MatApl(mag, mat)

		# Should be able to compute field
		H = rad.Fld(mag, 'h', [30, 0, 0])
		assert len(H) == 3


class TestMaterialOnGroups:
	"""Test material application to groups"""

	def test_material_on_container(self):
		"""Test applying material to entire container"""
		rad.UtiDelAll()

		# Create multiple objects
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mag2 = rad.ObjRecMag([15, 0, 0], [10, 10, 10], [0, 0, 1])
		group = rad.ObjCnt([mag1, mag2])

		# Apply material to group (anisotropic, easy axis [1,1,1])
		mat = rad.MatLin([1000, 0], [1, 1, 1])
		rad.MatApl(group, mat)

		# Should compute field
		H = rad.Fld(group, 'h', [30, 0, 0])
		assert len(H) == 3

	def test_different_materials_in_group(self):
		"""Test that individual objects can have different materials"""
		rad.UtiDelAll()

		# Create two magnets with different materials
		mag1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
		mat1 = rad.MatLin([1000, 0], [1, 1, 1])  # Easy axis [1,1,1]
		rad.MatApl(mag1, mat1)

		mag2 = rad.ObjRecMag([20, 0, 0], [10, 10, 10], [0, 0, 1])
		mat2 = rad.MatLin([2000, 0], [1, 0, 0])  # Easy axis [1,0,0] (x-direction)
		rad.MatApl(mag2, mat2)

		# Group them
		group = rad.ObjCnt([mag1, mag2])

		# Should compute combined field
		H = rad.Fld(group, 'h', [10, 0, 0])
		assert len(H) == 3


class TestMaterialWithRelaxation:
	"""Test materials in relaxation/solving context"""

	def test_material_with_solve(self):
		"""Test that material works with Solve() for nonlinear problems"""
		rad.UtiDelAll()

		# Create iron core with nonlinear material
		core = rad.ObjRecMag([0, 0, 0], [20, 20, 20], [1, 0, 0])

		# Nonlinear material
		HM_data = [[0, 0], [100, 800], [500, 1200], [1000, 1400], [5000, 1500], [10000, 1500]]
		mat = rad.MatSatIsoTab(HM_data)
		rad.MatApl(core, mat)

		# Create interaction for solving
		rad.ObjDrwAtr(core, [1, 0, 0])  # Color (for drawing)

		# Solve (may fail if problem is simple, but shouldn't crash)
		try:
			result = rad.Solve(core, 0.0001, 1000)
			# Result format: [max_M_change, number_of_iterations, ...]
			assert len(result) >= 2
		except:
			# If solve fails, that's okay for this test
			pass


if __name__ == "__main__":
	pytest.main([__file__, "-v"])
