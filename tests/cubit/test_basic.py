"""Basic tests for cubit_mesh_export module.

Note: These tests require Cubit to be installed and available.
They are primarily placeholder tests to establish the testing structure.
"""

import unittest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))

import cubit_mesh_export


class TestModuleImport(unittest.TestCase):
	"""Test that the module can be imported successfully."""

	def test_module_exists(self):
		"""Test that cubit_mesh_export module exists."""
		self.assertIsNotNone(cubit_mesh_export)

	def test_functions_exist(self):
		"""Test that all main export functions exist."""
		self.assertTrue(hasattr(cubit_mesh_export, 'export_exodus'))
		self.assertTrue(hasattr(cubit_mesh_export, 'export_Gmesh'))
		self.assertTrue(hasattr(cubit_mesh_export, 'export_gmesh'))
		self.assertTrue(hasattr(cubit_mesh_export, 'export_nastran'))
		self.assertTrue(hasattr(cubit_mesh_export, 'export_NGSolveCurvedMesh'))
		self.assertTrue(hasattr(cubit_mesh_export, 'export_meg'))
		self.assertTrue(hasattr(cubit_mesh_export, 'export_vtk'))
		self.assertTrue(hasattr(cubit_mesh_export, 'export_vtu'))

	def test_functions_callable(self):
		"""Test that export functions are callable."""
		self.assertTrue(callable(cubit_mesh_export.export_exodus))
		self.assertTrue(callable(cubit_mesh_export.export_Gmesh))
		self.assertTrue(callable(cubit_mesh_export.export_gmesh))
		self.assertTrue(callable(cubit_mesh_export.export_nastran))
		self.assertTrue(callable(cubit_mesh_export.export_NGSolveCurvedMesh))
		self.assertTrue(callable(cubit_mesh_export.export_meg))
		self.assertTrue(callable(cubit_mesh_export.export_vtk))
		self.assertTrue(callable(cubit_mesh_export.export_vtu))

	def test_aliases_point_to_same_function(self):
		"""Test that snake_case aliases point to the same functions."""
		self.assertIs(cubit_mesh_export.export_gmesh, cubit_mesh_export.export_Gmesh)
		self.assertIs(cubit_mesh_export.export_nastran, cubit_mesh_export.export_Nastran)

	def test_export_gmesh_signature(self):
		"""Test export_Gmesh has correct parameters."""
		import inspect
		sig = inspect.signature(cubit_mesh_export.export_Gmesh)
		params = list(sig.parameters.keys())
		self.assertIn('cubit', params)
		self.assertIn('FileName', params)
		self.assertIn('version', params)
		self.assertIn('DIM', params)
		self.assertEqual(sig.parameters['version'].default, "2.2")
		self.assertEqual(sig.parameters['DIM'].default, "auto")


class TestFunctionSignatures(unittest.TestCase):
	"""Test function signatures and default parameters."""

	def test_export_nastran_defaults(self):
		"""Test export_nastran has correct default parameters."""
		import inspect
		sig = inspect.signature(cubit_mesh_export.export_nastran)
		self.assertEqual(sig.parameters['DIM'].default, "3D")
		self.assertEqual(sig.parameters['PYRAM'].default, True)

	def test_export_meg_defaults(self):
		"""Test export_meg has correct default parameters."""
		import inspect
		sig = inspect.signature(cubit_mesh_export.export_meg)
		self.assertEqual(sig.parameters['DIM'].default, 'T')
		self.assertEqual(sig.parameters['MGR2'].default, None)

	def test_export_vtk_signature(self):
		"""Test export_vtk has correct parameters."""
		import inspect
		sig = inspect.signature(cubit_mesh_export.export_vtk)
		params = list(sig.parameters.keys())
		self.assertIn('cubit', params)
		self.assertIn('FileName', params)
		# Note: ORDER parameter was removed in v1.5.1 (auto-detection)

	def test_export_vtu_defaults(self):
		"""Test export_vtu has correct default parameters."""
		import inspect
		sig = inspect.signature(cubit_mesh_export.export_vtu)
		self.assertEqual(sig.parameters['binary'].default, False)

	def test_export_NGSolveCurvedMesh_signature(self):
		"""Test export_NGSolveCurvedMesh has correct parameters."""
		import inspect
		sig = inspect.signature(cubit_mesh_export.export_NGSolveCurvedMesh)
		params = list(sig.parameters.keys())
		self.assertIn('cubit', params)
		self.assertIn('order', params)

	def test_export_exodus_defaults(self):
		"""Test export_exodus has correct default parameters."""
		import inspect
		sig = inspect.signature(cubit_mesh_export.export_exodus)
		self.assertEqual(sig.parameters['overwrite'].default, True)
		self.assertEqual(sig.parameters['large_model'].default, False)


if __name__ == '__main__':
	unittest.main()
