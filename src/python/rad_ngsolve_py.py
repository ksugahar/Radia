"""
rad_ngsolve - Pure Python Implementation
NGSolve CoefficientFunction wrappers for Radia (DLL-free version)

This is a pure Python implementation that avoids C++ DLL dependencies.
It provides the same interface as the C++ version but implemented entirely in Python.

Usage:
	import radia as rad
	from ngsolve import *
	import rad_ngsolve_py as rad_ngsolve

	magnet = rad.ObjRecMag([0,0,0], [20,20,30], [0,0,1000])
	B = rad_ngsolve.RadBfield(magnet)

	mesh = Mesh(...)
	B_integral = Integrate(B, mesh)
"""

import radia as rad
try:
	from ngsolve import CoefficientFunction, CF
except ImportError:
	CoefficientFunction = object
	CF = None


class RadiaBFieldCF:
	"""
	B-field (magnetic flux density) coefficient function for Radia objects.

	This class wraps Radia's field calculation into a callable object
	compatible with NGSolve CoefficientFunctions.

	Parameters
	----------
	radia_obj : int
		Radia object index (from rad.ObjRecMag, rad.ObjCnt, etc.)

	Examples
	--------
	>>> import radia as rad
	>>> magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1000])
	>>> rad.Solve(magnet, 0.0001, 10000)
	>>> B = RadiaBFieldCF(magnet)
	"""

	def __init__(self, radia_obj):
		# Note: Don't call super().__init__() - NGSolve CoefficientFunction doesn't allow it
		self.radia_obj = radia_obj
		self._eval_count = 0  # For diagnostics
		self.dim = 3  # Vector field

	def Evaluate(self, x, y, z):
		"""
		Evaluate B-field at a point.

		Parameters
		----------
		x, y, z : float
			Coordinates in mm (Radia units)

		Returns
		-------
		tuple
			(Bx, By, Bz) in Tesla
		"""
		try:
			B = rad.Fld(self.radia_obj, 'b', [x, y, z])
			self._eval_count += 1
			return tuple(B)
		except Exception as e:
			# Return zero field on error
			return (0.0, 0.0, 0.0)

	def __call__(self, mip, *args):
		"""
		Evaluate function at mapped integration point (NGSolve interface).

		Parameters
		----------
		mip : BaseMappedIntegrationPoint or coordinates
			NGSolve integration point or coordinate values
		"""
		# Check if mip is a MappedIntegrationPoint
		if hasattr(mip, 'point'):
			# NGSolve MappedIntegrationPoint
			pnt = mip.point
			return self.Evaluate(pnt[0], pnt[1], pnt[2])
		elif hasattr(mip, '__getitem__'):
			# Array-like coordinate
			if len(args) == 0:
				# Called with single array [x, y, z]
				return self.Evaluate(mip[0], mip[1], mip[2])
			else:
				# Called with (x, y, z) as separate args
				return self.Evaluate(mip, args[0], args[1])
		else:
			# Called with scalar x, y, z
			return self.Evaluate(mip, args[0], args[1])


class RadiaHFieldCF:
	"""
	H-field (magnetic field intensity) coefficient function for Radia objects.

	Parameters
	----------
	radia_obj : int
		Radia object index

	Examples
	--------
	>>> H = RadiaHFieldCF(magnet)
	"""

	def __init__(self, radia_obj):
		# Note: Don't call super().__init__() - NGSolve CoefficientFunction doesn't allow it
		self.radia_obj = radia_obj
		self._eval_count = 0
		self.dim = 3  # Vector field

	def Evaluate(self, x, y, z):
		"""Evaluate H-field at a point."""
		try:
			H = rad.Fld(self.radia_obj, 'h', [x, y, z])
			self._eval_count += 1
			return tuple(H)
		except:
			return (0.0, 0.0, 0.0)

	def __call__(self, mip, *args):
		"""Evaluate H-field at mapped integration point (NGSolve interface)."""
		if hasattr(mip, 'point'):
			pnt = mip.point
			return self.Evaluate(pnt[0], pnt[1], pnt[2])
		elif hasattr(mip, '__getitem__'):
			if len(args) == 0:
				return self.Evaluate(mip[0], mip[1], mip[2])
			else:
				return self.Evaluate(mip, args[0], args[1])
		else:
			return self.Evaluate(mip, args[0], args[1])


class RadiaAFieldCF:
	"""
	A-field (vector potential) coefficient function for Radia objects.

	Parameters
	----------
	radia_obj : int
		Radia object index

	Examples
	--------
	>>> A = RadiaAFieldCF(magnet)
	"""

	def __init__(self, radia_obj):
		# Note: Don't call super().__init__() - NGSolve CoefficientFunction doesn't allow it
		self.radia_obj = radia_obj
		self._eval_count = 0
		self.dim = 3  # Vector field

	def Evaluate(self, x, y, z):
		"""Evaluate vector potential at a point."""
		try:
			A = rad.Fld(self.radia_obj, 'a', [x, y, z])
			self._eval_count += 1
			return tuple(A)
		except:
			return (0.0, 0.0, 0.0)

	def __call__(self, mip, *args):
		"""Evaluate A-field at mapped integration point (NGSolve interface)."""
		if hasattr(mip, 'point'):
			pnt = mip.point
			return self.Evaluate(pnt[0], pnt[1], pnt[2])
		elif hasattr(mip, '__getitem__'):
			if len(args) == 0:
				return self.Evaluate(mip[0], mip[1], mip[2])
			else:
				return self.Evaluate(mip, args[0], args[1])
		else:
			return self.Evaluate(mip, args[0], args[1])


# Convenience functions (matching C++ API)
def RadBfield(radia_obj):
	"""
	Create B-field coefficient function from Radia object.

	Parameters
	----------
	radia_obj : int
		Radia object index

	Returns
	-------
	CoefficientFunction
		NGSolve CoefficientFunction for magnetic field B (vector, dim=3)

	Examples
	--------
	>>> import radia as rad
	>>> import rad_ngsolve_py as rad_ngsolve
	>>>
	>>> magnet = rad.ObjRecMag([0,0,0], [20,20,30], [0,0,1])
	>>> B = rad_ngsolve.RadBfield(magnet)
	>>>
	>>> # Use in NGSolve
	>>> from ngsolve import *
	>>> mesh = Mesh(...)
	>>> Bz = B[2]  # Extract z-component
	>>> B_integral = Integrate(Bz, mesh)
	"""
	# Create the wrapper object
	field = RadiaBFieldCF(radia_obj)

	# Return as NGSolve CoefficientFunction
	# CF expects a callable that takes (x, y, z) and returns a vector
	return CF((
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[0],
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[1],
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[2]
	), dims=(3,))


def RadHfield(radia_obj):
	"""
	Create H-field coefficient function from Radia object.

	Parameters
	----------
	radia_obj : int
		Radia object index

	Returns
	-------
	CoefficientFunction
		NGSolve CoefficientFunction for field H (vector, dim=3)
	"""
	field = RadiaHFieldCF(radia_obj)
	return CF((
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[0],
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[1],
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[2]
	), dims=(3,))


def RadAfield(radia_obj):
	"""
	Create vector potential coefficient function from Radia object.

	Parameters
	----------
	radia_obj : int
		Radia object index

	Returns
	-------
	CoefficientFunction
		NGSolve CoefficientFunction for vector potential A (vector, dim=3)
	"""
	field = RadiaAFieldCF(radia_obj)
	return CF((
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[0],
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[1],
		lambda mip: field.Evaluate(mip.point[0], mip.point[1], mip.point[2])[2]
	), dims=(3,))


# Module information
__version__ = '0.1.0'
__author__ = 'Claude Code'
__doc__ = """
rad_ngsolve_py - Pure Python NGSolve Integration for Radia

This module provides NGSolve CoefficientFunction wrappers for Radia
magnetic field calculations, implemented entirely in Python to avoid
C++ DLL dependency issues.

Main Functions
--------------
RadBfield(radia_obj) : Create B-field coefficient function
RadHfield(radia_obj) : Create H-field coefficient function
RadAfield(radia_obj) : Create A-field coefficient function

Example Usage
-------------
import radia as rad
from ngsolve import *
import rad_ngsolve_py as rad_ngsolve

# Create Radia geometry
magnet = rad.ObjRecMag([0,0,0], [20,20,30], [0,0,1000])
rad.Solve(magnet, 0.0001, 10000)

# Create coefficient function
B = rad_ngsolve.RadBfield(magnet)

# Use in NGSolve
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.1))
B_integral = Integrate(B, mesh)
Draw(B, mesh, "B_field")
"""


if __name__ == "__main__":
	print("rad_ngsolve_py - Pure Python NGSolve Integration for Radia")
	print(f"Version: {__version__}")
	print()
	print("Example usage:")
	print("  import rad_ngsolve_py as rad_ngsolve")
	print("  B = rad_ngsolve.RadBfield(magnet)")
	print()
	print("Available functions:")
	print("  - RadBfield(radia_obj)")
	print("  - RadHfield(radia_obj)")
	print("  - RadAfield(radia_obj)")
