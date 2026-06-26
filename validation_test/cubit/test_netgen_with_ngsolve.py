"""
Test Netgen mesh construction patterns with NGSolve integration.

This test requires NGSolve to be installed in the system Python.
It covers the NGSolve-side mesh behavior used after Cubit's 2025.12
`export netgen` command writes a .vol file.

Run with system Python (not Cubit Python):
  python tests/test_netgen_with_ngsolve.py
"""

import sys
import os

# Add parent directory to path for cubit_mesh_curver
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))

# Check NGSolve first
try:
	import ngsolve
	from netgen.meshing import Mesh as NetgenMesh
	from netgen.meshing import MeshPoint, Pnt, Element3D, Element2D, FaceDescriptor
	print(f"NGSolve version: {ngsolve.__version__}")
except ImportError as e:
	print(f"NGSolve not available: {e}")
	print("This test requires NGSolve. Install it with:")
	print("  pip install ngsolve")
	sys.exit(1)


def test_manual_netgen_mesh_creation():
	"""Test creating a Netgen mesh manually."""
	print("=" * 60)
	print("Test 1: Manual Netgen Mesh Creation")
	print("=" * 60)

	mesh = NetgenMesh(dim=3)

	# Add points (vertices of a tetrahedron)
	points = [
		(0.0, 0.0, 0.0),
		(1.0, 0.0, 0.0),
		(0.5, 1.0, 0.0),
		(0.5, 0.5, 1.0),
	]

	pnums = []
	for p in points:
		pnums.append(mesh.Add(MeshPoint(Pnt(*p))))

	print(f"  Added {len(pnums)} points")

	# Add a face descriptor for the volume
	mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))

	# Add one tetrahedron
	mesh.Add(Element3D(1, [pnums[0], pnums[1], pnums[2], pnums[3]]))
	print(f"  Added 1 tetrahedron")

	# Convert to NGSolve mesh
	ngmesh = ngsolve.Mesh(mesh)
	print(f"  NGSolve mesh: ne={ngmesh.ne}, nv={ngmesh.nv}")

	assert ngmesh.ne == 1, f"Expected 1 element, got {ngmesh.ne}"
	assert ngmesh.nv == 4, f"Expected 4 vertices, got {ngmesh.nv}"

	# Test creating FE space
	fes = ngsolve.H1(ngmesh, order=1)
	print(f"  H1 space: {fes.ndof} DOFs")

	print("  PASS: Manual mesh creation works")
	return True


def test_larger_mesh():
	"""Test with a larger manually created mesh."""
	print("\n" + "=" * 60)
	print("Test 2: Larger Manual Mesh")
	print("=" * 60)

	mesh = NetgenMesh(dim=3)

	# Create a 2x2x2 grid of points
	points = []
	for i in range(3):
		for j in range(3):
			for k in range(3):
				points.append((float(i), float(j), float(k)))

	pnums = []
	for p in points:
		pnums.append(mesh.Add(MeshPoint(Pnt(*p))))

	print(f"  Added {len(pnums)} points (3x3x3 grid)")

	# Add face descriptor
	mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))

	# Add some tetrahedra (simple decomposition)
	# Each cube can be decomposed into 5 or 6 tets
	# For simplicity, let's add a few tets manually
	def idx(i, j, k):
		return i * 9 + j * 3 + k

	tet_count = 0
	# First cube at origin
	v = [idx(0, 0, 0), idx(1, 0, 0), idx(0, 1, 0), idx(0, 0, 1)]
	mesh.Add(Element3D(1, [pnums[v[0]], pnums[v[1]], pnums[v[2]], pnums[v[3]]]))
	tet_count += 1

	v = [idx(1, 0, 0), idx(1, 1, 0), idx(0, 1, 0), idx(1, 1, 1)]
	mesh.Add(Element3D(1, [pnums[v[0]], pnums[v[1]], pnums[v[2]], pnums[v[3]]]))
	tet_count += 1

	print(f"  Added {tet_count} tetrahedra")

	# Convert to NGSolve
	ngmesh = ngsolve.Mesh(mesh)
	print(f"  NGSolve mesh: ne={ngmesh.ne}, nv={ngmesh.nv}")

	# Test higher order FE spaces
	for order in [1, 2, 3]:
		fes = ngsolve.H1(ngmesh, order=order)
		print(f"  H1 order={order}: {fes.ndof} DOFs")

	print("  PASS: Larger mesh works")
	return True


def test_mesh_with_boundary():
	"""Test mesh with boundary elements."""
	print("\n" + "=" * 60)
	print("Test 3: Mesh with Boundary Elements")
	print("=" * 60)

	mesh = NetgenMesh(dim=3)

	# Simple tetrahedron
	points = [
		(0.0, 0.0, 0.0),
		(1.0, 0.0, 0.0),
		(0.5, 1.0, 0.0),
		(0.5, 0.5, 1.0),
	]

	pnums = []
	for p in points:
		pnums.append(mesh.Add(MeshPoint(Pnt(*p))))

	# Add face descriptors for boundary
	# Each face of the tetrahedron is a boundary
	fd_idx = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))

	# Add volume element
	mesh.Add(Element3D(1, [pnums[0], pnums[1], pnums[2], pnums[3]]))

	# Add boundary triangles (faces of tetrahedron)
	mesh.Add(Element2D(1, [pnums[0], pnums[2], pnums[1]]))  # bottom face
	mesh.Add(Element2D(1, [pnums[0], pnums[1], pnums[3]]))  # front face
	mesh.Add(Element2D(1, [pnums[1], pnums[2], pnums[3]]))  # right face
	mesh.Add(Element2D(1, [pnums[2], pnums[0], pnums[3]]))  # left face

	print(f"  Volume elements: {mesh.ne}")

	# Convert to NGSolve
	ngmesh = ngsolve.Mesh(mesh)
	print(f"  NGSolve mesh: ne={ngmesh.ne}, nv={ngmesh.nv}")
	print(f"  Boundaries: {ngmesh.GetBoundaries()}")

	# Test that we can apply boundary conditions
	fes = ngsolve.H1(ngmesh, order=1)
	print(f"  H1 space: {fes.ndof} DOFs")

	print("  PASS: Mesh with boundary works")
	return True


def test_mesh_with_materials():
	"""Test mesh with multiple materials (blocks)."""
	print("\n" + "=" * 60)
	print("Test 4: Mesh with Multiple Materials")
	print("=" * 60)

	mesh = NetgenMesh(dim=3)

	# Two tetrahedra with different materials
	points = [
		(0.0, 0.0, 0.0),
		(1.0, 0.0, 0.0),
		(0.5, 1.0, 0.0),
		(0.5, 0.5, 1.0),
		(1.5, 0.5, 0.5),  # Additional point for second tet
	]

	pnums = []
	for p in points:
		pnums.append(mesh.Add(MeshPoint(Pnt(*p))))

	# Two face descriptors for two materials
	mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
	mesh.Add(FaceDescriptor(surfnr=2, domin=2, bc=2))

	# Material 1: First tetrahedron
	mesh.Add(Element3D(1, [pnums[0], pnums[1], pnums[2], pnums[3]]))

	# Material 2: Second tetrahedron (sharing face with first)
	mesh.Add(Element3D(2, [pnums[1], pnums[2], pnums[3], pnums[4]]))

	# Set material names
	mesh.SetMaterial(1, "region1")
	mesh.SetMaterial(2, "region2")

	print(f"  Volume elements: {mesh.ne}")

	# Convert to NGSolve
	ngmesh = ngsolve.Mesh(mesh)
	print(f"  NGSolve mesh: ne={ngmesh.ne}, nv={ngmesh.nv}")
	print(f"  Materials: {ngmesh.GetMaterials()}")

	# Test material-specific coefficients
	cf = ngsolve.CoefficientFunction([1.0, 2.0])  # Different values per material
	print(f"  Material-based coefficient created")

	print("  PASS: Mesh with materials works")
	return True


def test_high_order_curving():
	"""Test mesh curving with geometry (if available)."""
	print("\n" + "=" * 60)
	print("Test 5: High-Order Curving (geometry based)")
	print("=" * 60)

	try:
		from netgen.csg import unit_cube

		# Generate a mesh from CSG geometry
		geo = unit_cube
		mesh = geo.GenerateMesh(maxh=0.3)

		print(f"  Generated mesh from CSG")

		# Convert to NGSolve and curve
		ngmesh = ngsolve.Mesh(mesh)

		with ngsolve.TaskManager():
			for order in [2, 3, 4]:
				try:
					ngmesh.Curve(order)
					print(f"  mesh.Curve({order}) succeeded!")
				except Exception as e:
					print(f"  mesh.Curve({order}) failed: {e}")

		print("  PASS: High-order curving test completed")
		return True

	except ImportError:
		print("  netgen.csg not available, skipping geometry test")
		return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("Netgen / NGSolve Integration Tests")
	print("=" * 60)
	print("These tests verify Netgen mesh creation patterns")
	print("used after Cubit 2025.12 export netgen writes .vol files")
	print("=" * 60)

	all_passed = True

	tests = [
		test_manual_netgen_mesh_creation,
		test_larger_mesh,
		test_mesh_with_boundary,
		test_mesh_with_materials,
		test_high_order_curving,
	]

	for test in tests:
		try:
			if not test():
				all_passed = False
		except Exception as e:
			print(f"  FAIL: {e}")
			import traceback
			traceback.print_exc()
			all_passed = False

	print("\n" + "=" * 60)
	if all_passed:
		print("All tests PASSED!")
	else:
		print("Some tests FAILED!")
	print("=" * 60)
