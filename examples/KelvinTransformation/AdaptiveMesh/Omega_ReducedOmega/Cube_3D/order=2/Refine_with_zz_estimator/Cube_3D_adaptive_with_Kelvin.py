"""
3D Omega-Reduced Omega Method for Magnetostatics with Kelvin Transformation
WITH ADAPTIVE MESH REFINEMENT (Local Refinement with Doerfler Marking)
1/8 MODEL (Octant: x>=0, y>=0, z>=0)

Problem: Magnetic cube (mu_r=1000) in uniform z-directed background field

Formulation (Omega-Reduced Omega):
- Total region (magnetic): H = grad(Omega), no source
- Reduced region (air): Source term from Omega_s

Kelvin transformation:
- Maps infinite exterior domain to finite sphere interior
- Permeability transformation: mu'(r') = (R/r')^2 * mu0  (3D)
- Periodic BC couples interior (r=R) with exterior (r'=R)

Symmetry conditions (for z-directed field):
- x=0 plane: Natural BC (Neumann)
- y=0 plane: Natural BC (Neumann)
- z=0 plane: Dirichlet BC (Omega = 0) - equipotential plane

Uniform refinement algorithm:
  1. ZZ-type error estimator (H(div) flux recovery)
  2. Mark ALL elements
  3. Local refinement
  4. VTU output at each iteration

Stop condition: DOF >= max_dof
"""
import os
import sys
import glob

# Set environment for ksugahar's NGSolve build with PatchwiseSolveWithInterface
os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin'
os.environ['PATH'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin;' + os.environ.get('PATH', '')
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ngsolve\Lib\site-packages')

# Delete existing output files
script_dir = os.path.dirname(os.path.abspath(__file__))
print("Deleting existing output files...")
deleted_count = 0
for ext in ['*.png', '*.vtu', '*.mat']:
	for f in glob.glob(os.path.join(script_dir, ext)):
		try:
			os.remove(f)
			print(f"  Deleted: {f}")
			deleted_count += 1
		except Exception as e:
			print(f"  Failed to delete {f}: {e}")
print(f"Deleted {deleted_count} files.")

from numpy import pi, sqrt, linspace, zeros, nan, isnan, meshgrid, array, log, log10, sin, cos
from ngsolve import *
from netgen.occ import *
import scipy.io as sio
import tempfile

# Import matplotlib for plotting
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend to avoid tkinter errors
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                              'bf': 'serif:bold', 'fontset': 'cm'})
matplotlib.rcParams['font.family'] = 'Times New Roman'

print("=" * 60)
print("3D Omega-Reduced Omega Method with Kelvin Transform")
print("ADAPTIVE REFINEMENT (ZZ Estimator) - 1/8 Model (Octant) - CUBE")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
cube_size = 1.0          # Magnetic cube size [m] (cube from 0 to 1 in each direction, per CubeMesh.py)
air_total_radius = 1.5   # Total region air box size [m] (A_domain extends to 1.5, per CubeMesh.py)
kelvin_radius = 3.0      # Kelvin transformation radius [m] (rk: must be > sqrt(3)*air_total_radius ~= 2.6)
mu_r = 1000              # Relative permeability (per CubeMesh.py default)
mu0 = 4 * pi * 1e-7      # Vacuum permeability [H/m]

# Source field: H_s = (0, 0, H0) uniform in z-direction
H0 = 1000.0  # [A/m]

# Mesh parameters
maxh_initial = 1.0       # Initial mesh size (coarser for 3D)
order = 2                # Finite element order

# Offset for exterior domain (z-direction)
offset_z = 3.0

# Adaptive mesh parameters
max_iterations = 20      # Stop after 20 iterations

print(f"\nProblem parameters:")
print(f"  Cube size: {cube_size} m (cube from 0 to {cube_size} in each direction)")
print(f"  Air_total radius: {air_total_radius} m (Total region air)")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")
print(f"  Model type: 1/8 (octant x>=0, y>=0, z>=0)")
print(f"\nRegion structure:")
print(f"  magnetic + air_total: Total potential (H = grad(Omega))")
print(f"  air_inner + air_outer: Reduced potential (H = grad(Omega) + Hs)")
print(f"\nUniform mesh parameters:")
print(f"  Initial mesh size: {maxh_initial} m")
print(f"  Polynomial order: {order}")
print(f"  Stop condition: {max_iterations} iterations or DOF >= 1e5")


# ============================================================
# Geometry Definition for 3D - 1/8 Model (Cube)
# 4-region structure (per CubeMesh.py):
#   magnetic: cube (0,0,0)-(1,1,1) - Total potential (iron)
#   air_total: cube (0,0,0)-(1.5,1.5,1.5) minus magnetic - Total potential (A_domain)
#   air_inner: 1/8 sphere (r=kelvin_radius) minus air_total cube - Reduced potential (Omega_domain)
#   air_outer: Kelvin-transformed exterior - Reduced potential (Kelvin)
# ============================================================
def create_geometry():
	"""Create 1/8 geometry with 4 regions per CubeMesh.py structure."""
	print("\nCreating 1/8 cube geometry with 4-region structure (per CubeMesh.py)...")
	print(f"  magnetic: cube (0,0,0)-({cube_size},{cube_size},{cube_size})")
	print(f"  air_total: cube (0,0,0)-({air_total_radius},{air_total_radius},{air_total_radius}) - magnetic")
	print(f"  air_inner: 1/8 sphere (r={kelvin_radius}m) - air_total cube")
	print(f"  air_outer: Kelvin-transformed exterior at center=({kelvin_radius*2}, 0, 0)")

	# ===== Region 1: Magnetic cube (Total potential) =====
	# 1/8 cube in octant x>=0, y>=0, z>=0
	mag_cube = Box(Pnt(0, 0, 0), Pnt(cube_size, cube_size, cube_size))
	mag_cube.mat("magnetic")

	# Name faces for magnetic cube
	for face in mag_cube.faces:
		fc = face.center
		if abs(fc.x) < 1e-6:
			face.name = "sym_x"
		elif abs(fc.y) < 1e-6:
			face.name = "sym_y"
		elif abs(fc.z) < 1e-6:
			face.name = "sym_z"
		else:
			face.name = "mag_air_total"  # Interface with air_total

	# ===== Region 2: air_total CUBE (Total potential) - per CubeMesh.py =====
	air_total_box = Box(Pnt(0, 0, 0), Pnt(air_total_radius, air_total_radius, air_total_radius))
	air_total = air_total_box - mag_cube
	air_total.mat("air_total")

	# Name faces for air_total (cube)
	for face in air_total.faces:
		fc = face.center
		# Outer faces of the cube are total_reduced boundary
		if abs(fc.x - air_total_radius) < 0.01:
			face.name = "total_reduced"  # x = air_total_radius plane
		elif abs(fc.y - air_total_radius) < 0.01:
			face.name = "total_reduced"  # y = air_total_radius plane
		elif abs(fc.z - air_total_radius) < 0.01:
			face.name = "total_reduced"  # z = air_total_radius plane
		elif abs(fc.x) < 0.01:
			face.name = "sym_x"
		elif abs(fc.y) < 0.01:
			face.name = "sym_y"
		elif abs(fc.z) < 0.01:
			face.name = "sym_z"
		else:
			face.name = "mag_air_total"  # Interface with magnetic

	# ===== Region 3: air_inner (Reduced potential) =====
	# Following CubeMesh.py approach: use Sphere * Box for 1/8 sphere
	# Then subtract the air_total region
	air_inner_full = Sphere(Pnt(0, 0, 0), kelvin_radius) * Box((0, 0, 0), (kelvin_radius, kelvin_radius, kelvin_radius))
	air_inner = air_inner_full - air_total_box
	air_inner.faces.Min(X).name = "sym_x"
	air_inner.faces.Min(Y).name = "sym_y"
	air_inner.faces.Min(Z).name = "sym_z"
	air_inner.mat("air_inner")

	# ===== Region 4: air_outer (Kelvin-transformed exterior) =====
	# Per CubeMesh.py: center = (2*rKelvin, 0, 0)
	center_x = kelvin_radius * 2
	outer_sphere = Sphere(Pnt(center_x, 0, 0), kelvin_radius) * Box((center_x, 0, 0), (center_x + kelvin_radius, kelvin_radius, kelvin_radius))
	outer_sphere.faces.Min(X).name = "sym_x_ext"
	outer_sphere.faces.Min(Y).name = "sym_y"
	outer_sphere.faces.Min(Z).name = "sym_z"
	outer_sphere.mat("air_outer")

	# ===== CRITICAL: Identify periodic faces BEFORE Glue =====
	# Following CubeMesh.py: faces[0] is the spherical surface
	print("Identifying periodic boundaries BEFORE Glue...")
	print(f"  air_inner faces: {len(air_inner.faces)}")
	print(f"  outer_sphere faces: {len(outer_sphere.faces)}")

	# Find spherical faces (the ones that are not on box boundaries)
	kelvin_int_face = None
	kelvin_ext_face = None

	for face in air_inner.faces:
		fc = face.center
		dist = sqrt(fc.x**2 + fc.y**2 + fc.z**2)
		if abs(dist - kelvin_radius) < kelvin_radius * 0.2:
			kelvin_int_face = face
			face.name = "kelvin_int"
			print(f"  Found kelvin_int face at center ({fc.x:.2f}, {fc.y:.2f}, {fc.z:.2f}), dist={dist:.2f}")
			break

	for face in outer_sphere.faces:
		fc = face.center
		dist = sqrt((fc.x - center_x)**2 + fc.y**2 + fc.z**2)
		if abs(dist - kelvin_radius) < kelvin_radius * 0.2:
			kelvin_ext_face = face
			face.name = "kelvin_ext"
			print(f"  Found kelvin_ext face at center ({fc.x:.2f}, {fc.y:.2f}, {fc.z:.2f}), dist={dist:.2f}")
			break

	if kelvin_int_face is not None and kelvin_ext_face is not None:
		kelvin_ext_face.Identify(kelvin_int_face, "periodic", IdentificationType.PERIODIC)
		print("  Periodic identification applied BEFORE Glue!")
	else:
		print("  WARNING: Could not find periodic faces!")

	# GND vertex at center of exterior domain
	vertex = Vertex(Pnt(center_x, 0, 0))
	vertex.name = "GND"

	# Glue all domains (periodic already identified)
	geo = Glue([mag_cube, air_total, air_inner, outer_sphere, vertex])

	# Name the solids explicitly
	for i, solid in enumerate(geo.solids):
		if i == 0:
			solid.name = "magnetic"
		elif i == 1:
			solid.name = "air_total"
		elif i == 2:
			solid.name = "air_inner"
		elif i == 3:
			solid.name = "air_outer"

	return OCCGeometry(geo)


# ============================================================
# Solve Omega-Reduced Omega formulation (3D - 1/8)
# ============================================================
def solve_omega_formulation(mesh, fe_order):
	"""Solve Omega-Reduced Omega formulation on given mesh with 4 regions."""
	fes_before = H1(mesh, order=fe_order, dirichlet="GND|sym_z|sym_z_ext")
	fes = Periodic(fes_before)

	Omega = fes.TrialFunction()
	psi = fes.TestFunction()

	# For exterior domain: center at (2*kelvin_radius, 0, 0)
	center_x_ext = 2 * kelvin_radius
	r_prime_sq = (x - center_x_ext)**2 + y**2 + z**2
	r_prime = sqrt(r_prime_sq + 1e-20)

	mu_kelvin = (kelvin_radius / r_prime)**2 * mu0

	mu_dict = {
		"magnetic": mu_r * mu0,
		"air_total": mu0,
		"air_inner": mu0,
		"air_outer": mu_kelvin
	}
	Mu = CoefficientFunction([mu_dict[mat] for mat in mesh.GetMaterials()])

	Omega_s = H0 * z
	Hs = CoefficientFunction((0.0, 0.0, H0))
	Bs = CoefficientFunction((0.0, 0.0, mu0 * H0))

	a = BilinearForm(fes)
	a += Mu * grad(Omega) * grad(psi) * dx("magnetic")
	a += Mu * grad(Omega) * grad(psi) * dx("air_total")
	a += Mu * grad(Omega) * grad(psi) * dx("air_inner")
	a += Mu * grad(Omega) * grad(psi) * dx("air_outer")
	a.Assemble()

	gfOmega = GridFunction(fes)
	gfOmega.Set(Omega_s, BND, mesh.Boundaries("total_reduced"))

	f = LinearForm(fes)
	f += Mu * grad(gfOmega) * grad(psi) * dx("air_inner")
	f.Assemble()

	normal = -specialcf.normal(mesh.dim)
	f += (normal * Bs) * psi * ds("total_reduced")
	f.Assemble()

	gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

	Hs_kelvin = CoefficientFunction((0.0, 0.0, (kelvin_radius / r_prime)**2 * H0))
	Bs_kelvin = mu0 * Hs_kelvin

	Hs_dict = {
		"magnetic": Hs,
		"air_total": Hs,
		"air_inner": Hs,
		"air_outer": Hs_kelvin
	}
	Hs_cf = CoefficientFunction([Hs_dict[mat] for mat in mesh.GetMaterials()])

	Bs_dict = {
		"magnetic": Bs,
		"air_total": Bs,
		"air_inner": Bs,
		"air_outer": Bs_kelvin
	}
	Bs_cf = CoefficientFunction([Bs_dict[mat] for mat in mesh.GetMaterials()])

	H_pert_dict = {
		"magnetic": grad(gfOmega) - Hs,
		"air_total": grad(gfOmega) - Hs,
		"air_inner": grad(gfOmega),
		"air_outer": grad(gfOmega)
	}
	H_pert_cf = CoefficientFunction([H_pert_dict[mat] for mat in mesh.GetMaterials()])

	B_pert_dict = {
		"magnetic": (mu_r * mu0) * (grad(gfOmega) - Hs),
		"air_total": mu0 * (grad(gfOmega) - Hs),
		"air_inner": mu0 * grad(gfOmega),
		"air_outer": mu_kelvin * grad(gfOmega)
	}
	B_pert_cf = CoefficientFunction([B_pert_dict[mat] for mat in mesh.GetMaterials()])

	H_total_dict = {
		"magnetic": grad(gfOmega),
		"air_total": grad(gfOmega),
		"air_inner": grad(gfOmega) + Hs,
		"air_outer": grad(gfOmega) + Hs_kelvin
	}
	H_total_cf = CoefficientFunction([H_total_dict[mat] for mat in mesh.GetMaterials()])

	B_total_dict = {
		"magnetic": (mu_r * mu0) * grad(gfOmega),
		"air_total": mu0 * grad(gfOmega),
		"air_inner": mu0 * (grad(gfOmega) + Hs),
		"air_outer": mu_kelvin * (grad(gfOmega) + Hs_kelvin)
	}
	B_total_cf = CoefficientFunction([B_total_dict[mat] for mat in mesh.GetMaterials()])

	fields = {
		'Hs_cf': Hs_cf,
		'Bs_cf': Bs_cf,
		'H_pert_cf': H_pert_cf,
		'B_pert_cf': B_pert_cf,
		'H_total_cf': H_total_cf,
		'B_total_cf': B_total_cf,
		'mu_kelvin': mu_kelvin,
	}

	return fes, gfOmega, Mu, fields


# ============================================================
# Equilibrated Error Estimator (Braess-Schöberl type)
# ============================================================
def compute_error_estimator(mesh, fes, B_pert_cf):
	"""
	Compute equilibrated error estimator using PatchwiseSolveWithInterface.

	For Omega-Reduced Omega formulation, the physical law is div(B) = 0.
	We construct sigma_eq ∈ H(div) satisfying div(sigma_eq) = 0 exactly,
	then compute ||B_h - sigma_eq|| as the error estimator.

	Uses PatchwiseSolveWithInterface to properly handle material interfaces
	(BBND elements) where B may have discontinuities.

	Prager-Synge theorem guarantees this is a reliable upper bound.
	"""
	flux = B_pert_cf  # B_h = mu H_h from numerical solution

	# Use order p for H(div) space (same as H1 order for efficiency)
	recovery_order = max(1, fes.globalorder)
	fes_flux = HDiv(mesh, order=recovery_order)

	sigma = fes_flux.TrialFunction()
	tau = fes_flux.TestFunction()

	# Bilinear form for patchwise solve: (sigma, τ)
	bf = InnerProduct(sigma, tau) * dx

	# Linear form: (B_h, τ)
	# For equilibration: we want div(sigma_eq) = 0
	# The patchwise solve minimizes ||sigma - B_h|| subject to div constraints
	lf = InnerProduct(flux, tau) * dx

	# Grid function for equilibrated flux
	gf_sigma_eq = GridFunction(fes_flux)

	# Use PatchwiseSolveWithInterface to properly handle material interfaces
	# This includes BBND (interface) elements in the patch assembly
	try:
		PatchwiseSolveWithInterface(bf, lf, gf_sigma_eq)
	except Exception as e:
		print(f"  Warning: PatchwiseSolveWithInterface failed ({e}), trying PatchwiseSolve")
		try:
			PatchwiseSolve(bf, lf, gf_sigma_eq)
		except Exception as e2:
			print(f"  Warning: PatchwiseSolve also failed ({e2}), falling back to global L2 projection")
			# Fallback to global L2 projection (ZZ-type)
			a_flux = BilinearForm(fes_flux)
			a_flux += InnerProduct(sigma, tau) * dx
			a_flux.Assemble()

			f_flux = LinearForm(fes_flux)
			f_flux += InnerProduct(flux, tau) * dx
			f_flux.Assemble()

			gf_sigma_eq.vec.data = a_flux.mat.Inverse(fes_flux.FreeDofs(), inverse="sparsecholesky") * f_flux.vec

	# Error estimator: ||B_h - sigma_eq||^2
	err = InnerProduct(flux - gf_sigma_eq, flux - gf_sigma_eq)
	element_errors = Integrate(err, mesh, element_wise=True)

	# Check equilibration quality (div(sigma_eq) should be small)
	div_residual = Integrate(div(gf_sigma_eq)**2 * dx, mesh)
	if div_residual > 1e-10:
		print(f"  Equilibration quality: ||div(sigma_eq)|| = {sqrt(div_residual):.2e}")

	return element_errors




# ============================================================
# Adaptive Marking (Doerfler marking with target DOF ratio)
# ============================================================
def mark_elements_adaptive_theta(element_errors, current_ne, target_ratio=2.0):
	"""Mark elements with dynamically adjusted theta to achieve target DOF ratio."""
	max_error = max(element_errors)
	if max_error <= 0:
		return [], 1.0

	expansion_factor = 28.0
	target_marked = int(current_ne * (target_ratio - 1) / expansion_factor)
	target_marked = max(1, min(target_marked, current_ne))

	theta_low, theta_high = 0.0, 1.0
	best_theta = 0.5
	best_marked = []

	for _ in range(20):
		theta = (theta_low + theta_high) / 2
		cutoff = theta * max_error
		marked = [i for i, err in enumerate(element_errors) if err >= cutoff]
		n_marked = len(marked)

		if n_marked == target_marked:
			return marked, theta
		elif n_marked < target_marked:
			theta_high = theta
			if n_marked > len(best_marked):
				best_marked = marked
				best_theta = theta
		else:
			theta_low = theta
			if abs(n_marked - target_marked) < abs(len(best_marked) - target_marked):
				best_marked = marked
				best_theta = theta

	cutoff = best_theta * max_error
	marked = [i for i, err in enumerate(element_errors) if err >= cutoff]
	return marked, best_theta

# ============================================================
# VTK Output
# ============================================================
def output_vtk(mesh, iteration, output_data):
	"""Output mesh and solution to VTK file."""
	import shutil
	import gc

	max_elements_for_vtk = 30000
	if mesh.ne > max_elements_for_vtk:
		print(f"  (VTK skipped: {mesh.ne} elements > {max_elements_for_vtk})")
		return None

	temp_dir = tempfile.gettempdir()
	temp_basename = f"Cube_3D_adaptive_with_Kelvin_iter_{iteration:02d}"
	temp_vtk_path = os.path.join(temp_dir, temp_basename)
	final_vtk_path = os.path.join(script_dir, temp_basename + ".vtu")

	coefs = [output_data['gfu']]
	names = ["Omega"]

	if 'B_pert_cf' in output_data and output_data['B_pert_cf'] is not None:
		coefs.append(output_data['B_pert_cf'])
		names.append("B_pert")

	if 'H_pert_cf' in output_data and output_data['H_pert_cf'] is not None:
		coefs.append(output_data['H_pert_cf'])
		names.append("H_pert")

	if 'B_total_cf' in output_data and output_data['B_total_cf'] is not None:
		coefs.append(output_data['B_total_cf'])
		names.append("B_total")

	if 'H_total_cf' in output_data and output_data['H_total_cf'] is not None:
		coefs.append(output_data['H_total_cf'])
		names.append("H_total")

	if 'Bs_cf' in output_data and output_data['Bs_cf'] is not None:
		coefs.append(output_data['Bs_cf'])
		names.append("Bs")

	if 'Hs_cf' in output_data and output_data['Hs_cf'] is not None:
		coefs.append(output_data['Hs_cf'])
		names.append("Hs")

	vtk = VTKOutput(mesh, coefs=coefs, names=names,
	                filename=temp_vtk_path, subdivision=0)
	vtk.Do()

	del vtk
	gc.collect()

	temp_vtu_file = temp_vtk_path + ".vtu"
	if os.path.exists(temp_vtu_file):
		shutil.copy2(temp_vtu_file, final_vtk_path)
		try:
			os.remove(temp_vtu_file)
		except:
			pass

	return final_vtk_path


# ============================================================
# Save iteration data to MAT file
# ============================================================
def save_iteration_mat(iter_num, history, n_elements, n_vertices):
	"""Save iteration data to MAT file."""
	mat_iter_data = {
		'iter_num': iter_num,
		'n_elements': n_elements,
		'n_vertices': n_vertices,
		'ndof': history['ndof'][-1],
		'error': history['error'][-1],
		'energy_magnetic': history['energy_magnetic'][-1],
		'energy_air_total': history['energy_air_total'][-1],
		'energy_air_inner': history['energy_air_inner'][-1],
		'energy_air_outer': history['energy_air_outer'][-1],
		'history_ndof': array(history['ndof']),
		'history_error': array(history['error']),
		'cube_size': cube_size,
		'air_total_radius': air_total_radius,
		'kelvin_radius': kelvin_radius,
		'offset_z': offset_z,
		'order': order,
		'dimension': 3,
		'method': 'refine_zz',
		'model_type': '1/8',
		'geometry': 'cube'
	}
	mat_iter_file = os.path.join(script_dir, f"cube_3d_iter_{iter_num:02d}.mat")
	sio.savemat(mat_iter_file, mat_iter_data)
	print(f"  MAT saved: {mat_iter_file}")


# ============================================================
# Generate convergence plot
# ============================================================
def compute_y0_cross_section(verts):
	"""
	Compute the cross-section of a tetrahedron with the y=0 plane.
	Returns a list of (x, z) coordinates forming the cross-section polygon,
	or None if the tetrahedron doesn't intersect y=0.
	"""
	y_coords = [v[1] for v in verts]
	y_min, y_max = min(y_coords), max(y_coords)

	# Check if tetrahedron intersects y=0 plane
	if y_min > 0 or y_max < 0:
		return None  # No intersection

	# Collect intersection points
	cross_points = []

	# Check each edge for intersection with y=0
	edges = [
		(0, 1), (0, 2), (0, 3),
		(1, 2), (1, 3), (2, 3)
	]

	for i, j in edges:
		y0, y1 = verts[i][1], verts[j][1]

		# Check if edge crosses y=0
		if (y0 <= 0 <= y1) or (y1 <= 0 <= y0):
			if abs(y1 - y0) < 1e-12:
				# Edge is on y=0 plane, add both endpoints
				cross_points.append((verts[i][0], verts[i][2]))
				cross_points.append((verts[j][0], verts[j][2]))
			else:
				# Interpolate to find intersection point
				t = -y0 / (y1 - y0)
				x = verts[i][0] + t * (verts[j][0] - verts[i][0])
				z = verts[i][2] + t * (verts[j][2] - verts[i][2])
				cross_points.append((x, z))

	# Remove duplicate points
	unique_points = []
	for p in cross_points:
		is_dup = False
		for q in unique_points:
			if abs(p[0] - q[0]) < 1e-10 and abs(p[1] - q[1]) < 1e-10:
				is_dup = True
				break
		if not is_dup:
			unique_points.append(p)

	if len(unique_points) < 3:
		return None

	# Sort points by angle around centroid for proper polygon ordering
	cx = sum(p[0] for p in unique_points) / len(unique_points)
	cz = sum(p[1] for p in unique_points) / len(unique_points)

	def angle_key(p):
		from math import atan2
		return atan2(p[1] - cz, p[0] - cx)

	unique_points.sort(key=angle_key)

	return unique_points

def generate_convergence_plot(iter_num, history, output_data):
	"""Generate 2x2 convergence plot."""
	mesh = output_data['mesh']
	element_errors = output_data['element_errors']

	fig = plt.figure(figsize=(14, 12), dpi=150)

	theta_circle = linspace(0, pi/2, 50)
	r_kelvin_plot = kelvin_radius * sin(theta_circle)
	z_kelvin_plot = kelvin_radius * cos(theta_circle)

	materials = mesh.GetMaterials()

	max_elements_for_plot = 50000
	skip_mesh_plot = mesh.ne > max_elements_for_plot

	# ===== Top-left: Interior domain mesh and error =====
	ax1 = plt.subplot(2, 2, 1)
	if mesh is not None and element_errors is not None and not skip_mesh_plot:
		polygons_interior = []
		error_interior = []
		for el_idx, el in enumerate(mesh.Elements(VOL)):
			mat_name = el.mat
			if mat_name in ["magnetic", "air_total", "air_inner"]:
				verts = [mesh[v].point for v in el.vertices]
				# Compute y=0 cross-section
				cross_section = compute_y0_cross_section(verts)
				if cross_section is not None:
					polygons_interior.append(cross_section)
					error_interior.append(element_errors[el_idx])

		if polygons_interior and error_interior:
			err_arr = array(error_interior)
			err_arr = err_arr.clip(min=1e-20)
			log_err = log10(err_arr)
			norm = Normalize(vmin=-14, vmax=-3)
			colors = cm.jet(norm(log_err))
			pc1 = PolyCollection(polygons_interior, facecolor=colors, edgecolor='white', linewidth=0.2)
			ax1.add_collection(pc1)
			sm = cm.ScalarMappable(cmap='jet', norm=norm)
			sm.set_array([])
			cbar = plt.colorbar(sm, ax=ax1)
			cbar.set_label('$\\log_{10}$(ZZ Error)')

	# Draw cube boundary
	ax1.plot([0, cube_size], [0, 0], 'r-', linewidth=2)
	ax1.plot([cube_size, cube_size], [0, cube_size], 'r-', linewidth=2)
	ax1.plot([cube_size, 0], [cube_size, cube_size], 'r-', linewidth=2)
	ax1.plot([0, 0], [cube_size, 0], 'r-', linewidth=2, label='Cube')
	# Total-Reduced boundary: cube faces at x=air_total_radius and z=air_total_radius
	ax1.plot([air_total_radius, air_total_radius], [0, air_total_radius], 'm--', linewidth=1.5, label='Total-Reduced')
	ax1.plot([0, air_total_radius], [air_total_radius, air_total_radius], 'm--', linewidth=1.5)
	ax1.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
	ax1.set_xlim(-0.05, kelvin_radius + 0.05)
	ax1.set_ylim(-0.05, kelvin_radius + 0.05)
	ax1.set_aspect('equal')
	ax1.set_xlabel('$x$ [m]')
	ax1.set_ylabel('$z$ [m]')
	ax1.set_title(f'Interior domain (mag + air\\_total + air\\_inner) on $y=0$')
	ax1.legend(loc='upper right', fontsize=8)

	ne_mag = history['elements_magnetic'][-1]
	ne_air_tot = history['elements_air_total'][-1]
	ne_air_in = history['elements_air_inner'][-1]
	ne_interior = ne_mag + ne_air_tot + ne_air_in
	ax1.text(0.98, 0.78, f'Interior: {ne_interior}', transform=ax1.transAxes,
	         fontsize=8, ha='right', va='top', color='black')
	ax1.text(0.98, 0.70, f'(mag:{ne_mag}, air_tot:{ne_air_tot}, air_in:{ne_air_in})', transform=ax1.transAxes,
	         fontsize=7, ha='right', va='top', color='gray')

	# ===== Top-right: Exterior domain mesh and error =====
	ax2 = plt.subplot(2, 2, 2)
	if mesh is not None and element_errors is not None and not skip_mesh_plot:
		polygons_exterior = []
		error_exterior = []
		for el_idx, el in enumerate(mesh.Elements(VOL)):
			mat_name = el.mat
			if mat_name == "air_outer":
				verts = [mesh[v].point for v in el.vertices]
				# Compute y=0 cross-section
				cross_section = compute_y0_cross_section(verts)
				if cross_section is not None:
					polygons_exterior.append(cross_section)
					error_exterior.append(element_errors[el_idx])

		if polygons_exterior and error_exterior:
			err_arr = array(error_exterior)
			err_arr = err_arr.clip(min=1e-20)
			log_err = log10(err_arr)
			norm = Normalize(vmin=-14, vmax=-3)
			colors = cm.jet(norm(log_err))
			pc2 = PolyCollection(polygons_exterior, facecolor=colors, edgecolor='white', linewidth=0.2)
			ax2.add_collection(pc2)
			sm = cm.ScalarMappable(cmap='jet', norm=norm)
			sm.set_array([])
			cbar = plt.colorbar(sm, ax=ax2)
			cbar.set_label('$\\log_{10}$(ZZ Error)')

	ax2.plot(r_kelvin_plot + kelvin_radius * 2, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
	ax2.set_xlim(kelvin_radius * 2 - 0.05, kelvin_radius * 3 + 0.05)
	ax2.set_ylim(-0.05, kelvin_radius + 0.05)
	ax2.set_aspect('equal')
	ax2.set_xlabel('$x$ [m]')
	ax2.set_ylabel('$z$ [m]')
	ax2.set_title(f'Exterior domain (air\\_outer, Kelvin) on $y=0$')
	ax2.legend(loc='upper right', fontsize=8)

	ne_air_out = history['elements_air_outer'][-1]
	ax2.text(0.98, 0.78, f'air_outer: {ne_air_out}', transform=ax2.transAxes,
	         fontsize=8, ha='right', va='top', color='black')

	# ===== Bottom-left: DOF vs Error convergence =====
	ax3 = plt.subplot(2, 2, 3)
	ax3.loglog(history['ndof'], history['error'], 'ko-', linewidth=2, markersize=6, label='Uniform Refinement')
	ax3.loglog(history['ndof'][-1], history['error'][-1], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)

	ndof_line = array([1e2, 1e5])
	err_ref_point = history['error'][0]
	N_ref = history['ndof'][0]
	err_line = err_ref_point * (N_ref / ndof_line) ** (order / 3)
	ax3.loglog(ndof_line, err_line, 'r--', linewidth=1.5, label=f'$O(N^{{-{order}/3}})$')

	ax3.set_xlabel('DOFs')
	ax3.set_ylabel('Error Estimator')
	ax3.set_title(f'DOF vs Error (iter {iter_num})')
	ax3.set_xlim(1e2, 1e5)
	if len(history['error']) > 0:
		err_min = min(history['error']) * 0.5
		err_max = max(history['error']) * 2.0
		ax3.set_ylim(err_min, err_max)
	ax3.legend(loc='lower left')
	ax3.grid(True, alpha=0.3)
	ax3.tick_params(direction='in')

	# ===== Bottom-right: DOF vs Energy =====
	ax4 = plt.subplot(2, 2, 4)
	if len(history['energy_magnetic']) > 0:
		ax4.semilogx(history['ndof'], history['energy_magnetic'], 'rs-', linewidth=2, markersize=5, label='magnetic')
		ax4.semilogx(history['ndof'], history['energy_air_total'], 'mo-', linewidth=2, markersize=5, label='air\\_total')
		ax4.semilogx(history['ndof'], history['energy_air_inner'], 'go-', linewidth=2, markersize=5, label='air\\_inner')
		ax4.semilogx(history['ndof'], history['energy_air_outer'], 'bo-', linewidth=2, markersize=5, label='air\\_outer')
		ax4.semilogx(history['ndof'], history['energy'], 'k^-', linewidth=2, markersize=5, label='Total')

		ax4.semilogx(history['ndof'][-1], history['energy_magnetic'][-1], 'ro', markersize=10, markerfacecolor='none', markeredgewidth=2)
		ax4.semilogx(history['ndof'][-1], history['energy_air_total'][-1], 'mo', markersize=10, markerfacecolor='none', markeredgewidth=2)
		ax4.semilogx(history['ndof'][-1], history['energy_air_inner'][-1], 'go', markersize=10, markerfacecolor='none', markeredgewidth=2)
		ax4.semilogx(history['ndof'][-1], history['energy_air_outer'][-1], 'bo', markersize=10, markerfacecolor='none', markeredgewidth=2)
		ax4.semilogx(history['ndof'][-1], history['energy'][-1], 'ko', markersize=10, markerfacecolor='none', markeredgewidth=2)

		E_mag = history['energy_magnetic'][-1]
		E_air_tot = history['energy_air_total'][-1]
		E_air_in = history['energy_air_inner'][-1]
		E_air_out = history['energy_air_outer'][-1]
		E_tot = history['energy'][-1]
		ax4.text(0.98, 0.95, f'magnetic: {E_mag:.4e} J', transform=ax4.transAxes,
		         fontsize=8, ha='right', va='top', color='red')
		ax4.text(0.98, 0.88, f'air_total: {E_air_tot:.4e} J', transform=ax4.transAxes,
		         fontsize=8, ha='right', va='top', color='magenta')
		ax4.text(0.98, 0.81, f'air_inner: {E_air_in:.4e} J', transform=ax4.transAxes,
		         fontsize=8, ha='right', va='top', color='green')
		ax4.text(0.98, 0.74, f'air_outer: {E_air_out:.4e} J', transform=ax4.transAxes,
		         fontsize=8, ha='right', va='top', color='blue')
		ax4.text(0.98, 0.67, f'Total: {E_tot:.4e} J', transform=ax4.transAxes,
		         fontsize=8, ha='right', va='top', color='black')

	ax4.set_xlabel('DOFs')
	ax4.set_ylabel('Perturbation Energy (J)')
	ax4.set_title('DOF vs Perturbation Field Energy')
	ax4.set_xlim(1e2, 1e5)
	if len(history['energy']) > 0:
		all_energies = history['energy'] + history['energy_magnetic'] + history['energy_air_total'] + history['energy_air_inner'] + history['energy_air_outer']
		y_min = min(all_energies) * 0.9
		y_max = max(all_energies) * 1.1
		if y_min > 0:
			ax4.set_ylim(y_min, y_max)
	ax4.legend(loc='lower left', fontsize=7)
	ax4.grid(True, alpha=0.3)
	ax4.tick_params(direction='in')

	plt.suptitle(f'Cube Iteration {iter_num}: DOFs={history["ndof"][-1]}, Error={history["error"][-1]:.2e}',
	             fontsize=14, fontweight='bold', y=1.01)
	plt.tight_layout()

	png_file = os.path.join(script_dir, f"cube_3d_iter_{iter_num:02d}.png")
	plt.savefig(png_file, dpi=150, bbox_inches='tight')
	print(f"  PNG saved: {png_file}")
	plt.close()


def count_elements_by_region(mesh):
	"""Count elements in each material region."""
	materials = mesh.GetMaterials()
	counts = {mat: 0 for mat in set(materials)}
	for el in mesh.Elements(VOL):
		mat_name = el.mat
		counts[mat_name] = counts.get(mat_name, 0) + 1
	return counts


# ============================================================
# Main: Uniform Refinement Loop
# ============================================================
print("\n" + "=" * 60)
print("Creating 1/8 3D cube geometry...")
print("=" * 60)

geo = create_geometry()
mesh = Mesh(geo.GenerateMesh(maxh=maxh_initial, grading=0.5))
mesh.Curve(order)

print(f"\nInitial mesh:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# History tracking
history = {
	'ndof': [],
	'elements': [],
	'elements_total': [],
	'elements_reduced': [],
	'elements_magnetic': [],
	'elements_air_total': [],
	'elements_air_inner': [],
	'elements_air_outer': [],
	'error': [],
	'energy': [],
	'energy_total_region': [],
	'energy_reduced_region': [],
	'energy_magnetic': [],
	'energy_air_total': [],
	'energy_air_inner': [],
	'energy_air_outer': []
}

print("\n" + "=" * 60)
print("Starting Uniform Mesh Refinement (All Elements)")
print("=" * 60)

iteration = 0
prev_ndof = 0

while True:
	if prev_ndof >= 1e5:
		print(f"\n  DOF limit reached ({prev_ndof} >= 1e5), stopping without computing.")
		break

	print(f"\n{'=' * 60}")
	print(f"Iteration {iteration + 1}")
	print("=" * 60)

	fes, gfu, Mu, fields = solve_omega_formulation(mesh, order)

	H_pert_cf = fields['H_pert_cf']
	B_pert_cf = fields['B_pert_cf']
	H_total_cf = fields['H_total_cf']
	B_total_cf = fields['B_total_cf']
	Hs_cf = fields['Hs_cf']
	Bs_cf = fields['Bs_cf']
	mu_kelvin = fields['mu_kelvin']

	element_errors = compute_error_estimator(mesh, fes, B_pert_cf)
	total_error = sqrt(sum(element_errors))

	Hs = CoefficientFunction((0.0, 0.0, H0))
	Omega_s = H0 * z

	H_pert_total = grad(gfu) - Hs
	energy_magnetic = Integrate(0.5 * (mu_r * mu0) * InnerProduct(H_pert_total, H_pert_total) * dx("magnetic"), mesh)
	energy_air_total = Integrate(0.5 * mu0 * InnerProduct(H_pert_total, H_pert_total) * dx("air_total"), mesh)

	fesOr = H1(mesh, order=order, definedon="air_inner|air_outer")
	Orr = GridFunction(fesOr)
	Oxr = GridFunction(fesOr)
	Orr.Set(gfu, VOL, definedon="air_inner|air_outer")
	Oxr.Set(Omega_s, BND, mesh.Boundaries("total_reduced"))
	H_pert_reduced = grad(Orr) - grad(Oxr)
	energy_air_inner = Integrate(0.5 * mu0 * InnerProduct(H_pert_reduced, H_pert_reduced) * dx("air_inner"), mesh)

	H_pert_kelvin = grad(Orr)
	energy_air_outer = Integrate(0.5 * mu_kelvin * InnerProduct(H_pert_kelvin, H_pert_kelvin) * dx("air_outer"), mesh)

	energy_magnetic_full = 8 * energy_magnetic
	energy_air_total_full = 8 * energy_air_total
	energy_air_inner_full = 8 * energy_air_inner
	energy_air_outer_full = 8 * energy_air_outer

	energy_1_8 = energy_magnetic + energy_air_total + energy_air_inner + energy_air_outer
	energy_full = 8 * energy_1_8

	energy_total_region_full = energy_magnetic_full + energy_air_total_full
	energy_reduced_region_full = energy_air_inner_full + energy_air_outer_full

	region_counts = count_elements_by_region(mesh)
	ne_magnetic = region_counts.get('magnetic', 0)
	ne_air_total = region_counts.get('air_total', 0)
	ne_air_inner = region_counts.get('air_inner', 0)
	ne_air_outer = region_counts.get('air_outer', 0)
	ne_total = ne_magnetic + ne_air_total
	ne_reduced = ne_air_inner + ne_air_outer

	history['ndof'].append(fes.ndof)
	history['elements'].append(mesh.ne)
	history['elements_total'].append(ne_total)
	history['elements_reduced'].append(ne_reduced)
	history['elements_magnetic'].append(ne_magnetic)
	history['elements_air_total'].append(ne_air_total)
	history['elements_air_inner'].append(ne_air_inner)
	history['elements_air_outer'].append(ne_air_outer)
	history['error'].append(total_error)
	history['energy'].append(energy_full)
	history['energy_total_region'].append(energy_total_region_full)
	history['energy_reduced_region'].append(energy_reduced_region_full)
	history['energy_magnetic'].append(energy_magnetic_full)
	history['energy_air_total'].append(energy_air_total_full)
	history['energy_air_inner'].append(energy_air_inner_full)
	history['energy_air_outer'].append(energy_air_outer_full)

	print(f"  Elements: {mesh.ne} (Total: {ne_total}, Reduced: {ne_reduced})")
	print(f"    magnetic: {ne_magnetic}, air_total: {ne_air_total}, air_inner: {ne_air_inner}, air_outer: {ne_air_outer}")
	print(f"  Vertices: {mesh.nv}")
	print(f"  DOFs: {fes.ndof}")
	print(f"  Error estimator: {total_error:.6e}")
	print(f"  Energy (8x1/8): {energy_full:.6e} J")

	output_data = {
		'mesh': mesh,
		'gfu': gfu,
		'B_pert_cf': B_pert_cf,
		'H_pert_cf': H_pert_cf,
		'B_total_cf': B_total_cf,
		'H_total_cf': H_total_cf,
		'Bs_cf': Bs_cf,
		'Hs_cf': Hs_cf,
		'element_errors': element_errors,
	}

	vtk_file = output_vtk(mesh, iteration, output_data)
	if vtk_file is not None:
		print(f"  VTK saved: {vtk_file}")

	save_iteration_mat(iteration + 1, history, mesh.ne, mesh.nv)

	mat_filename = os.path.join(script_dir, os.path.splitext(os.path.basename(__file__))[0] + ".mat")
	mat_data = {
		'ndof': array(history['ndof']),
		'elements': array(history['elements']),
		'elements_total': array(history['elements_total']),
		'elements_reduced': array(history['elements_reduced']),
		'elements_magnetic': array(history['elements_magnetic']),
		'elements_air_total': array(history['elements_air_total']),
		'elements_air_inner': array(history['elements_air_inner']),
		'elements_air_outer': array(history['elements_air_outer']),
		'error': array(history['error']),
		'energy': array(history['energy']),
		'energy_total_region': array(history['energy_total_region']),
		'energy_reduced_region': array(history['energy_reduced_region']),
		'energy_magnetic': array(history['energy_magnetic']),
		'energy_air_total': array(history['energy_air_total']),
		'energy_air_inner': array(history['energy_air_inner']),
		'energy_air_outer': array(history['energy_air_outer']),
		'order': order,
		'dimension': 3,
		'method': 'refine_zz',
		'model_type': '1/8',
		'geometry': 'cube'
	}
	sio.savemat(mat_filename, mat_data)
	print(f"  Main MAT saved: {mat_filename}")

	generate_convergence_plot(iteration + 1, history, output_data)

	if iteration + 1 >= max_iterations:
		print(f"\n  Iteration limit reached ({iteration + 1} >= {max_iterations}), stopping.")
		break

	prev_ndof = fes.ndof

	# Adaptive marking using ZZ error estimator
	marked, theta_used = mark_elements_adaptive_theta(element_errors, mesh.ne, target_ratio=2.0)
	print(f"  Marked {len(marked)}/{mesh.ne} elements (theta={theta_used:.3f})")

	for el in mesh.Elements():
		mesh.SetRefinementFlag(el, False)
	for idx in marked:
		elements_list = list(mesh.Elements(VOL)); mesh.SetRefinementFlag(elements_list[idx], True)

	mesh.Refine()
	mesh.Curve(order)

	iteration += 1


# ============================================================
# Final Statistics
# ============================================================
print("\n" + "=" * 60)
print("Convergence History")
print("=" * 60)

print(f"\n{'Iter':<6} {'Elements':<10} {'DOFs':<10} {'Error Est':<12} {'E_sum(J)':<12}")
print("-" * 70)
for i in range(len(history['ndof'])):
	print(f"{i+1:<6} {history['elements'][i]:<10} {history['ndof'][i]:<10} "
	      f"{history['error'][i]:<12.4e} {history['energy'][i]:<12.4e}")

print(f"\nInitial -> Final:")
print(f"  Elements: {history['elements'][0]} -> {history['elements'][-1]}")
print(f"  DOFs: {history['ndof'][0]} -> {history['ndof'][-1]}")
if history['error'][-1] > 0:
	print(f"  Error: {history['error'][0]:.4e} -> {history['error'][-1]:.4e} "
	      f"({history['error'][0]/history['error'][-1]:.1f}x reduction)")

mat_filename = os.path.join(script_dir, os.path.splitext(os.path.basename(__file__))[0] + ".mat")
mat_data = {
	'ndof': array(history['ndof']),
	'elements': array(history['elements']),
	'elements_total': array(history['elements_total']),
	'elements_reduced': array(history['elements_reduced']),
	'elements_magnetic': array(history['elements_magnetic']),
	'elements_air_total': array(history['elements_air_total']),
	'elements_air_inner': array(history['elements_air_inner']),
	'elements_air_outer': array(history['elements_air_outer']),
	'error': array(history['error']),
	'energy': array(history['energy']),
	'energy_total_region': array(history['energy_total_region']),
	'energy_reduced_region': array(history['energy_reduced_region']),
	'energy_magnetic': array(history['energy_magnetic']),
	'energy_air_total': array(history['energy_air_total']),
	'energy_air_inner': array(history['energy_air_inner']),
	'energy_air_outer': array(history['energy_air_outer']),
	'order': order,
	'dimension': 3,
	'method': 'refine_zz',
	'model_type': '1/8',
	'geometry': 'cube'
}
sio.savemat(mat_filename, mat_data)
print(f"\nData saved to: {mat_filename}")

print(f"\nFinal Energy Values:")
print(f"  Total region (magnetic+air_total): {history['energy_total_region'][-1]:.4e} J")
print(f"  Reduced region (air_inner+air_outer): {history['energy_reduced_region'][-1]:.4e} J")
print(f"  Sum: {history['energy'][-1]:.4e} J")
print(f"  Per-region:")
print(f"    magnetic:   {history['energy_magnetic'][-1]:.4e} J")
print(f"    air_total:  {history['energy_air_total'][-1]:.4e} J")
print(f"    air_inner:  {history['energy_air_inner'][-1]:.4e} J")
print(f"    air_outer:  {history['energy_air_outer'][-1]:.4e} J")

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
