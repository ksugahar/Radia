"""
3D Omega-Reduced Omega Method for Magnetostatics with Kelvin Transformation
WITH MAXH-BASED MESH REFINEMENT (Remeshing with decreasing maxh)
1/8 MODEL (Octant: x>=0, y>=0, z>=0)

Problem: Magnetic sphere (mu_r=100) in uniform z-directed background field

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

Maxh-based algorithm:
  1. Generate mesh with current maxh
  2. Solve and compute error estimator
  3. Reduce maxh by factor ~0.794 (to double DOF in 3D)
  4. Regenerate mesh and repeat

Stop condition: DOF >= max_dof
"""
import os
import glob

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
from ngsolve import TaskManager
from netgen.occ import *
from radia.kelvin_source import kelvin_mu_factor_3d_cf, build_material_cf
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
print("MAXH-BASED REFINEMENT - 1/8 Model (Octant) - SPHERE")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
magnetic_radius = 0.25   # Magnetic sphere radius [m]
air_total_radius = 0.5   # Total region air sphere radius [m] (surrounds magnetic)
kelvin_radius = 1.0      # Kelvin transformation radius [m]
mu_r = 100               # Relative permeability
mu0 = 4 * pi * 1e-7      # Vacuum permeability [H/m]

# Source field: H_s = (0, 0, H0) uniform in z-direction
H0 = 1000.0  # [A/m]

# Mesh parameters
maxh_initial = 0.3       # Initial mesh size (same as other methods)
order = 4                # Finite element order

# Offset for exterior domain (z-direction)
offset_z = 3.0

# Maxh refinement parameters
max_iterations = 20      # Stop after 20 iterations
maxh_factor = 0.794      # Factor to reduce maxh each iteration (2^(-1/3) for ~2x DOF in 3D)

print(f"\nProblem parameters:")
print(f"  Magnetic sphere radius: {magnetic_radius} m")
print(f"  Air_total radius: {air_total_radius} m (Total region air)")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")
print(f"  Model type: 1/8 (octant x>=0, y>=0, z>=0)")
print(f"\nRegion structure:")
print(f"  magnetic + air_total: Total potential (H = grad(Omega))")
print(f"  air_inner + air_outer: Reduced potential (H = grad(Omega) + Hs)")
print(f"\nMaxh-based mesh parameters:")
print(f"  Initial mesh size: {maxh_initial} m")
print(f"  Polynomial order: {order}")
print(f"  Stop condition: {max_iterations} iterations or DOF >= 1e5")
print(f"  Maxh reduction factor: {maxh_factor} per iteration")


# ============================================================
# Geometry Definition for 3D - 1/8 Model (Sphere)
# 4-region structure:
#   magnetic: sphere (Total potential)
#   air_total: sphere around magnetic (Total potential)
#   air_inner: between air_total and Kelvin boundary (Reduced potential)
#   air_outer: Kelvin-transformed exterior (Reduced potential)
# ============================================================
def create_geometry():
	"""Create 1/8 geometry with 4 regions and periodic boundary conditions."""
	print("\nCreating 1/8 sphere geometry with 4-region structure...")
	print(f"  magnetic: sphere (r={magnetic_radius}m)")
	print(f"  air_total: sphere (r={air_total_radius}m) around magnetic")
	print(f"  air_inner: shell ({air_total_radius}m to {kelvin_radius}m)")
	print(f"  air_outer: Kelvin-transformed exterior")

	# Create cutting boxes for 1/8 symmetry
	z_max_cut = offset_z + kelvin_radius * 2
	cut_x = Box(Pnt(-kelvin_radius*2, -kelvin_radius*2, -kelvin_radius*2),
	            Pnt(0, kelvin_radius*2, z_max_cut))
	cut_y = Box(Pnt(-kelvin_radius*2, -kelvin_radius*2, -kelvin_radius*2),
	            Pnt(kelvin_radius*2, 0, z_max_cut))
	cut_z = Box(Pnt(-kelvin_radius*2, -kelvin_radius*2, -kelvin_radius*2),
	            Pnt(kelvin_radius*2, kelvin_radius*2, 0))

	# ===== Region 1: Magnetic sphere (Total potential) =====
	mag_sphere_full = Sphere(Pnt(0, 0, 0), magnetic_radius)
	mag_sphere = mag_sphere_full - cut_x - cut_y - cut_z
	mag_sphere.mat("magnetic")

	for face in mag_sphere.faces:
		fc = face.center
		dist = sqrt(fc.x**2 + fc.y**2 + fc.z**2)
		if abs(fc.x) < 0.05:
			face.name = "sym_x"
		elif abs(fc.y) < 0.05:
			face.name = "sym_y"
		elif abs(fc.z) < 0.05:
			face.name = "sym_z"
		elif dist > magnetic_radius * 0.8:
			face.name = "mag_air_total"

	# ===== Region 2: air_total sphere (Total potential) =====
	air_total_sphere_full = Sphere(Pnt(0, 0, 0), air_total_radius)
	air_total_sphere = air_total_sphere_full - cut_x - cut_y - cut_z
	air_total = air_total_sphere - mag_sphere_full
	air_total.mat("air_total")

	for face in air_total.faces:
		fc = face.center
		dist = sqrt(fc.x**2 + fc.y**2 + fc.z**2)
		if dist > air_total_radius * 0.8:
			face.name = "total_reduced"
		elif abs(fc.x) < 0.05:
			face.name = "sym_x"
		elif abs(fc.y) < 0.05:
			face.name = "sym_y"
		elif abs(fc.z) < 0.05:
			face.name = "sym_z"
		else:
			face.name = "mag_air_total"

	# ===== Region 3: air_inner shell (Reduced potential) =====
	kelvin_sphere_full = Sphere(Pnt(0, 0, 0), kelvin_radius)
	kelvin_sphere = kelvin_sphere_full - cut_x - cut_y - cut_z
	air_inner = kelvin_sphere - air_total_sphere
	air_inner.mat("air_inner")

	for face in air_inner.faces:
		fc = face.center
		dist = sqrt(fc.x**2 + fc.y**2 + fc.z**2)
		if abs(dist - kelvin_radius) < kelvin_radius * 0.2:
			face.name = "kelvin_int"
		elif dist < air_total_radius * 1.2 and dist > air_total_radius * 0.5:
			face.name = "total_reduced"
		elif abs(fc.x) < 0.05:
			face.name = "sym_x"
		elif abs(fc.y) < 0.05:
			face.name = "sym_y"
		elif abs(fc.z) < 0.05:
			face.name = "sym_z"

	# ===== Region 4: air_outer (Kelvin-transformed exterior) =====
	outer_sphere_full = Sphere(Pnt(0, 0, offset_z), kelvin_radius)
	half_x_ext = HalfSpace(Pnt(0, 0, 0), Vec(-1, 0, 0))
	half_y_ext = HalfSpace(Pnt(0, 0, 0), Vec(0, -1, 0))
	half_z_ext = HalfSpace(Pnt(0, 0, offset_z), Vec(0, 0, -1))
	outer_sphere = outer_sphere_full * half_x_ext * half_y_ext * half_z_ext
	outer_sphere.mat("air_outer")

	for face in outer_sphere.faces:
		fc = face.center
		if abs(fc.x) < 0.1:
			face.name = "sym_x"
		elif abs(fc.y) < 0.1:
			face.name = "sym_y"
		elif abs(fc.z - offset_z) < 0.1:
			face.name = "sym_z_ext"
		else:
			face.name = "kelvin_ext"

	vertex = Vertex(Pnt(0, 0, offset_z))
	vertex.name = "GND"

	geo = Glue([mag_sphere, air_total, air_inner, outer_sphere, vertex])

	for i, solid in enumerate(geo.solids):
		if i == 0:
			solid.name = "magnetic"
		elif i == 1:
			solid.name = "air_total"
		elif i == 2:
			solid.name = "air_inner"
		elif i == 3:
			solid.name = "air_outer"

	print("Identifying periodic boundaries...")

	kelvin_int_faces = []
	kelvin_ext_faces = []

	for solid in geo.solids:
		for face in solid.faces:
			if face.name == "kelvin_int":
				kelvin_int_faces.append(face)
				print(f"  Found kelvin_int face in solid '{solid.name}'")
			elif face.name == "kelvin_ext":
				kelvin_ext_faces.append(face)
				print(f"  Found kelvin_ext face in solid '{solid.name}'")

	print(f"  Total kelvin_int faces: {len(kelvin_int_faces)}")
	print(f"  Total kelvin_ext faces: {len(kelvin_ext_faces)}")

	if len(kelvin_int_faces) > 0 and len(kelvin_ext_faces) > 0:
		for int_face in kelvin_int_faces:
			for ext_face in kelvin_ext_faces:
				int_face.Identify(ext_face, "periodic", IdentificationType.PERIODIC)
		print(f"  Periodic identification applied: {len(kelvin_int_faces)} x {len(kelvin_ext_faces)} face pairs")
	else:
		print("  WARNING: Could not find periodic faces!")

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

	r_prime_sq = x**2 + y**2 + (z - offset_z)**2
	r_prime = sqrt(r_prime_sq + 1e-20)

	# Kelvin material modulation via centralized helper (Nagamine CEFC 2026).
	mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=(0.0, 0.0, offset_z),
	                                           R=kelvin_radius)
	mu_kelvin = mu0 * mu_kelvin_factor  # alias for downstream energy integration
	Mu = build_material_cf(
		mesh, mu0, mu_kelvin_factor,
		outer_keyword="air_outer",
		overrides={"magnetic": mu_r * mu0},
	)

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
# ZZ Error Estimator
# ============================================================
def compute_error_estimator(mesh, fes, H_pert_cf, mu0, mu_r):
	"""Compute ZZ-type error estimator using H(div) flux recovery."""
	B_bounded_dict = {
		"magnetic": (mu_r * mu0) * H_pert_cf,
		"air_total": mu0 * H_pert_cf,
		"air_inner": mu0 * H_pert_cf,
		"air_outer": mu0 * H_pert_cf
	}
	flux = CoefficientFunction([B_bounded_dict[mat] for mat in mesh.GetMaterials()])

	recovery_order = max(1, fes.globalorder - 1)
	fes_flux = HDiv(mesh, order=recovery_order)
	gf_flux = GridFunction(fes_flux)

	sigma = fes_flux.TrialFunction()
	tau = fes_flux.TestFunction()
	a_flux = BilinearForm(fes_flux)
	a_flux += InnerProduct(sigma, tau) * dx
	a_flux.Assemble()

	f_flux = LinearForm(fes_flux)
	f_flux += InnerProduct(flux, tau) * dx
	f_flux.Assemble()

	gf_flux.vec.data = a_flux.mat.Inverse(fes_flux.FreeDofs(), inverse="sparsecholesky") * f_flux.vec

	err = InnerProduct(flux - gf_flux, flux - gf_flux)
	element_errors = Integrate(err, mesh, element_wise=True)
	return element_errors


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
	temp_basename = f"Sphere_3D_maxh_with_Kelvin_iter_{iteration:02d}"
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
def save_iteration_mat(iter_num, history, n_elements, n_vertices, current_maxh):
	"""Save iteration data to MAT file."""
	mat_iter_data = {
		'iter_num': iter_num,
		'n_elements': n_elements,
		'n_vertices': n_vertices,
		'ndof': history['ndof'][-1],
		'error': history['error'][-1],
		'maxh': current_maxh,
		'energy_magnetic': history['energy_magnetic'][-1],
		'energy_air_total': history['energy_air_total'][-1],
		'energy_air_inner': history['energy_air_inner'][-1],
		'energy_air_outer': history['energy_air_outer'][-1],
		'history_ndof': array(history['ndof']),
		'history_error': array(history['error']),
		'history_maxh': array(history['maxh']),
		'magnetic_radius': magnetic_radius,
		'air_total_radius': air_total_radius,
		'kelvin_radius': kelvin_radius,
		'offset_z': offset_z,
		'order': order,
		'dimension': 3,
		'method': 'maxh',
		'model_type': '1/8',
		'geometry': 'sphere'
	}
	mat_iter_file = os.path.join(script_dir, f"sphere_3d_iter_{iter_num:02d}.mat")
	sio.savemat(mat_iter_file, mat_iter_data)
	print(f"  MAT saved: {mat_iter_file}")


# ============================================================
# Generate convergence plot
# ============================================================
def generate_convergence_plot(iter_num, history, output_data, current_maxh):
	"""Generate 2x2 convergence plot."""
	mesh = output_data['mesh']
	element_errors = output_data['element_errors']

	fig = plt.figure(figsize=(14, 12), dpi=150)

	theta_circle = linspace(0, pi/2, 50)
	r_kelvin_plot = kelvin_radius * sin(theta_circle)
	z_kelvin_plot = kelvin_radius * cos(theta_circle)

	materials = mesh.GetMaterials()
	y_tol = current_maxh * 0.1

	max_elements_for_plot = 50000
	skip_mesh_plot = mesh.ne > max_elements_for_plot

	# ===== Top-left: Interior domain mesh and error =====
	ax1 = plt.subplot(2, 2, 1)
	if mesh is not None and element_errors is not None and not skip_mesh_plot:
		triangles_interior = []
		error_interior = []
		for el_idx, el in enumerate(mesh.Elements(VOL)):
			mat_name = el.mat
			if mat_name in ["magnetic", "air_total", "air_inner"]:
				verts = [mesh[v].point for v in el.vertices]
				y_coords = [v[1] for v in verts]
				if min(y_coords) < y_tol:
					face_verts = [v for v in verts if abs(v[1]) < y_tol]
					if len(face_verts) >= 3:
						xz_coords = [(v[0], v[2]) for v in face_verts[:3]]
						triangles_interior.append(xz_coords)
						error_interior.append(element_errors[el_idx])

		if triangles_interior and error_interior:
			err_arr = array(error_interior)
			err_arr = err_arr.clip(min=1e-20)
			log_err = log10(err_arr)
			norm = Normalize(vmin=-14, vmax=-3)
			colors = cm.jet(norm(log_err))
			pc1 = PolyCollection(triangles_interior, facecolor=colors, edgecolor='white', linewidth=0.2)
			ax1.add_collection(pc1)
			sm = cm.ScalarMappable(cmap='jet', norm=norm)
			sm.set_array([])
			cbar = plt.colorbar(sm, ax=ax1)
			cbar.set_label('$\\log_{10}$(ZZ Error)')

	r_mag_plot = magnetic_radius * sin(theta_circle)
	z_mag_plot = magnetic_radius * cos(theta_circle)
	ax1.plot(r_mag_plot, z_mag_plot, 'r-', linewidth=2, label='Magnetic')
	r_air_total_plot = air_total_radius * sin(theta_circle)
	z_air_total_plot = air_total_radius * cos(theta_circle)
	ax1.plot(r_air_total_plot, z_air_total_plot, 'm--', linewidth=1.5, label='Total-Reduced')
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
		triangles_exterior = []
		error_exterior = []
		for el_idx, el in enumerate(mesh.Elements(VOL)):
			mat_name = el.mat
			if mat_name == "air_outer":
				verts = [mesh[v].point for v in el.vertices]
				y_coords = [v[1] for v in verts]
				if min(y_coords) < y_tol:
					face_verts = [v for v in verts if abs(v[1]) < y_tol]
					if len(face_verts) >= 3:
						xz_coords = [(v[0], v[2]) for v in face_verts[:3]]
						triangles_exterior.append(xz_coords)
						error_exterior.append(element_errors[el_idx])

		if triangles_exterior and error_exterior:
			err_arr = array(error_exterior)
			err_arr = err_arr.clip(min=1e-20)
			log_err = log10(err_arr)
			norm = Normalize(vmin=-14, vmax=-3)
			colors = cm.jet(norm(log_err))
			pc2 = PolyCollection(triangles_exterior, facecolor=colors, edgecolor='white', linewidth=0.2)
			ax2.add_collection(pc2)
			sm = cm.ScalarMappable(cmap='jet', norm=norm)
			sm.set_array([])
			cbar = plt.colorbar(sm, ax=ax2)
			cbar.set_label('$\\log_{10}$(ZZ Error)')

	ax2.plot(r_kelvin_plot, z_kelvin_plot + offset_z, 'g--', linewidth=1.5, label='Kelvin boundary')
	ax2.set_xlim(-0.05, kelvin_radius + 0.05)
	ax2.set_ylim(offset_z - 0.05, offset_z + kelvin_radius + 0.05)
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
	ax3.loglog(history['ndof'], history['error'], 'ko-', linewidth=2, markersize=6, label='Maxh-based')
	ax3.loglog(history['ndof'][-1], history['error'][-1], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)

	ndof_line = array([1e2, 1e5])
	err_ref_point = history['error'][0]
	N_ref = history['ndof'][0]
	err_line = err_ref_point * (N_ref / ndof_line) ** (order / 3)
	ax3.loglog(ndof_line, err_line, 'r--', linewidth=1.5, label=f'$O(N^{{-{order}/3}})$')

	ax3.set_xlabel('DOFs')
	ax3.set_ylabel('Error Estimator')
	ax3.set_title(f'DOF vs Error (iter {iter_num}, maxh={current_maxh:.4f})')
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

	plt.suptitle(f'Sphere Iteration {iter_num}: DOFs={history["ndof"][-1]}, maxh={current_maxh:.4f}, Error={history["error"][-1]:.2e}',
	             fontsize=14, fontweight='bold', y=1.01)
	plt.tight_layout()

	png_file = os.path.join(script_dir, f"sphere_3d_iter_{iter_num:02d}.png")
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
# Main: Maxh-based Refinement Loop
# ============================================================
print("\n" + "=" * 60)
print("Creating 1/8 3D sphere geometry...")
print("=" * 60)

geo = create_geometry()

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
	'energy_air_outer': [],
	'maxh': []
}

print("\n" + "=" * 60)
print("Starting Maxh-based Mesh Refinement")
print("=" * 60)

iteration = 0
current_maxh = maxh_initial
prev_ndof = 0

with TaskManager():
	while True:
		if prev_ndof >= 1e5:
			print(f"\n  DOF limit reached ({prev_ndof} >= 1e5), stopping without computing.")
			break

		print(f"\n{'=' * 60}")
		print(f"Iteration {iteration + 1} (maxh = {current_maxh:.6f})")
		print("=" * 60)

		# Generate new mesh with current maxh
		mesh = Mesh(geo.GenerateMesh(maxh=current_maxh, grading=0.5))
		mesh.Curve(order)

		print(f"  Mesh generated with maxh={current_maxh:.6f}")
		print(f"  Elements: {mesh.ne}")
		print(f"  Vertices: {mesh.nv}")

		# Solve
		fes, gfu, Mu, fields = solve_omega_formulation(mesh, order)

		H_pert_cf = fields['H_pert_cf']
		B_pert_cf = fields['B_pert_cf']
		H_total_cf = fields['H_total_cf']
		B_total_cf = fields['B_total_cf']
		Hs_cf = fields['Hs_cf']
		Bs_cf = fields['Bs_cf']
		mu_kelvin = fields['mu_kelvin']

		# Compute error estimator
		element_errors = compute_error_estimator(mesh, fes, H_pert_cf, mu0, mu_r)
		total_error = sqrt(sum(element_errors))

		# Energy calculation
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

		# Record history
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
		history['maxh'].append(current_maxh)

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

		save_iteration_mat(iteration + 1, history, mesh.ne, mesh.nv, current_maxh)

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
			'maxh': array(history['maxh']),
			'order': order,
			'dimension': 3,
			'method': 'maxh',
			'model_type': '1/8',
			'geometry': 'sphere'
		}
		sio.savemat(mat_filename, mat_data)
		print(f"  Main MAT saved: {mat_filename}")

		generate_convergence_plot(iteration + 1, history, output_data, current_maxh)

		if iteration + 1 >= max_iterations:
			print(f"\n  Iteration limit reached ({iteration + 1} >= {max_iterations}), stopping.")
			break

		prev_ndof = fes.ndof

		# Reduce maxh for next iteration
		current_maxh *= maxh_factor

		iteration += 1


	# ============================================================
	# Final Statistics
	# ============================================================
	print("\n" + "=" * 60)
	print("Convergence History")
	print("=" * 60)

	print(f"\n{'Iter':<6} {'maxh':<10} {'Elements':<10} {'DOFs':<10} {'Error Est':<12} {'E_sum(J)':<12}")
	print("-" * 70)
	for i in range(len(history['ndof'])):
		print(f"{i+1:<6} {history['maxh'][i]:<10.6f} {history['elements'][i]:<10} {history['ndof'][i]:<10} "
		      f"{history['error'][i]:<12.4e} {history['energy'][i]:<12.4e}")

	print(f"\nInitial -> Final:")
	print(f"  maxh: {history['maxh'][0]:.6f} -> {history['maxh'][-1]:.6f}")
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
		'maxh': array(history['maxh']),
		'order': order,
		'dimension': 3,
		'method': 'maxh',
		'model_type': '1/8',
		'geometry': 'sphere'
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
