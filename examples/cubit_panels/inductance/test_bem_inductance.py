"""
BEM inductance extraction test - energy method with LaplaceSL.

Creates a torus with gap in Cubit, meshes (tet or hex), exports to NGSolve,
runs LaplaceSL BEM to extract self-inductance via energy method:
    L = mu_0 * J^T * SL * J   (for I=1 toroidal current)

Verified against Neumann formula: L = mu_0 * R * (ln(8R/a) - 2).

Usage:
    python test_bem_inductance.py                    # tet mesh, order=1
    python test_bem_inductance.py --mesh hex         # hex sweep mesh, order=1
    python test_bem_inductance.py --mesh tet --order 2  # tet mesh, order=2

Requires: NGSolve (imported FIRST), Cubit, netgen fork with CalcSurfacesOfNode
"""

import argparse
import math
import os
import sys
import time

# NGSolve MUST be imported BEFORE cubit (DLL conflict avoidance)
import numpy as np
from ngsolve import (Mesh, HDivSurface, TaskManager, ds,
                     Integrate, CF, BND, GridFunction, Norm, sqrt, x, y, z)
from ngsolve.bem import LaplaceSL

# Setup paths
radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src', 'radia')
if os.path.abspath(radia_src) not in sys.path:
	sys.path.insert(0, os.path.abspath(radia_src))

cubit_path = os.environ.get("CUBIT_PATH", r"C:\Program Files\Coreform Cubit 2025.3\bin")
if cubit_path not in sys.path:
	sys.path.append(cubit_path)

for p in list(sys.path):
	if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
		sys.path.remove(p)

import cubit
import io, contextlib
with contextlib.redirect_stderr(io.StringIO()):
	cubit.init(["cubit", "-nojournal", "-batch", "-noinit"])

for p in list(sys.path):
	if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
		sys.path.remove(p)

import cubit_mesh_export

MU_0 = 4.0 * math.pi * 1e-7

# ============================================================
# Parameters
# ============================================================
R = 0.05    # Major radius [m]
a = 0.005   # Minor radius [m]
gap = 0.005 # Gap width [m]


def create_and_mesh(mesh_type="tet", interval=8):
	"""Create torus with gap, mesh."""
	cubit.cmd('reset')
	cubit.cmd(f'create torus major radius {R} minor radius {a}')
	cubit.cmd(f'brick x {3*a} y {gap} z {3*a}')
	cubit.cmd(f'move Volume 2 x {R} y 0 z 0 include_merged')
	cubit.cmd('subtract volume 2 from volume 1')

	# Mesh
	vol_ids = list(cubit.get_entities("volume"))
	if mesh_type == "hex":
		cubit.cmd('volume all scheme sweep')
	else:
		cubit.cmd('volume all scheme tetmesh')
	cubit.cmd(f'curve all interval {interval}')
	cubit.cmd('mesh volume all')

	ne_hex = sum(len(cubit.get_volume_hexes(vid)) for vid in vol_ids)
	ne_tet = sum(len(cubit.get_volume_tets(vid)) for vid in vol_ids)
	print(f"  Cubit mesh: {ne_hex} hexes, {ne_tet} tets")

	# Register blocks: volume + surface elements (tri AND quad)
	cubit.cmd("set duplicate block elements on")
	for bid in list(cubit.get_block_id_list()):
		cubit.cmd(f"delete block {bid}")
	for i, vid in enumerate(vol_ids):
		cubit.cmd(f"block {i + 1} add volume {vid}")
	bnd_id = len(vol_ids) + 1
	cubit.cmd(f"block {bnd_id} add tri all")
	cubit.cmd(f"block {bnd_id} add face all")
	cubit.cmd(f'block {bnd_id} name "conductor"')

	n_tri = len(cubit.get_block_tris(bnd_id))
	n_quad = len(cubit.get_block_faces(bnd_id))
	print(f"  Surface elements: {n_tri} tris, {n_quad} quads")

	return vol_ids


def run_bem(order=1):
	"""Export to NGSolve, run LaplaceSL BEM with energy method."""
	t0 = time.perf_counter()
	mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
		cubit, order=order, surface_only=True,
		split_quads=True)
	t_export = time.perf_counter() - t0
	print(f"  export_NGSolveCurvedMesh(order={order}): {t_export:.2f}s")

	print(f"  NGSolve mesh: dim={mesh.dim}, ne={mesh.ne}, nse={mesh.GetNE(BND)}, nv={mesh.nv}")

	# Surface area verification
	area = float(Integrate(CF(1), mesh, VOL_or_BND=BND))
	area_exact = 4 * math.pi**2 * R * a  # Full torus (gap is small)
	print(f"  Surface area: {area*1e4:.4f} cm^2 (exact ~{area_exact*1e4:.4f} cm^2)")

	# HDivSurface
	fes = HDivSurface(mesh, order=0)
	u, v = fes.TnT()
	ndof = fes.ndof
	print(f"  HDivSurface: ndof={ndof}")

	# Toroidal current for I=1: J = e_phi / (2*pi*a)
	# e_phi = (-y, x, 0) / sqrt(x^2 + y^2)
	r_cf = sqrt(x*x + y*y)
	J_toroidal = CF((-y/r_cf, x/r_cf, 0)) / (2 * math.pi * a)

	# Project onto HDivSurface
	gf_J = GridFunction(fes)
	gf_J.Set(J_toroidal, definedon=mesh.Boundaries(".*"), dual=True)
	J_vec = gf_J.vec.FV().NumPy().copy()

	# LaplaceSL operator (use_fmm=False for reproducibility and faster dense extraction)
	# TaskManager must be used for BOTH setup and extraction, or neither.
	# See: Joachim Schoeberl feedback (2026-03-22)
	t0 = time.perf_counter()
	with TaskManager():
		L_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
	t_op = time.perf_counter() - t0
	print(f"  LaplaceSL operator (use_fmm=False): {t_op:.2f}s")

	# Dense SL matrix extraction using ToDense() (optimized for BEM operators)
	print(f"  Extracting {ndof}x{ndof} dense matrix via ToDense()...")
	t0 = time.perf_counter()
	with TaskManager():
		SL = L_op.mat.ToDense().NumPy()
	t_extract = time.perf_counter() - t0
	print(f"  Dense extraction: {t_extract:.1f}s")

	# Energy method: L = mu_0 * J^T * SL * J  (for I=1)
	L_total = MU_0 * float(J_vec @ SL @ J_vec)

	# Neumann formula (external inductance, thin torus R >> a)
	L_neumann = MU_0 * R * (math.log(8 * R / a) - 2)

	# Diagnostics
	sym_err = np.linalg.norm(SL - SL.T) / np.linalg.norm(SL)

	print(f"\n{'='*50}")
	print(f"Results (order={order}):")
	print(f"  L (BEM energy): {L_total*1e9:.4f} nH")
	print(f"  L (Neumann):    {L_neumann*1e9:.4f} nH")
	print(f"  Error:          {(L_total - L_neumann) / L_neumann * 100:.2f}%")
	print(f"  SL symmetry:    ||SL-SL^T||/||SL|| = {sym_err:.2e}")
	print(f"  DOFs:           {ndof}")
	print(f"  Total time:     {t_export + t_op + t_extract:.1f}s")
	print(f"{'='*50}")

	# GMSH export: surface current distribution |J|
	try:
		from gmsh_post_export import GmshPostExport

		elem_J = Integrate(Norm(gf_J), mesh, VOL_or_BND=BND, element_wise=True)
		elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)

		num_nodes = mesh.nv
		node_sum = np.zeros(num_nodes)
		node_cnt = np.zeros(num_nodes)
		for el in mesh.Elements(BND):
			avg_J = abs(elem_J[el.nr]) / max(abs(elem_A[el.nr]), 1e-30)
			for vtx in el.vertices:
				node_sum[vtx.nr] += avg_J
				node_cnt[vtx.nr] += 1
		node_J = np.where(node_cnt > 0, node_sum / node_cnt, 0.0)

		gmsh_file = os.path.join(
			os.path.dirname(os.path.abspath(__file__)),
			"inductance_J.msh"
		)
		post = GmshPostExport(mesh, boundary=True)
		post.add_field("|J|", node_J, ncomp=1)
		post.write(gmsh_file)
		print(f"  GMSH output: {gmsh_file}")
	except Exception as ex:
		print(f"  GMSH export failed: {ex}")

	return {"L_nH": L_total * 1e9, "L_neumann_nH": L_neumann * 1e9,
	        "error_pct": (L_total - L_neumann) / L_neumann * 100,
	        "ndof": ndof, "order": order}


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="BEM inductance test (energy method)")
	parser.add_argument("--mesh", choices=["tet", "hex"], default="tet",
	                    help="Mesh type: tet (tri surface) or hex (quad surface)")
	parser.add_argument("--order", type=int, default=1, help="Curve order (1-5)")
	parser.add_argument("--interval", type=int, default=8, help="Curve interval for hex mesh")
	args = parser.parse_args()

	print(f"BEM Inductance Test (energy method): mesh={args.mesh}, order={args.order}")
	print(f"Torus: R={R*1e3:.1f}mm, a={a*1e3:.1f}mm, gap={gap*1e3:.1f}mm")

	vol_ids = create_and_mesh(args.mesh, args.interval)
	result = run_bem(args.order)
