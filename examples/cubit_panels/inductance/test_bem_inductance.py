"""
BEM inductance extraction test - Hodge decomposition with LaplaceSL.

Creates a torus in Cubit, exports surface mesh to NGSolve,
solves EFIE via Hodge decomposition to extract self-inductance.

Method:
  1. Hodge decomposition: harmonic = ker([D; C^T M_J])
  2. LaplaceSL projected onto harmonic subspace
  3. Generalized eigenvalue -> toroidal mode -> L
  4. GMSH v2.2 export with |J| field

Usage:
    python test_bem_inductance.py
    python test_bem_inductance.py --order 2

Requires: NGSolve (imported FIRST), Cubit, scipy, netgen fork
"""

import argparse
import math
import os
import sys
import time

import numpy as np
from scipy.linalg import null_space, eigh
from ngsolve import *
from ngsolve.bem import LaplaceSL

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
R = 0.05
a = 0.005


def extract_rect(mat, nrow, ncol):
	M = np.zeros((nrow, ncol))
	ei = mat.CreateRowVector(); col = mat.CreateColVector()
	for j in range(ncol):
		ei[:] = 0; ei[j] = 1.0; mat.Mult(ei, col)
		for i in range(nrow): M[i, j] = col[i]
	return M


def extract_dense(mat, n):
	M = np.zeros((n, n))
	ei = mat.CreateColVector(); col = mat.CreateColVector()
	for j in range(n):
		ei[:] = 0; ei[j] = 1.0; mat.Mult(ei, col)
		for i in range(n): M[i, j] = col[i]
	return M


def create_and_mesh():
	"""Create full torus in Cubit, tet mesh."""
	cubit.cmd('reset')
	cubit.cmd(f'create torus major radius {R} minor radius {a}')
	vid = list(cubit.get_entities("volume"))[0]
	cubit.cmd(f'volume {vid} scheme tetmesh')
	cubit.cmd(f'mesh volume {vid}')
	cubit.cmd("set duplicate block elements on")
	for bid in list(cubit.get_block_id_list()):
		cubit.cmd(f"delete block {bid}")
	cubit.cmd(f'block 1 add volume {vid}')
	cubit.cmd(f'block 2 add tri all')
	cubit.cmd(f'block 2 name "conductor"')
	ne = len(cubit.get_volume_tets(vid))
	n_tri = len(cubit.get_entities("tri"))
	print(f"  Cubit: {ne} tets, {n_tri} surface tris")


def run_hodge(order=2):
	"""Export mesh, solve via Hodge decomposition, export GMSH."""
	t0 = time.perf_counter()
	mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
		cubit, order=order, surface_only=True, split_quads=True)
	t_export = time.perf_counter() - t0

	nse = mesh.GetNE(BND)
	nv = mesh.nv
	ne = sum(1 for e in mesh.edges)
	euler = nv - ne + nse
	area = float(Integrate(CF(1), mesh, VOL_or_BND=BND))
	area_exact = 4 * math.pi**2 * R * a

	print(f"  Mesh: nse={nse}, nv={nv}, ne={ne}, Euler={euler}")
	print(f"  Area: {area*1e4:.2f} cm^2 (exact={area_exact*1e4:.2f}, err={(area-area_exact)/area_exact*100:+.2f}%)")
	print(f"  export: {t_export:.1f}s")

	# Hodge decomposition
	fes_J = HDivSurface(mesh, order=0)
	fes_L2 = SurfaceL2(mesh, order=0)
	n_J, n_f = fes_J.ndof, fes_L2.ndof

	u_J, q = fes_J.TrialFunction(), fes_L2.TestFunction()
	bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
	bf_D += div(u_J.Trace()) * q * ds
	bf_D.Assemble()
	D = extract_rect(bf_D.mat, n_f, n_J)

	C = np.zeros((n_J, nv))
	for i, e in enumerate(mesh.edges):
		vv = list(e.vertices)
		if i < n_J:
			C[i, vv[0].nr] = -1
			C[i, vv[1].nr] = +1

	u, v = fes_J.TnT()
	bf_M = BilinearForm(fes_J)
	bf_M += InnerProduct(u.Trace(), v.Trace()) * ds
	bf_M.Assemble()
	M_J = extract_dense(bf_M.mat, n_J)

	c_h_all = null_space(np.vstack([D, C.T @ M_J]), rcond=1e-10)
	n_harm = c_h_all.shape[1]
	print(f"  Hodge: n_J={n_J}, n_harm={n_harm}")

	# LaplaceSL
	t0 = time.perf_counter()
	jt, jv = fes_J.TnT()
	with TaskManager():
		V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
		SL = V_op.mat.ToDense().NumPy()
	t_sl = time.perf_counter() - t0
	print(f"  LaplaceSL + ToDense: {t_sl:.1f}s")

	# Eigenvalue on harmonic subspace
	SL_harm = c_h_all.T @ SL @ c_h_all
	M_harm = c_h_all.T @ M_J @ c_h_all
	eigvals_h, eigvecs_h = eigh(SL_harm, M_harm)

	L_neumann = MU_0 * R * (math.log(8 * R / a) - 2)

	print(f"\n{'='*50}")
	print(f"Results (order={order}):")
	for k in range(min(n_harm, 3)):
		L_k = MU_0 * eigvals_h[k] * R / a
		err = (L_k - L_neumann) / L_neumann * 100
		label = "poloidal" if k == 0 else "TOROIDAL" if k == 1 else f"mode{k}"
		print(f"  {label:10s}: L = {L_k*1e9:.4f} nH  ({err:+.2f}%)")
	print(f"  Neumann:    L = {L_neumann*1e9:.4f} nH")
	print(f"{'='*50}")

	# GMSH export: toroidal current distribution
	if n_harm >= 2:
		c_tor = c_h_all @ eigvecs_h[:, 1]
		L_tor = MU_0 * eigvals_h[1] * R / a
	else:
		c_tor = c_h_all @ eigvecs_h[:, 0]
		L_tor = MU_0 * eigvals_h[0] * R / a

	gf_J = GridFunction(fes_J)
	gf_J.vec.FV().NumPy()[:] = c_tor

	elem_J = Integrate(Norm(gf_J), mesh, VOL_or_BND=BND, element_wise=True)
	elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
	ns, nc = np.zeros(nv), np.zeros(nv)
	for el in mesh.Elements(BND):
		aj = abs(elem_J[el.nr]) / max(abs(elem_A[el.nr]), 1e-30)
		for vtx in el.vertices:
			ns[vtx.nr] += aj
			nc[vtx.nr] += 1
	node_J = np.where(nc > 0, ns / nc, 0.0)

	try:
		from gmsh_post_export import GmshPostExport
		gmsh_file = os.path.join(
			os.path.dirname(os.path.abspath(__file__)),
			f"inductance_J_p{order}.msh")
		post = GmshPostExport(mesh, boundary=True)
		post.add_field("|J_toroidal|", node_J, ncomp=1)
		post.write(gmsh_file)
	except Exception as ex:
		print(f"  GMSH export failed: {ex}")

	return {"L_nH": L_tor * 1e9, "L_neumann_nH": L_neumann * 1e9,
	        "n_J": n_J, "n_harm": n_harm, "order": order}


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="BEM inductance (Hodge decomposition)")
	parser.add_argument("--order", type=int, default=2, help="Curve order (1-2)")
	args = parser.parse_args()

	print(f"BEM Inductance (Hodge): R={R*1e3:.0f}mm, a={a*1e3:.0f}mm, order={args.order}")
	create_and_mesh()
	result = run_hodge(args.order)
