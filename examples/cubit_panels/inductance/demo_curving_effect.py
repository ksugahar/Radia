"""
Demonstrate effect of mesh curving (p=1 vs p=2) on BEM inductance accuracy.

Creates both tri (tet mesh) and quad (hex sweep) surface meshes from Cubit,
then compares p=1 (flat) vs p=2 (curved) for each.

Key finding: p=2 dramatically improves geometric accuracy (surface area),
which translates to improved BEM inductance accuracy on coarse meshes.

Usage:
    python demo_curving_effect.py

Requires: NGSolve (imported FIRST), Cubit, netgen fork with CalcSurfacesOfNode
"""

import math
import os
import sys
import time

# NGSolve MUST be imported BEFORE cubit
import numpy as np
from ngsolve import (Mesh, HDivSurface, TaskManager, ds,
                     Integrate, CF, BND, GridFunction, Norm, sqrt, x, y, z)
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
R = 0.05    # Major radius [m]
a = 0.005   # Minor radius [m]
gap = 0.005 # Gap width [m]
L_neumann = MU_0 * R * (math.log(8 * R / a) - 2)
area_exact = 4 * math.pi**2 * R * a


def create_tet_torus(interval):
	"""Create torus with gap, tet mesh (tri surface)."""
	cubit.cmd('reset')
	cubit.cmd(f'create torus major radius {R} minor radius {a}')
	cubit.cmd(f'brick x {3*a} y {gap} z {3*a}')
	cubit.cmd(f'move Volume 2 x {R} y 0 z 0 include_merged')
	cubit.cmd('subtract volume 2 from volume 1')

	vid = list(cubit.get_entities("volume"))[0]
	cubit.cmd(f'volume {vid} scheme tetmesh')
	cubit.cmd(f'curve all interval {interval}')
	cubit.cmd(f'mesh volume {vid}')

	cubit.cmd("set duplicate block elements on")
	for bid in list(cubit.get_block_id_list()):
		cubit.cmd(f"delete block {bid}")
	cubit.cmd(f'block 1 add volume {vid}')
	cubit.cmd(f'block 2 add tri all')
	cubit.cmd(f'block 2 name "conductor"')

	ne = len(cubit.get_volume_tets(vid))
	n_tri = len(cubit.get_entities("tri"))
	return ne, n_tri, "tri"


def create_hex_torus(interval_src=2, interval_sweep=12):
	"""Create torus by sweep, hex mesh (quad surface)."""
	cubit.cmd('reset')
	cubit.cmd(f'create surface circle radius {a} yplane')
	cubit.cmd(f'move surface 1 x {R}')
	cubit.cmd(f'sweep surface 1 zaxis angle 350')

	vid = list(cubit.get_entities("volume"))[0]
	surfaces = cubit.get_relatives("volume", vid, "surface")
	planar = [sid for sid in surfaces if cubit.surface(sid).is_planar()]
	cubit.cmd(f'volume {vid} scheme sweep source surface {planar[0]} target surface {planar[1]}')

	# Source face mesh density
	cubit.cmd(f'surface {planar[0]} size {a / interval_src}')
	cubit.cmd(f'mesh surface {planar[0]}')

	# Sweep direction intervals
	for sid in surfaces:
		if not cubit.surface(sid).is_planar():
			curves = cubit.get_relatives("surface", sid, "curve")
			for cid in curves:
				length = cubit.get_curve_length(cid)
				if length > 0.01:
					cubit.cmd(f'curve {cid} interval {interval_sweep}')

	cubit.cmd(f'mesh volume {vid}')

	cubit.cmd("set duplicate block elements on")
	for bid in list(cubit.get_block_id_list()):
		cubit.cmd(f"delete block {bid}")
	cubit.cmd(f'block 1 add volume {vid}')
	cubit.cmd(f'block 2 add face all')
	cubit.cmd(f'block 2 name "conductor"')

	ne = len(cubit.get_volume_hexes(vid))
	n_quad = len(cubit.get_entities("quad"))
	return ne, n_quad, "quad"


def run_bem(order):
	"""Run BEM inductance extraction at given curve order."""
	mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
		cubit, order=order, surface_only=True, split_quads=True)

	area = float(Integrate(CF(1), mesh, VOL_or_BND=BND))
	area_err = (area - area_exact) / area_exact * 100

	fes = HDivSurface(mesh, order=0)
	u, v = fes.TnT()
	ndof = fes.ndof

	r_cf = sqrt(x * x + y * y)
	J_tor = CF((-y / r_cf, x / r_cf, 0)) / (2 * math.pi * a)
	gf_J = GridFunction(fes)
	gf_J.Set(J_tor, definedon=mesh.Boundaries(".*"), dual=True)
	J_vec = gf_J.vec.FV().NumPy().copy()

	t0 = time.perf_counter()
	with TaskManager():
		L_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
		SL = L_op.mat.ToDense().NumPy()
	t_bem = time.perf_counter() - t0

	L_total = MU_0 * float(J_vec @ SL @ J_vec)
	L_err = (L_total - L_neumann) / L_neumann * 100

	# GMSH export
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
	except Exception:
		node_J = None

	return {
		"nse": mesh.GetNE(BND), "ndof": ndof,
		"area_err": area_err, "L_nH": L_total * 1e9,
		"L_err": L_err, "t_bem": t_bem,
		"mesh": mesh, "node_J": node_J,
	}


if __name__ == "__main__":
	print(f"Neumann formula: L = {L_neumann*1e9:.2f} nH")
	print(f"Exact area: {area_exact*1e4:.2f} cm^2")
	print()

	results = []

	# === Tri (tet mesh) ===
	print("=" * 70)
	print("  TRI surface (from tet mesh)")
	print("=" * 70)
	for interval in [3, 5, 8]:
		ne, n_surf, surf_type = create_tet_torus(interval)
		for order in [1, 2]:
			r = run_bem(order)
			label = f"tri int={interval}"
			print(f"  {label:15s} p={order}  nse={r['nse']:5d}  ndof={r['ndof']:5d}"
			      f"  area_err={r['area_err']:+6.2f}%  L={r['L_nH']:8.2f}nH  L_err={r['L_err']:+6.2f}%")
			results.append({"type": "tri", "interval": interval, "order": order, **r})
		print()

	# === Quad (hex sweep mesh) ===
	print("=" * 70)
	print("  QUAD surface (from hex sweep mesh)")
	print("=" * 70)
	for isrc, isw in [(2, 12), (3, 18)]:
		ne, n_surf, surf_type = create_hex_torus(isrc, isw)
		for order in [1, 2]:
			r = run_bem(order)
			label = f"quad {isrc}x{isw}"
			print(f"  {label:15s} p={order}  nse={r['nse']:5d}  ndof={r['ndof']:5d}"
			      f"  area_err={r['area_err']:+6.2f}%  L={r['L_nH']:8.2f}nH  L_err={r['L_err']:+6.2f}%")
			results.append({"type": "quad", "src": isrc, "sw": isw, "order": order, **r})
		print()

	# === Summary ===
	print("=" * 70)
	print("  SUMMARY: p=1 vs p=2 area improvement")
	print("=" * 70)
	for i in range(0, len(results), 2):
		p1, p2 = results[i], results[i + 1]
		if p1["area_err"] != 0:
			ratio = abs(p1["area_err"]) / max(abs(p2["area_err"]), 1e-10)
		else:
			ratio = 1.0
		typ = p1.get("type", "?")
		print(f"  {typ:5s}  p=1: area_err={p1['area_err']:+6.2f}%, L_err={p1['L_err']:+6.2f}%"
		      f"  ->  p=2: area_err={p2['area_err']:+6.2f}%, L_err={p2['L_err']:+6.2f}%"
		      f"  (area {ratio:.0f}x better)")
