"""
spiral_peec_mna_sweep.py

Parameter sweep: panel mode and eps_eff to find best match with measurement.
Uses proper PEEC MNA formulation: Y = A * inv(Z_branch) * A^T + jw * C_eff
"""

import sys
import numpy as np
import time
import os
import pandas as pd

import radia  # noqa: F401
sys.path.insert(0, 'S:/Radia/01_GitHub/src/radia')
from peec_matrices import PEECBuilder

MU_0 = 4.0 * np.pi * 1e-7
EPS_0 = 8.854187817e-12

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================
# Read geometry (same as spiral_peec_mna.py)
# ============================================================
df = pd.read_excel(os.path.join(SCRIPT_DIR, 'provided_by_sato', 'point.xlsx'), header=0)
df3 = df.iloc[:, :3].apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan).dropna()
x_mm = df3.iloc[:, 0].to_numpy()
y_mm = df3.iloc[:, 1].to_numpy()
z_mm = df3.iloc[:, 2].to_numpy()

keep = [0]
for i in range(1, len(x_mm)):
	dist = np.sqrt((x_mm[i]-x_mm[i-1])**2 + (y_mm[i]-y_mm[i-1])**2 + (z_mm[i]-z_mm[i-1])**2)
	if dist > 1e-10:
		keep.append(i)

pts_m = np.column_stack([x_mm[keep], y_mm[keep], z_mm[keep]]) * 1e-3
n_pts = len(pts_m)
L_wire = np.sum(np.linalg.norm(np.diff(pts_m, axis=0), axis=1))

w = 1.750e-3
pitch = 2.000e-3
sigma = 5.8e7
R_meas_dc = 1.05
h = (1.0/sigma) * L_wire / (w * R_meas_dc)
n_seg = n_pts - 1

print(f"Geometry: {n_pts} pts, {n_seg} segs, wire={L_wire*1e3:.1f} mm, h={h*1e6:.1f} um")


def make_dowell_Zs_func(n_loop, R_dc_per_edge, h, sigma, w, pitch, p=21):
	eta = w / pitch
	def Zs_func(freq):
		if freq <= 0:
			return np.zeros(n_loop, dtype=complex)
		omega = 2.0 * np.pi * freq
		delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
		Delta = (h / delta) * np.sqrt(eta)
		D2 = 2.0 * Delta
		M1 = (np.sinh(D2) + np.sin(D2)) / (np.cosh(D2) - np.cos(D2))
		D1 = (np.sinh(Delta) - np.sin(Delta)) / (np.cosh(Delta) + np.cos(Delta))
		Fr = Delta * (M1 + (2.0 / 3.0) * (p**2 - 1) * D1)
		Fr = max(Fr, 1.0)
		return R_dc_per_edge * (Fr - 1.0) + 0j
	return Zs_func


def run_mna(panel_mode, eps_eff, freqs):
	"""Build PEEC model and run MNA sweep. Returns Z_port array."""
	panels_per_seg = {0: 6, 1: 2, 2: 4, 3: 2, 4: 1}[panel_mode]

	builder = PEECBuilder()
	node_ids = []
	for p_pt in pts_m:
		nid = builder.add_node_at(p_pt[0], p_pt[1], p_pt[2], area=pitch*w)
		node_ids.append(nid)
	for i in range(n_seg):
		builder.add_connected_segment(
			node_ids[i], node_ids[i+1],
			w, h, sigma, cross_section_type=0, nwinc=1, nhinc=1
		)
	builder.add_port(node_ids[0], node_ids[-1])
	n_panels = builder.generate_face_panels(mode=panel_mode, eps_r=1.0)

	topo = builder.build_topology(include_star=True)
	n_loop = topo['n_loop']
	n_star = topo['n_star']
	n_nodes = topo['n_nodes']
	seg_nodes = np.array(topo['segment_nodes'])

	L_mat = np.array(topo['L']).reshape(n_loop, n_loop)
	R_dc = np.array(topo['R'])
	P_mat = np.array(topo['P']).reshape(n_star, n_star)

	# Gathering matrix G
	G = np.zeros((n_nodes, n_star))
	for s in range(n_seg):
		nf, nt = seg_nodes[s, 0], seg_nodes[s, 1]
		for k in range(panels_per_seg):
			p = s * panels_per_seg + k
			G[nf, p] = 0.5
			G[nt, p] = 0.5

	# C_eff with dielectric scaling
	P_eff = P_mat / eps_eff
	X = np.linalg.solve(P_eff, G.T)
	C_eff = G @ X

	# Incidence matrix
	A_inc = np.zeros((n_nodes, n_loop))
	for s in range(n_loop):
		A_inc[seg_nodes[s, 0], s] = +1.0
		A_inc[seg_nodes[s, 1], s] = -1.0

	# Dowell Zs
	Zs_func = make_dowell_Zs_func(n_loop, R_dc, h, sigma, w, pitch, p=21)

	# Ground and reduce
	port = topo['ports'][0]
	ground_node = port[1]
	reduced_nodes = [n for n in range(n_nodes) if n != ground_node]
	node_to_reduced = {ni: ri for ri, ni in enumerate(reduced_nodes)}

	A_red = A_inc[reduced_nodes, :]
	C_red = C_eff[np.ix_(reduced_nodes, reduced_nodes)]
	n_reduced = len(reduced_nodes)

	I_port = np.zeros(n_reduced, dtype=complex)
	I_port[node_to_reduced[port[0]]] = 1.0

	# Sweep
	Z_port = np.zeros(len(freqs), dtype=complex)
	for fi, freq in enumerate(freqs):
		omega = 2.0 * np.pi * freq
		Zs = Zs_func(freq)
		Z_branch = np.diag(R_dc + Zs) + 1j * omega * L_mat
		Y_branch = np.linalg.inv(Z_branch)
		Y_node = A_red @ Y_branch @ A_red.T + 1j * omega * C_red
		V_node = np.linalg.solve(Y_node, I_port)
		Z_port[fi] = V_node[node_to_reduced[port[0]]]

	return Z_port, np.sum(L_mat), np.sum(C_eff)


def find_srf(freqs, Z_port):
	"""Find SRF where Im(Z) crosses zero."""
	for i in range(len(freqs) - 1):
		if np.imag(Z_port[i]) > 0 and np.imag(Z_port[i+1]) < 0:
			f1, f2 = freqs[i], freqs[i+1]
			z1, z2 = np.imag(Z_port[i]), np.imag(Z_port[i+1])
			return f1 - z1 * (f2 - f1) / (z2 - z1)
	return None


# ============================================================
# Parameter sweep
# ============================================================
freqs = np.logspace(2, 7.7, 50)

configs = [
	(3, 1.0,  "mode=3 sides, vacuum"),
	(3, 2.7,  "mode=3 sides, eps_eff=2.7"),
	(3, 4.4,  "mode=3 sides, eps_eff=4.4"),
	(4, 1.0,  "mode=4 top, vacuum"),
	(4, 2.7,  "mode=4 top, eps_eff=2.7"),
	(4, 4.4,  "mode=4 top, eps_eff=4.4"),
]

results = []
print(f"\n{'Config':<30s}  {'C_eff[pF]':>9s}  {'SRF[MHz]':>9s}  {'L_dc[uH]':>9s}")
print("-" * 65)

for mode, eps_eff, label in configs:
	t0 = time.perf_counter()
	Z_port, L_dc, C_total = run_mna(mode, eps_eff, freqs)
	dt = time.perf_counter() - t0
	srf = find_srf(freqs, Z_port)
	srf_str = f"{srf/1e6:.2f}" if srf else "none"
	print(f"{label:<30s}  {C_total*1e12:9.2f}  {srf_str:>9s}  {L_dc*1e6:9.2f}  ({dt:.0f}s)")
	results.append((mode, eps_eff, label, Z_port, L_dc, C_total, srf))

# ============================================================
# Plot comparison
# ============================================================
try:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import openpyxl

	# Load measurement
	wb = openpyxl.load_workbook(os.path.join(SCRIPT_DIR, 'provided_by_sato', '14tAir_\u5468\u6ce2\u6570\u7279\u6027\u6e2c\u5b9a\u5024.xlsx'), data_only=True)
	ws = wb['Sheet1']
	m_freqs, m_Rs, m_Ls = [], [], []
	for row in ws.iter_rows(min_row=4, values_only=True):
		if row[1] is not None:
			try:
				m_freqs.append(float(row[1]))
				m_Rs.append(float(row[4]))
				m_Ls.append(float(row[5]))
			except (ValueError, TypeError):
				pass
	m_freqs = np.array(m_freqs)
	m_Rs = np.array(m_Rs)
	m_Ls = np.array(m_Ls)

	fig, axes = plt.subplots(1, 2, figsize=(14, 6))
	colors = ['gray', 'blue', 'darkblue', 'lightgreen', 'green', 'darkgreen']

	for i, (mode, eps_eff, label, Z_port, L_dc, C_total, srf) in enumerate(results):
		L_eff = np.imag(Z_port) / (2*np.pi*freqs)
		R_eff = np.real(Z_port)
		srf_str = f" SRF={srf/1e6:.1f}M" if srf else ""
		axes[0].semilogx(freqs, L_eff * 1e6, color=colors[i], linewidth=1.5,
						 label=f"{label}{srf_str}")
		axes[1].loglog(freqs, np.abs(R_eff), color=colors[i], linewidth=1.5,
					   label=f"{label}")

	# Measurement
	mask_L = m_Ls > 0
	axes[0].semilogx(m_freqs[mask_L], m_Ls[mask_L] * 1e6, 'r-', linewidth=2, label='Measurement')
	axes[1].loglog(m_freqs, m_Rs, 'r-', linewidth=2, label='Measurement')

	axes[0].set_xlabel('Frequency [Hz]')
	axes[0].set_ylabel('Inductance [uH]')
	axes[0].set_title('L(f) -- PEEC MNA parameter sweep')
	axes[0].legend(fontsize=7)
	axes[0].grid(True, alpha=0.3)
	axes[0].set_ylim(-50, 100)

	axes[1].set_xlabel('Frequency [Hz]')
	axes[1].set_ylabel('|R(f)| [Ohm]')
	axes[1].set_title('R(f) -- PEEC MNA parameter sweep')
	axes[1].legend(fontsize=7)
	axes[1].grid(True, alpha=0.3)

	plt.tight_layout()
	out_path = os.path.join(SCRIPT_DIR, 'spiral_peec_mna_sweep.png')
	plt.savefig(out_path, dpi=150)
	plt.close()
	print(f"\nFigure saved: {out_path}")
except Exception as e:
	print(f"Plot error: {e}")

print("\nDone.")
