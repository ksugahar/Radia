"""
spiral_peec_mna_compare.py

Detailed L(f) and R(f) comparison: PEEC MNA vs measurement.
Best model: mode=3 (side panels), eps_eff=2.7 (half-space image).
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
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


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


# ============================================================
# Read geometry
# ============================================================
df = pd.read_excel(os.path.join(SCRIPT_DIR, 'provided_by_sato', 'point.xlsx'), header=0)
df3 = df.iloc[:, :3].apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan).dropna()
x_mm, y_mm, z_mm = df3.iloc[:, 0].to_numpy(), df3.iloc[:, 1].to_numpy(), df3.iloc[:, 2].to_numpy()

keep = [0]
for i in range(1, len(x_mm)):
	if np.sqrt((x_mm[i]-x_mm[i-1])**2 + (y_mm[i]-y_mm[i-1])**2 + (z_mm[i]-z_mm[i-1])**2) > 1e-10:
		keep.append(i)

pts_m = np.column_stack([x_mm[keep], y_mm[keep], z_mm[keep]]) * 1e-3
n_pts = len(pts_m)
n_seg = n_pts - 1
L_wire = np.sum(np.linalg.norm(np.diff(pts_m, axis=0), axis=1))

w = 1.750e-3
pitch = 2.000e-3
sigma = 5.8e7
R_meas_dc = 1.05
h = (1.0/sigma) * L_wire / (w * R_meas_dc)

print(f"Geometry: {n_pts} pts, {n_seg} segs, wire={L_wire*1e3:.1f} mm, h={h*1e6:.1f} um")

# ============================================================
# Build PEEC model: mode=3, eps_eff=2.7
# ============================================================
PANEL_MODE = 3
EPS_R = 4.4
eps_eff = (1.0 + EPS_R) / 2.0
panels_per_seg = 2  # mode=3

print(f"Panel mode={PANEL_MODE}, eps_eff={eps_eff:.2f}")

builder = PEECBuilder()
node_ids = []
for p_pt in pts_m:
	nid = builder.add_node_at(p_pt[0], p_pt[1], p_pt[2], area=pitch*w)
	node_ids.append(nid)
for i in range(n_seg):
	builder.add_connected_segment(node_ids[i], node_ids[i+1], w, h, sigma,
								  cross_section_type=0, nwinc=1, nhinc=1)
builder.add_port(node_ids[0], node_ids[-1])
builder.generate_face_panels(mode=PANEL_MODE, eps_r=1.0)

t0 = time.perf_counter()
topo = builder.build_topology(include_star=True)
dt_build = time.perf_counter() - t0

n_loop = topo['n_loop']
n_star = topo['n_star']
n_nodes = topo['n_nodes']
seg_nodes = np.array(topo['segment_nodes'])

L_mat = np.array(topo['L']).reshape(n_loop, n_loop)
R_dc = np.array(topo['R'])
P_mat = np.array(topo['P']).reshape(n_star, n_star)
L_dc = np.sum(L_mat)
print(f"n_loop={n_loop}, n_star={n_star}, L_DC={L_dc*1e6:.2f} uH, build={dt_build:.1f}s")

# Gathering matrix G
G = np.zeros((n_nodes, n_star))
for s in range(n_seg):
	nf, nt = seg_nodes[s, 0], seg_nodes[s, 1]
	for k in range(panels_per_seg):
		G[nf, s*panels_per_seg+k] = 0.5
		G[nt, s*panels_per_seg+k] = 0.5

# C_eff
P_eff = P_mat / eps_eff
X = np.linalg.solve(P_eff, G.T)
C_eff = G @ X
print(f"C_eff_total = {np.sum(C_eff)*1e12:.2f} pF")

# Incidence matrix
A_inc = np.zeros((n_nodes, n_loop))
for s in range(n_loop):
	A_inc[seg_nodes[s, 0], s] = +1.0
	A_inc[seg_nodes[s, 1], s] = -1.0

# Ground and reduce
port = topo['ports'][0]
ground_node = port[1]
reduced_nodes = [n for n in range(n_nodes) if n != ground_node]
node_to_reduced = {ni: ri for ri, ni in enumerate(reduced_nodes)}
A_red = A_inc[reduced_nodes, :]
C_red = C_eff[np.ix_(reduced_nodes, reduced_nodes)]
I_port = np.zeros(len(reduced_nodes), dtype=complex)
I_port[node_to_reduced[port[0]]] = 1.0

Zs_func = make_dowell_Zs_func(n_loop, R_dc, h, sigma, w, pitch, p=21)

# ============================================================
# Fine frequency sweep (100 points, log-spaced)
# ============================================================
freqs = np.logspace(2, 7.2, 100)  # 100 Hz to ~16 MHz
print(f"\nMNA sweep ({len(freqs)} points)...")
t0 = time.perf_counter()

Z_port = np.zeros(len(freqs), dtype=complex)
for fi, freq in enumerate(freqs):
	omega = 2.0 * np.pi * freq
	Zs = Zs_func(freq)
	Z_branch = np.diag(R_dc + Zs) + 1j * omega * L_mat
	Y_branch = np.linalg.inv(Z_branch)
	Y_node = A_red @ Y_branch @ A_red.T + 1j * omega * C_red
	V_node = np.linalg.solve(Y_node, I_port)
	Z_port[fi] = V_node[node_to_reduced[port[0]]]

dt_sweep = time.perf_counter() - t0
print(f"Sweep: {dt_sweep:.0f} s ({dt_sweep/len(freqs):.1f} s/freq)")

L_eff = np.imag(Z_port) / (2*np.pi*freqs)
R_eff = np.real(Z_port)

# SRF
for i in range(len(freqs)-1):
	if np.imag(Z_port[i]) > 0 and np.imag(Z_port[i+1]) < 0:
		f1, f2 = freqs[i], freqs[i+1]
		z1, z2 = np.imag(Z_port[i]), np.imag(Z_port[i+1])
		srf = f1 - z1*(f2-f1)/(z2-z1)
		print(f"SRF = {srf/1e6:.2f} MHz")
		break

# ============================================================
# Load measurement
# ============================================================
import openpyxl
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

# ============================================================
# Numerical comparison at specific frequencies
# ============================================================
print(f"\n{'freq':>10s}  {'L_peec':>8s}  {'L_meas':>8s}  {'err_L':>7s}  {'R_peec':>8s}  {'R_meas':>8s}  {'err_R':>7s}")
for f_target in [1e3, 1e4, 1e5, 3e5, 5e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6]:
	ip = np.argmin(np.abs(freqs - f_target))
	im = np.argmin(np.abs(m_freqs - f_target))
	Lp = L_eff[ip]*1e6
	Lm = m_Ls[im]*1e6
	Rp = R_eff[ip]
	Rm = m_Rs[im]
	err_L = (Lp - Lm) / Lm * 100 if abs(Lm) > 0.01 else float('nan')
	err_R = (Rp - Rm) / Rm * 100 if Rm > 0.01 else float('nan')
	print(f"{freqs[ip]:10.0f}  {Lp:8.2f}  {Lm:8.2f}  {err_L:6.1f}%  {Rp:8.3f}  {Rm:8.3f}  {err_R:6.1f}%")

# ============================================================
# Plot
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# --- (a) L(f) full range ---
ax = axes[0, 0]
ax.semilogx(freqs, L_eff*1e6, 'b-', linewidth=2, label=f'PEEC MNA (mode={PANEL_MODE}, eps_eff={eps_eff:.1f})')
mask_L = m_Ls > 0
ax.semilogx(m_freqs[mask_L], m_Ls[mask_L]*1e6, 'r-', linewidth=1.5, alpha=0.7, label='Measurement')
ax.axhline(y=L_dc*1e6, color='gray', ls='--', alpha=0.5, label=f'L_DC={L_dc*1e6:.1f} uH')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Inductance [uH]')
ax.set_title('(a) L(f) -- full range')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(-50, 200)

# --- (b) L(f) zoom to pre-resonance ---
ax = axes[0, 1]
mask_f = freqs < 5e6
ax.semilogx(freqs[mask_f], L_eff[mask_f]*1e6, 'b-', linewidth=2, label='PEEC MNA')
mask_mf = (m_freqs < 5e6) & (m_Ls > 0)
ax.semilogx(m_freqs[mask_mf], m_Ls[mask_mf]*1e6, 'r-', linewidth=1.5, alpha=0.7, label='Measurement')
ax.axhline(y=L_dc*1e6, color='gray', ls='--', alpha=0.5)
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Inductance [uH]')
ax.set_title('(b) L(f) -- zoom (f < 5 MHz)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(15, 55)

# --- (c) R(f) full range ---
ax = axes[1, 0]
ax.loglog(freqs, np.abs(R_eff), 'b-', linewidth=2, label='PEEC MNA + Dowell p=21')
ax.loglog(m_freqs, m_Rs, 'r-', linewidth=1.5, alpha=0.7, label='Measurement')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('|R(f)| [Ohm]')
ax.set_title('(c) R(f) -- full range')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# --- (d) R(f) zoom to low frequency ---
ax = axes[1, 1]
mask_f2 = freqs < 3e6
ax.loglog(freqs[mask_f2], np.abs(R_eff[mask_f2]), 'b-', linewidth=2, label='PEEC MNA + Dowell p=21')
mask_mf2 = m_freqs < 3e6
ax.loglog(m_freqs[mask_mf2], m_Rs[mask_mf2], 'r-', linewidth=1.5, alpha=0.7, label='Measurement')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('|R(f)| [Ohm]')
ax.set_title('(d) R(f) -- zoom (f < 3 MHz)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.suptitle(f'PEEC MNA vs Measurement: 14-turn spiral on FR4\n'
			 f'mode={PANEL_MODE}, eps_eff={eps_eff:.1f}, Dowell p=21, L_DC={L_dc*1e6:.1f} uH',
			 fontsize=12)
plt.tight_layout()
out_path = os.path.join(SCRIPT_DIR, 'spiral_peec_mna_compare.png')
plt.savefig(out_path, dpi=150)
plt.close()
print(f"\nFigure saved: {out_path}")
print("Done.")
