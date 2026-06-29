"""
spiral_peec_mna_pfit.py

Find optimal Dowell p parameter by fitting to measurement R(f).
Also test a two-parameter model: Dowell low-frequency + saturation at high frequency.
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


def make_dowell_Zs(n_loop, R_dc, h, sigma, w, pitch, p):
	eta = w / pitch
	def Zs_func(freq):
		if freq <= 0:
			return np.zeros(n_loop, dtype=complex)
		omega = 2*np.pi*freq
		delta = np.sqrt(2/(omega*MU_0*sigma))
		Delta = (h/delta)*np.sqrt(eta)
		D2 = 2*Delta
		M1 = (np.sinh(D2)+np.sin(D2))/(np.cosh(D2)-np.cos(D2))
		D1 = (np.sinh(Delta)-np.sin(Delta))/(np.cosh(Delta)+np.cos(Delta))
		Fr = Delta*(M1 + (2/3)*(p**2-1)*D1)
		Fr = max(Fr, 1.0)
		return R_dc*(Fr - 1.0) + 0j
	return Zs_func


def make_proximity_Zs(n_loop, R_dc, h, sigma, w, pitch, K_prox):
	"""Physical proximity model for thin conductors (h << delta).

	R_prox(f) = K_prox * R_dc * (omega * mu0 * sigma * h^2 / 2)

	For h << delta: proximity loss ~ omega^2 * h^2 (eddy currents from external B)
	For h ~ delta:  transition to ~ omega^0.5 (skin-effect limited penetration)

	Uses Dowell's D1 function to capture the h/delta transition naturally.
	K_prox is fitted to measurement (encodes the geometric proximity field strength).
	"""
	eta = w / pitch

	def Zs_func(freq):
		if freq <= 0:
			return np.zeros(n_loop, dtype=complex)
		omega = 2*np.pi*freq
		delta = np.sqrt(2/(omega*MU_0*sigma))
		Delta = (h/delta)*np.sqrt(eta)

		# Dowell D1: captures the h/delta transition naturally
		# D1 ~ Delta^3/6 for small Delta, D1 ~ 1 for large Delta
		D1 = (np.sinh(Delta)-np.sin(Delta))/(np.cosh(Delta)+np.cos(Delta))

		# Proximity Fr: K_prox plays the role of (2/3)(p^2-1)
		Fr_prox = K_prox * Delta * D1

		# Skin effect: M1 term (negligible for thin conductors but included)
		D2 = 2*Delta
		M1 = (np.sinh(D2)+np.sin(D2))/(np.cosh(D2)-np.cos(D2))
		Fr_skin = Delta * M1

		Fr = Fr_skin + Fr_prox
		Fr = max(Fr, 1.0)
		return R_dc * (Fr - 1.0) + 0j

	return Zs_func


# ============================================================
# Read geometry, build PEEC
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

builder = PEECBuilder()
node_ids = []
for p_pt in pts_m:
	node_ids.append(builder.add_node_at(p_pt[0], p_pt[1], p_pt[2], area=pitch*w))
for i in range(n_seg):
	builder.add_connected_segment(node_ids[i], node_ids[i+1], w, h, sigma,
								  cross_section_type=0, nwinc=1, nhinc=1)
builder.add_port(node_ids[0], node_ids[-1])
builder.generate_face_panels(mode=3, eps_r=1.0)

topo = builder.build_topology(include_star=True)
n_loop, n_star, n_nodes = topo['n_loop'], topo['n_star'], topo['n_nodes']
seg_nodes = np.array(topo['segment_nodes'])
L_mat = np.array(topo['L']).reshape(n_loop, n_loop)
R_dc = np.array(topo['R'])
P_mat = np.array(topo['P']).reshape(n_star, n_star)

# C_eff
eps_eff = (1.0 + 4.4) / 2.0
G = np.zeros((n_nodes, n_star))
for s in range(n_seg):
	nf, nt = seg_nodes[s, 0], seg_nodes[s, 1]
	G[nf, 2*s] = G[nt, 2*s] = 0.5
	G[nf, 2*s+1] = G[nt, 2*s+1] = 0.5
X = np.linalg.solve(P_mat / eps_eff, G.T)
C_eff = G @ X

# MNA setup
A_inc = np.zeros((n_nodes, n_loop))
for s in range(n_loop):
	A_inc[seg_nodes[s, 0], s] = +1.0
	A_inc[seg_nodes[s, 1], s] = -1.0

port = topo['ports'][0]
ground_node = port[1]
reduced_nodes = [n for n in range(n_nodes) if n != ground_node]
node_to_reduced = {ni: ri for ri, ni in enumerate(reduced_nodes)}
A_red = A_inc[reduced_nodes, :]
C_red = C_eff[np.ix_(reduced_nodes, reduced_nodes)]
I_port = np.zeros(len(reduced_nodes), dtype=complex)
I_port[node_to_reduced[port[0]]] = 1.0


def mna_sweep(Zs_func, freqs):
	Z_port = np.zeros(len(freqs), dtype=complex)
	for fi, freq in enumerate(freqs):
		omega = 2*np.pi*freq
		Zs = Zs_func(freq)
		Z_branch = np.diag(R_dc + Zs) + 1j*omega*L_mat
		Y_branch = np.linalg.inv(Z_branch)
		Y_node = A_red @ Y_branch @ A_red.T + 1j*omega*C_red
		V_node = np.linalg.solve(Y_node, I_port)
		Z_port[fi] = V_node[node_to_reduced[port[0]]]
	return Z_port


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
# Scan Dowell p values for best R(f) fit
# ============================================================
freqs_fit = np.array([1e5, 3e5, 5e5, 1e6, 2e6])
R_meas_fit = np.array([m_Rs[np.argmin(np.abs(m_freqs-f))] for f in freqs_fit])

# Quick scan with coarse freq grid
freqs_coarse = np.logspace(4, 6.5, 30)
print("Scanning Dowell p for best R(f) fit (100kHz - 2MHz)...")
best_err = 1e10
best_p = 0

for p_test in np.arange(3, 25, 0.5):
	Zs_func = make_dowell_Zs(n_loop, R_dc, h, sigma, w, pitch, p=p_test)
	Z_port = mna_sweep(Zs_func, freqs_fit)
	R_peec = np.real(Z_port)
	err = np.mean(((R_peec - R_meas_fit)/R_meas_fit)**2)
	if err < best_err:
		best_err = err
		best_p = p_test

print(f"Best Dowell p = {best_p:.1f} (RMSE = {np.sqrt(best_err)*100:.1f}%)")

# Now scan K_prox for the physical proximity model
print("\nScanning K_prox for proximity model...")
best_err2 = 1e10
best_K = 0

for K_test in np.arange(10, 500, 5):
	Zs_func = make_proximity_Zs(n_loop, R_dc, h, sigma, w, pitch, K_prox=K_test)
	Z_port = mna_sweep(Zs_func, freqs_fit)
	R_peec = np.real(Z_port)
	err = np.mean(((R_peec - R_meas_fit)/R_meas_fit)**2)
	if err < best_err2:
		best_err2 = err
		best_K = K_test

print(f"Best K_prox = {best_K:.0f} (RMSE = {np.sqrt(best_err2)*100:.1f}%)")

# Equivalent p for K_prox: K = (2/3)(p^2-1) -> p = sqrt(3K/2+1)
p_equiv = np.sqrt(3*best_K/2 + 1)
print(f"Equivalent Dowell p = {p_equiv:.1f}")

# ============================================================
# Final comparison with fine freq grid
# ============================================================
freqs = np.logspace(2, 7.2, 80)
print(f"\nFinal comparison (80 freq points)...")

models = [
	(f"Dowell p={best_p:.0f} (fit)", make_dowell_Zs(n_loop, R_dc, h, sigma, w, pitch, p=best_p)),
	("Dowell p=21 (old)",  make_dowell_Zs(n_loop, R_dc, h, sigma, w, pitch, p=21)),
	(f"Proximity K={best_K:.0f}", make_proximity_Zs(n_loop, R_dc, h, sigma, w, pitch, K_prox=best_K)),
]

all_results = []
for label, Zs_func in models:
	t0 = time.perf_counter()
	Z_port = mna_sweep(Zs_func, freqs)
	dt = time.perf_counter() - t0
	srf = None
	for i in range(len(freqs)-1):
		if np.imag(Z_port[i]) > 0 and np.imag(Z_port[i+1]) < 0:
			f1, f2 = freqs[i], freqs[i+1]
			z1, z2 = np.imag(Z_port[i]), np.imag(Z_port[i+1])
			srf = f1 - z1*(f2-f1)/(z2-z1)
			break
	srf_s = f"{srf/1e6:.2f}M" if srf else "none"
	print(f"  {label:<25s}: SRF={srf_s}, {dt:.0f}s")
	all_results.append((label, Z_port, srf))

# Table
print(f"\n{'freq':>8s}  {'Meas':>7s}", end="")
for label, _, _ in all_results:
	print(f"  {label[:16]:>16s}", end="")
print()
for f_target in [1e4, 1e5, 3e5, 5e5, 1e6, 2e6, 3e6, 5e6]:
	im = np.argmin(np.abs(m_freqs - f_target))
	ip = np.argmin(np.abs(freqs - f_target))
	print(f"{f_target/1e6:7.1f}M  {m_Rs[im]:7.2f}", end="")
	for label, Z_port, _ in all_results:
		R = np.real(Z_port[ip])
		err = (R - m_Rs[im])/m_Rs[im]*100
		print(f"  {R:7.2f} ({err:+5.0f}%)", end="")
	print()

# ============================================================
# Plot
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
colors = ['blue', 'gray', 'green']

# L(f)
ax = axes[0]
for i, (label, Z_port, srf) in enumerate(all_results):
	L_eff = np.imag(Z_port)/(2*np.pi*freqs)
	srf_s = f" SRF={srf/1e6:.1f}M" if srf else ""
	ax.semilogx(freqs, L_eff*1e6, color=colors[i], linewidth=1.5, label=f"{label}{srf_s}")
mask_L = m_Ls > 0
ax.semilogx(m_freqs[mask_L], m_Ls[mask_L]*1e6, 'r-', linewidth=2, label='Measurement')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('L [uH]')
ax.set_title('L(f)')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_ylim(-50, 200)

# R(f)
ax = axes[1]
for i, (label, Z_port, srf) in enumerate(all_results):
	R_eff = np.real(Z_port)
	ax.loglog(freqs, np.abs(R_eff), color=colors[i], linewidth=1.5, label=label)
ax.loglog(m_freqs, m_Rs, 'r-', linewidth=2, label='Measurement')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('R [Ohm]')
ax.set_title('R(f)')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

plt.suptitle(f'PEEC MNA: Optimized R(f) models vs measurement', fontsize=12)
plt.tight_layout()
out_path = os.path.join(SCRIPT_DIR, 'spiral_peec_mna_pfit.png')
plt.savefig(out_path, dpi=150)
plt.close()
print(f"\nFigure saved: {out_path}")
