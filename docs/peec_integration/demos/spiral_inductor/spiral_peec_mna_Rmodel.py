"""
spiral_peec_mna_Rmodel.py

Compare different frequency-dependent resistance models with measurement.
All use the same proper MNA formulation: Y = A inv(Z) A^T + jw C_eff

Models:
  1. Dowell p=21 (current)
  2. Dowell p=7  (reduced proximity)
  3. Bessel SIBC (skin effect only, no proximity)
  4. Bessel SIBC + sqrt(f) proximity (empirical)
"""

import numpy as np
import time
import os
import pandas as pd
from scipy.special import jv as besselj  # Bessel functions

from radia.peec_matrices import PEECBuilder

MU_0 = 4.0 * np.pi * 1e-7
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# Surface impedance models
# ============================================================

def make_dowell_Zs(n_loop, R_dc, h, sigma, w, pitch, p):
	"""Dowell winding loss: Fr = Delta*(M1 + (2/3)(p^2-1)*D1)"""
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


def make_bessel_Zs(n_loop, R_dc, h, sigma, w):
	"""Bessel SIBC for rectangular conductor (skin effect only).

	Z_s/R_dc = (1+j)*xi * coth((1+j)*xi)  where xi = h/(2*delta)

	This is the exact solution for a conducting slab of thickness h
	with current flowing along the length.
	"""
	def Zs_func(freq):
		if freq <= 0:
			return np.zeros(n_loop, dtype=complex)
		omega = 2*np.pi*freq
		delta = np.sqrt(2/(omega*MU_0*sigma))
		xi = h / (2*delta)
		# Z_ac/R_dc = (1+j)*xi * coth((1+j)*xi)
		z = (1+1j)*xi
		# coth(z) = cosh(z)/sinh(z), handle overflow
		if abs(z) > 500:
			ratio = 1.0  # coth -> 1 for large z
		else:
			ratio = np.cosh(z)/np.sinh(z)
		Zac_over_Rdc = z * ratio
		# Fr = Re(Zac/Rdc), and Zs = Rdc*(Zac/Rdc - 1)
		return R_dc * (Zac_over_Rdc - 1.0)
	return Zs_func


def make_bessel_prox_Zs(n_loop, R_dc, h, sigma, w, pitch, alpha=0.8):
	"""Bessel SIBC + empirical sqrt(f) proximity correction.

	R(f) = R_dc * [Fr_skin + alpha * (n_adj) * G_prox(Delta)]

	where G_prox(Delta) = Delta^2/(Delta^2+1) is a saturation function
	that rises as Delta^2 at low freq and saturates to 1 at high freq.
	n_adj = effective number of adjacent turns.
	"""
	eta = w / pitch

	def Zs_func(freq):
		if freq <= 0:
			return np.zeros(n_loop, dtype=complex)
		omega = 2*np.pi*freq
		delta = np.sqrt(2/(omega*MU_0*sigma))

		# Skin effect (Bessel SIBC)
		xi = h / (2*delta)
		z = (1+1j)*xi
		if abs(z) > 500:
			ratio = 1.0
		else:
			ratio = np.cosh(z)/np.sinh(z)
		Fr_skin = z * ratio

		# Proximity correction
		# Use Dowell's D1 function but with moderate p
		Delta = (h/delta)*np.sqrt(eta)
		D1 = (np.sinh(Delta)-np.sin(Delta))/(np.cosh(Delta)+np.cos(Delta))
		# Saturating proximity: rises as Delta^2 for small Delta, saturates
		Fr_prox = alpha * Delta * D1

		return R_dc * (Fr_skin - 1.0 + Fr_prox)
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

print(f"Geometry: {n_pts} pts, {n_seg} segs, h={h*1e6:.1f} um")

# ============================================================
# Build PEEC model (shared across all R models)
# ============================================================
PANEL_MODE = 3
EPS_R = 4.4
eps_eff = (1.0 + EPS_R) / 2.0

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

topo = builder.build_topology(include_star=True)
n_loop = topo['n_loop']
n_star = topo['n_star']
n_nodes = topo['n_nodes']
seg_nodes = np.array(topo['segment_nodes'])

L_mat = np.array(topo['L']).reshape(n_loop, n_loop)
R_dc = np.array(topo['R'])
P_mat = np.array(topo['P']).reshape(n_star, n_star)

# Gathering, C_eff, incidence (same for all models)
G = np.zeros((n_nodes, n_star))
for s in range(n_seg):
	nf, nt = seg_nodes[s, 0], seg_nodes[s, 1]
	G[nf, 2*s] = 0.5; G[nt, 2*s] = 0.5
	G[nf, 2*s+1] = 0.5; G[nt, 2*s+1] = 0.5

X = np.linalg.solve(P_mat / eps_eff, G.T)
C_eff = G @ X

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
# Models to compare
# ============================================================
freqs = np.logspace(2, 7.2, 80)

models = [
	("Dowell p=21",  make_dowell_Zs(n_loop, R_dc, h, sigma, w, pitch, p=21)),
	("Dowell p=7",   make_dowell_Zs(n_loop, R_dc, h, sigma, w, pitch, p=7)),
	("Bessel SIBC",  make_bessel_Zs(n_loop, R_dc, h, sigma, w)),
	("Bessel+prox a=0.8", make_bessel_prox_Zs(n_loop, R_dc, h, sigma, w, pitch, alpha=0.8)),
]

results = []
for label, Zs_func in models:
	t0 = time.perf_counter()
	Z_port = mna_sweep(Zs_func, freqs)
	dt = time.perf_counter() - t0
	# SRF
	srf = None
	for i in range(len(freqs)-1):
		if np.imag(Z_port[i]) > 0 and np.imag(Z_port[i+1]) < 0:
			f1, f2 = freqs[i], freqs[i+1]
			z1, z2 = np.imag(Z_port[i]), np.imag(Z_port[i+1])
			srf = f1 - z1*(f2-f1)/(z2-z1)
			break
	srf_s = f"{srf/1e6:.2f} MHz" if srf else "none"
	print(f"  {label:<25s}: SRF={srf_s}, {dt:.0f}s")
	results.append((label, Z_port, srf))

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
# Numerical comparison table
# ============================================================
print(f"\nR(f) comparison:")
print(f"{'freq':>8s}  {'Meas':>7s}", end="")
for label, _, _ in results:
	print(f"  {label[:12]:>12s}", end="")
print()

for f_target in [1e4, 1e5, 3e5, 5e5, 1e6, 2e6, 3e6, 5e6]:
	im = np.argmin(np.abs(m_freqs - f_target))
	ip = np.argmin(np.abs(freqs - f_target))
	print(f"{f_target/1e6:7.1f}M  {m_Rs[im]:7.2f}", end="")
	for label, Z_port, _ in results:
		R = np.real(Z_port[ip])
		print(f"  {R:12.2f}", end="")
	print()

# ============================================================
# Plot
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
colors = ['blue', 'green', 'orange', 'purple']

# (a) L(f) zoom
ax = axes[0, 0]
for i, (label, Z_port, srf) in enumerate(results):
	L_eff = np.imag(Z_port)/(2*np.pi*freqs)
	ax.semilogx(freqs, L_eff*1e6, color=colors[i], linewidth=1.5, label=label)
mask_L = (m_Ls > 0) & (m_freqs < 8e6)
ax.semilogx(m_freqs[mask_L], m_Ls[mask_L]*1e6, 'r-', linewidth=2, label='Measurement')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('L [uH]')
ax.set_title('(a) L(f)')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_ylim(15, 80)

# (b) L(f) near SRF
ax = axes[0, 1]
for i, (label, Z_port, srf) in enumerate(results):
	L_eff = np.imag(Z_port)/(2*np.pi*freqs)
	ax.semilogx(freqs, L_eff*1e6, color=colors[i], linewidth=1.5, label=label)
mask_L2 = (m_Ls > 0)
ax.semilogx(m_freqs[mask_L2], m_Ls[mask_L2]*1e6, 'r-', linewidth=2, label='Measurement')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('L [uH]')
ax.set_title('(b) L(f) -- full range')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_ylim(-100, 400)

# (c) R(f) zoom < 3 MHz
ax = axes[1, 0]
for i, (label, Z_port, srf) in enumerate(results):
	R_eff = np.real(Z_port)
	mask = freqs < 3e6
	ax.loglog(freqs[mask], np.abs(R_eff[mask]), color=colors[i], linewidth=1.5, label=label)
mask_m = m_freqs < 3e6
ax.loglog(m_freqs[mask_m], m_Rs[mask_m], 'r-', linewidth=2, label='Measurement')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('R [Ohm]')
ax.set_title('(c) R(f) zoom (f < 3 MHz)')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# (d) R(f) full range
ax = axes[1, 1]
for i, (label, Z_port, srf) in enumerate(results):
	R_eff = np.real(Z_port)
	ax.loglog(freqs, np.abs(R_eff), color=colors[i], linewidth=1.5, label=label)
ax.loglog(m_freqs, m_Rs, 'r-', linewidth=2, label='Measurement')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('R [Ohm]')
ax.set_title('(d) R(f) full range')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

plt.suptitle('R(f) model comparison: PEEC MNA + different skin/proximity models', fontsize=12)
plt.tight_layout()
out_path = os.path.join(SCRIPT_DIR, 'spiral_peec_mna_Rmodel.png')
plt.savefig(out_path, dpi=150)
plt.close()
print(f"\nFigure saved: {out_path}")
