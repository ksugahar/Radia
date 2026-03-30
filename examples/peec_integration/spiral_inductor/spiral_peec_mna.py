"""
spiral_peec_mna.py

Proper PEEC MNA frequency sweep for the spiral inductor.

Replaces the incorrect M_LS Schur complement approach with the standard
PEEC Modified Nodal Analysis:

    Y_node = A x Z_branch⁻¹ x A^T + jomega x C_eff

where:
    C_eff = G x P⁻¹ x G^T  (panel-to-node gathered capacitance)

This correctly produces self-resonance because:
    - Low freq:  Z ~= jomegaL  (inductive)
    - High freq: Z ~= -j/(omegaC)  (capacitive)
    - SRF:       Im(Z) = 0  (resonance)

The previous Schur complement Z_eff = Z_LL - Z_LS Z_SS⁻¹ Z_SL gives:
    Z_eff = R + jomegaL + jomega^3 x M_LS P⁻¹ M_LS^T
which is ALWAYS inductive — no resonance possible.
"""

import sys
import numpy as np
import time
import os
import pandas as pd

# Import radia first to set up MKL DLL directories on Windows
import radia  # noqa: F401 -- triggers DLL path setup in __init__.py

sys.path.insert(0, 'S:/Radia/01_GitHub/src/radia')
from peec_matrices import PEECBuilder

MU_0 = 4.0 * np.pi * 1e-7
EPS_0 = 8.854187817e-12


def make_dowell_Zs_func(n_loop, R_dc_per_edge, h, sigma, w, pitch, p=21):
	"""Full Dowell winding loss with skin + proximity effect."""
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
# Read point.xlsx geometry
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
df = pd.read_excel(os.path.join(SCRIPT_DIR, 'provided_by_sato', 'point.xlsx'), header=0)
df3 = df.iloc[:, :3].apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan).dropna()
x_mm = df3.iloc[:, 0].to_numpy()
y_mm = df3.iloc[:, 1].to_numpy()
z_mm = df3.iloc[:, 2].to_numpy()

# Remove duplicate consecutive points
keep = [0]
for i in range(1, len(x_mm)):
	dist = np.sqrt((x_mm[i]-x_mm[i-1])**2 + (y_mm[i]-y_mm[i-1])**2 + (z_mm[i]-z_mm[i-1])**2)
	if dist > 1e-10:
		keep.append(i)

pts_mm = np.column_stack([x_mm[keep], y_mm[keep], z_mm[keep]])
pts_m = pts_mm * 1e-3  # mm -> meters
n_pts = len(pts_m)
L_wire = np.sum(np.linalg.norm(np.diff(pts_m, axis=0), axis=1))

print(f"Geometry: {n_pts} points, {n_pts-1} segments, wire={L_wire*1e3:.1f} mm")

# ============================================================
# Parameters
# ============================================================
w = 1.750e-3       # track width [m]
pitch = 2.000e-3   # center-to-center [m]
sigma = 5.8e7      # copper
rho = 1.0 / sigma
R_meas_dc = 1.05
h = rho * L_wire / (w * R_meas_dc)  # derived thickness
print(f"h = {h*1e6:.2f} um (derived from R_DC={R_meas_dc} Ohm)")

# ============================================================
# Build PEEC model
# ============================================================
PANEL_MODE = 3       # 3=left+right sides, 4=top only
EPS_R = 4.4          # FR4 relative permittivity

print(f"\nBuilding PEEC model (panel mode={PANEL_MODE}, eps_r={EPS_R})...")

builder = PEECBuilder()
node_ids = []
for p_pt in pts_m:
	nid = builder.add_node_at(p_pt[0], p_pt[1], p_pt[2], area=pitch*w)
	node_ids.append(nid)

n_seg = n_pts - 1
for i in range(n_seg):
	builder.add_connected_segment(
		node_ids[i], node_ids[i+1],
		w, h, sigma,
		cross_section_type=0, nwinc=1, nhinc=1
	)
builder.add_port(node_ids[0], node_ids[-1])

# Generate face panels (eps_r=1.0 for free-space P, we handle dielectric separately)
n_panels = builder.generate_face_panels(mode=PANEL_MODE, eps_r=1.0)
print(f"Generated {n_panels} face panels (mode={PANEL_MODE})")

# Panels per segment (depends on mode)
panels_per_seg = {0: 6, 1: 2, 2: 4, 3: 2, 4: 1}[PANEL_MODE]

t0 = time.perf_counter()
topo = builder.build_topology(include_star=True)
dt_build = time.perf_counter() - t0

n_loop = topo['n_loop']
n_star = topo['n_star']
n_nodes = topo['n_nodes']
print(f"n_loop={n_loop}, n_star={n_star}, n_nodes={n_nodes}, build={dt_build:.1f} s")

# ============================================================
# Extract matrices
# ============================================================
L_mat = np.array(topo['L']).reshape(n_loop, n_loop)
R_dc = np.array(topo['R'])
P_mat = np.array(topo['P']).reshape(n_star, n_star)
seg_nodes = np.array(topo['segment_nodes'])  # (n_loop, 2)

L_dc = np.sum(L_mat)
print(f"L_DC = {L_dc*1e6:.2f} uH")

# ============================================================
# Build gathering matrix G (n_nodes x n_panels)
# ============================================================
# Each panel's charge is split 50/50 between the two endpoint nodes
# of its parent segment.
print("\nBuilding gathering matrix G...")
G = np.zeros((n_nodes, n_star))
for s in range(n_seg):
	nf = seg_nodes[s, 0]
	nt = seg_nodes[s, 1]
	for k in range(panels_per_seg):
		p = s * panels_per_seg + k
		G[nf, p] = 0.5
		G[nt, p] = 0.5

print(f"G shape: {G.shape}, non-zeros: {np.count_nonzero(G)}")

# ============================================================
# Apply dielectric: half-space image method
# ============================================================
# For charges on a dielectric half-space (FR4 below, air above):
#   P_eff = P_free x (1 - alpha)  for self/mutual at same height
#   where alpha = -(eps_r - 1) / (eps_r + 1)
#
# For EMQS quasi-static image:
#   G_image(r, r') = alpha / |r - r'_image|
#   where r'_image is reflected through substrate surface
#
# For a flat spiral (all conductors at same z ~= 0):
#   Image distance ~= original distance (since conductor thickness << lateral dims)
#   So P_image ~= P_free, and P_total = P_free x (1 + alpha)
#
# Effective permittivity for charges on interface: eps_eff = (1 + eps_r) / 2
# This is equivalent to P_total = P_free / eps_eff

eps_eff = (1.0 + EPS_R) / 2.0
print(f"Dielectric model: half-space image, eps_eff = (1+{EPS_R})/2 = {eps_eff:.2f}")

P_eff = P_mat / eps_eff

# ============================================================
# Compute C_eff = G x P_eff⁻¹ x G^T (n_nodes x n_nodes)
# ============================================================
print("Computing C_eff = G * inv(P) * G^T ...")
t0 = time.perf_counter()

# Solve P_eff x X = G^T for X  (more stable than explicit inversion)
X = np.linalg.solve(P_eff, G.T)  # (n_star x n_nodes)
C_eff = G @ X  # (n_nodes x n_nodes)

dt_ceff = time.perf_counter() - t0
print(f"C_eff computed in {dt_ceff:.1f} s, shape={C_eff.shape}")

C_total_eff = np.sum(C_eff)
print(f"C_total_eff = {C_total_eff*1e12:.2f} pF")

# Lumped SRF estimate
if C_total_eff > 0 and L_dc > 0:
	f_srf_lumped = 1.0 / (2*np.pi*np.sqrt(L_dc * C_total_eff))
	print(f"Lumped SRF estimate: {f_srf_lumped/1e6:.1f} MHz")

C_target = 1.0 / ((2*np.pi*5e6)**2 * L_dc)
print(f"C needed for SRF=5MHz: {C_target*1e12:.1f} pF")

# ============================================================
# Build incidence matrix A (n_nodes x n_loop)
# ============================================================
A_inc = np.zeros((n_nodes, n_loop))
for s in range(n_loop):
	nf = seg_nodes[s, 0]
	nt = seg_nodes[s, 1]
	A_inc[nf, s] = +1.0  # current leaves node_from
	A_inc[nt, s] = -1.0  # current enters node_to

# ============================================================
# Proper MNA frequency sweep
# ============================================================
Zs_func = make_dowell_Zs_func(
	n_loop, R_dc, h, sigma, w, pitch, p=21
)

# Port info
port = topo['ports'][0]
node_pos = port[0]
node_neg = port[1]

# Ground the negative terminal node
# Build reduced system by removing the ground node row/col
ground_node = node_neg
all_nodes = list(range(n_nodes))
reduced_nodes = [n for n in all_nodes if n != ground_node]
n_reduced = len(reduced_nodes)

# Map from original node index to reduced index
node_to_reduced = {}
for ri, ni in enumerate(reduced_nodes):
	node_to_reduced[ni] = ri

# Reduced matrices
A_red = A_inc[reduced_nodes, :]  # (n_reduced x n_loop)
C_red = C_eff[np.ix_(reduced_nodes, reduced_nodes)]  # (n_reduced x n_reduced)

# Port excitation: inject 1A at node_pos (in reduced indexing)
I_port = np.zeros(n_reduced, dtype=complex)
I_port[node_to_reduced[node_pos]] = 1.0

# Frequency sweep
freqs = np.logspace(2, 7.7, 50)
print(f"\nMNA Frequency sweep ({len(freqs)} points)...")
t0 = time.perf_counter()

Z_port = np.zeros(len(freqs), dtype=complex)
for fi, freq in enumerate(freqs):
	omega = 2.0 * np.pi * freq

	# Surface impedance
	Zs = Zs_func(freq)

	# Z_branch = diag(R_dc + Zs) + jomega L
	Z_branch = np.diag(R_dc + Zs) + 1j * omega * L_mat

	# Y_branch = Z_branch^{-1}
	Y_branch = np.linalg.inv(Z_branch)

	# Y_node = A^T x Y_branch x A + jomega x C_eff  (reduced system)
	# Note: A_red is (n_reduced x n_loop), so A_red @ Y_branch @ A_red.T is correct
	Y_node = A_red @ Y_branch @ A_red.T + 1j * omega * C_red

	# Solve Y_node x V = I_port
	V_node = np.linalg.solve(Y_node, I_port)

	# Z_port = V_pos (since I = 1A, ground at neg)
	Z_port[fi] = V_node[node_to_reduced[node_pos]]

dt_sweep = time.perf_counter() - t0
print(f"Sweep time: {dt_sweep:.1f} s ({dt_sweep/len(freqs):.2f} s/freq)")

L_eff = np.imag(Z_port) / (2*np.pi*freqs)
R_eff = np.real(Z_port)

# Find SRF (where Im(Z) crosses zero from positive to negative)
srf_freqs = []
for i in range(len(freqs) - 1):
	if np.imag(Z_port[i]) > 0 and np.imag(Z_port[i+1]) < 0:
		# Linear interpolation
		f1, f2 = freqs[i], freqs[i+1]
		z1, z2 = np.imag(Z_port[i]), np.imag(Z_port[i+1])
		f_srf = f1 - z1 * (f2 - f1) / (z2 - z1)
		srf_freqs.append(f_srf)

if srf_freqs:
	print(f"\n*** SRF detected: {srf_freqs[0]/1e6:.2f} MHz ***")
else:
	print("\nNo SRF detected in frequency range")

print(f"\n  {'freq':>10s}  {'L [uH]':>8s}  {'R [Ohm]':>10s}  {'Im(Z)':>10s}")
for f_target in [100, 1e3, 1e4, 1e5, 5e5, 1e6, 2e6, 5e6, 8e6, 1e7, 2e7, 5e7]:
	idx = np.argmin(np.abs(freqs - f_target))
	print(f"  {freqs[idx]:10.0f}  {L_eff[idx]*1e6:8.3f}  {R_eff[idx]:10.4f}  {np.imag(Z_port[idx]):10.1f}")

# ============================================================
# Plot comparison with measurement
# ============================================================
try:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import openpyxl

	# Load measurement
	wb = openpyxl.load_workbook(os.path.join(SCRIPT_DIR, 'provided_by_sato', '14tAir_周波数特性測定値.xlsx'), data_only=True)
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

	fig, axes = plt.subplots(1, 2, figsize=(14, 5))

	srf_str = f', SRF={srf_freqs[0]/1e6:.1f}MHz' if srf_freqs else ''

	# (a) Inductance
	label_peec = f'PEEC MNA (mode={PANEL_MODE}, eps_eff={eps_eff:.1f})'
	axes[0].semilogx(freqs, L_eff * 1e6, 'b-', linewidth=2, label=label_peec)
	mask_L = m_Ls > 0
	axes[0].semilogx(m_freqs[mask_L], m_Ls[mask_L] * 1e6, 'r-', alpha=0.7, linewidth=1.5, label='Measurement')
	axes[0].axhline(y=L_dc*1e6, color='gray', linestyle='--', alpha=0.5, label=f'L_DC={L_dc*1e6:.1f} uH')
	if srf_freqs:
		axes[0].axvline(x=srf_freqs[0], color='green', linestyle=':', alpha=0.7, label=f'SRF={srf_freqs[0]/1e6:.1f} MHz')
	axes[0].set_xlabel('Frequency [Hz]')
	axes[0].set_ylabel('Inductance [uH]')
	axes[0].set_title(f'L(f) -- Proper PEEC MNA{srf_str}')
	axes[0].legend()
	axes[0].grid(True, alpha=0.3)
	axes[0].set_ylim(bottom=-50, top=100)

	# (b) Resistance
	axes[1].loglog(freqs, np.abs(R_eff), 'b-', linewidth=2, label=f'{label_peec} + Dowell p=21')
	axes[1].loglog(m_freqs, m_Rs, 'r-', alpha=0.7, linewidth=1.5, label='Measurement')
	if srf_freqs:
		axes[1].axvline(x=srf_freqs[0], color='green', linestyle=':', alpha=0.7, label=f'SRF={srf_freqs[0]/1e6:.1f} MHz')
	axes[1].set_xlabel('Frequency [Hz]')
	axes[1].set_ylabel('|R(f)| [Ohm]')
	axes[1].set_title('R(f) -- Proper PEEC MNA + Dowell proximity')
	axes[1].legend()
	axes[1].grid(True, alpha=0.3)

	plt.tight_layout()
	out_path = os.path.join(SCRIPT_DIR, 'spiral_peec_mna.png')
	plt.savefig(out_path, dpi=150)
	plt.close()
	print(f"\nFigure saved: {out_path}")
except Exception as e:
	print(f"Plot error: {e}")

print("\nDone.")
