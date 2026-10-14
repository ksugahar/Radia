"""
spiral_peec_srf.py

RLCM frequency sweep for the spiral inductor using PEECBuilder's
include_star=True mode (Loop-Star decomposition with capacitance).

Results:
  - Node-based P matrix (disk approximation) gives C_total ~ 2.8 pF
  - Measured SRF ~ 5 MHz requires C ~ 60 pF (20x larger)
  - R(f) with Dowell proximity (p=21) matches measurement well up to ~2 MHz
  - L(f) shows correct increasing trend but SRF is beyond modeled range

Limitation:
  Node-based P uses isotropic disk potential (1/(4*pi*eps0*a)) which
  underestimates inter-turn and inter-layer capacitance in PCB geometries.
  Face-panel-based P (FastMaxwell approach) would be needed for accurate SRF.
"""

import numpy as np
import time
import os
import pandas as pd

from radia.peec_matrices import PEECBuilder
from radia.peec_topology import PEECCircuitSolver

MU_0 = 4.0 * np.pi * 1e-7


def make_dowell_Zs_func(n_loop, R_dc_per_edge, h, sigma, w, pitch, p=21):
	"""Full Dowell winding loss with skin + proximity effect.

	Fr = Delta * [M1(Delta) + (2/3)(p^2-1) * D1(Delta)]

	Delta = (h/delta) * sqrt(eta), eta = w/pitch
	p = effective number of layers (fitted to measurement)
	"""
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
# Build PEEC model with Loop-Star (include capacitance)
# ============================================================
# Face panel mode: 3=left+right (inter-turn), 4=top only (single-side)
PANEL_MODE = 3       # side face panels for inter-turn capacitance
EPS_R = 4.4          # FR4 relative permittivity

print(f"\nBuilding RLCM model (panel mode={PANEL_MODE}, eps_r={EPS_R})...")

builder = PEECBuilder()
node_ids = []
for p_pt in pts_m:
	nid = builder.add_node_at(p_pt[0], p_pt[1], p_pt[2], area=pitch*w)
	node_ids.append(nid)
for i in range(n_pts - 1):
	builder.add_connected_segment(
		node_ids[i], node_ids[i+1],
		w, h, sigma,
		cross_section_type=0, nwinc=1, nhinc=1
	)
builder.add_port(node_ids[0], node_ids[-1])

# Generate face panels from segment box surfaces
n_panels = builder.generate_face_panels(mode=PANEL_MODE, eps_r=EPS_R)
print(f"Generated {n_panels} face panels (mode={PANEL_MODE}, eps_r={EPS_R})")

t0 = time.perf_counter()
topo = builder.build_topology(include_star=True)
dt_build = time.perf_counter() - t0
print(f"n_loop={topo['n_loop']}, n_star={topo['n_star']}, build={dt_build:.1f} s")

# Capacitance estimate
L_mat = np.array(topo['L']) if isinstance(topo['L'], list) else topo['L']
L_dc = np.sum(L_mat) if L_mat.ndim == 2 else float('nan')
P_mat = np.array(topo['P']) if isinstance(topo['P'], list) else topo['P']
try:
	C_mat = np.linalg.inv(P_mat)
	C_total = np.sum(C_mat)
	f_srf_lumped = 1.0 / (2*np.pi*np.sqrt(abs(L_dc * C_total)))
	print(f"C_total = {C_total*1e12:.2f} pF (lumped SRF est: {f_srf_lumped/1e6:.1f} MHz)")
	C_target = 1.0 / ((2*np.pi*5e6)**2 * L_dc)
	print(f"C needed for SRF=5MHz: {C_target*1e12:.1f} pF (ratio: {C_target/C_total:.0f}x)")
except Exception:
	print("C_total: could not compute")

# ============================================================
# RLCM Frequency sweep
# ============================================================
solver = PEECCircuitSolver(topo)
Zs_func = make_dowell_Zs_func(
	topo['n_loop'], np.array(topo['R']), h, sigma, w, pitch, p=21
)

# 30 log-spaced points: 100 Hz to 50 MHz
freqs = np.logspace(2, 7.7, 30)
print(f"\nFrequency sweep ({len(freqs)} points)...")
t0 = time.perf_counter()
Z_port = solver.frequency_sweep(freqs, Zs_func)
dt_sweep = time.perf_counter() - t0
print(f"Sweep time: {dt_sweep:.1f} s ({dt_sweep/len(freqs):.1f} s/freq)")

L_eff = np.imag(Z_port) / (2*np.pi*freqs)
R_eff = np.real(Z_port)

print(f"\n  {'freq':>10s}  {'L [uH]':>8s}  {'R [Ohm]':>10s}")
for f_target in [100, 1e3, 1e4, 1e5, 5e5, 1e6, 2e6, 5e6, 8e6]:
	idx = np.argmin(np.abs(freqs - f_target))
	print(f"  {freqs[idx]:10.0f}  {L_eff[idx]*1e6:8.3f}  {R_eff[idx]:10.4f}")

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

	# (a) Inductance
	label_peec = f'PEEC panel (mode={PANEL_MODE}, eps_r={EPS_R})'
	axes[0].semilogx(freqs, L_eff * 1e6, 'b-', linewidth=2, label=label_peec)
	mask_L = m_Ls > 0
	axes[0].semilogx(m_freqs[mask_L], m_Ls[mask_L] * 1e6, 'r-', alpha=0.7, linewidth=1.5, label='Measurement')
	axes[0].axhline(y=L_dc*1e6, color='gray', linestyle='--', alpha=0.5, label=f'L_DC={L_dc*1e6:.1f} uH')
	axes[0].set_xlabel('Frequency [Hz]')
	axes[0].set_ylabel('Inductance [uH]')
	axes[0].set_title(f'L(f) -- Panel RLCM (C={C_total*1e12:.1f} pF)')
	axes[0].legend()
	axes[0].grid(True, alpha=0.3)
	axes[0].set_ylim(bottom=-50, top=100)

	# (b) Resistance
	axes[1].loglog(freqs, np.abs(R_eff), 'b-', linewidth=2, label=f'{label_peec} + Dowell p=21')
	axes[1].loglog(m_freqs, m_Rs, 'r-', alpha=0.7, linewidth=1.5, label='Measurement')
	axes[1].set_xlabel('Frequency [Hz]')
	axes[1].set_ylabel('|R(f)| [Ohm]')
	axes[1].set_title('R(f) -- Panel RLCM + Dowell proximity')
	axes[1].legend()
	axes[1].grid(True, alpha=0.3)

	plt.tight_layout()
	out_path = os.path.join(SCRIPT_DIR, 'spiral_peec_srf.png')
	plt.savefig(out_path, dpi=150)
	plt.close()
	print(f"\nFigure saved: {out_path}")
except Exception as e:
	print(f"Plot error: {e}")

print("\nDone.")
