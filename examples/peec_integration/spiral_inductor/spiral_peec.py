"""
spiral_peec.py

14ターン矩形スパイラルインダクタの PEEC ベースラインモデル。

point.xlsx の中心線データを読み込み、L(f) と R(f) を計算する。
静電容量モデルを含まないため SRF は再現しない (ベースライン比較用)。

導体厚さ: 測定 R_DC = 1.05 Ohm から逆算 (30.1 um)
抵抗モデル: 表皮効果のみ (coth モデル)

測定データ: 14tAir_周波数特性測定値.xlsx
  - L_DC = 17.39 uH, R_DC = 1.05 Ohm
  - 周波数範囲: 100 Hz -- 8 MHz
"""

import sys
import numpy as np
import time
import os
import pandas as pd

# Import radia first to set up MKL DLL directories on Windows
import radia  # noqa: F401 -- triggers DLL path setup in __init__.py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
from peec_matrices import PEECBuilder

MU_0 = 4.0 * np.pi * 1e-7

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# Read point.xlsx geometry
# ============================================================

def load_geometry():
    """Load centerline waypoints from point.xlsx."""
    df = pd.read_excel(
        os.path.join(SCRIPT_DIR, 'provided_by_sato', 'point.xlsx'), header=0)
    df3 = df.iloc[:, :3].apply(
        pd.to_numeric, errors='coerce').replace(
        [np.inf, -np.inf], np.nan).dropna()
    x_mm = df3.iloc[:, 0].to_numpy()
    y_mm = df3.iloc[:, 1].to_numpy()
    z_mm = df3.iloc[:, 2].to_numpy()

    # Remove duplicate consecutive points
    keep = [0]
    for i in range(1, len(x_mm)):
        dist = np.sqrt((x_mm[i]-x_mm[i-1])**2 +
                       (y_mm[i]-y_mm[i-1])**2 +
                       (z_mm[i]-z_mm[i-1])**2)
        if dist > 1e-10:
            keep.append(i)

    pts_mm = np.column_stack([x_mm[keep], y_mm[keep], z_mm[keep]])
    pts_m = pts_mm * 1e-3  # mm -> meters
    return pts_m


# ============================================================
# Skin effect model
# ============================================================

def skin_effect_resistance(R_dc, h, sigma, freqs):
    """Compute frequency-dependent resistance with skin effect (coth model).

    R_ac/R_dc = (h/2) / delta * coth((h/2) / delta)
    where delta = sqrt(2 / (omega * mu_0 * sigma))
    """
    t = h / 2.0  # half-thickness
    omega = 2.0 * np.pi * freqs
    R_ac = np.zeros_like(freqs)
    for i, om in enumerate(omega):
        if om > 0:
            delta = np.sqrt(2.0 / (om * MU_0 * sigma))
            x = t / delta
            if x < 0.01:
                R_ac[i] = R_dc
            else:
                R_ac[i] = R_dc * x / np.tanh(x)
        else:
            R_ac[i] = R_dc
    return R_ac


# ============================================================
# Load measurement data
# ============================================================

def load_measurements():
    """Load measurement data from Excel file."""
    try:
        import openpyxl
    except ImportError:
        print("  (openpyxl not available)")
        return None

    xlsx_path = os.path.join(SCRIPT_DIR,
                             'provided_by_sato', '14tAir_周波数特性測定値.xlsx')
    if not os.path.exists(xlsx_path):
        return None

    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb['Sheet1']

    freqs, Rs, Ls = [], [], []
    for row in ws.iter_rows(min_row=4, values_only=True):
        if row[1] is not None:
            try:
                freqs.append(float(row[1]))
                Rs.append(float(row[4]))
                Ls.append(float(row[5]))
            except (ValueError, TypeError):
                pass

    return {'freq': np.array(freqs), 'Rs': np.array(Rs), 'Ls': np.array(Ls)}


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 60)
    print("14-Turn Rectangular Spiral Inductor - PEEC Baseline Model")
    print("  Geometry: point.xlsx centerline (1653 nodes)")
    print("  Resistance: skin effect only (coth model)")
    print("=" * 60)

    # --- Load point.xlsx geometry ---
    pts_m = load_geometry()
    n_pts = len(pts_m)
    L_wire = np.sum(np.linalg.norm(np.diff(pts_m, axis=0), axis=1))
    print(f"\nGeometry: {n_pts} nodes, {n_pts-1} segments")
    print(f"  Wire length: {L_wire*1e3:.1f} mm")

    # --- Parameters ---
    w = 1.750e-3       # track width [m]
    pitch = 2.000e-3   # center-to-center [m]
    sigma = 5.8e7      # copper conductivity [S/m]
    rho = 1.0 / sigma
    R_meas_dc = 1.05   # measured R_DC [Ohm]
    L_meas = 17.388e-6 # measured L_DC [H]

    # Derive conductor thickness from measured R_DC
    h = rho * L_wire / (w * R_meas_dc)
    print(f"  Conductor thickness: h = {h*1e6:.1f} um "
          f"(derived from R_DC = {R_meas_dc} Ohm)")
    print(f"  Track width: {w*1e3:.3f} mm, Pitch: {pitch*1e3:.3f} mm")

    # --- Build PEEC model ---
    print(f"\nBuilding PEEC model (L matrix only, no capacitance)...")
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

    t0 = time.perf_counter()
    topo = builder.build_topology(include_star=False)
    dt_build = time.perf_counter() - t0

    n_loop = topo['n_loop']
    print(f"  n_loop = {n_loop}, build time = {dt_build:.1f} s")

    # --- Extract L matrix and compute L_DC ---
    L_mat = np.array(topo['L']).reshape(n_loop, n_loop)
    R_dc_vec = np.array(topo['R'])
    R_dc_total = np.sum(R_dc_vec)
    L_dc = np.sum(L_mat)

    # Inductance breakdown
    L_self = np.trace(L_mat)
    L_positive = np.sum(L_mat[L_mat > 0]) - L_self
    L_negative = np.sum(L_mat[L_mat < 0])

    L_err = (L_dc - L_meas) / L_meas * 100
    R_err = (R_dc_total - R_meas_dc) / R_meas_dc * 100

    print(f"\n--- Results ---")
    print(f"  L_DC  = {L_dc*1e6:.3f} uH  (meas: {L_meas*1e6:.3f} uH, "
          f"err: {L_err:+.1f}%)")
    print(f"  R_DC  = {R_dc_total:.4f} Ohm (meas: {R_meas_dc:.4f} Ohm, "
          f"err: {R_err:+.1f}%)")
    print(f"\n  Inductance breakdown:")
    print(f"    Self (diagonal):      {L_self*1e6:.3f} uH")
    print(f"    Same-direction (+): {L_positive*1e6:+.3f} uH")
    print(f"    Anti-parallel  (-): {L_negative*1e6:+.3f} uH")
    print(f"    Total:              {L_dc*1e6:.3f} uH")

    # --- Frequency sweep ---
    print(f"\n--- Frequency sweep (skin effect only) ---")
    freqs = np.logspace(2, np.log10(8e6), 200)
    L_eff = np.full_like(freqs, L_dc)
    R_eff = skin_effect_resistance(R_dc_total, h, sigma, freqs)

    print(f"  {'freq':>10s}  {'L [uH]':>8s}  {'R [Ohm]':>8s}  {'R/R_dc':>6s}")
    for f_target in [100, 1e3, 1e4, 1e5, 1e6, 5e6]:
        idx = np.argmin(np.abs(freqs - f_target))
        print(f"  {freqs[idx]:10.0f}  {L_eff[idx]*1e6:8.3f}  "
              f"{R_eff[idx]:8.4f}  {R_eff[idx]/R_dc_total:6.2f}")

    # --- Load measurement data ---
    meas_data = load_measurements()

    # --- Plot comparison ---
    plot_comparison(freqs, L_eff, R_eff, meas_data, L_dc, R_dc_total)

    print("\nDone.")


def plot_comparison(freqs, L_eff, R_eff, meas_data, L_dc, R_dc):
    """Plot PEEC vs measurement comparison."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not available)")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # (a) Inductance
    axes[0].semilogx(freqs, L_eff * 1e6, 'b-', linewidth=2,
                     label=f'PEEC point.xlsx ({L_dc*1e6:.2f} uH)')
    if meas_data is not None:
        mask = meas_data['Ls'] > 0
        axes[0].semilogx(meas_data['freq'][mask], meas_data['Ls'][mask] * 1e6,
                         'r-', alpha=0.7, linewidth=1.5,
                         label='Measurement (17.39 uH)')
    axes[0].set_xlabel('Frequency [Hz]')
    axes[0].set_ylabel('Inductance [uH]')
    axes[0].set_title('L(f) -- 14-turn spiral (point.xlsx, skin effect only)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(bottom=0, top=25)

    # (b) Resistance
    axes[1].loglog(freqs, R_eff, 'b-', linewidth=2,
                   label=f'PEEC skin (R_dc={R_dc:.3f} Ohm)')
    if meas_data is not None:
        axes[1].loglog(meas_data['freq'], meas_data['Rs'],
                       'r-', alpha=0.7, linewidth=1.5, label='Measurement')
    axes[1].set_xlabel('Frequency [Hz]')
    axes[1].set_ylabel('Resistance [Ohm]')
    axes[1].set_title('R(f) -- skin effect only (coth model)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = os.path.join(SCRIPT_DIR, 'spiral_peec.png')
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"\n  Figure saved: {out_path}")


if __name__ == '__main__':
    main()
