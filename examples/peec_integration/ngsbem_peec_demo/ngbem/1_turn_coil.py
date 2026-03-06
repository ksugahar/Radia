"""
1_turn_coil.py

ngbem 単体でコイルのインダクタンスを計算するデモ。

Method:
  LaplaceSL (vector single-layer potential) on HDivSurface function space
  + Hodge decomposition で harmonic (loop) current mode を抽出し、
  BEM 二次形式から L を算出する。

Part 1: 単巻円形ループ
  - BEM 計算 + 収束確認
  - 解析式 L = mu_0 R [ln(8R/w) - 0.5] との比較

Part 2: N ターンソレノイド（空芯コイル）
  - 自己インダクタンス: BEM (Part 1 の結果)
  - 相互インダクタンス: Neumann 公式（完全楕円積分）
  - Wheeler 近似式 / 長岡係数との比較

Part 3: 導体遮蔽 (SIBC)
  - ShieldBEMSIBC (ngbem_eddy.py) によるアルミ板遮蔽
  - フィラメントモデルでコイル-遮蔽板カップリング
  - 周波数依存の L_eff, R_eff

Requires: ngsolve, ngsolve-ngsbem, scipy, numpy
"""

import sys
import os
import numpy as np
import time
from scipy.linalg import null_space
from scipy.special import ellipk, ellipe

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from ngbem_eddy import ShieldBEMSIBC

MU_0 = 4.0 * np.pi * 1e-7   # [H/m]


# ===========================================================
# Matrix extraction utilities
# ===========================================================

def _extract_dense(mat, n):
    """NGSolve BaseMatrix -> dense numpy array (n x n)."""
    ei = mat.CreateColVector()
    col = mat.CreateColVector()
    ei[:] = 0; ei[0] = 1.0
    mat.Mult(ei, col)
    is_complex = isinstance(col[0], complex)
    M = np.zeros((n, n), dtype=complex if is_complex else float)
    for j in range(n):
        M[j, 0] = col[j]
    for i in range(1, n):
        ei[:] = 0; ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(n):
            M[j, i] = col[j]
    return M


def _extract_rect(mat, nrow, ncol):
    """NGSolve BaseMatrix -> dense numpy array (nrow x ncol)."""
    ei = mat.CreateRowVector()
    col = mat.CreateColVector()
    M = np.zeros((nrow, ncol))
    for i in range(ncol):
        ei[:] = 0; ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(nrow):
            M[j, i] = col[j]
    return M


# ===========================================================
# Geometry: circular ring (annulus) surface mesh
# ===========================================================

def create_circular_ring_mesh(R_center, trace_width, maxh,
                               z=0.0, label="conductor"):
    """Create annular ring surface mesh using Netgen OCC.

    Args:
        R_center:    Center radius of the ring trace [m]
        trace_width: Radial width of the ring [m]
        maxh:        Maximum mesh element size [m]
        z:           Z-coordinate of the ring plane [m]
        label:       Boundary label for BEM

    Returns:
        NGSolve Mesh (surface triangulation in 3D)
    """
    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import Mesh

    R_out = R_center + trace_width / 2
    R_in  = R_center - trace_width / 2

    wp_o = WorkPlane(Axes(Pnt(0, 0, z), Dir(0, 0, 1), Dir(1, 0, 0)))
    outer = wp_o.Circle(R_out).Face()

    wp_i = WorkPlane(Axes(Pnt(0, 0, z), Dir(0, 0, 1), Dir(1, 0, 0)))
    inner = wp_i.Circle(R_in).Face()

    ring = outer - inner
    ring.faces.name = label

    geo = OCCGeometry(ring)
    return Mesh(geo.GenerateMesh(maxh=maxh))


# ===========================================================
# Analytical formulas
# ===========================================================

def analytical_L_circular_loop(R, w):
    """External inductance of circular loop (flat tape, width w).

    L = mu_0 * R * [ln(8R/w) - 1/2]

    GMD of a line segment of length w from itself: ln(GMD) = ln(w) - 3/2
    => L = mu_0 R [ln(8R/(w*exp(-3/2))) - 2]
         = mu_0 R [ln(8R/w) - 1/2]
    """
    return MU_0 * R * (np.log(8.0 * R / w) - 0.5)


def neumann_M_coaxial(R1, R2, d):
    """Mutual inductance of two coaxial circular loops (Neumann formula).

    M = mu_0 * sqrt(R1*R2) * [(2/k - k)*K(k^2) - (2/k)*E(k^2)]
    where k^2 = 4*R1*R2 / ((R1+R2)^2 + d^2)
    """
    k2 = 4.0 * R1 * R2 / ((R1 + R2)**2 + d**2)
    k = np.sqrt(k2)
    K_val = ellipk(k2)
    E_val = ellipe(k2)
    return MU_0 * np.sqrt(R1 * R2) * ((2.0/k - k) * K_val - (2.0/k) * E_val)


def wheeler_L_solenoid(N, R, length):
    """Wheeler's approximation for single-layer air-core solenoid.

    L = mu_0 * N^2 * pi * R^2 / (length + 0.9 * R)
    """
    return MU_0 * N**2 * np.pi * R**2 / (length + 0.9 * R)


# ===========================================================
# Conductor shielding (SIBC) helpers
# ===========================================================

def discretize_ring_coil(R_center, N_seg=60, z=0.0):
    """Discretize circular ring centerline into N filament segments.

    Args:
        R_center: Coil center radius [m]
        N_seg:    Number of segments
        z:        Z-coordinate of the coil plane [m]

    Returns:
        topo_dict: {'segment_centers', 'segment_directions', 'segment_lengths'}
    """
    theta = np.linspace(0, 2 * np.pi, N_seg + 1)
    x = R_center * np.cos(theta)
    y = R_center * np.sin(theta)

    centers    = np.zeros((N_seg, 3))
    directions = np.zeros((N_seg, 3))
    lengths    = np.zeros(N_seg)

    for i in range(N_seg):
        p0 = np.array([x[i], y[i], z])
        p1 = np.array([x[i + 1], y[i + 1], z])
        diff = p1 - p0
        seg_len = np.linalg.norm(diff)
        centers[i]    = 0.5 * (p0 + p1)
        directions[i] = diff / seg_len
        lengths[i]    = seg_len

    return {
        'segment_centers':    centers,
        'segment_directions': directions,
        'segment_lengths':    lengths,
    }


def create_shield_plate_mesh(plate_size, plate_thickness, gap, maxh):
    """Create a rectangular conductor plate mesh (3D volume) below z=0.

    Args:
        plate_size:      Side length of square plate [m]
        plate_thickness: Plate thickness [m]
        gap:             Air gap from coil plane (z=0) to plate top [m]
        maxh:            Maximum mesh element size [m]

    Returns:
        NGSolve Mesh (3D volume mesh with boundary surface)
    """
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh

    half  = plate_size / 2
    z_top = -gap
    z_bot = z_top - plate_thickness

    block = Box(Pnt(-half, -half, z_bot), Pnt(half, half, z_top))
    block.solids.name = "conductor"
    block.faces.name  = "surface"

    geo = OCCGeometry(block)
    return Mesh(geo.GenerateMesh(maxh=maxh))


# ===========================================================
# BEM inductance computation (ngbem standalone)
# ===========================================================

def compute_loop_inductance(mesh, R_center, trace_width,
                             sigma, thickness, label="conductor"):
    """Compute loop self-inductance using LaplaceSL BEM.

    Steps:
      1. HDivSurface (RWG edge basis) + SurfaceL2 (face basis) 構築
      2. Divergence 行列 D, 接続行列 C 構築
      3. Hodge 分解で harmonic (循環) mode c_h を抽出
      4. LaplaceSL on HDivSurface で vector single layer V_A を組立
      5. L/R ratio (正規化非依存) から L を算出

    Args:
        mesh:        NGSolve surface mesh (annular ring)
        R_center:    Center radius [m]
        trace_width: Ring width [m]
        sigma:       Conductivity [S/m]
        thickness:   Conductor thickness [m]
        label:       Boundary label

    Returns:
        L:    Self-inductance [H]
        info: Diagnostic dictionary
    """
    from ngsolve import (HDivSurface, SurfaceL2, BilinearForm,
                         ds, BND, TaskManager, InnerProduct)
    from ngsolve import div as ng_div
    from ngsolve.bem import LaplaceSL

    # --- Function spaces ---
    fes_J  = HDivSurface(mesh, order=0)    # edge (RWG) basis -> loop DOFs
    fes_L2 = SurfaceL2(mesh, order=0)      # face (piecewise const) -> star DOFs
    n_J = fes_J.ndof
    n_v = mesh.nv
    n_f = fes_L2.ndof

    # --- D: divergence matrix (n_face x n_edge) ---
    u_J  = fes_J.TrialFunction()
    q_L2 = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += ng_div(u_J.Trace()) * q_L2 * ds
    bf_D.Assemble()
    D = _extract_rect(bf_D.mat, n_f, n_J)

    # --- C: signed vertex-edge incidence (n_edge x n_vertex) ---
    C = np.zeros((n_J, n_v))
    for e_idx, edge in enumerate(mesh.edges):
        verts = list(edge.vertices)
        C[e_idx, verts[0].nr] = -1
        C[e_idx, verts[1].nr] = +1

    # --- M_J: edge mass matrix ---
    u2, v2 = fes_J.TnT()
    bf_M = BilinearForm(fes_J)
    bf_M += InnerProduct(u2.Trace(), v2.Trace()) * ds
    bf_M.Assemble()
    M_J = np.real(_extract_dense(bf_M.mat, n_J))

    # --- Harmonic mode: kernel of [D; C^T M_J] ---
    #   - D c_h = 0        (divergence-free: no charge accumulation)
    #   - C^T M_J c_h = 0  (curl-free in graph sense: pure circulation)
    #   For genus-1 surface (ring), dim(ker) = 1
    constraint = np.vstack([D, C.T @ M_J])
    null_h = null_space(constraint, rcond=1e-10)

    if null_h.shape[1] == 0:
        raise RuntimeError(
            "No harmonic mode found. Geometry must have a hole (genus >= 1).")
    if null_h.shape[1] > 1:
        print(f"  Warning: {null_h.shape[1]} harmonic modes found "
              f"(expected 1 for single ring)")

    c_h = null_h[:, 0]
    energy = c_h @ M_J @ c_h

    # --- V_A: LaplaceSL on HDivSurface (vector single layer) ---
    #   V_A[i,j] = int int J_i(r) . J_j(r') / (4 pi |r-r'|) dS dS'
    #   NOTE: .Trace() is required for correct BEM integration on
    #         surface meshes with boundary edges.
    j_trial = fes_J.TrialFunction()
    j_test  = fes_J.TestFunction()
    with TaskManager():
        V_op = LaplaceSL(
            j_trial.Trace() * ds(label)
        ) * j_test.Trace() * ds(label)

    V_A = np.real(_extract_dense(V_op.mat, n_J))
    V_A_proj = c_h @ V_A @ c_h

    # --- Inductance via L/R ratio (normalization-independent) ---
    #   L/R = mu_0 * V_A_proj / (R_sheet * energy)
    #   L   = (L/R) * R_loop
    R_sheet   = 1.0 / (sigma * thickness)
    perimeter = 2.0 * np.pi * R_center
    R_loop    = R_sheet * perimeter / trace_width

    LR_ratio = MU_0 * V_A_proj / (R_sheet * energy)
    L = LR_ratio * R_loop

    info = {
        'n_edges':    n_J,
        'n_faces':    n_f,
        'n_vertices': n_v,
        'n_harmonic': null_h.shape[1],
        'div_check':  np.linalg.norm(D @ c_h),
        'curl_check': np.linalg.norm(C.T @ M_J @ c_h),
        'energy':     energy,
        'V_A_proj':   V_A_proj,
        'R_sheet':    R_sheet,
        'R_loop':     R_loop,
    }
    return L, info


# ===========================================================
# Main
# ===========================================================

def main():
    # ---- Parameters ----
    R = 10e-3           # 10 mm coil radius (center of trace)
    w = 1e-3            # 1 mm trace width
    t = 35e-6           # 35 um copper thickness
    sigma = 5.8e7       # Cu conductivity [S/m]

    print("=" * 60)
    print(" Air-core coil inductance -- ngbem standalone")
    print("=" * 60)

    # ========================================================
    # Part 1: Single circular loop
    # ========================================================
    print("\n--- Part 1: Single circular loop ---")
    print(f"  R = {R*1e3:.1f} mm,  w = {w*1e3:.1f} mm,  t = {t*1e6:.0f} um Cu")

    L_ana = analytical_L_circular_loop(R, w)
    print(f"\n  Analytical: L = mu_0*R*[ln(8R/w) - 0.5]"
          f" = {L_ana*1e9:.2f} nH")

    # BEM convergence study
    print(f"\n  BEM convergence:")
    print(f"  {'maxh':>8s}  {'elem':>6s}  {'edges':>6s}  "
          f"{'L [nH]':>10s}  {'err [%]':>8s}  {'time':>8s}")
    print("  " + "-" * 58)

    L_best = None
    for maxh in [2.0e-3, 1.5e-3, 1.0e-3]:
        t0 = time.perf_counter()
        mesh = create_circular_ring_mesh(R, w, maxh, label="conductor")
        L_bem, info = compute_loop_inductance(
            mesh, R, w, sigma, t, label="conductor")
        dt = time.perf_counter() - t0

        err = abs(L_bem - L_ana) / L_ana * 100
        print(f"  {maxh*1e3:7.1f}mm  {info['n_faces']:6d}  "
              f"{info['n_edges']:6d}  {L_bem*1e9:10.2f}  "
              f"{err:8.1f}  {dt*1e3:6.0f} ms")
        L_best = L_bem

    print(f"\n  Diagnostics (finest mesh):")
    print(f"    ||D*c_h||      = {info['div_check']:.1e}  (should be ~0)")
    print(f"    ||C^T*M*c_h||  = {info['curl_check']:.1e}  (should be ~0)")
    print(f"    harmonic modes  = {info['n_harmonic']}  (expected 1)")
    print(f"\n  Result:")
    print(f"    BEM:        L = {L_best*1e9:.2f} nH")
    print(f"    Analytical: L = {L_ana*1e9:.2f} nH")
    print(f"    DC resistance R = {info['R_loop']*1e3:.2f} mOhm")
    print(f"    L/R = {L_best / info['R_loop'] * 1e6:.2f} us")

    # ========================================================
    # Part 2: N-turn air-core solenoid
    # ========================================================
    N     = 10
    pitch = 2e-3               # 2 mm turn pitch
    coil_length = (N - 1) * pitch

    print(f"\n\n--- Part 2: {N}-turn air-core solenoid ---")
    print(f"  N = {N},  pitch = {pitch*1e3:.1f} mm,  "
          f"length = {coil_length*1e3:.1f} mm,  R = {R*1e3:.1f} mm")

    # Self-inductance per turn (from BEM)
    L_self = L_best
    print(f"\n  Self-inductance per turn (BEM): {L_self*1e9:.2f} nH")

    # Mutual inductance table (Neumann formula)
    print(f"\n  Mutual inductance (Neumann formula):")
    for k in range(1, min(5, N)):
        M_k = neumann_M_coaxial(R, R, k * pitch)
        print(f"    M(dn={k}): {M_k*1e9:.2f} nH  "
              f"(coupling = {M_k / L_self * 100:.1f}%)")

    # Total inductance: L = sum_ij M_ij
    L_total = 0.0
    for i in range(N):
        L_total += L_self
        for j in range(i + 1, N):
            L_total += 2 * neumann_M_coaxial(R, R, (j - i) * pitch)

    # Wheeler reference
    L_wheeler = wheeler_L_solenoid(N, R, coil_length)
    L_inf     = MU_0 * N**2 * np.pi * R**2 / coil_length

    print(f"\n  Total inductance:")
    print(f"    BEM + Neumann:  {L_total*1e6:.4f} uH")
    print(f"    Wheeler approx: {L_wheeler*1e6:.4f} uH")
    print(f"    Infinite sol.:  {L_inf*1e6:.4f} uH")
    print(f"    Error (vs Wheeler): "
          f"{abs(L_total - L_wheeler) / L_wheeler * 100:.1f}%")

    # Nagaoka coefficient
    K_N = L_total / L_inf
    print(f"\n  Nagaoka coefficient K_N = {K_N:.4f}")
    print(f"    (K_N = 1.0 for infinite solenoid, < 1 for finite)")

    # ========================================================
    # Part 3: Conductor shielding (BEM + SIBC)
    # ========================================================
    print(f"\n\n--- Part 3: Conductor shielding (SIBC) ---")

    # Coil discretized into filament segments
    N_seg = 60
    topo_dict = discretize_ring_coil(R, N_seg=N_seg, z=0.0)
    print(f"  Coil: {N_seg} filament segments, R = {R*1e3:.1f} mm")

    # Aluminum shield plate below coil
    plate_size = 30e-3         # 30 mm square
    plate_thickness = 2e-3     # 2 mm thick
    gap = 0.5e-3               # 0.5 mm air gap
    sigma_al = 3.7e7           # aluminum conductivity [S/m]

    print(f"  Shield: {plate_size*1e3:.0f}mm Al plate, "
          f"t={plate_thickness*1e3:.1f}mm, gap={gap*1e3:.1f}mm")
    print(f"  sigma = {sigma_al:.1e} S/m")

    # Create shield mesh and solver
    t0 = time.perf_counter()
    plate_mesh = create_shield_plate_mesh(
        plate_size, plate_thickness, gap, maxh=3e-3)
    shield = ShieldBEMSIBC(plate_mesh, sigma=sigma_al, mu_r=1.0)
    shield.assemble(intorder=4)
    dt_asm = time.perf_counter() - t0
    print(f"  Assembly: {dt_asm*1e3:.0f} ms")

    # Frequency sweep
    freqs_sibc = np.logspace(2, 6, 25)   # 100 Hz -- 1 MHz
    L_eff_arr = np.zeros(len(freqs_sibc))
    R_eff_arr = np.zeros(len(freqs_sibc))

    t0 = time.perf_counter()
    for k, f in enumerate(freqs_sibc):
        omega_k = 2.0 * np.pi * f
        Delta_Z = shield.compute_impedance_matrix(f, topo_dict)
        # Single-turn: all segments in series, same current
        Delta_Z_total = np.sum(Delta_Z)
        L_eff_arr[k] = L_best + np.imag(Delta_Z_total) / omega_k
        R_eff_arr[k] = info['R_loop'] + np.real(Delta_Z_total)
    dt_sweep = time.perf_counter() - t0
    print(f"  Freq sweep ({len(freqs_sibc)} pts): {dt_sweep:.1f} s")

    # Physical checks
    dL = L_eff_arr - L_best
    dR = R_eff_arr - info['R_loop']
    print(f"\n  Physical checks:")
    print(f"    dL < 0 (Lenz's law):  min dL = {np.min(dL)*1e9:+.2f} nH  "
          f"{'OK' if np.all(dL <= 1e-15) else 'WARN'}")
    print(f"    dR > 0 (eddy loss):   min dR = {np.min(dR)*1e3:+.4f} mOhm  "
          f"{'OK' if np.all(dR >= -1e-15) else 'WARN'}")

    # Summary table
    print(f"\n  {'freq':>10s}  {'L_eff [nH]':>10s}  {'dL [nH]':>10s}  "
          f"{'R_eff [mOhm]':>12s}  {'dR [mOhm]':>10s}")
    print("  " + "-" * 60)
    # Print ~6 representative points
    idx_show = np.linspace(0, len(freqs_sibc) - 1, 6, dtype=int)
    for idx in idx_show:
        f = freqs_sibc[idx]
        print(f"  {f:10.0f}  {L_eff_arr[idx]*1e9:10.2f}  "
              f"{dL[idx]*1e9:+10.2f}  {R_eff_arr[idx]*1e3:12.3f}  "
              f"{dR[idx]*1e3:+10.4f}")

    # Skin depth at representative frequencies
    print(f"\n  Skin depth (aluminum):")
    for f_ref in [1e3, 1e4, 1e5, 1e6]:
        delta_s = np.sqrt(2.0 / (2*np.pi*f_ref * MU_0 * sigma_al))
        print(f"    f={f_ref:.0e}: delta = {delta_s*1e3:.2f} mm"
              f"  (plate t={plate_thickness*1e3:.1f}mm)")

    # ========================================================
    # Optional: Plot
    # ========================================================
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # (a) Mutual inductance vs turn separation
        dn_arr = np.arange(1, N)
        M_vals = np.array([neumann_M_coaxial(R, R, k * pitch)
                           for k in dn_arr])
        axes[0, 0].stem(dn_arr, M_vals * 1e9)
        axes[0, 0].set_xlabel('Turn separation $\\Delta n$')
        axes[0, 0].set_ylabel('Mutual inductance [nH]')
        axes[0, 0].set_title('Mutual inductance vs turn separation')
        axes[0, 0].grid(True, alpha=0.3)

        # (b) Cumulative inductance vs number of turns
        L_cumul = []
        for n in range(1, N + 1):
            L_n = n * L_self
            for i in range(n):
                for j in range(i + 1, n):
                    L_n += 2 * neumann_M_coaxial(R, R, (j - i) * pitch)
            L_cumul.append(L_n)

        n_arr = np.arange(1, N + 1)
        axes[0, 1].plot(n_arr, np.array(L_cumul) * 1e6,
                     'o-', label='BEM + Neumann')
        axes[0, 1].plot(n_arr, [L_self * n**2 * 1e6 for n in n_arr],
                     '--', alpha=0.5, label='$N^2 \\times L_{self}$')
        L_wh = [wheeler_L_solenoid(n, R, max((n - 1) * pitch, pitch))
                for n in n_arr]
        axes[0, 1].plot(n_arr, np.array(L_wh) * 1e6,
                     's--', alpha=0.5, label='Wheeler')
        axes[0, 1].set_xlabel('Number of turns $N$')
        axes[0, 1].set_ylabel('Total inductance [$\\mu$H]')
        axes[0, 1].set_title('Solenoid inductance vs turns')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)

        # (c) SIBC: L_eff vs frequency
        axes[1, 0].semilogx(freqs_sibc, L_eff_arr * 1e9, 'b-o',
                             markersize=3, label='L_eff (with shield)')
        axes[1, 0].axhline(L_best * 1e9, color='gray', linestyle=':',
                            alpha=0.5, label='L_air')
        axes[1, 0].set_xlabel('Frequency [Hz]')
        axes[1, 0].set_ylabel('Inductance [nH]')
        axes[1, 0].set_title('Conductor shielding: L_eff(f)')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)

        # (d) SIBC: R_eff vs frequency
        axes[1, 1].loglog(freqs_sibc, R_eff_arr * 1e3, 'r-o',
                           markersize=3, label='R_eff (with shield)')
        axes[1, 1].axhline(info['R_loop'] * 1e3, color='gray',
                            linestyle=':', alpha=0.5, label='R_dc')
        axes[1, 1].set_xlabel('Frequency [Hz]')
        axes[1, 1].set_ylabel('Resistance [mOhm]')
        axes[1, 1].set_title('Conductor shielding: R_eff(f)')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = os.path.join(os.path.dirname(__file__),
                                '1_turn_coil.png')
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"\n  Figure saved: {out_path}")
    except ImportError:
        print("\n  (matplotlib not available, skipping plot)")

    print("\nDone.")


if __name__ == '__main__':
    main()
