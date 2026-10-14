"""Transverse gradient (Gx) coil design via the stream function method.

The HARD 2D companion to ``demo_coil_design_gz.py``.  A transverse gradient
(target Bz = Gx*x) is NOT axisymmetric, so the surface current on the cylinder
(radius a) cannot be reduced to azimuthal rings.  The stream function lives on
the full 2D cylinder surface, psi(phi, z), and its equal-current contours form
a non-trivial "fingerprint" pattern that varies with phi -- the classic
fingerprint/saddle coil of MRI gradient design.

Pipeline (same machinery as the Gz demo, but 2D):

    target Bz = Gx*x  in a DSV
        ->  least-norm surface current  (ACA+)+TSVD pseudo-inverse
        ->  stream function psi(phi, z) on the cylinder surface grid
        ->  equal-current contours of psi  =  discrete wires
        ->  verify Bz over the DSV vs Gx*x

Basis
-----
At each surface node (phi_p, z_q) we place a small quad current LOOP on the
cylinder around that node: corners at (phi_p +- dphi/2, z_q +- dz/2) mapped to
3D (a*cos(phi), a*sin(phi), z).  A unit loop has unit "stream value"; a smooth
field of overlapping loops with per-node weight psi(phi,z) realises a surface
current K = n_hat x grad(psi) (interior segments of adjacent loops cancel,
leaving only the net circulation == the contour wires).  The matrix entry
A(i,j) = Bz at observation point i from unit loop j.

Field kernel
------------
A(i,j) and the final verification both use an EXACT numpy Biot-Savart segment
formula (the HACApK ``_seg_Hz`` form, copied from tests/test_stream_function.py)
-- NOT ``rad.ObjFlmCur``, which has a confirmed ~10x-wrong Bz bug at certain
loop orientations (phi=90 deg).  We work in Hz units (scale-free) for the solve
and convert to Tesla (B = MU0 * H) only in the verification.

Run:  python demo_coil_design_gx.py [--nphi 24] [--nz 40] [--nlevels 12]
"""
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

from radia.stream_function import aca_tsvd, pseudo_inverse_solve
from radia.biot_savart import h_segments_batch, MU0


# --------------------------------------------------------------------------
# Exact Biot-Savart segment formula (HACApK form) -- copied verbatim from
# tests/test_stream_function.py.  _seg_Hz returns Hz from one straight segment;
# we use it directly to build PLAIN (non-mirrored) quad loops on the cylinder.
# (The test's _loop_Hz additionally mirrors through z, which we do NOT want
# here.)  Verified to agree with radia.biot_savart.h_segments_batch to machine
# precision for a CCW-ordered quad loop, so the solve matrix and the
# verification field are bit-consistent.
# --------------------------------------------------------------------------
def _seg_Hz(O, P1, P2):
    eps = 1.0e-15
    d = P2 - P1
    R21 = float(d @ d)
    o1 = O - P1
    o2 = O - P2
    OR1 = np.sqrt(o1 @ o1)
    OR2 = np.sqrt(o2 @ o2)
    O121 = float(o1 @ d)
    O221 = float(o2 @ d)
    Rc12 = O121 / OR1
    Rc21 = -O221 / OR2
    L = o1 - O121 * d / R21
    L1 = float(L @ L) + eps
    f = (Rc12 + Rc21) / (4.0 * np.pi * L1 * R21)
    return (d[0] * L[1] - d[1] * L[0]) * f  # Hz


def _loop_Hz(obs, corners):
    """Hz at obs from a plain quad current loop (corners (4,3), CCW), unit I."""
    hz = 0.0
    for c in range(4):
        n = (c + 1) % 4
        hz += _seg_Hz(obs, corners[c], corners[n])
    return hz


# --------------------------------------------------------------------------
# Cylinder surface geometry
# --------------------------------------------------------------------------
def cyl_point(a, phi, z):
    return np.array([a * np.cos(phi), a * np.sin(phi), z])


def loop_corners(a, phi, z, dphi, dz):
    """4 corners of a quad loop on the cylinder around (phi, z), CCW.

    Order (phi-, z-), (phi+, z-), (phi+, z+), (phi-, z+).  On the OUTWARD-
    normal side of the cylinder this traverses CCW as seen from outside, which
    matches the +z reference loop sign used in the agreement check.
    """
    return np.array([
        cyl_point(a, phi - dphi / 2, z - dz / 2),
        cyl_point(a, phi + dphi / 2, z - dz / 2),
        cyl_point(a, phi + dphi / 2, z + dz / 2),
        cyl_point(a, phi - dphi / 2, z + dz / 2),
    ])


def make_dsv(dsv, n):
    """3D grid of observation points inside a sphere of radius dsv (the DSV)."""
    g = np.linspace(-dsv, dsv, n)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    r = np.linalg.norm(pts, axis=1)
    return pts[r <= dsv + 1e-12]


# --------------------------------------------------------------------------
# Contour psi(phi, z) -> wires (each wire = list of 3D segments, current = dI)
# --------------------------------------------------------------------------
def contour_wires(psi_zphi, phi_grid, z_grid, a, n_levels):
    """Equal-current contours of psi on the (phi, z) grid.

    psi_zphi has shape (n_z, n_phi).  We pad one extra phi column (periodic
    wrap) so contours can cross the phi seam, then run matplotlib's contour.
    Returns (wires, dI) where wires is a list of segment-lists; each segment is
    a (p1, p2) pair of 3D points, and each wire carries current dI.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lo, hi = float(psi_zphi.min()), float(psi_zphi.max())
    if hi <= lo:
        return [], 0.0
    dI = (hi - lo) / n_levels
    levels = lo + (np.arange(n_levels) + 0.5) * dI

    # Periodic wrap in phi: append first column at phi + 2pi.
    phi_ext = np.concatenate([phi_grid, [phi_grid[0] + 2.0 * np.pi]])
    psi_ext = np.column_stack([psi_zphi, psi_zphi[:, :1]])  # (n_z, n_phi+1)

    PHI, Z = np.meshgrid(phi_ext, z_grid)  # both (n_z, n_phi+1)

    fig = plt.figure()
    cs = plt.contour(PHI, Z, psi_ext, levels=levels)
    plt.close(fig)

    wires = []
    # Matplotlib API differs across versions: prefer allsegs.
    allsegs = getattr(cs, "allsegs", None)
    if allsegs is None:
        return [], dI
    for segs_at_level in allsegs:
        for poly in segs_at_level:
            if len(poly) < 2:
                continue
            pts3d = np.array([cyl_point(a, p[0], p[1]) for p in poly])
            segs = [(pts3d[k], pts3d[k + 1]) for k in range(len(pts3d) - 1)]
            if segs:
                wires.append(segs)
    return wires, dI


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nphi", type=int, default=24, help="surface nodes around phi")
    ap.add_argument("--nz", type=int, default=40, help="surface nodes along z")
    ap.add_argument("--nlevels", type=int, default=12, help="equal-current contour levels")
    ap.add_argument("--k", type=int, default=0, help="TSVD modes (0 = auto)")
    ap.add_argument("--ndsv", type=int, default=7, help="DSV grid points per axis")
    args = ap.parse_args()

    a = 0.15           # cylinder radius [m]
    L = 0.50           # cylinder length [m]
    Gx = 1.0           # target transverse gradient dBz/dx [arb units]
    dsv = 0.10         # DSV radius [m] (target region |r| <= dsv)

    n_phi, n_z = args.nphi, args.nz

    # --- surface node grid (phi, z) --------------------------------------
    phi_grid = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    z_grid = np.linspace(-L / 2, L / 2, n_z)
    dphi = 2.0 * np.pi / n_phi
    dz = L / (n_z - 1)

    # Basis loops: one per node, index j = q * n_phi + p  (z-major)
    corners_list = []
    for q, zq in enumerate(z_grid):
        for p, pp in enumerate(phi_grid):
            corners_list.append(loop_corners(a, pp, zq, dphi, dz))
    N = len(corners_list)

    # --- observation points (DSV) + transverse-gradient target ----------
    obs = make_dsv(dsv, args.ndsv)
    M = obs.shape[0]
    B_target = Gx * obs[:, 0]    # Bz = Gx * x  (Hz units; scale-free solve)

    # --- matrix entry A(i,j) = Hz at obs[i] from unit loop j -------------
    def entry(i, j):
        return _loop_Hz(obs[i], corners_list[j])

    modes = min(M, N) if args.k <= 0 else min(args.k, M, N)
    res = aca_tsvd(M, N, entry, modes=modes, kmax=min(M, N),
                   aca_eps=1.0e-8)
    k_use = res.modes if args.k <= 0 else min(args.k, res.modes)
    psi = pseudo_inverse_solve(res, B_target, k_mode=k_use)   # (N,) stream values

    # continuous-current residual (sanity: A psi should match Gx*x)
    A_psi = np.array([sum(entry(i, j) * psi[j] for j in range(N))
                      for i in range(M)])
    cont_res = float(np.linalg.norm(A_psi - B_target) /
                     (np.linalg.norm(B_target) + 1e-30))

    print("=== Cylindrical Gx (transverse) gradient coil"
          " (stream function + (ACA+)+TSVD) ===")
    print(f"cylinder a={a*1e3:.0f} mm, L={L*1e3:.0f} mm | surface grid"
          f" n_phi={n_phi} x n_z={n_z} -> N={N} loop basis")
    print(f"DSV radius={dsv*1e3:.0f} mm, M={M} target points, Gx={Gx} (Bz=Gx*x)")
    print(f"ACA+ k_aca = {res.k_aca} (of min(M,N)={min(M, N)}),  "
          f"TSVD modes used = {k_use}")
    print(f"continuous surface current: ||A psi - Gx x||/||Gx x|| = {cont_res:.3e}")

    # --- contour psi(phi, z) -> equal-current wires ----------------------
    psi_zphi = psi.reshape(n_z, n_phi)         # (n_z, n_phi) fingerprint
    wires, dI = contour_wires(psi_zphi, phi_grid, z_grid, a, args.nlevels)
    n_wire = len(wires)
    print(f"discretised into {n_wire} equal-current contour wires"
          f" (I/wire = {dI:.4g} per unit psi)")

    # --- verify the DISCRETE wire coil over the DSV ----------------------
    # B = MU0 * sum_wires h_segments_batch(wire, obs, current=dI).
    # The sign of dI per wire follows the contour orientation; matplotlib
    # orients each closed loop consistently with the field's gradient, so a
    # single global current dI reproduces the target up to an overall sign
    # which we fix by least-squares scaling below.
    Bz_wires = np.zeros(M)
    for segs in wires:
        H = h_segments_batch(segs, obs, current=dI)
        Bz_wires += MU0 * H[:, 2]

    # Fix the (scale-free) overall gain: the solve was done in Hz units, so
    # compare shapes via the best-fit scalar gain g minimizing ||g*Bz - Gx x||.
    denom = float(np.dot(Bz_wires, Bz_wires))
    g = float(np.dot(Bz_wires, B_target) / denom) if denom > 0 else 0.0
    Bz_fit = g * Bz_wires
    rms_err = float(np.linalg.norm(Bz_fit - B_target) /
                    (np.linalg.norm(B_target) + 1e-30))
    print(f"discrete wire coil: best-fit gain g = {g:.4g}")
    print(f"Gx field-match relative RMS error over DSV = {rms_err:.3e}")

    # single-stroke connection
    single_stroke = False
    print("single-stroke connection: not done (separate contour loops)")

    # --- plot: psi(phi,z) fingerprint + Bz-vs-x scatter ------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

        PHI, Z = np.meshgrid(np.degrees(phi_grid), z_grid * 1e3)
        cf = ax1.contourf(PHI, Z, psi_zphi, levels=20, cmap="RdBu_r")
        lo, hi = float(psi_zphi.min()), float(psi_zphi.max())
        if hi > lo:
            lv = lo + (np.arange(args.nlevels) + 0.5) * (hi - lo) / args.nlevels
            ax1.contour(PHI, Z, psi_zphi, levels=lv, colors="k", linewidths=0.7)
        fig.colorbar(cf, ax=ax1, label="psi (stream value)")
        ax1.set_xlabel(r"$\phi$ (deg)")
        ax1.set_ylabel("$z$ (mm)")
        ax1.text(0.02, 0.96, "(a)", transform=ax1.transAxes, fontsize=10)

        ax2.scatter(obs[:, 0] * 1e3, Bz_fit, s=14, alpha=0.6,
                    label="discrete coil Bz (gain-fit)")
        xs = np.linspace(obs[:, 0].min(), obs[:, 0].max(), 50)
        ax2.plot(xs * 1e3, Gx * xs, "k--", label="target Gx*x")
        ax2.set_xlabel("$x$ (mm)")
        ax2.set_ylabel("$B_z$ (arb.)")
        ax2.text(0.02, 0.96, "(b)", transform=ax2.transAxes, fontsize=10)
        ax2.legend()
        ax2.grid(alpha=0.3)

        out = HERE / "demo_coil_design_gx.png"
        fig.tight_layout()
        fig.savefig(out, dpi=110)
        print(f"Saved plot: {out.name}")
    except ImportError:
        print("(matplotlib not installed; skipped plot)")


if __name__ == "__main__":
    main()
