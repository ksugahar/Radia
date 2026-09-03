"""FE-direct order-3 stream-function coil design on an ARBITRARY (curved) former.

This is the generalisation of the planar / cylindrical stream-function (SF)
coil design to a coil surface of ANY shape -- demonstrated on a SPHERE, the
case the pre-fixed basis-loop representation CANNOT do (no structured (phi, z)
grid on a sphere).  The motivation is manufacturability: design the SF directly
ON the real former (sphere / conformal bobbin / 3D-printed freeform) so the
single-stroke wire lies on that surface and is windable / printable.  The NMR
shim coil is the driving application (a coil around the magnet that cancels the
measured field inhomogeneity).

Pipeline:
  [1] mesh the SOLID former, put H1 order-3 psi on its BOUNDARY (surface FEM
      done robustly via `definedon=mesh.Boundaries`, `ds`, `.Trace()`).
  [2] general surface kernel: surface current K = n_hat x grad_s psi (reduces
      to the planar z_hat x grad psi).  Per target j,
        A[j,i] = (mu0/4pi) int_S [(n x grad_s phi_i) x (r_t - r')]_z
                 / |r_t - r'|^3 dS
  [3] H1 min-seminorm solve (smoothest psi exactly hitting the target field).
  [4] sample psi at boundary vertices -> griddata to a (phi, theta) grid ->
      contour -> single-stroke (latitude-ring spiral for axisymmetric targets)
      -> embed on the sphere.
  [5] manufacturability: real inter-turn spacing (>= conductor width).
  [6] optional single-current SHEET-METAL distortion (radial dr + tangential)
      to cancel the single-stroke connector field.

Validated (R=150 mm coil, DSV r=80 mm):
  uniform Bz : cres 3e-15, single-stroke 0.24 %
  Z2 shim    : cres 5e-14, single-stroke 4.3 %, spacing 10.5 mm,
               + sheet-metal distort -> 0.36 % at ~2 mm bend (one current)

Run: python demo_sphere_fe_direct.py --target z2 --distort
"""
import argparse
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import griddata
from netgen.occ import Sphere, Pnt, OCCGeometry
from ngsolve import (Mesh, H1, LinearForm, BilinearForm, grad, specialcf, CF,
                     x, y, z, ds, Cross, TaskManager, BND)

HERE = Path(__file__).resolve().parent
from demo_sf_to_peec_gx import make_dsv, contour_polylines_phi_z, bz_fast

MU0_4PI = 1.0e-7


# --------------------------------------------------------------------------
# [1][2] surface mesh + general Biot-Savart kernel  K = n_hat x grad_s psi
# --------------------------------------------------------------------------
def build_surface(radius, maxh, order):
    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), radius))
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        mesh.Curve(order)            # curve order == field order (isoparametric)
        fes = H1(mesh, order=order, definedon=mesh.Boundaries(".*"))
    fi = np.where(np.array(fes.FreeDofs()))[0]
    return mesh, fes, fi


def assemble_A(fes, obs):
    """A[j,i] = Bz at obs[j] from the i-th surface H1 basis function's current."""
    n = specialcf.normal(3)
    pos = CF((x, y, z))
    v = fes.TestFunction()
    Kv = Cross(n, grad(v).Trace())   # surface current of the test function
    A = np.zeros((len(obs), fes.ndof))
    with TaskManager():
        for j, tj in enumerate(obs):
            d = CF((float(tj[0]), float(tj[1]), float(tj[2]))) - pos
            r3 = (d[0] ** 2 + d[1] ** 2 + d[2] ** 2) ** 1.5
            Lf = LinearForm(fes)
            Lf += MU0_4PI * Cross(Kv, d)[2] / r3 * ds
            Lf.Assemble()
            A[j, :] = Lf.vec.FV().NumPy()
    return A


def h1_min_seminorm(fes, fi, A, B):
    """psi = argmin psi^T S psi  s.t. A psi = B  (S = surface H1 stiffness)."""
    u, v = fes.TnT()
    au = BilinearForm(fes, symmetric=True)
    au += grad(u).Trace() * grad(v).Trace() * ds + 1.0e-8 * u * v * ds
    with TaskManager():
        au.Assemble()
    S = np.array(au.mat.ToDense())[np.ix_(fi, fi)]
    Af = A[:, fi]
    SiAt = np.linalg.solve(S, Af.T)
    psi = np.zeros(fes.ndof)
    psi[fi] = SiAt @ np.linalg.solve(Af @ SiAt, B)
    return psi


# --------------------------------------------------------------------------
# [4] sample psi (boundary vertices) -> regular (phi, theta) grid
# --------------------------------------------------------------------------
def sample_to_grid(psi, mesh, radius, nph, nth, th_pad=0.18):
    bv = sorted({vv.nr for el in mesh.Elements(BND) for vv in el.vertices})
    pts = np.array([mesh.vertices[i].point for i in bv])
    THv = np.arccos(np.clip(pts[:, 2] / radius, -1, 1))
    PHv = np.arctan2(pts[:, 1], pts[:, 0]) % (2 * np.pi)
    phg = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    thg = np.linspace(th_pad, np.pi - th_pad, nth)
    PG, TG = np.meshgrid(phg, thg)
    P3 = np.concatenate([PHv - 2 * np.pi, PHv, PHv + 2 * np.pi])  # periodic phi
    grid = griddata((P3, np.tile(THv, 3)), np.tile(psi[bv], 3), (PG, TG),
                    method="cubic")
    grid = np.nan_to_num(grid, nan=float(np.nanmean(psi[bv])))
    return grid, phg, thg


def _dedupe(a, tol=1e-6):
    keep = np.concatenate([[True],
                           np.linalg.norm(np.diff(a, axis=0), axis=1) > tol])
    return a[keep]


def sphere_embed(phc, thc, radius):
    return np.column_stack([radius * np.sin(thc) * np.cos(phc),
                            radius * np.sin(thc) * np.sin(phc),
                            radius * np.cos(thc)])


def single_stroke(grid, phg, thg, radius, nlevels):
    """Spiral single-stroke for an axisymmetric (latitude-ring) coil: order the
    contours by theta, open each at the phi nearest the previous chain end
    (minimise the connector), embed on the sphere.  Returns (path, phi, theta,
    contour_id)."""
    polys, _, _ = contour_polylines_phi_z(grid, phg, thg, nlevels)
    ps = sorted([p for p in polys if len(p) >= 3], key=lambda p: np.mean(p[:, 1]))
    chain, cid, prev = [], [], None
    for k, p in enumerate(ps):
        if prev is None:
            i = int(np.argmin(p[:, 0]))
        else:
            i = int(np.argmin(np.abs(((p[:, 0] - prev + np.pi) % (2 * np.pi))
                                     - np.pi)))
        opened = _dedupe(np.vstack([p[i:], p[:i]]))
        chain.append(opened)
        cid.append(np.full(len(opened), k))
        prev = opened[-1, 0]
    chain = np.vstack(chain)
    cid = np.concatenate(cid)
    ph, th = chain[:, 0], chain[:, 1]
    path = sphere_embed(ph, th, radius)
    return path, ph, th, cid, len(ps)


def inter_turn_spacing(path, ph, cid, seam=0.4, sub=3):
    """Min distance between DIFFERENT turns, measured away from the phi=0 seam
    (the closed-circle endpoints + seam connection points are not real
    clearances)."""
    away = (ph > seam) & (ph < 2 * np.pi - seam)
    P, C = path[away][::sub], cid[away][::sub]
    if len(P) < 2:
        return float("nan")
    d2 = np.sum((P[:, None, :] - P[None, :, :]) ** 2, axis=2)
    d2[C[:, None] == C[None, :]] = np.inf
    return float(np.sqrt(d2.min())) * 1e3


# --------------------------------------------------------------------------
# [6] single-current SHEET-METAL distortion on the sphere
# --------------------------------------------------------------------------
def coil_distort_sphere(ph, th, radius, obs_fit, B_fit, obs_eval, B_eval,
                        n_grid=6, n_iter=8, lam_disp_rel=0.1, step=0.8,
                        fd=1.0e-4):
    """Bend the wire (radial dr + azimuthal + polar) with ONE current to cancel
    the single-stroke connector field.  Smooth (phi, theta) control grid
    (phi-periodic), Gauss-Newton + displacement-Tikhonov penalty."""
    nfc, ntc = n_grid, n_grid
    dfc = 2 * np.pi / nfc
    tcg = np.linspace(th.min(), th.max(), ntc)
    dtc = tcg[1] - tcg[0]
    nc = nfc * ntc

    def interp(C, phc, thc):
        ff = (phc % (2 * np.pi)) / dfc
        i0 = np.floor(ff).astype(int) % nfc
        i1 = (i0 + 1) % nfc
        af = ff - np.floor(ff)
        fz = np.clip((thc - tcg[0]) / dtc, 0, ntc - 1.001)
        j0 = fz.astype(int)
        az = fz - j0
        return (C[i0, j0] * (1 - af) * (1 - az) + C[i1, j0] * af * (1 - az)
                + C[i0, j0 + 1] * (1 - af) * az + C[i1, j0 + 1] * af * az)

    sin_th = np.sin(np.clip(th, 0.1, np.pi - 0.1))

    def deform(dofs):
        dr = interp(dofs[0:nc].reshape(nfc, ntc), ph, th)
        dsf = interp(dofs[nc:2 * nc].reshape(nfc, ntc), ph, th)
        dst = interp(dofs[2 * nc:3 * nc].reshape(nfc, ntc), ph, th)
        return sphere_embed(ph + dsf / (radius * sin_th), th + dst / radius,
                            radius + dr)

    def fitI(Bu):
        return float(np.dot(Bu, B_fit) / np.dot(Bu, Bu))

    den = np.linalg.norm(B_eval) + 1e-30

    def rms(w):
        Bf = bz_fast(w, 1.0, obs_fit)
        return 100.0 * float(np.linalg.norm(fitI(Bf) * bz_fast(w, 1.0, obs_eval)
                                            - B_eval) / den)

    ndof = 3 * nc
    dofs = np.zeros(ndof)
    p0 = deform(dofs)
    best = (rms(p0), p0, 0.0)
    for _ in range(int(n_iter)):
        w = deform(dofs)
        Bf = bz_fast(w, 1.0, obs_fit)
        base = fitI(Bf) * Bf
        r = B_fit - base
        J = np.empty((len(obs_fit), ndof))
        for d in range(ndof):
            d2 = dofs.copy()
            d2[d] += fd
            Bu = bz_fast(deform(d2), 1.0, obs_fit)
            J[:, d] = (fitI(Bu) * Bu - base) / fd
        JTJ = J.T @ J
        md = float(np.mean(np.diag(JTJ))) + 1e-30
        dofs += step * np.linalg.solve(
            JTJ + (0.03 + lam_disp_rel) * md * np.eye(ndof),
            J.T @ r - lam_disp_rel * md * dofs)
        w = deform(dofs)
        m = rms(w)
        if m < best[0]:
            best = (m, w.copy(), float(np.max(np.abs(dofs))) * 1e3)
    return best


# --------------------------------------------------------------------------
def target_field(obs, name):
    if name == "uniform":                     # l=1: uniform Bz
        return np.ones(len(obs))
    if name == "z2":                          # l=2,m=0 shim: 2z^2 - x^2 - y^2
        f = 2 * obs[:, 2] ** 2 - obs[:, 0] ** 2 - obs[:, 1] ** 2
        return f / np.max(np.abs(f))
    raise ValueError(f"unknown target {name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, default=0.15, help="coil sphere [m]")
    ap.add_argument("--dsv", type=float, default=0.08, help="DSV radius [m]")
    ap.add_argument("--order", type=int, default=3, help="H1 + curve order")
    ap.add_argument("--maxh", type=float, default=0.03)
    ap.add_argument("--target", choices=["uniform", "z2"], default="z2")
    ap.add_argument("--nlevels", type=int, default=14)
    ap.add_argument("--distort", action="store_true",
                    help="single-current sheet-metal distortion to cut the "
                         "single-stroke RMS")
    ap.add_argument("--distort-penalty", type=float, default=0.1)
    ap.add_argument("--cond-width", type=float, default=2.0,
                    help="conductor width [mm] for the manufacturability gate")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    print(f"=== FE-direct order-{args.order} SF coil on a SPHERE former "
          f"(R={args.radius*1e3:.0f}mm, DSV r={args.dsv*1e3:.0f}mm, "
          f"target={args.target}) ===")
    t0 = time.time()
    mesh, fes, fi = build_surface(args.radius, args.maxh, args.order)
    print(f"[1] sphere surface: H1 order={args.order} (isoparametric), "
          f"surface dofs={len(fi)}  [{time.time()-t0:.1f}s]")

    obs = make_dsv(args.dsv, 7)
    obs_e = make_dsv(args.dsv, 9)
    Bt = target_field(obs, args.target)
    Bte = target_field(obs_e, args.target)

    A = assemble_A(fes, obs)
    psi = h1_min_seminorm(fes, fi, A, Bt)
    cres = float(np.linalg.norm(A @ psi - Bt) / np.linalg.norm(Bt))
    print(f"[2][3] general surface kernel + H1 min-seminorm: "
          f"continuous-SF residual = {cres:.2e}")

    grid, phg, thg = sample_to_grid(psi, mesh, args.radius, 96, 60)
    path, ph, th, cid, ncirc = single_stroke(grid, phg, thg, args.radius,
                                             args.nlevels)
    Bu = bz_fast(path, 1.0, obs_e)
    Iw = float(np.dot(Bu, Bte) / np.dot(Bu, Bu))
    rms = 100.0 * float(np.linalg.norm(Iw * Bu - Bte) / np.linalg.norm(Bte))
    print(f"[4] single-stroke: {ncirc} contours, {len(path)} pts, "
          f"RMS over DSV = {rms:.2f}%")

    pitch = inter_turn_spacing(path, ph, cid)
    gate = "OK" if pitch >= args.cond_width else "TIGHT (needs spacing distort)"
    print(f"[5] manufacturability: real inter-turn spacing = {pitch:.1f} mm "
          f"(conductor width {args.cond_width:.1f} mm) -> {gate}")

    dist_path = None
    if args.distort:
        rms_d, dist_path, disp = coil_distort_sphere(
            ph, th, args.radius, obs, Bt, obs_e, Bte,
            lam_disp_rel=args.distort_penalty)
        print(f"[6] sheet-metal distortion (1 current): RMS {rms:.2f}% -> "
              f"{rms_d:.2f}%  (max bend {disp:.1f} mm, no extra feeds)")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            plt.rcParams["pdf.fonttype"] = 42
            from mpl_toolkits.mplot3d.art3d import Line3DCollection
            P = (dist_path if dist_path is not None else path) * 1e3
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111, projection="3d")
            uu, vv = np.mgrid[0:2 * np.pi:24j, 0:np.pi:16j]
            r3 = args.radius * 1e3
            ax.plot_wireframe(r3 * np.cos(uu) * np.sin(vv),
                              r3 * np.sin(uu) * np.sin(vv), r3 * np.cos(vv),
                              color="0.85", lw=0.3)
            seg = np.concatenate([P[:-1, None], P[1:, None]], axis=1)
            lc = Line3DCollection(seg, cmap="viridis", lw=1.0)
            lc.set_array(np.linspace(0, 1, len(seg)))
            ax.add_collection3d(lc)
            for s in "xyz":
                getattr(ax, f"set_{s}lim")(-r3 * 1.1, r3 * 1.1)
            out = HERE / "demo_sphere_fe_direct.png"
            fig.tight_layout()
            fig.savefig(out, dpi=120)
            print(f"Saved plot: {out.name}")
        except ImportError:
            print("(matplotlib not installed; skipped plot)")


if __name__ == "__main__":
    main()
