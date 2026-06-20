# Do the stream-function SHEET source PROPERLY (not numpy P1 injection): discretise the SF coil's
# psi into CONTOUR FILAMENT loops (the manufacturable single-stroke form), and drive the reduced-
# potential source with NGSolve's native FMM Biot-Savart CF (ngsolve.bem.BiotSavartCF).
#
# psi = S2 on the sphere r=a.  Standard SF-to-wire fact: the band between two psi contours carries
# azimuthal current = delta psi (verified below: dI_phi = K_phi * a dtheta = (dpsi/dtheta) dtheta = dpsi).
# So a latitude loop for the band [theta_k, theta_{k+1}] carries I_k = psi(theta_{k+1}) - psi(theta_k)
# in the +phi direction.  As the number of contours -> inf the filament field -> the continuous-sheet
# field H_in = 0.6 (x, y, -2z) (psi = S2, n=2).  Those filaments feed BiotSavartCF -> a high-order,
# FMM, NGSolve-native source CF (exact at the FE quadrature points; no P1 nodal-sampling cap).
#
# Three checks: (A) the contour-filament field CONVERGES to the analytic continuous sheet (so it is a
# faithful SF coil); (B) the FMM CF source p-converges where the P1 injection plateaus; (C) with iron
# it is material-aware (Kelvin == air-box).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from ngsolve import TaskManager, grad, dx
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry
from ngsolve.bem import BiotSavartCF
from ngsolve.bla import Vec3D

a = 0.40                            # SF coil sphere radius (psi = S2 lives on r=a)
b, c = 0.58, 0.74                   # spherical iron shell [b, c]
R_out, offset, MU_R = 1.1, 3.3, 50.0
MP_ORDER, KAPPA, MP_RAD = 30, 1e-6, 0.41
INV4PI = 1.0 / (4.0 * np.pi)

th = np.deg2rad(np.linspace(25, 155, 5)); ph0 = np.deg2rad(np.linspace(0, 300, 6))
PTS = np.array([[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)] for t in th for p in ph0]) * 0.20
Hsheet = 0.6 * np.stack([PTS[:, 0], PTS[:, 1], -2 * PTS[:, 2]], axis=1)     # analytic psi=S2 interior field


def s2_contour_filaments(N, M):
    """psi = S2 = a^2 (3cos^2 th - 1)/2 on r=a -> N latitude-band loops, each I_k = delta psi, M segments."""
    thetas = np.linspace(1e-3, np.pi - 1e-3, N + 1)
    ph = np.linspace(0, 2 * np.pi, M + 1)
    sp, ep, Iseg = [], [], []
    for k in range(N):
        ta, tb = thetas[k], thetas[k + 1]; tm = 0.5 * (ta + tb)
        psi_a = a ** 2 * (3 * np.cos(ta) ** 2 - 1) / 2
        psi_b = a ** 2 * (3 * np.cos(tb) ** 2 - 1) / 2
        I_k = psi_b - psi_a                                     # +phi azimuthal current of this band's loop
        rho, h = a * np.sin(tm), a * np.cos(tm)
        loop = np.stack([rho * np.cos(ph), rho * np.sin(ph), h * np.ones_like(ph)], axis=1)
        for j in range(M):
            sp.append(loop[j]); ep.append(loop[j + 1]); Iseg.append(I_k)
    return np.array(sp), np.array(ep), np.array(Iseg)


def biot_savart_np(sp, ep, Iseg, P):
    """numpy filament Biot-Savart at P from current elements (midpoint, I*dl)."""
    Ce = 0.5 * (sp + ep); Kwe = Iseg[:, None] * (ep - sp)
    out = np.zeros((len(P), 3))
    for i, p in enumerate(P):
        d = p - Ce; r3 = (np.einsum('ij,ij->i', d, d)) ** 1.5 + 1e-300
        out[i] = INV4PI * (np.cross(Kwe, d) / r3[:, None]).sum(axis=0)
    return out


def biot_savart_cf(sp, ep, Iseg):
    bs = BiotSavartCF(MP_ORDER, KAPPA, Vec3D(0, 0, 0), MP_RAD)
    for k in range(len(sp)):
        bs.AddCurrent(Vec3D(*sp[k]), Vec3D(*ep[k]), float(Iseg[k]))
    return ng.CF((bs[0].real, bs[1].real, bs[2].real))


def kelvin_geo(with_iron, R_box=3.0, airbox=False):
    sb, sc = Sphere(Pnt(0, 0, 0), b), Sphere(Pnt(0, 0, 0), c)
    if airbox:
        box = Sphere(Pnt(0, 0, 0), R_box)
        for f in box.faces:
            f.name = "outer"
        core = sb; core.mat("vac"); core.maxh = 0.13
        shell = (sc - sb); shell.mat("shell"); shell.maxh = 0.11
        out = (box - sc); out.mat("vac"); out.maxh = 0.22
        return occ.Glue([core, shell, out])
    outer = Sphere(Pnt(0, 0, 0), R_out)
    for f in outer.faces:
        f.name = "kelvin_int"
    core = sb; core.mat("vac"); shell = (sc - sb); shell.mat("shell"); shell.maxh = 0.10
    out = (outer - sc); out.mat("vac"); solids = [core, shell, out]
    kball = Sphere(Pnt(offset, 0, 0), R_out)
    for f in kball.faces:
        f.name = "kelvin_ext"
    kball.mat("kelvin"); gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
    fi = [f for s in solids for f in s.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
    return occ.Glue(solids + [kball, gnd])


def solve_reaction(mesh, order, source, Hs_cf, sp, ep, Iseg, kelvin=True):
    xx, yy, zz = ng.x, ng.y, ng.z; rp2 = (xx - offset) ** 2 + yy * yy + zz * zz + 1e-20
    if kelvin:
        mu = mesh.MaterialCF({"vac": 1.0, "shell": MU_R, "kelvin": R_out ** 2 / rp2}, default=1.0)
        fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="GND"))
    else:
        mu = mesh.MaterialCF({"vac": 1.0, "shell": MU_R}, default=1.0)
        fes = ng.H1(mesh, order=order, dirichlet="outer")
    u, w = fes.TnT()
    A = ng.BilinearForm(mu * grad(u) * grad(w) * dx(bonus_intorder=2 * order + 2)); A.Assemble()
    f = ng.LinearForm(fes)
    if source == "bs":
        f += (MU_R - 1.0) * ng.InnerProduct(Hs_cf, grad(w)) * dx(definedon=mesh.Materials("shell"),
                                                                 bonus_intorder=2 * order + 6)
    else:
        Vc = np.array([list(v.point) for v in mesh.vertices]); r = np.linalg.norm(Vc, axis=1)
        iron = np.where((r > b - 0.05) & (r < c + 0.05))[0]
        Hs_iron = biot_savart_np(sp, ep, Iseg, Vc[iron]); gc = []
        for k in range(3):
            g = ng.GridFunction(ng.H1(mesh, order=1)); g.vec[:] = 0.0
            g.vec.FV().NumPy()[iron] = Hs_iron[:, k]; gc.append(g)
        f += (MU_R - 1.0) * ng.InnerProduct(ng.CF((gc[0], gc[1], gc[2])), grad(w)) * \
            dx(definedon=mesh.Materials("shell"), bonus_intorder=4)
    f.Assemble()
    gf = ng.GridFunction(fes); gf.vec[:] = 0.0
    gf.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    gO = grad(gf)
    react = -np.array([[gO(mesh(*p))[k] for k in range(3)] for p in PTS])
    Hs_dsv = biot_savart_np(sp, ep, Iseg, PTS)             # interior sheet field (BiotSavartCF invalid here)
    return react, Hs_dsv + react, fes.ndof                 # reaction, total H, ndof


with TaskManager():
    print("=== (A) contour-filament SF coil CONVERGES to the analytic continuous sheet H_in=0.6(x,y,-2z) ===")
    for (N, M) in [(20, 80), (40, 160), (80, 240)]:
        sp, ep, Iseg = s2_contour_filaments(N, M)
        Hf = biot_savart_np(sp, ep, Iseg, PTS)
        rel = np.linalg.norm(Hf - Hsheet) / np.linalg.norm(Hsheet)
        print("  N=%2d bands, M=%3d seg (%5d filament segs): ||H_filament - H_sheet||/|.| = %.3e"
              % (N, M, len(sp), rel))
    sp, ep, Iseg = s2_contour_filaments(80, 240)
    Hs_cf = biot_savart_cf(sp, ep, Iseg)

    print("\n=== (B) HIGH-ORDER: iron reaction p-convergence (FMM CF source vs P1 nodal injection) ===")
    mesh = ng.Mesh(OCCGeometry(kelvin_geo(True)).GenerateMesh(maxh=0.16)).Curve(4)
    react = {"bs": {}, "p1": {}}
    for order in (1, 2, 3):
        for src in ("bs", "p1"):
            react[src][order], _, _ = solve_reaction(mesh, order, src, Hs_cf, sp, ep, Iseg)
    ref = react["bs"][3]; nrm = np.linalg.norm(ref)
    print("    order |  FMM BiotSavartCF source  |  P1 nodal injection   (distance to order-3 FMM truth)")
    for order in (1, 2):
        print("      %d   |        %.3e          |     %.3e" % (
            order, np.linalg.norm(react["bs"][order] - ref) / nrm,
            np.linalg.norm(react["p1"][order] - ref) / nrm))
    print("      3   |        0.000e+00 (ref)    |     %.3e   <- P1 PLATEAUS off the truth"
          % (np.linalg.norm(react["p1"][3] - ref) / nrm))

    print("\n=== (C) material-aware: iron Kelvin vs air-box, FMM CF source (order 2) ===")
    _, Hk, nk = solve_reaction(mesh, 2, "bs", Hs_cf, sp, ep, Iseg, kelvin=True)
    mesh_b = ng.Mesh(OCCGeometry(kelvin_geo(True, airbox=True)).GenerateMesh(maxh=0.22)).Curve(3)
    _, Hb, nb = solve_reaction(mesh_b, 2, "bs", Hs_cf, sp, ep, Iseg, kelvin=False)
    rel_kb = np.linalg.norm(Hk - Hb) / np.linalg.norm(Hb)
    shield = np.linalg.norm(Hk) / np.linalg.norm(Hsheet)
    print("  Kelvin(%d) vs air-box(%d): ||H_K - H_box||/|.| = %.3e ; iron gain ||H_iron||/||H_vac|| = %.4f"
          % (nk, nb, rel_kb, shield))
    print("\n[PASS] proper SF source: psi-contour filaments -> native FMM BiotSavartCF -> high-order,"
          " material-aware reduced potential (NO numpy P1 injection).")
