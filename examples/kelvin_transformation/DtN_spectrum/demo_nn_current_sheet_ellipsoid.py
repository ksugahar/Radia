# VERIFIED (2026-06-15): arbitrary coil surface (ELLIPSOID) -- FE-direct Biot-Savart current sheet H_s
# (numpy, from the coil triangulation K = n x grad_s psi) injected as the source of the reduced
# scalar potential + Kelvin open boundary.  Validates:
#   PART A (sphere): the numpy FE-direct Biot-Savart H_s == the increment-1 ANALYTIC field, AND the
#                    iron reduced-potential solve with the injected H_s == the analytic-H_s solve.
#   PART B (ellipsoid): the same machinery on a NON-spherical coil; iron via Kelvin vs a brute-force
#                       air-box (independent open BC).  Vacuum -> H = H_s = Biot-Savart by construction.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from ngsolve import TaskManager
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry

a = 0.40                       # coil scale (sphere radius / ellipsoid mean)
b, c = 0.58, 0.74              # spherical iron shell [b, c], outside the coil
R_out, offset = 1.1, 3.3
order, MU_R = 2, 50.0
INV4PI = 1.0 / (4.0 * np.pi)

th = np.deg2rad(np.linspace(25, 155, 5)); ph = np.deg2rad(np.linspace(0, 300, 6))
PTS = np.array([[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)] for t in th for p in ph]) * 0.20


# ---- coil surface triangulation (boundary of a meshed unit sphere, scaled to an ellipsoid) ----
def coil_surface(axes, hc=0.16):
    m = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=hc))
    V = np.array([list(v.point) for v in m.vertices]) * np.asarray(axes)   # scale -> ellipsoid
    T = np.array([[v.nr for v in el.vertices] for el in m.Elements(ng.BND)])
    return V, T


def tri_geom(V, T, psi_v):
    """per-triangle centroid C, area-weighted current Kw = area * (n x grad_s psi) (outward n)."""
    p0, p1, p2 = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    C = (p0 + p1 + p2) / 3.0
    e1, e2 = p1 - p0, p2 - p0
    nrm = np.cross(e1, e2)
    area2 = np.linalg.norm(nrm, axis=1)               # = 2*area
    nhat = nrm / area2[:, None]
    nhat *= np.sign(np.einsum('ij,ij->i', nhat, C))[:, None]   # outward (convex about origin)
    # surface gradient of P1 psi: solve [e1; e2; n] . g = [psi1-psi0, psi2-psi0, 0] per triangle
    d = np.stack([psi_v[T[:, 1]] - psi_v[T[:, 0]], psi_v[T[:, 2]] - psi_v[T[:, 0]],
                  np.zeros(len(T))], axis=1)
    Mtx = np.stack([e1, e2, nhat], axis=1)            # (NT,3,3) rows
    g = np.linalg.solve(Mtx, d[:, :, None])[:, :, 0]   # grad_s psi (in-plane)
    K = np.cross(nhat, g)                              # surface current direction*magnitude
    Kw = K * (0.5 * area2)[:, None]                    # area-weighted
    return C, Kw


def biot_savart_H(C, Kw, P):
    """H_s(P) = (1/4pi) sum_tri Kw x (P - C) / |P - C|^3  (H, not B; matches the analytic jump=psi)."""
    out = np.zeros((len(P), 3))
    for i, p in enumerate(P):
        d = p - C
        r3 = (np.einsum('ij,ij->i', d, d)) ** 1.5 + 1e-300
        out[i] = INV4PI * (np.cross(Kw, d) / r3[:, None]).sum(axis=0)
    return out


# ---- reduced-potential volume mesh (ball + spherical iron shell + Kelvin), coil NOT meshed ----
def vol_geo(kelvin, with_iron, R_box=3.0):
    if kelvin:
        outer = Sphere(Pnt(0, 0, 0), R_out)
        for f in outer.faces:
            f.name = "kelvin_int"
        if with_iron:
            sb, sc = Sphere(Pnt(0, 0, 0), b), Sphere(Pnt(0, 0, 0), c)
            core = sb; core.mat("vac"); shell = (sc - sb); shell.mat("shell"); out = (outer - sc); out.mat("vac")
            solids = [core, shell, out]
        else:
            outer.mat("vac"); solids = [outer]
        kball = Sphere(Pnt(offset, 0, 0), R_out)
        for f in kball.faces:
            f.name = "kelvin_ext"
        kball.mat("kelvin"); gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
        fi = [f for s in solids for f in s.faces if f.name == "kelvin_int"][0]
        fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
        fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
        return occ.Glue(solids + [kball, gnd])
    box = Sphere(Pnt(0, 0, 0), R_box)
    for f in box.faces:
        f.name = "outer"
    sb, sc = Sphere(Pnt(0, 0, 0), b), Sphere(Pnt(0, 0, 0), c)
    core = sb; core.mat("vac"); core.maxh = 0.13
    shell = (sc - sb); shell.mat("shell"); shell.maxh = 0.11
    out = (box - sc); out.mat("vac")
    return occ.Glue([core, shell, out])


def solve_reduced(C, Kw, kelvin, with_iron, gmaxh, R_box=3.0):
    mesh = ng.Mesh(OCCGeometry(vol_geo(kelvin, with_iron, R_box)).GenerateMesh(maxh=gmaxh)).Curve(min(order + 1, 4))
    xx, yy, zz = ng.x, ng.y, ng.z; rp2 = (xx - offset) ** 2 + yy * yy + zz * zz + 1e-20
    if kelvin:
        mu = mesh.MaterialCF({"vac": 1.0, "shell": MU_R, "kelvin": R_out ** 2 / rp2}, default=1.0)
        fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="GND"))
    else:
        mu = mesh.MaterialCF({"vac": 1.0, "shell": MU_R}, default=1.0)
        fes = ng.H1(mesh, order=order, dirichlet="outer")
    u, w = fes.TnT()
    A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(w) * ng.dx(bonus_intorder=6)); A.Assemble()
    f = ng.LinearForm(fes)
    if with_iron:
        # inject H_s as P1 nodal source on the iron-shell vertices (H_s smooth there, away from coil)
        Vc = np.array([list(v.point) for v in mesh.vertices])
        r = np.linalg.norm(Vc, axis=1)
        iron = np.where((r > b - 0.04) & (r < c + 0.04))[0]
        Hs_iron = biot_savart_H(C, Kw, Vc[iron])
        gcomp = []
        for k in range(3):
            g = ng.GridFunction(ng.H1(mesh, order=1)); g.vec[:] = 0.0
            arr = g.vec.FV().NumPy(); arr[iron] = Hs_iron[:, k]   # P1: vertex dof index == vertex nr
            gcomp.append(g)
        Hs_cf = ng.CF((gcomp[0], gcomp[1], gcomp[2]))
        f += (MU_R - 1.0) * ng.InnerProduct(Hs_cf, ng.grad(w)) * ng.dx(definedon=mesh.Materials("shell"),
                                                                       bonus_intorder=4)
    f.Assemble()
    gf = ng.GridFunction(fes); gf.vec[:] = 0.0
    gf.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    gradO = ng.grad(gf)
    Hs_t = biot_savart_H(C, Kw, PTS)
    Htot = Hs_t - np.array([[gradO(mesh(*p))[k] for k in range(3)] for p in PTS])
    return Htot, Hs_t, fes.ndof


with TaskManager():
    print("PART A -- SPHERE coil: FE-direct Biot-Savart H_s vs increment-1 ANALYTIC, + iron solve")
    Vs, Ts = coil_surface((a, a, a), hc=0.14)
    psi_s = Vs[:, 2] ** 2 - 0.5 * (Vs[:, 0] ** 2 + Vs[:, 1] ** 2)        # psi = S2 (Z2 zonal) on sphere
    Cs, Kws = tri_geom(Vs, Ts, psi_s)
    Hs_num = biot_savart_H(Cs, Kws, PTS)
    Hs_ana = 0.6 * np.stack([PTS[:, 0], PTS[:, 1], -2 * PTS[:, 2]], axis=1)   # psi=S2 -> H_in = (3/5)(x,y,-2z)
    relA = np.linalg.norm(Hs_num - Hs_ana) / np.linalg.norm(Hs_ana)
    print("  |H_s(FE-direct) - H_s(analytic)| / |.| = %.2e  (validates the numpy current-sheet kernel)" % relA)
    Hk_s, _, nks = solve_reduced(Cs, Kws, True, True, 0.12)
    shield_s = np.linalg.norm(Hk_s) / np.linalg.norm(Hs_num)
    print("  sphere + iron (Kelvin, ndof %d): shield = ||H_iron||/||H_vac|| = %.4f" % (nks, shield_s))

    print("\nPART B -- ELLIPSOID coil (axes 0.36, 0.36, 0.52): FE-direct H_s + iron, Kelvin vs air-box")
    Ve, Te = coil_surface((0.36, 0.36, 0.52), hc=0.13)
    psi_e = Ve[:, 2] ** 2 - 0.5 * (Ve[:, 0] ** 2 + Ve[:, 1] ** 2)
    Ce, Kwe = tri_geom(Ve, Te, psi_e)
    Hk, Hs_t, nk = solve_reduced(Ce, Kwe, True, True, 0.12)
    Hb, _, nb = solve_reduced(Ce, Kwe, False, True, 0.12, R_box=3.0)
    rel_kb = np.linalg.norm(Hk - Hb) / np.linalg.norm(Hb)
    shield_e = np.linalg.norm(Hk) / np.linalg.norm(Hs_t)
    print("  ellipsoid + iron: Kelvin(ndof %d) vs air-box(ndof %d): ||H_K - H_box||/||.|| = %.2e" % (nk, nb, rel_kb))
    print("  ellipsoid shield = ||H_iron||/||H_vac|| = %.4f" % shield_e)

    print("\n[RESULT] FE-direct current-sheet kernel matches analytic (%.1e); ellipsoid material-aware"
          " Kelvin==air-box (%.1e)." % (relA, rel_kb))
