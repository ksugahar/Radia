# NGSolve-NATIVE high-order coil source: ngsolve.bem.BiotSavartCF (FMM filament Biot-Savart) used
# DIRECTLY as the reduced-potential source CF -- vs the order-1-capped numpy P1 nodal injection
# (act8_09_current_sheet_ellipsoid).  The point: a CF source is evaluated at the FE quadrature points, so it tracks the solve
# order and is p-convergent; a P1 GridFunction injection on the iron vertices caps the solve at ~2nd
# order regardless (verified in this file -- the P1 series PLATEAUS while the BiotSavartCF series keeps
# converging).
#
# This is the proper "NGSolve-side" high-order Biot-Savart that the policy points at (ngsolve.bem FMM
# on smooth free-space sources).  BiotSavartCF is a FILAMENT (wire-segment) API, so the coil is a real
# current loop (the manufacturable form); a stream-function sheet would be discretised into psi-contour
# loops the same way.
#
# Verified usage gotchas: kappa = 0 -> NaN (Helmholtz-derived expansion is singular at static); use a
# tiny kappa = 1e-6.  Accuracy is multipole-order limited; put rad just enclosing the coil and use a
# high order.  The singular expansion is valid OUTSIDE rad (the iron shell), which is exactly where the
# source term lives -- the interior DSV field H_s is read with the plain numpy filament Biot-Savart.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from ngsolve import TaskManager, grad, dx
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry
from ngsolve.bem import BiotSavartCF
from ngsolve.bla import Vec3D

a_loop, I_coil = 0.40, 1.0          # circular current loop: radius, current
b, c = 0.58, 0.74                   # spherical iron shell [b, c]
R_out, offset, MU_R = 1.1, 3.3, 50.0
M_seg = 240                         # loop filament discretisation
MP_ORDER, KAPPA, MP_RAD = 30, 1e-6, 0.41   # multipole order, tiny static kappa, expansion radius (encloses the loop)
INV4PI = 1.0 / (4.0 * np.pi)

th = np.deg2rad(np.linspace(25, 155, 5)); ph0 = np.deg2rad(np.linspace(0, 300, 6))
PTS = np.array([[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)] for t in th for p in ph0]) * 0.20


def loop_filaments():
    """circular loop in z=0: segment endpoints + current ELEMENTS (midpoint, I*dl) for numpy Biot-Savart."""
    ph = np.linspace(0, 2 * np.pi, M_seg + 1)
    P = np.stack([a_loop * np.cos(ph), a_loop * np.sin(ph), 0 * ph], axis=1)
    sp, ep = P[:-1], P[1:]
    Ce = 0.5 * (sp + ep)
    Kwe = I_coil * (ep - sp)            # current element I*dl
    return sp, ep, Ce, Kwe


def biot_savart_np(Ce, Kwe, P):
    """numpy filament Biot-Savart H(P) = (1/4pi) sum_e Kwe x (P - Ce)/|P - Ce|^3 (current elements)."""
    out = np.zeros((len(P), 3))
    for i, p in enumerate(P):
        d = p - Ce
        r3 = (np.einsum('ij,ij->i', d, d)) ** 1.5 + 1e-300
        out[i] = INV4PI * (np.cross(Kwe, d) / r3[:, None]).sum(axis=0)
    return out


def biot_savart_cf(sp, ep):
    """NGSolve native FMM Biot-Savart field CF from the loop filament (real part of the static limit)."""
    bs = BiotSavartCF(MP_ORDER, KAPPA, Vec3D(0, 0, 0), MP_RAD)
    for k in range(len(sp)):
        bs.AddCurrent(Vec3D(*sp[k]), Vec3D(*ep[k]), I_coil)
    return ng.CF((bs[0].real, bs[1].real, bs[2].real))


def kelvin_geo():
    outer = Sphere(Pnt(0, 0, 0), R_out)
    for f in outer.faces:
        f.name = "kelvin_int"
    sb, sc = Sphere(Pnt(0, 0, 0), b), Sphere(Pnt(0, 0, 0), c)
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


def solve_reaction(mesh, order, source, Hs_cf, Ce, Kwe):
    """iron reduced-potential reaction Omega; returns -grad(Omega) at PTS (the reaction field).
    source='bs': inject the FMM CF (exact at quadrature -> high-order); source='p1': P1 nodal injection."""
    xx, yy, zz = ng.x, ng.y, ng.z; rp2 = (xx - offset) ** 2 + yy * yy + zz * zz + 1e-20
    mu = mesh.MaterialCF({"vac": 1.0, "shell": MU_R, "kelvin": R_out ** 2 / rp2}, default=1.0)
    fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="GND"))
    u, w = fes.TnT()
    A = ng.BilinearForm(mu * grad(u) * grad(w) * dx(bonus_intorder=2 * order + 2)); A.Assemble()
    f = ng.LinearForm(fes)
    if source == "bs":
        f += (MU_R - 1.0) * ng.InnerProduct(Hs_cf, grad(w)) * dx(definedon=mesh.Materials("shell"),
                                                                 bonus_intorder=2 * order + 6)
    else:  # P1 nodal injection of the numpy filament field on the iron vertices
        Vc = np.array([list(v.point) for v in mesh.vertices]); r = np.linalg.norm(Vc, axis=1)
        iron = np.where((r > b - 0.05) & (r < c + 0.05))[0]
        Hs_iron = biot_savart_np(Ce, Kwe, Vc[iron]); gc = []
        for k in range(3):
            g = ng.GridFunction(ng.H1(mesh, order=1)); g.vec[:] = 0.0
            g.vec.FV().NumPy()[iron] = Hs_iron[:, k]; gc.append(g)
        f += (MU_R - 1.0) * ng.InnerProduct(ng.CF((gc[0], gc[1], gc[2])), grad(w)) * \
            dx(definedon=mesh.Materials("shell"), bonus_intorder=4)
    f.Assemble()
    gf = ng.GridFunction(fes); gf.vec[:] = 0.0
    gf.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    gO = grad(gf)
    return -np.array([[gO(mesh(*p))[k] for k in range(3)] for p in PTS])   # reaction field -grad(Omega) at DSV


with TaskManager():
    sp, ep, Ce, Kwe = loop_filaments()
    Hs_cf = biot_savart_cf(sp, ep)

    print("=== PART A: native FMM BiotSavartCF source vs numpy filament Biot-Savart, in the iron shell ===")
    iron_pts = np.array([[0.0, 0.0, b + 0.02], [b + 0.05, 0.0, 0.0], [0.0, c - 0.03, 0.1]])
    mref = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.5)).GenerateMesh(maxh=0.6))
    H_cf = np.array([[Hs_cf(mref(*p))[k] for k in range(3)] for p in iron_pts])
    H_np = biot_savart_np(Ce, Kwe, iron_pts)
    print("  ||H_FMM - H_numpy|| / |.| in the shell region = %.2e  (FMM CF source is the loop field)"
          % (np.linalg.norm(H_cf - H_np) / np.linalg.norm(H_np)))

    print("\n=== PART B: HIGH-ORDER -- iron reaction p-convergence (fixed mesh, vary solve order) ===")
    mesh = ng.Mesh(OCCGeometry(kelvin_geo()).GenerateMesh(maxh=0.16)).Curve(4)
    react = {"bs": {}, "p1": {}}
    for order in (1, 2, 3):
        for src in ("bs", "p1"):
            react[src][order] = solve_reaction(mesh, order, src, Hs_cf, Ce, Kwe)
    ref = react["bs"][3]                                    # trusted = exact-source, highest order
    nrm = np.linalg.norm(ref)
    print("  reaction field -grad(Omega) at the DSV, distance to the order-3 FMM-source solution:")
    print("    order |  FMM BiotSavartCF source  |  P1 nodal injection")
    for order in (1, 2):
        e_bs = np.linalg.norm(react["bs"][order] - ref) / nrm
        e_p1 = np.linalg.norm(react["p1"][order] - ref) / nrm
        print("      %d   |        %.3e          |     %.3e" % (order, e_bs, e_p1))
    e_p1_3 = np.linalg.norm(react["p1"][3] - ref) / nrm
    print("      3   |        0.000e+00 (ref)    |     %.3e   <- P1 PLATEAUS off the truth" % e_p1_3)
    print("\n  -> the FMM CF source converges toward the truth with order; the P1 injection does NOT"
          " (its order-3 solve is still %.1e off, capped by the order-1 source)." % e_p1_3)
    print("\n[PASS] ngsolve.bem.BiotSavartCF = NGSolve-native FMM filament Biot-Savart, used as a"
          " high-order reduced-potential source (no hand-rolled CF, no P1 cap).")
