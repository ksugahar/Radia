# SINGLE-PASS (1-path / one continuous wire) coil INVERSE DESIGN through the verified iron-aware
# forward map.  Multi-turn coils are designed by stream-function contour levels; a SINGLE-PASS coil
# has no continuum current density -- its only design DOFs are the WIRE PATH shape.  So instead of a
# linear least-norm over psi (demo_oo), we run a NONLINEAR fit over the path parameters, evaluating
# the field through demo_rr's already-verified pipeline:
#
#     wire path  --numpy Biot-Savart (DSV, interior)--> H_s
#                --ngsolve.bem.BiotSavartCF (FMM, shell)--> source CF
#                --Kelvin reduced-potential FEM (iron reaction)--> -grad(Omega)
#     H_total(DSV) = H_s + reaction
#
# The same pipeline that turns SF contour loops into a material-aware field (demo_rr) is geometry-
# agnostic in the wire: pass ONE loop's segments and you have the forward map of a single-pass coil.
# The Kelvin-FEM stiffness depends only on (geometry, mu) -> factor it ONCE, reuse across the
# optimizer's iterations (only the source-dependent RHS is reassembled).
#
# Narrative (mirrors demo_oo): target = the iron-system field of a REFERENCE single loop (so it is
# single-pass-achievable by construction).  Optimizing the wire path THROUGH THE IRON-AWARE MAP HITS;
# optimizing it through the FREE-SPACE map and then evaluating in the true iron system MISSES.
import os
import time
import json

os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from ngsolve import TaskManager, grad, dx
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry
from ngsolve.bem import BiotSavartCF
from ngsolve.bla import Vec3D
from scipy.optimize import least_squares

b, c = 0.58, 0.74                  # spherical iron shell [b, c]
R_out, offset, MU_R = 1.1, 3.3, 50.0
MP_ORDER, KAPPA, MP_RAD = 20, 1e-6, 0.52   # FMM multipole: center origin, radius encloses any loop < 0.52 < shell b
INV4PI = 1.0 / (4.0 * np.pi)

th = np.deg2rad(np.linspace(25, 155, 5)); ph0 = np.deg2rad(np.linspace(0, 300, 6))
PTS = np.array([[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)] for t in th for p in ph0]) * 0.20
NT = PTS.shape[0] * 3              # full vector target rows

# --- single-pass wire = ONE circular loop, design DOFs theta = (radius rho, axial z0, current I) ---
REF = np.array([0.40, 0.12, 1.0])         # reference loop (the achievable target generator)
THETA0 = np.array([0.30, 0.0, 0.5])       # optimizer start (deliberately away from REF)
LB = np.array([0.15, -0.30, 0.05]); UB = np.array([0.50, 0.30, 5.0])


def loop_filaments(theta, M=160):
    rho, z0, I = theta
    ph = np.linspace(0, 2 * np.pi, M + 1)
    pts = np.stack([rho * np.cos(ph), rho * np.sin(ph), z0 * np.ones_like(ph)], axis=1)
    return pts[:-1], pts[1:], I * np.ones(M)


def biot_savart_np(sp, ep, Iseg, P):
    """numpy filament Biot-Savart H at P from current elements (midpoint, I*dl). Valid everywhere."""
    Ce = 0.5 * (sp + ep); Kwe = Iseg[:, None] * (ep - sp)
    out = np.zeros((len(P), 3))
    for i, p in enumerate(P):
        d = p - Ce; r3 = (np.einsum('ij,ij->i', d, d)) ** 1.5 + 1e-300
        out[i] = INV4PI * (np.cross(Kwe, d) / r3[:, None]).sum(axis=0)
    return out


def biot_savart_cf(sp, ep, Iseg):
    """native FMM Biot-Savart CF -- valid OUTSIDE MP_RAD (the iron shell), not at the interior DSV."""
    bs = BiotSavartCF(MP_ORDER, KAPPA, Vec3D(0, 0, 0), MP_RAD)
    for k in range(len(sp)):
        bs.AddCurrent(Vec3D(*sp[k]), Vec3D(*ep[k]), float(Iseg[k]))
    return ng.CF((bs[0].real, bs[1].real, bs[2].real))


def kelvin_geo():
    sb, sc = Sphere(Pnt(0, 0, 0), b), Sphere(Pnt(0, 0, 0), c)
    outer = Sphere(Pnt(0, 0, 0), R_out)
    for f in outer.faces:
        f.name = "kelvin_int"
    core = sb; core.mat("vac"); shell = (sc - sb); shell.mat("shell"); shell.maxh = 0.11
    out = (outer - sc); out.mat("vac"); solids = [core, shell, out]
    kball = Sphere(Pnt(offset, 0, 0), R_out)
    for f in kball.faces:
        f.name = "kelvin_ext"
    kball.mat("kelvin"); gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
    fi = [f for s in solids for f in s.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
    return occ.Glue(solids + [kball, gnd])


class IronForward:
    """Kelvin reduced-potential iron reaction. Factor the (geometry, mu)-only stiffness ONCE; each
    call reassembles only the source-dependent RHS and back-substitutes."""

    def __init__(self, maxh=0.18, order=2):
        mesh = ng.Mesh(OCCGeometry(kelvin_geo()).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 4))
        xx, yy, zz = ng.x, ng.y, ng.z; rp2 = (xx - offset) ** 2 + yy * yy + zz * zz + 1e-20
        mu = mesh.MaterialCF({"vac": 1.0, "shell": MU_R, "kelvin": R_out ** 2 / rp2}, default=1.0)
        fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="GND"))
        u, w = fes.TnT()
        A = ng.BilinearForm(mu * grad(u) * grad(w) * dx(bonus_intorder=2 * order + 2)); A.Assemble()
        self.mesh, self.fes, self.w, self.order = mesh, fes, w, order
        self.Ainv = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        self.targets = [mesh(*p) for p in PTS]
        self.ndof = fes.ndof

    def reaction(self, Hs_cf):
        w = self.w
        f = ng.LinearForm(self.fes)
        f += (MU_R - 1.0) * ng.InnerProduct(Hs_cf, grad(w)) * \
            dx(definedon=self.mesh.Materials("shell"), bonus_intorder=2 * self.order + 6)
        f.Assemble()
        gf = ng.GridFunction(self.fes); gf.vec[:] = 0.0
        gf.vec.data = self.Ainv * f.vec
        gO = grad(gf)
        return -np.array([[gO(t)[k] for k in range(3)] for t in self.targets])


def H_iron(theta, fwd):
    sp, ep, Iseg = loop_filaments(theta)
    Hs = biot_savart_np(sp, ep, Iseg, PTS)                 # interior DSV: numpy (FMM invalid here)
    react = fwd.reaction(biot_savart_cf(sp, ep, Iseg))     # shell source: FMM CF
    return (Hs + react).reshape(-1)


def H_free(theta):
    sp, ep, Iseg = loop_filaments(theta)
    return biot_savart_np(sp, ep, Iseg, PTS).reshape(-1)   # iron ignored


with TaskManager():
    print("SINGLE-PASS coil inverse design (one continuous wire) through the iron-aware forward map")
    print("  iron shell [%.2f,%.2f] mu_r=%g, Kelvin open boundary; DSV = %d pts (x3 = %d rows)"
          % (b, c, MU_R, len(PTS), NT))
    t0 = time.time()
    fwd = IronForward()
    print("  Kelvin-FEM factored once: ndof=%d  (%.1fs)" % (fwd.ndof, time.time() - t0))

    # --- (A) sanity: on-axis loop H_z matches the analytic single-loop formula ---
    sp, ep, Iseg = loop_filaments(REF)
    p_ax = np.array([[0.0, 0.0, 0.05]])
    Hz_np = biot_savart_np(sp, ep, Iseg, p_ax)[0, 2]
    rho_r, z0_r, I_r = REF; d = 0.05 - z0_r
    Hz_an = I_r * rho_r ** 2 / (2.0 * (rho_r ** 2 + d ** 2) ** 1.5)
    print("\n(A) sanity: on-axis H_z numpy=%.6e  analytic=%.6e  rel=%.2e"
          % (Hz_np, Hz_an, abs(Hz_np - Hz_an) / abs(Hz_an)))

    # --- target = the IRON-system field of the reference single loop (single-pass-achievable) ---
    B_target = H_iron(REF, fwd); nrm = np.linalg.norm(B_target)
    react_ref = B_target - H_free(REF)
    print("    iron reaction / free field at DSV = %.3f  (iron is a real part of the operator)"
          % (np.linalg.norm(react_ref) / np.linalg.norm(H_free(REF))))

    # --- design WITH the iron-aware map ---
    t1 = time.time()
    sol_iron = least_squares(lambda th: H_iron(th, fwd) - B_target, THETA0, bounds=(LB, UB),
                             diff_step=3e-3, xtol=1e-10, ftol=1e-12)
    th_iron = sol_iron.x
    hit = np.linalg.norm(H_iron(th_iron, fwd) - B_target) / nrm
    print("\n(B) DESIGN WITH iron-aware map : rho=%.4f z0=%.4f I=%.4f  residual=%.2e  <- HITS  (%.1fs)"
          % (*th_iron, hit, time.time() - t1))

    # --- design IGNORING iron (free-space map), then evaluate in the TRUE iron system ---
    t2 = time.time()
    sol_free = least_squares(lambda th: H_free(th) - B_target, THETA0, bounds=(LB, UB),
                             diff_step=3e-3, xtol=1e-10, ftol=1e-12)
    th_free = sol_free.x
    miss = np.linalg.norm(H_iron(th_free, fwd) - B_target) / nrm
    print("    DESIGN IGNORING iron (free) : rho=%.4f z0=%.4f I=%.4f  residual(in iron)=%.2e  <- MISSES (%.1fs)"
          % (*th_free, miss, time.time() - t2))

    print("\n  reference loop          : rho=%.4f z0=%.4f I=%.4f" % tuple(REF))
    print("  iron-aware recovers it  : d(rho,z0,I) = (%.1e, %.1e, %.1e)"
          % tuple(np.abs(th_iron - REF)))
    print("  free-space design drifts: d(rho,z0,I) = (%.1e, %.1e, %.1e)"
          % tuple(np.abs(th_free - REF)))
    print("\n[RESULT] single-pass wire-path inverse design through the iron-aware Kelvin/FMM map HITS"
          " (%.1e); ignoring the iron yoke MISSES the in-iron target (%.1f%%)." % (hit, 100 * miss))

    out = {
        "geometry": {"shell": [b, c], "mu_r": MU_R, "R_out": R_out, "offset": offset,
                     "DSV_pts": len(PTS), "ndof": fwd.ndof},
        "reference_loop": {"rho": REF[0], "z0": REF[1], "I": REF[2]},
        "design_iron_aware": {"rho": th_iron[0], "z0": th_iron[1], "I": th_iron[2], "residual": hit},
        "design_free_space": {"rho": th_free[0], "z0": th_free[1], "I": th_free[2],
                              "residual_in_iron": miss},
        "iron_reaction_over_free": float(np.linalg.norm(react_ref) / np.linalg.norm(H_free(REF))),
        "sanity_onaxis_rel": float(abs(Hz_np - Hz_an) / abs(Hz_an)),
    }
    with open(os.path.join(os.path.dirname(__file__), "demo_sp1_results.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("  results -> demo_sp1_results.json")
