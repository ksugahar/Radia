# SINGLE-PASS coil WITH FEED LEADS -- the case where the leads carry the FULL current and their field
# is NOT negligible (induction-heating-style coils).  In a multi-turn coil the working part has N*I
# while each lead has only I, so the lead field is ~1/N and is dropped.  In a SINGLE-PASS coil the lead
# carries the same I as the working conductor -> the lead field is a FIRST-ORDER part of the design,
# not a perturbation.
#
# A single-pass coil is an OPEN conductor (terminal to terminal); current continuity (div J = 0) means
# the circuit must be closed.  We model the coil as ONE closed polyline:
#     terminal_B --lead_in--> End_B --arc(+phi)--> End_A --lead_out--> terminal_A --chord--> terminal_B
# The two radial leads sit at azimuth +-gap/2 and carry opposite radial currents -> as the gap -> 0 the
# lead pair cancels; a finite gap leaves a first-order lead field (the busbar chord closes it far out).
#
# The verified iron-aware forward map (act8_13_sf_contour_filament_fmm_source/act7_12_kelvin_pml_complementary: numpy Biot-Savart at the interior DSV + FMM
# BiotSavartCF source in the shell + Kelvin reduced-potential iron reaction) is GEOMETRY-AGNOSTIC in
# the wire, so the leads enter the map simply by appending their segments -- and the iron reaction is
# applied to the lead field automatically.
#
# Shown: (A) the lead field at the DSV is first-order and grows with lead separation; (B) a design that
# INCLUDES the leads in the iron-aware model HITS the true target, while a design that OMITS the leads
# (treats the coil as a bare loop) MISSES when evaluated against the true leads+iron field.
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

b, c = 0.58, 0.74
R_out, offset, MU_R = 1.1, 3.3, 50.0
MP_ORDER, KAPPA, MP_RAD = 14, 1e-6, 0.52       # FMM: center origin, radius encloses leads (R_term<0.52<shell b)
INV4PI = 1.0 / (4.0 * np.pi)
R_TERM = 0.50                                   # lead terminal radius (busbar), inside the shell
GAP_DEG = 20.0                                  # design lead separation (azimuthal gap of the break)

# TARGET = a WORKPIECE patch right next to the coil (the induction-heating geometry), NOT a far central
# DSV.  This is where the leads matter: a probe shows lead/arc field is 2% at a far central DSV (r=0.20)
# but 33% / 49% / 74% on a workpiece patch at standoff 0.10 / 0.06 / 0.03 from the conductor.  The leads
# carry the full current and pass close to the workpiece -> their field is first-order AND dominant.
_phw = np.linspace(-np.deg2rad(40), np.deg2rad(40), 9)
_rw = np.linspace(0.30, 0.50, 5)
_zwp = 0.06                                       # workpiece plane (standoff 0.06 below the coil plane z0=0.12)
PTS = np.array([[r * np.cos(a), r * np.sin(a), _zwp] for a in _phw for r in _rw])
NT = PTS.shape[0] * 3

REF = np.array([0.40, 0.12, 1.0])               # reference working loop (rho, z0, I)
THETA0 = np.array([0.30, 0.0, 0.5])
LB = np.array([0.15, -0.30, 0.05]); UB = np.array([0.50, 0.30, 5.0])


def _arc(rho, z0, a0, a1, n):
    ph = np.linspace(a0, a1, n + 1)
    return np.stack([rho * np.cos(ph), rho * np.sin(ph), z0 * np.ones_like(ph)], axis=1)


def single_pass_filaments(theta, gap, R_term=R_TERM, M_arc=180, M_lead=16, M_chord=6):
    """closed polyline: lead_in + working arc(long way) + lead_out + terminal chord, all carrying I."""
    rho, z0, I = theta; g = gap
    tB = np.array([R_term * np.cos(g / 2), R_term * np.sin(g / 2), z0])     # terminal B (+gap/2)
    eB = np.array([rho * np.cos(g / 2), rho * np.sin(g / 2), z0])          # arc end B
    eA = np.array([rho * np.cos(g / 2), -rho * np.sin(g / 2), z0])         # arc end A (-gap/2)
    tA = np.array([R_term * np.cos(g / 2), -R_term * np.sin(g / 2), z0])    # terminal A
    pieces = [
        np.linspace(tB, eB, M_lead + 1),                                   # lead_in (radial)
        _arc(rho, z0, g / 2, 2 * np.pi - g / 2, M_arc),                     # working arc (long way, +phi)
        np.linspace(eA, tA, M_lead + 1),                                   # lead_out (radial)
        np.linspace(tA, tB, M_chord + 1),                                  # terminal chord (busbar return)
    ]
    pts = np.vstack([p if i == 0 else p[1:] for i, p in enumerate(pieces)])
    return pts[:-1], pts[1:], I * np.ones(len(pts) - 1)


def leads_only_filaments(theta, gap, R_term=R_TERM, M_lead=16, M_chord=6):
    """just the lead pair + chord (no working arc) -- to measure the lead field in isolation."""
    rho, z0, I = theta; g = gap
    tB = np.array([R_term * np.cos(g / 2), R_term * np.sin(g / 2), z0])
    eB = np.array([rho * np.cos(g / 2), rho * np.sin(g / 2), z0])
    eA = np.array([rho * np.cos(g / 2), -rho * np.sin(g / 2), z0])
    tA = np.array([R_term * np.cos(g / 2), -R_term * np.sin(g / 2), z0])
    pieces = [np.linspace(tB, eB, M_lead + 1), np.linspace(eA, tA, M_lead + 1),
              np.linspace(tA, tB, M_chord + 1)]
    pts = np.vstack([p if i == 0 else p[1:] for i, p in enumerate(pieces)])
    # NOTE: open chain (eB..eA missing) -> only meaningful as a FIELD probe at the DSV, not a circuit
    return pts[:-1], pts[1:], I * np.ones(len(pts) - 1)


def bare_loop_filaments(theta, M=200):
    """full closed loop (the 'ignore the leads' model)."""
    rho, z0, I = theta
    ph = np.linspace(0, 2 * np.pi, M + 1)
    pts = np.stack([rho * np.cos(ph), rho * np.sin(ph), z0 * np.ones_like(ph)], axis=1)
    return pts[:-1], pts[1:], I * np.ones(M)


def biot_savart_np(sp, ep, Iseg, P):
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


def kelvin_geo():
    sb, sc = Sphere(Pnt(0, 0, 0), b), Sphere(Pnt(0, 0, 0), c)
    outer = Sphere(Pnt(0, 0, 0), R_out)
    for f in outer.faces:
        f.name = "kelvin_int"
    core = sb; core.mat("vac"); shell = (sc - sb); shell.mat("shell"); shell.maxh = 0.12
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
    def __init__(self, maxh=0.20, order=2):
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
            dx(definedon=self.mesh.Materials("shell"), bonus_intorder=2 * self.order + 2)
        f.Assemble()
        gf = ng.GridFunction(self.fes); gf.vec[:] = 0.0
        gf.vec.data = self.Ainv * f.vec
        gO = grad(gf)
        return -np.array([[gO(t)[k] for k in range(3)] for t in self.targets])


def H_iron(sp, ep, Iseg, fwd):
    Hs = biot_savart_np(sp, ep, Iseg, PTS)
    react = fwd.reaction(biot_savart_cf(sp, ep, Iseg))
    return (Hs + react).reshape(-1)


with TaskManager():
    print("SINGLE-PASS coil WITH FEED LEADS -- leads carry the full current (induction-heating-style)")
    print("  iron shell [%.2f,%.2f] mu_r=%g, Kelvin; lead terminal R=%.2f, gap=%g deg; workpiece patch=%d pts"
          % (b, c, MU_R, R_TERM, GAP_DEG, len(PTS)))
    t0 = time.time(); fwd = IronForward()
    print("  Kelvin-FEM factored once: ndof=%d  (%.1fs)\n" % (fwd.ndof, time.time() - t0))

    # --- (A) the lead field is FIRST-ORDER: sweep the lead separation (azimuthal gap) ---
    print("(A) lead field vs working-arc field AT THE WORKPIECE PATCH (free space), sweeping the lead gap:")
    print("    gap[deg] | ||H_leads|| / ||H_arc||  (first-order AND large -- target is next to the coil)")
    for gd in (5.0, 10.0, 20.0, 40.0):
        g = np.deg2rad(gd)
        sp_a, ep_a, I_a = single_pass_filaments(REF, g)
        H_full = biot_savart_np(sp_a, ep_a, I_a, PTS)
        sp_l, ep_l, I_l = leads_only_filaments(REF, g)
        H_lead = biot_savart_np(sp_l, ep_l, I_l, PTS)
        H_arc = H_full - H_lead
        print("     %5.0f   |   %.3f" % (gd, np.linalg.norm(H_lead) / np.linalg.norm(H_arc)))

    # --- target = the TRUE iron-system field of the single-pass coil (working arc + leads) ---
    gap = np.deg2rad(GAP_DEG)
    sp_r, ep_r, I_r = single_pass_filaments(REF, gap)
    B_target = H_iron(sp_r, ep_r, I_r, fwd); nrm = np.linalg.norm(B_target)

    # --- design A: model INCLUDES the leads (full single-pass coil) in the iron-aware map ---
    def resid_withleads(theta):
        sp, ep, Iseg = single_pass_filaments(theta, gap)
        return H_iron(sp, ep, Iseg, fwd) - B_target

    t1 = time.time()
    sa = least_squares(resid_withleads, THETA0, bounds=(LB, UB), diff_step=3e-3, xtol=1e-7, ftol=1e-9)
    hit = np.linalg.norm(resid_withleads(sa.x)) / nrm
    print("\n(B) DESIGN INCLUDING leads in the iron-aware model:")
    print("    rho=%.4f z0=%.4f I=%.4f  residual=%.2e  <- HITS  (%.1fs)" % (*sa.x, hit, time.time() - t1))

    # --- design B: model OMITS the leads (bare loop) -> evaluate against the TRUE (leads+iron) target ---
    def resid_noleads(theta):
        sp, ep, Iseg = bare_loop_filaments(theta)
        return H_iron(sp, ep, Iseg, fwd) - B_target

    t2 = time.time()
    sb_ = least_squares(resid_noleads, THETA0, bounds=(LB, UB), diff_step=3e-3, xtol=1e-7, ftol=1e-9)
    sp_b, ep_b, I_b = single_pass_filaments(sb_.x, gap)     # build the REAL coil at the no-lead optimum
    miss = np.linalg.norm(H_iron(sp_b, ep_b, I_b, fwd) - B_target) / nrm
    print("    DESIGN OMITTING leads (bare-loop model):")
    print("    rho=%.4f z0=%.4f I=%.4f  residual(true coil)=%.2e  <- MISSES  (%.1fs)"
          % (*sb_.x, miss, time.time() - t2))

    print("\n  reference single-pass coil : rho=%.4f z0=%.4f I=%.4f (gap=%g deg)" % (*REF, GAP_DEG))
    print("  with-leads design recovers : d(rho,z0,I) = (%.1e, %.1e, %.1e)" % tuple(np.abs(sa.x - REF)))
    print("  no-leads design drifts     : d(rho,z0,I) = (%.1e, %.1e, %.1e)" % tuple(np.abs(sb_.x - REF)))
    print("\n[RESULT] single-pass coil with feed leads: the lead field is first-order; INCLUDING the leads"
          " in the iron-aware map HITS (%.1e), OMITTING them MISSES the true coil field (%.1f%%)."
          % (hit, 100 * miss))

    out = {
        "geometry": {"shell": [b, c], "mu_r": MU_R, "R_term": R_TERM, "gap_deg": GAP_DEG,
                     "DSV_pts": len(PTS), "ndof": fwd.ndof},
        "reference": {"rho": REF[0], "z0": REF[1], "I": REF[2]},
        "design_with_leads": {"rho": sa.x[0], "z0": sa.x[1], "I": sa.x[2], "residual": hit},
        "design_omit_leads": {"rho": sb_.x[0], "z0": sb_.x[1], "I": sb_.x[2], "residual_true": miss},
    }
    with open(os.path.join(os.path.dirname(__file__), "demo_sp2_results.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("  results -> demo_sp2_results.json")
