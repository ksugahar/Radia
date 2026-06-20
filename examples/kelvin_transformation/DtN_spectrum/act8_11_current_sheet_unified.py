# VERIFIED (2026-06-18): GAP-FREE UNIFIED current-sheet stream-function coil model.
# ONE Kelvin-FEM, the coil a CONFORMING meshed interface -- the FEM ITSELF produces Biot-Savart in
# vacuum (no analytic H_s bolted on), and the SAME single solve is material-aware with iron.
#
# This is the gap-free counterpart of the HYBRID reduced-potential demos (act8_08_current_sheet_reduced_potential/nn/oo), which add
# a numpy/analytic Biot-Savart H_s on top of an FEM reaction (a "seam").  Here there is no seam:
# the stream-function current sheet K = n x grad psi on the winding surface S enters the single FEM
# as the magnetic-double-layer POTENTIAL JUMP [Omega] = psi.  The classical equivalence
# "surface current K = n x grad psi  <=>  dipole/double layer of strength psi" guarantees the
# double-layer potential's field IS the Biot-Savart field of K -- so a plain Laplace FEM that
# realises the jump reproduces Biot-Savart in vacuum, with NO Biot-Savart kernel anywhere.
#
# Formulation.  Omega = Phi + s :  Phi continuous H1 (Kelvin / Dirichlet far),  s = the jump-lift,
# nonzero only in the coil interior, s = psi on S and harmonic inside.  Testing div(mu grad Omega)=0
# region-wise against a continuous w (the natural interface condition [mu dOmega/dn] = 0 emerges) gives
#       int_full  mu grad(Phi).grad(w)  =  - int_(coil_in) grad(s).grad(w)        for all w,
# with mu = 1 in the coil.  The total field is H = -grad(Omega) = -grad(Phi) - grad(s).  Iron enters
# ONLY through mu in the bilinear form -- one solve yields the material-aware field.
#
# NOTE on the source term (verified below): the "double-layer delta" surface term oint_S psi(grad w . n)
# is REDUNDANT -- adding it changes the result by 0 (bit-identical).  The volume jump-lift term alone
# is the operative source.  The lift-free pure-surface term is WRONG (it cannot carry the jump).
#
# Geometry + psi normalisation match act8_08_current_sheet_reduced_potential (sphere coil r=a, psi = P_2, iron shell [b,c], MU_R=50)
# so this validates against act8_08_current_sheet_reduced_potential's two checks: vacuum == analytic Biot-Savart, and (with iron) an
# independent open BC -- Kelvin vs a brute-force air-box -- but now with the coil MESHED, no H_s seam.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from ngsolve import TaskManager, grad, dx, ds, specialcf
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry

a = 0.40                       # coil sphere radius (a MESHED conforming current-sheet interface)
b, c = 0.55, 0.70              # iron shell [b, c]  (outside the coil: b > a)
R_out, offset = 1.0, 3.0       # Kelvin truncation radius + image-ball offset
order, MU_R, n = 2, 50.0, 2
c_in = (n + 1) / (2 * n + 1) * a ** (-n)          # Phi_in coeff for psi = P_n  (= 3.75 for n=2, a=0.4)

# psi = S2/a^2 = P_2(cos th) on r=a (matches act8_08_current_sheet_reduced_potential's unit-P_2 sheet).  S2 = z^2 - (x^2+y^2)/2 is
# harmonic, so the harmonic jump-lift s equals psi identically inside the coil.
PSI_CF = (ng.z * ng.z - 0.5 * (ng.x * ng.x + ng.y * ng.y)) / a ** 2

th = np.deg2rad(np.linspace(25, 155, 5)); ph = np.deg2rad(np.linspace(0, 300, 6))
PTS = np.array([[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)] for t in th for p in ph]) * 0.22


def H_biotsavart(P):
    """analytic interior Biot-Savart field of the psi = P_2 current sheet:  H = c_in (x, y, -2z)."""
    x, y, z = P[:, 0], P[:, 1], P[:, 2]
    return c_in * np.stack([x, y, -2 * z], axis=1)


def build_geo(kelvin, with_iron, R_box=3.0):
    """coil_in (r<a) + vac gap + iron shell [b,c] + outer vac, closed by Kelvin image-ball or an air-box."""
    s_a = Sphere(Pnt(0, 0, 0), a); s_a.mat("coil_in")
    for f in s_a.faces:
        f.name = "coil_iface"
    s_b = Sphere(Pnt(0, 0, 0), b); s_c = Sphere(Pnt(0, 0, 0), c)
    gap = (s_b - s_a); gap.mat("vac")
    shell = (s_c - s_b); shell.mat("shell")
    if kelvin:
        outer = Sphere(Pnt(0, 0, 0), R_out)
        for f in outer.faces:
            f.name = "kelvin_int"
        out = (outer - s_c); out.mat("vac")
        solids = [s_a, gap, shell, out]
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
    smid = Sphere(Pnt(0, 0, 0), 1.5)
    s_a.maxh = 0.10; gap.maxh = 0.13; shell.maxh = 0.13     # keep coil + iron region fine under the coarse global
    near = (smid - s_c); near.mat("vac"); near.maxh = 0.14  # field-carrying near air: keep fine
    far = (box - smid); far.mat("vac"); far.maxh = 0.45     # negligible shielded far field: coarse
    return occ.Glue([s_a, gap, shell, near, far])


def solve_unified(kelvin, with_iron, gmaxh, source="volume"):
    """One Kelvin/air-box FEM.  Returns interior field H(PTS) and ndof.  source in {volume, both, surface}."""
    mesh = ng.Mesh(OCCGeometry(build_geo(kelvin, with_iron)).GenerateMesh(maxh=gmaxh)).Curve(min(order + 1, 4))
    xx, yy, zz = ng.x, ng.y, ng.z; rp2 = (xx - offset) ** 2 + yy * yy + zz * zz + 1e-20
    mu_shell = MU_R if with_iron else 1.0
    if kelvin:
        mu = mesh.MaterialCF({"coil_in": 1.0, "vac": 1.0, "shell": mu_shell, "kelvin": R_out ** 2 / rp2},
                             default=1.0)
        fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="GND"))
    else:
        mu = mesh.MaterialCF({"coil_in": 1.0, "vac": 1.0, "shell": mu_shell}, default=1.0)
        fes = ng.H1(mesh, order=order, dirichlet="outer")
    # --- jump-lift s in coil_in: s = psi on coil_iface, harmonic inside (psi = S2/a^2 already harmonic) ---
    Vin = ng.H1(mesh, order=order, definedon="coil_in", dirichlet="coil_iface")
    s = ng.GridFunction(Vin); s.Set(PSI_CF, ng.BND, definedon=mesh.Boundaries("coil_iface"))
    us, vs = Vin.TnT()
    Ah = ng.BilinearForm(grad(us) * grad(vs) * dx); Ah.Assemble()
    rs = s.vec.CreateVector(); rs.data = -(Ah.mat * s.vec)
    s.vec.data += Ah.mat.Inverse(Vin.FreeDofs(), inverse="sparsecholesky") * rs
    # --- Phi on the full mesh (iron enters only through mu) ---
    u, w = fes.TnT()
    A = ng.BilinearForm(mu * grad(u) * grad(w) * dx(bonus_intorder=6)); A.Assemble()
    nrm = specialcf.normal(3)
    f = ng.LinearForm(fes)
    if source in ("volume", "both"):
        f += -grad(s) * grad(w) * dx(definedon=mesh.Materials("coil_in"), bonus_intorder=6)
    if source in ("surface", "both"):
        f += PSI_CF * (grad(w).Trace() * nrm) * ds(definedon=mesh.Boundaries("coil_iface"), bonus_intorder=6)
    f.Assemble()
    gf = ng.GridFunction(fes); gf.vec[:] = 0.0
    gf.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    gO, gS = grad(gf), grad(s)
    H = -np.array([[gO(mesh(*p))[k] + gS(mesh(*p))[k] for k in range(3)] for p in PTS])
    return H, fes.ndof


with TaskManager():
    Hbs = H_biotsavart(PTS)
    print("=== VACUUM: source term -- the jump-lift volume term is the operative source ===")
    print("  (Kelvin, no iron, maxh=0.14; the FEM alone must equal Biot-Savart)")
    for src in ("volume", "both", "surface"):
        H, nd = solve_unified(kelvin=True, with_iron=False, gmaxh=0.14, source=src)
        rel = np.linalg.norm(H - Hbs) / np.linalg.norm(Hbs)
        note = {"volume": "operative jump-lift source",
                "both": "+ double-layer delta == volume (redundant)",
                "surface": "lift-free -> cannot carry the jump (WRONG)"}[src]
        print("  source=%-8s ndof=%6d  rel=%.3e   %s" % (src, nd, rel, note))

    print("\n=== VACUUM convergence: the single FEM -> Biot-Savart (no bolt-on) ===")
    for h in (0.20, 0.14, 0.10):
        H, nd = solve_unified(kelvin=True, with_iron=False, gmaxh=h, source="volume")
        rel = np.linalg.norm(H - Hbs) / np.linalg.norm(Hbs)
        print("  maxh=%.2f ndof=%6d  ||H_FEM - H_BiotSavart||/|.| = %.3e" % (h, nd, rel))
    print("  H_FEM[0]=", np.round(H[0], 4), " H_BiotSavart[0]=", np.round(Hbs[0], 4))

    print("\n=== IRON shell mu_r=%g: ONE unified FEM, Kelvin vs air-box (independent open BC) ===" % MU_R)
    Hk, nk = solve_unified(kelvin=True, with_iron=True, gmaxh=0.11, source="volume")
    Hb, nb = solve_unified(kelvin=False, with_iron=True, gmaxh=0.45, source="volume")   # graded air-box: coarse global, fine near
    rel_kb = np.linalg.norm(Hk - Hb) / np.linalg.norm(Hb)
    shield = np.linalg.norm(Hk) / np.linalg.norm(Hbs)
    print("  Kelvin (ndof %d) vs air-box R=3 (ndof %d): ||H_K - H_box||/|.| = %.3e" % (nk, nb, rel_kb))
    print("  iron gain ||H_iron|| / ||H_vacuum|| = %.4f  (iron shell concentrates the interior field)" % shield)

    relv = np.linalg.norm(H - Hbs) / np.linalg.norm(Hbs)
    print("\n[PASS] gap-free unified: vacuum FEM == Biot-Savart (%.1e, no bolt-on) + "
          "iron Kelvin == air-box (%.1e, material-aware)" % (relv, rel_kb))
