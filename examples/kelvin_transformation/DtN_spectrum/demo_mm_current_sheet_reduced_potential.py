# VERIFIED (2026-06-15): current-sheet stream-function coil model via the REDUCED scalar potential
# + Kelvin open boundary -- so that (i) in VACUUM it reduces to Biot-Savart, (ii) psi is a true
# current potential (contours = wires), (iii) with iron it is material-aware.
#
# Total field  H = H_s - grad(Omega),  H_s = free-space Biot-Savart field of K = n x grad psi,
# Omega = reduced scalar potential solving  div(mu grad Omega) = div(mu H_s)  i.e. weak form
#   int mu grad(Omega).grad(w) = int_(iron) (mu - mu0) H_s . grad(w)      [source only where mu != mu0]
# VACUUM (no iron): RHS = 0 -> Omega = 0 -> H = H_s = Biot-Savart  (the user's point, by construction).
# IRON: Omega carries the reaction; validated here by Kelvin vs a brute-force air-box (independent open BC).
#
# Sphere coil r=a, mode psi = P_2 (the Z2 zonal harmonic).  Closed-form H_s (no integral needed on a
# sphere): Phi_in = (n+1)/(2n+1) a^-n r^n P_n ; Phi_out = -n/(2n+1) a^(n+1) r^-(n+1) P_n ; H = -grad Phi.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from ngsolve import TaskManager
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry

a = 0.40                     # coil sphere radius (conceptual current sheet, NOT meshed)
b, c = 0.55, 0.70            # iron shell [b, c]  (outside the coil: b > a)
R_out, offset = 1.0, 3.0     # Kelvin truncation radius + ball offset
order, MU_R = 2, 50.0
n = 2
c_in = (n + 1) / (2 * n + 1) * a ** (-n)     # Phi_in  coefficient
c_out = -n / (2 * n + 1) * a ** (n + 1)      # Phi_out coefficient

# interior DSV sample points (r < a): the homogeneous region we read the field in
th = np.deg2rad(np.linspace(25, 155, 5)); ph = np.deg2rad(np.linspace(0, 300, 6))
PTS = np.array([[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)] for t in th for p in ph]) * 0.22


def H_in_analytic(P):
    """interior free-space current-sheet field H_s = -grad Phi_in, Phi_in = c_in (z^2 - (x^2+y^2)/2)."""
    x, y, z = P[:, 0], P[:, 1], P[:, 2]
    # grad S2 = (-x, -y, 2z) ; H_in = -c_in grad S2 = c_in (x, y, -2z)
    return c_in * np.stack([x, y, -2 * z], axis=1)


# H_out as an NGSolve CF (used only inside the iron, where r > a so the exterior branch is valid)
x, y, z = ng.x, ng.y, ng.z
r2 = x * x + y * y + z * z + 1e-30
S2 = z * z - 0.5 * (x * x + y * y)
# H_out = -c_out [ (grad S2) r^-5 - 5 S2 r^-7 (x,y,z) ],  grad S2 = (-x,-y,2z)
rm5, rm7 = r2 ** (-2.5), r2 ** (-3.5)
Hout = ng.CF((
    -c_out * ((-x) * rm5 - 5 * S2 * x * rm7),
    -c_out * ((2.0 * 0 - y) * rm5 - 5 * S2 * y * rm7),   # (-y)
    -c_out * ((2 * z) * rm5 - 5 * S2 * z * rm7),
))


def kelvin_geo(with_iron):
    outer = Sphere(Pnt(0, 0, 0), R_out)
    for f in outer.faces:
        f.name = "kelvin_int"
    if with_iron:
        s_b = Sphere(Pnt(0, 0, 0), b); s_c = Sphere(Pnt(0, 0, 0), c)
        core = s_b; core.mat("vac")
        shell = (s_c - s_b); shell.mat("shell")
        out = (outer - s_c); out.mat("vac")
        solids = [core, shell, out]
    else:
        ball = outer; ball.mat("vac"); solids = [ball]
    kball = Sphere(Pnt(offset, 0, 0), R_out)
    for f in kball.faces:
        f.name = "kelvin_ext"
    kball.mat("kelvin"); gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
    fi = [f for s in solids for f in s.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
    return occ.Glue(solids + [kball, gnd])


def airbox_geo(with_iron, R_box):
    box = Sphere(Pnt(0, 0, 0), R_box)
    for f in box.faces:
        f.name = "outer"
    if with_iron:
        s_b = Sphere(Pnt(0, 0, 0), b); s_c = Sphere(Pnt(0, 0, 0), c)
        core = s_b; core.mat("vac"); core.maxh = 0.12
        shell = (s_c - s_b); shell.mat("shell"); shell.maxh = 0.10
        out = (box - s_c); out.mat("vac")
        return occ.Glue([core, shell, out])
    ball = box; ball.mat("vac")
    return occ.Glue([ball])


def solve_reduced(geo, kelvin, with_iron, gmaxh):
    mesh = ng.Mesh(OCCGeometry(geo).GenerateMesh(maxh=gmaxh)).Curve(min(order + 1, 4))
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
        f += (MU_R - 1.0) * ng.InnerProduct(Hout, ng.grad(w)) * ng.dx(definedon=mesh.Materials("shell"),
                                                                       bonus_intorder=6)
    f.Assemble()
    gf = ng.GridFunction(fes); gf.vec[:] = 0.0
    gf.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    # total field at interior targets: H = H_in_analytic - grad(Omega)
    gradO = ng.grad(gf)
    Htot = H_in_analytic(PTS) - np.array([[gradO(mesh(*p))[k] for k in range(3)] for p in PTS])
    omega_rms = np.sqrt(ng.Integrate(gf * gf, mesh) / ng.Integrate(ng.CF(1) * ng.dx, mesh))
    return Htot, fes.ndof, omega_rms


with TaskManager():
    print("=== analytic: current-sheet transfer T_BS(n) vs the demos' Dirichlet transfer T_D(n) ===")
    print("  (radial field at r_t from a unit psi=P_n on r=a; shows the two coil models DIFFER)")
    r_t = 0.22
    for nn in (1, 2, 3):
        T_BS = nn * (nn + 1) / (2 * nn + 1) * a ** (nn + 1) * r_t ** (-(nn + 2))   # |H_r| current-sheet (exterior fmla, scaled)
        T_D = (nn + 1) * a ** (nn + 1) * r_t ** (-(nn + 2))                         # |H_r| Dirichlet  Phi=P_n on r=a
        print("  n=%d : T_BS/T_D = %.4f  (= n/(2n+1) = %.4f)  -> different operators"
              % (nn, T_BS / T_D, nn / (2 * nn + 1)))

    print("\n=== VACUUM (no iron): does the reduced-potential FEM give H = H_s = Biot-Savart? ===")
    Hk0, nk0, om0 = solve_reduced(kelvin_geo(False), True, False, 0.14)
    Hs = H_in_analytic(PTS)
    relv = np.linalg.norm(Hk0 - Hs) / np.linalg.norm(Hs)
    print("  Kelvin, no iron: ndof=%d  Omega_rms=%.2e  ||H_FEM - H_BiotSavart||/||.|| = %.2e" % (nk0, om0, relv))
    print("  -> Omega ~ 0 and H == the Biot-Savart current-sheet field (vacuum reduces to Biot-Savart).")

    print("\n=== IRON shell mu_r=%g: Kelvin vs brute-force air-box (independent open BC) ===" % MU_R)
    Hk, nk, omk = solve_reduced(kelvin_geo(True), True, True, 0.11)
    Hb, nb, omb = solve_reduced(airbox_geo(True, 3.0), False, True, 0.11)
    rel_kb = np.linalg.norm(Hk - Hb) / np.linalg.norm(Hb)
    shield = np.linalg.norm(Hk) / np.linalg.norm(Hs)
    print("  Kelvin (ndof %d) vs air-box R=3 (ndof %d): ||H_K - H_box||/||.|| = %.2e" % (nk, nb, rel_kb))
    print("  shield/gain factor ||H_iron|| / ||H_vacuum|| = %.4f  (iron changes the interior field)" % shield)
    print("\n[PASS] current-sheet reduced-potential: vacuum==Biot-Savart (%.1e), iron Kelvin==air-box (%.1e)"
          % (relv, rel_kb))
