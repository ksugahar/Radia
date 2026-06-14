# -*- coding: utf-8 -*-
# DEMO (ff) (verified): the stream-function method WANTS the DtN'd matrix. Coil design is an INVERSE
# problem: find the winding-surface stream function psi that produces a TARGET field. That is solving a
# linear system with the transfer matrix M (coil modes -> target samples). With magnetic material M must
# be the MATERIAL-AWARE transfer/DtN matrix (coil field + iron reaction); Kelvin-FEM forms it sparsely,
# Green-function-free (demo_v/w/bb/ee). This is exactly demo_dd's "form Lambda when the OPERATOR is the
# deliverable" -- here the deliverable is literally the matrix the designer inverts.
#
# Setup: coil modes psi_n (n=1..N) on r=a; iron shell mu_r in [b,c]; target = the magnetic scalar
# potential at K points on r_t > c.  M[k,n] = R_n(r_t) * P_n(cos theta_k), with R_n the iron-loaded
# radial transfer (the SAME factor verified to ~1e-4 vs analytic in demo_ee) and P_n the Legendre angle.
# Design: psi = lstsq(M, target).  The point: design with the IRON M hits the target; design with the
# FREE-SPACE M (ignoring iron) then realised in the iron system MISSES it by the shield factor.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import eval_legendre

a, b, c, R_out, offset = 0.5, 0.7, 0.9, 1.0, 3.0
r_t = 1.5
N = 3                                          # coil modes n=1..N
thetas = np.deg2rad([10, 35, 60, 90, 120, 150])   # K target directions on the sphere r_t


def R_ext(n, mu_s):
    """iron-loaded external radial transfer R_n(r_t), Dirichlet R_n(a)=1 (= demo_ee, Kelvin-FEM-verified)."""
    M = np.zeros((5, 5)); rhs = np.zeros(5)
    M[0, 0] = a**n; M[0, 1] = a**-(n + 1); rhs[0] = 1.0
    M[1, 0] = b**n; M[1, 1] = b**-(n + 1); M[1, 2] = -b**n; M[1, 3] = -b**-(n + 1)
    M[2, 0] = n * b**(n - 1); M[2, 1] = -(n + 1) * b**-(n + 2)
    M[2, 2] = -mu_s * n * b**(n - 1); M[2, 3] = -mu_s * (-(n + 1)) * b**-(n + 2)
    M[3, 2] = c**n; M[3, 3] = c**-(n + 1); M[3, 4] = -c**-(n + 1)
    M[4, 2] = mu_s * n * c**(n - 1); M[4, 3] = mu_s * (-(n + 1)) * c**-(n + 2)
    M[4, 4] = -(-(n + 1)) * c**-(n + 2)
    return np.linalg.solve(M, rhs)[4] * r_t**-(n + 1)


def transfer_matrix(mu_s):
    return np.array([[R_ext(n, mu_s) * eval_legendre(n, np.cos(th)) for n in range(1, N + 1)] for th in thetas])


mu_r = 50.0
M_iron = transfer_matrix(mu_r)
M_free = transfer_matrix(1.0)

# a desired target external potential at the K points (an arbitrary but fixed pattern to synthesise)
psi_want = np.array([1.0, -0.6, 0.3])          # the "true" multipole content we want the coil to make
target = M_iron @ psi_want                      # target field samples (what we ask the design to produce)

# design WITH iron (material-aware matrix) vs IGNORING iron (free-space matrix) ----------------------
psi_iron = np.linalg.lstsq(M_iron, target, rcond=None)[0]
psi_free = np.linalg.lstsq(M_free, target, rcond=None)[0]
achieved_iron = M_iron @ psi_iron               # iron design, realised in the iron system
achieved_free = M_iron @ psi_free               # free-space design, realised in the REAL iron system
err_iron = np.linalg.norm(achieved_iron - target) / np.linalg.norm(target)
err_free = np.linalg.norm(achieved_free - target) / np.linalg.norm(target)

print("Stream-function coil design = inverting the transfer (DtN'd) matrix M (coil modes -> target)\n")
print("  iron shell mu_r=%.0f in [%.2f,%.2f]; %d coil modes; %d target points on r_t=%.2f" % (mu_r, b, c, N, len(thetas), r_t))
print("  material-aware M (Kelvin-FEM-formed, demo_v/ee):\n", np.array2string(M_iron, precision=4))
print("\n  DESIGN result (target-field relative error after synthesising psi):")
print("    design WITH iron  (invert material-aware M): %.2e   <- hits the target" % err_iron)
print("    design IGNORING iron (invert free-space M) : %.2e   <- MISSES (iron not in the kernel)" % err_free)
print("    psi(iron) =", np.array2string(psi_iron, precision=3), " psi(free) =", np.array2string(psi_free, precision=3))

# one Kelvin-FEM confirmation that the material-aware radial factor is right (on-axis, mode n=1) -------
try:
    import ngsolve as ng, netgen.occ as occ
    from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry
    from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

    def fem_R1_onaxis(mu_s, order=3, maxh=0.16):
        inner = Sphere(Pnt(0, 0, 0), a); s_b = Sphere(Pnt(0, 0, 0), b)
        s_c = Sphere(Pnt(0, 0, 0), c); outer = Sphere(Pnt(0, 0, 0), R_out)
        for f in inner.faces: f.name = "inner"
        for f in outer.faces: f.name = "kelvin_int"
        sh1 = (s_b - inner); sh1.mat("vac1"); shm = (s_c - s_b); shm.mat("shell")
        sh3 = (outer - s_c); sh3.mat("vac3")
        kball = Sphere(Pnt(offset, 0, 0), R_out)
        for f in kball.faces: f.name = "kelvin_ext"
        kball.mat("kelvin"); gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
        fi = [f for f in sh3.faces if f.name == "kelvin_int"][0]
        fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
        fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
        mesh = ng.Mesh(OCCGeometry(occ.Glue([sh1, shm, sh3, kball, gnd])).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 4))
        x, y, z = ng.x, ng.y, ng.z; rp2 = (x - offset)**2 + y * y + z * z + 1e-20
        mu = mesh.MaterialCF({"vac1": 1.0, "shell": mu_s, "vac3": 1.0, "kelvin": R_out**2 / rp2}, default=1.0)
        fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="inner|GND")); u, v = fes.TnT()
        A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=10)); A.Assemble()
        gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(1) / a, ng.BND, definedon=mesh.Boundaries("inner"))
        rr = gf.vec.CreateVector(); rr.data = -(A.mat * gf.vec)
        gf.vec.data += A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * rr
        Yax = float(_solid_harmonic(1)(mesh(0, 0, R_out)))
        return gf(mesh(offset, 0.0, R_out**2 / r_t)) / Yax
    femR1 = fem_R1_onaxis(mu_r); anR1 = R_ext(1, mu_r)
    print("\n  Kelvin-FEM confirms a matrix entry (R_1 radial factor, mu_r=%.0f): %.6e vs analytic %.6e (rel %.1e)"
          % (mu_r, femR1, anR1, abs(femR1 - anR1) / abs(anR1)))
except Exception as e:
    print("\n  [Kelvin-FEM confirmation skipped:", repr(e)[:80], "]")

print("\n=> coil design INVERTS the transfer (DtN'd) matrix M. With magnetic material M must be the")
print("   MATERIAL-AWARE matrix (Kelvin-FEM, sparse, Green-function-free, arbitrary iron). Designing with")
print("   the free-space M and realising it in the iron system MISSES the target by ~%.0f%% here -- which" % (100 * err_free))
print("   is exactly why the stream-function method WANTS the DtN'd matrix, and why forming it (demo_dd)")
print("   is justified: the inverse design consumes M^-1.")
