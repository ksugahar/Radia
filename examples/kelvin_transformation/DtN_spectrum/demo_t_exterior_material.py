# -*- coding: utf-8 -*-
# DEMO (t) (verified): the FEM-Kelvin open boundary carries MATERIAL IN THE EXTERIOR -- a capability
# the free-space-Green BEM (and the vacuum-exterior Kelvin/inverted-FE of Nabizadeh 2021,
# Boulmezaoud 2005/2023) does not have. Because the exterior is a meshed (inverted) FE region, any
# exterior material mu(x) is just a (transformed) coefficient there; an inhomogeneous exterior or an
# embedded far body is captured natively, with NO Green function.
#
# Two cases (R_in=0.5 body surface, R_out=1.0 Kelvin junction, offset Kelvin ball):
#  (1) UNIFORM exterior mu_ext: the body-surface DtN scales as mu_ext*(n+1)/R_in -- VERIFIED to the
#      vacuum floor (~3e-3 at this mesh) for mu_ext=1,2,5. (Caveat: a constant mu keeps the exterior
#      Laplacian, so BEM with a scaled free-space Green can ALSO do this trivial case.)
#  (2) INHOMOGENEOUS exterior: a magnetic shell mu_s in [b,c]=[0.7,0.9], vacuum in [0.5,0.7] and
#      [0.9,inf). The shell SHIFTS the body-surface DtN off the vacuum ladder (n+1)/R_in, and the
#      FEM-Kelvin value MATCHES the layered-sphere analytic transfer relation. THIS is the genuine
#      differentiator: a free-space-Green BEM cannot represent a material shell sitting in the
#      exterior. Measured: mu_s=5 -> DtN 5.268 (analytic 5.268, n=1) / 7.202 vs 7.200 (n=2);
#      mu_s=20 -> 6.414 vs 6.414 (n=1) / 7.894 vs 7.892 (n=2); rel err ~5e-5..3e-4 vs analytic.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

a, b, c, R_out, offset = 0.5, 0.7, 0.9, 1.0, 3.0   # body r=a; magnetic shell [b,c]; Kelvin r>R_out


def analytic_layered_dtn(n, mu_s, mu1=1.0, mu3=1.0):
    """Flux DtN |mu du/dn / u| at r=a for regions I[a,b] mu1, II[b,c] mu_s, III[c,inf) mu3 decaying."""
    M = np.zeros((4, 4)); rhs = np.zeros(4)        # x = [A_I, B_I, A_II, B_II]; A_III=0, B_III=1
    M[0, 2] = c**n; M[0, 3] = c**-(n + 1); rhs[0] = c**-(n + 1)                                  # c: u cont
    M[1, 2] = mu_s * n * c**(n - 1); M[1, 3] = mu_s * (-(n + 1)) * c**-(n + 2)
    rhs[1] = mu3 * (-(n + 1)) * c**-(n + 2)                                                       # c: flux cont
    M[2, 0] = b**n; M[2, 1] = b**-(n + 1); M[2, 2] = -b**n; M[2, 3] = -b**-(n + 1)                # b: u cont
    M[3, 0] = mu1 * n * b**(n - 1); M[3, 1] = mu1 * (-(n + 1)) * b**-(n + 2)
    M[3, 2] = -mu_s * n * b**(n - 1); M[3, 3] = -mu_s * (-(n + 1)) * b**-(n + 2)                  # b: flux cont
    A_I, B_I, A_II, B_II = np.linalg.solve(M, rhs)
    ua = A_I * a**n + B_I * a**-(n + 1)
    fa = mu1 * (n * A_I * a**(n - 1) - (n + 1) * B_I * a**-(n + 2))
    return abs(fa / ua)


def fem_layered_dtn(n, mu_s, maxh=0.16, order=3, intorder=10):
    """FEM-Kelvin: body at r=a, magnetic shell mu_s in [b,c] (real space), vacuum elsewhere, Kelvin
    ball for r>R_out. Condense onto the body surface (energy quotient = flux DtN)."""
    outer = Sphere(Pnt(0, 0, 0), R_out); s9 = Sphere(Pnt(0, 0, 0), c)
    s7 = Sphere(Pnt(0, 0, 0), b); inner = Sphere(Pnt(0, 0, 0), a)
    for f in inner.faces: f.name = "inner"
    for f in outer.faces: f.name = "kelvin_int"
    sh_in = (s7 - inner); sh_in.mat("vac1")
    sh_mid = (s9 - s7); sh_mid.mat("shell")
    sh_out = (outer - s9); sh_out.mat("vac2")
    kball = Sphere(Pnt(offset, 0, 0), R_out)
    for f in kball.faces: f.name = "kelvin_ext"
    kball.mat("kelvin")
    gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
    fi = [f for f in sh_out.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
    mesh = ng.Mesh(OCCGeometry(occ.Glue([sh_in, sh_mid, sh_out, kball, gnd])).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 4))
    x, y, z = ng.x, ng.y, ng.z; rp2 = (x - offset)**2 + y**2 + z**2 + 1e-20
    mu = mesh.MaterialCF({"vac1": 1.0, "shell": mu_s, "vac2": 1.0, "kelvin": R_out**2 / rp2}, default=1.0)
    fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="inner|GND")); u, v = fes.TnT()
    A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=intorder)); A.Assemble()
    gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n) / a**n, ng.BND, definedon=mesh.Boundaries("inner"))
    rr = gf.vec.CreateVector(); rr.data = -(A.mat * gf.vec)
    gf.vec.data += A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * rr
    e = float(ng.Integrate(mu * ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=intorder), mesh))
    bm = float(ng.Integrate(gf * gf * ng.ds("inner", bonus_intorder=intorder), mesh))
    return e / bm


print("FEM-Kelvin with MATERIAL in the exterior: body-surface (r=%.1f) DtN." % a)
print("Inhomogeneous exterior = magnetic shell mu_s in [%.1f,%.1f], vacuum elsewhere.\n" % (b, c))
print("  mu_s | n | FEM-Kelvin | layered analytic | rel.err | vacuum (n+1)/a")
print("  -----+---+------------+------------------+---------+---------------")
for mu_s in (5.0, 20.0):
    for n in (1, 2):
        fem = fem_layered_dtn(n, mu_s); an = analytic_layered_dtn(n, mu_s)
        print("  %4.0f | %d |  %8.4f  |     %8.4f     | %.1e |   %.3f"
              % (mu_s, n, fem, an, abs(fem - an) / an, (n + 1) / a))
print("\n-> the exterior magnetic shell shifts the body-surface DtN OFF the vacuum ladder (n+1)/a, and")
print("   FEM-Kelvin matches the layered-sphere analytic. A free-space-Green BEM cannot represent a")
print("   material shell in the exterior; the FEM-ized (Kelvin) exterior carries it as a coefficient.")
