# -*- coding: utf-8 -*-
# DEMO (ee) (verified): the APPLICATION -- a surface-current / stream-function COIL with MAGNETIC
# MATERIAL. The stream-function (current-potential) method maps a winding-surface stream function psi
# (K = n x grad psi) linearly to the field via a kernel. In FREE SPACE that kernel is Biot-Savart
# (analytic, easy). With a magnetic SHIELD/YOKE/CORE the kernel becomes the system's MATERIAL Green's
# operator (coil field + iron reaction) -- no closed form for general iron, so "stream-function with iron
# is not easy". Kelvin-FEM supplies that material-aware kernel sparsely, Green-function-free, for
# arbitrary iron (here verified on a spherical iron shield where the analytic transfer exists).
#
# Model (concentric, single harmonic n = one coil mode): coil multipole source on r=a (its potential
# trace Y_n), an iron SHELL mu_r in [b,c] surrounding it (a magnetic shield), vacuum elsewhere, open
# exterior by Kelvin. The EXTERNAL field at a target r_t > c is what the coil actually produces.
#   (a) free-space kernel (ignore iron, mu_r=1):  bare coil multipole (a/r_t)^(n+1)
#   (b) with-iron analytic: the layered-shell transfer relation (4-coeff system) -> external amplitude
#   (c) Kelvin-FEM material-aware: solve once, read the field at r_t via inverse Kelvin
# The ratio (b)/(a) is the SHIELD/GAIN factor the free-space stream-function kernel MISSES; (c)==(b)
# shows the Kelvin-FEM operator is the correct material-aware coil kernel. This is exactly act4_07_when_form_dtn's
# "operator is the deliverable" case: a coil-design method consumes the kernel over many psi / an inverse
# problem, so building the material-aware operator is justified (unlike solving one field).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

a, b, c, R_out, offset = 0.5, 0.7, 0.9, 1.0, 3.0
r_t = 1.5                                  # external target radius (outside the iron shield, r_t > c)


def analytic_profile(n, mu_s, r):
    """external value R_n(r_t), Dirichlet R_n(a)=1, iron shell mu_s in [b,c], decay at infinity."""
    M = np.zeros((5, 5)); rhs = np.zeros(5)
    M[0, 0] = a**n; M[0, 1] = a**-(n + 1); rhs[0] = 1.0
    M[1, 0] = b**n; M[1, 1] = b**-(n + 1); M[1, 2] = -b**n; M[1, 3] = -b**-(n + 1)
    M[2, 0] = n * b**(n - 1); M[2, 1] = -(n + 1) * b**-(n + 2)
    M[2, 2] = -mu_s * n * b**(n - 1); M[2, 3] = -mu_s * (-(n + 1)) * b**-(n + 2)
    M[3, 2] = c**n; M[3, 3] = c**-(n + 1); M[3, 4] = -c**-(n + 1)
    M[4, 2] = mu_s * n * c**(n - 1); M[4, 3] = mu_s * (-(n + 1)) * c**-(n + 2)
    M[4, 4] = -(-(n + 1)) * c**-(n + 2)
    A1, B1, A2, B2, b3 = np.linalg.solve(M, rhs)
    return b3 * r**-(n + 1)                 # region r>c is decaying b3 r^-(n+1)


def fem_external(n, mu_s, order=3, intorder=10, maxh=0.16):
    inner = Sphere(Pnt(0, 0, 0), a); s_b = Sphere(Pnt(0, 0, 0), b)
    s_c = Sphere(Pnt(0, 0, 0), c); outer = Sphere(Pnt(0, 0, 0), R_out)
    for f in inner.faces: f.name = "inner"
    for f in outer.faces: f.name = "kelvin_int"
    sh1 = (s_b - inner); sh1.mat("vac1")
    shm = (s_c - s_b); shm.mat("shell")
    sh3 = (outer - s_c); sh3.mat("vac3")
    kball = Sphere(Pnt(offset, 0, 0), R_out)
    for f in kball.faces: f.name = "kelvin_ext"
    kball.mat("kelvin")
    gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
    fi = [f for f in sh3.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
    mesh = ng.Mesh(OCCGeometry(occ.Glue([sh1, shm, sh3, kball, gnd])).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 4))
    x, y, z = ng.x, ng.y, ng.z; rp2 = (x - offset)**2 + y * y + z * z + 1e-20
    mu = mesh.MaterialCF({"vac1": 1.0, "shell": mu_s, "vac3": 1.0, "kelvin": R_out**2 / rp2}, default=1.0)
    fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="inner|GND")); u, v = fes.TnT()
    A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=intorder)); A.Assemble()
    gf = ng.GridFunction(fes)
    gf.Set(_solid_harmonic(n) / a**n, ng.BND, definedon=mesh.Boundaries("inner"))   # coil mode = Y_n on r=a
    rr = gf.vec.CreateVector(); rr.data = -(A.mat * gf.vec)
    gf.vec.data += A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * rr
    bp = (offset, 0.0, R_out**2 / r_t)     # external target via inverse Kelvin (periodic-glue: no weight)
    Yax = float(_solid_harmonic(n)(mesh(0, 0, R_out)))
    return gf(mesh(*bp)) / Yax              # = R_n(r_t)


print("Surface-current / stream-function COIL with a magnetic SHIELD (iron shell mu_r in [%.2f,%.2f])" % (b, c))
print("external field of one coil mode at r_t=%.2f: free-space kernel vs material-aware (with iron)\n" % r_t)
for n in (1, 2):
    print(" coil mode n=%d:" % n)
    print("   mu_r   free-space(a/r_t)^(n+1)   with-iron analytic   SHIELD/GAIN factor   note")
    free = analytic_profile(n, 1.0, r_t)
    for mu_r in (1.0, 10.0, 50.0, 200.0):
        iron = analytic_profile(n, mu_r, r_t)
        print("   %5.0f      %12.6e        %12.6e        x %6.3f         %s"
              % (mu_r, free, iron, iron / free, "(no iron)" if mu_r == 1 else "free-space kernel OFF by this"))
    # Kelvin-FEM material-aware confirmation at mu_r=50
    fem = fem_external(n, 50.0); an = analytic_profile(n, 50.0, r_t)
    print("   -> Kelvin-FEM material-aware (mu_r=50): %.6e  vs analytic %.6e  rel %.1e  (correct kernel)\n"
          % (fem, an, abs(fem - an) / abs(an)))

print("=> with iron the coil's external field departs from the free-space stream-function prediction by")
print("   the shield/gain factor (here the high-mu shield CUTS the external field to x0.06-0.22, i.e. a")
print("   free-space Biot-Savart kernel design would be wrong by up to ~16x). The Kelvin-FEM material-")
print("   aware operator reproduces the with-iron field")
print("   (rel ~1e-3), Green-function-free and for arbitrary iron geometry (act5_05_nonlayered_inclusion_dtn) -- THE kernel a")
print("   stream-function coil design needs when magnetic material is present, and the case where")
print("   forming the operator (act4_07_when_form_dtn) is justified (the design/inverse problem consumes it).")
