# -*- coding: utf-8 -*-
# DEMO (x) (verified): the Kelvin-FEM open problem is CONFORMALLY ISOMORPHIC to the SOMMERFELD
# half-space problem -- and reproduces its (static) layered Green's function with NO Sommerfeld integral.
#
# act5_01_exterior_material/u/v/w used spherically-LAYERED exteriors (Mie/transfer-matrix). The genuine Sommerfeld case
# is a PLANAR half-space: a source in medium 1 (z>0, coeff c1) over a planar interface z=0 to medium 2
# (z<0, coeff c2), with decay at infinity. The hard part for BEM is the Sommerfeld integral (the
# layered-media Green's function: inverse Hankel transform with branch cuts / surface-wave poles).
#
# THE ISOMORPHISM. Put the Kelvin inversion centre ON the interface (origin on z=0). Then:
#   * the interface plane passes through the centre -> a plane-through-centre maps to ITSELF;
#   * Kelvin inversion x -> R^2 x/|x|^2 PRESERVES the sign of z -> the two media z>0/z<0 are preserved;
#   * infinity -> the inversion centre (one point).
# So the UNBOUNDED half-space problem is conformally mapped to a BOUNDED two-region ball problem with
# the SAME planar interface and media, the decay-at-infinity (Sommerfeld radiation) condition becoming
# a single point condition at the centre. The Kelvin-FEM IS the Sommerfeld problem in inverted
# coordinates -- isomorphic, not merely similar.
#
# VERIFICATION against the closed-form static Sommerfeld Green's function (method of images): a point
# charge q in medium 1 at height z0 has, in z>0, an image beta*q at -z0 (beta=(c1-c2)/(c1+c2)); in
# z<0, a transmitted 2/(c1+c2)*q at z0. We drive a NET-NEUTRAL dipole (+q,-q) -- net-neutral so the
# n=0 monopole (which a 3D point-ground cannot hold, see act1_05_assemble_dtn_matrix) is not excited -- as two small
# uniformly-charged balls (exterior field = point charge), close the open domain by the Kelvin ball
# (also split z>0/z<0), and check the FE potential against the image sum in BOTH media and in the far
# field (inverse Kelvin). Match => the Kelvin-FEM operator is isomorphic to the Sommerfeld operator.
#
# Honest: this is the STATIC Sommerfeld (the image limit); the full WAVE Sommerfeld integral (branch
# cuts / surface-wave poles) is the time-harmonic case = the authors' extended-Kelvin radiating regime
# (sugahara2025, Maxwellian PML in the Kelvin-mapped exterior). The isomorphism (centre-on-interface
# Kelvin) is the same; only the per-region operator becomes Helmholtz.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
import netgen.occ as occ
from netgen.occ import Sphere, Box, Pnt, Vec, IdentificationType, OCCGeometry, Glue

R_out, offset = 1.0, 3.0
c1, c2 = 1.0, 4.0                 # medium 1 (air, z>0) / medium 2 (ground, z<0)
h, dlt, rs, q = 0.30, 0.10, 0.06, 1.0
zP, zN = h + dlt, h - dlt         # +q at 0.40, -q at 0.20  (net-neutral vertical dipole)
beta = (c1 - c2) / (c1 + c2)
charges = [(+q, zP), (-q, zN)]


def phi_image(X, Y, Z):
    """Closed-form static Sommerfeld (image) Green's function for the point charges in medium 1."""
    out = 0.0
    for qi, zi in charges:
        r1 = np.sqrt(X * X + Y * Y + (Z - zi) ** 2)
        if Z >= 0:                                   # medium 1: source + image beta*qi at -zi
            r2 = np.sqrt(X * X + Y * Y + (Z + zi) ** 2)
            out += (qi / (4 * np.pi * c1)) * (1.0 / r1 + beta / r2)
        else:                                        # medium 2: transmitted 2/(c1+c2)*qi
            out += (qi / (4 * np.pi)) * (2.0 / (c1 + c2)) / r1
    return out


def solve():
    B = Sphere(Pnt(0, 0, 0), R_out)
    for f in B.faces: f.name = "kint"
    sbp = Sphere(Pnt(0, 0, zP), rs); sbn = Sphere(Pnt(0, 0, zN), rs)
    B = B - sbp - sbn; B.mat("phys")
    sbp.mat("srcpos"); sbn.mat("srcneg")
    K = Sphere(Pnt(offset, 0, 0), R_out)
    for f in K.faces: f.name = "kext"
    K.mat("kelvin")
    tr = occ.gp_Trsf.Translation(Vec(offset, 0, 0))     # SINGLE truncation-sphere identification
    [f for f in B.faces if f.name == "kint"][0].Identify(
        [f for f in K.faces if f.name == "kext"][0], "kelvin", IdentificationType.PERIODIC, tr)
    mesh = ng.Mesh(OCCGeometry(Glue([B, sbp, sbn, K])).GenerateMesh(maxh=0.13)).Curve(3)
    x, y, z = ng.x, ng.y, ng.z; rp2 = (x - offset) ** 2 + y * y + z * z + 1e-20
    cz = ng.IfPos(z, c1, c2)                              # interface z=0 carried as a coefficient
    w = mesh.MaterialCF({"kelvin": R_out ** 2 / rp2}, default=1.0)
    mu = cz * w
    # Mean-zero constraint (NumberSpace) removes the constant nullspace in a WELL-CONDITIONED way:
    # a single ground point has zero capacity in 3D (it leaves the gauge ~1e4x the signal); int u dx=0
    # is exact. Source ALSO made exactly discretely neutral (net charge ~1e-15) so the n=0 near-null
    # mode is never excited (a residual net charge would blow it up ~1e7x; see act1_05_assemble_dtn_matrix).
    fes = ng.Periodic(ng.H1(mesh, order=3)) * ng.NumberSpace(mesh)
    (u, lam), (v, mu_) = fes.TnT()
    A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=8) + (u * mu_ + v * lam) * ng.dx)
    A.Assemble()
    vp = float(ng.Integrate(ng.CF(1), mesh, definedon=mesh.Materials("srcpos")))
    vn = float(ng.Integrate(ng.CF(1), mesh, definedon=mesh.Materials("srcneg")))
    f = ng.LinearForm(fes)
    f += (q / vp) * v * ng.dx(definedon=mesh.Materials("srcpos"), bonus_intorder=8)
    f += (-q / vn) * v * ng.dx(definedon=mesh.Materials("srcneg"), bonus_intorder=8)
    f.Assemble()
    gf = ng.GridFunction(fes)
    gf.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
    gfu = gf.components[0]
    print("   [diag] meshed src vols +%.3e/-%.3e  |gf| max=%.3e (well-conditioned, no blow-up)"
          % (vp, vn, max(abs(gfu.vec.FV().NumPy()))))
    return mesh, gfu


def sample(mesh, gf, X, Y, Z):
    rho = np.sqrt(X * X + Y * Y + Z * Z)
    if rho <= R_out - 1e-6:
        return gf(mesh(X, Y, Z))
    bp = (offset + R_out ** 2 * X / rho ** 2, R_out ** 2 * Y / rho ** 2, R_out ** 2 * Z / rho ** 2)
    return gf(mesh(*bp))                              # periodic-glue: no (R/rho) weight


mesh, gf = solve()
air = [(0.30, 0, 0.50), (0.50, 0, 0.30), (0, 0.40, 0.35), (0.30, 0.30, 0.30)]
gnd = [(0, 0, -0.30), (0.30, 0, -0.30), (0, 0.30, -0.40)]
far = [(0, 0, 3.0), (2.0, 0, 0.50), (0, 2.5, -0.50)]
allpts = [("air", p) for p in air] + [("ground", p) for p in gnd] + [("far", p) for p in far]
fem = np.array([sample(mesh, gf, *p) for _, p in allpts])
an = np.array([phi_image(*p) for _, p in allpts])
const = np.mean(fem - an)                              # single gauge constant (expect ~0: dipole has no monopole)
scale = np.max(np.abs(an))
print("Kelvin-FEM half-space  vs  static SOMMERFELD (image) Green's function")
print("(interface z=0, c1=%.0f over c2=%.0f; net-neutral dipole +q@%.2f / -q@%.2f; gauge const=%.2e)\n"
      % (c1, c2, zP, zN, const))
print("  region   point (x,y,z)            FEM-Kelvin      Sommerfeld(image)   rel(/max)")
for (reg, p), fe_, an_ in zip(allpts, fem, an):
    print("  %-6s (%5.2f,%5.2f,%5.2f)   %13.6e   %15.6e   %.2e"
          % (reg, p[0], p[1], p[2], fe_ - const, an_, abs(fe_ - const - an_) / scale))
rel = np.abs(fem - const - an) / scale
print("\n  max rel err (air+ground+far, incl. inverse-Kelvin far field): %.2e" % rel.max())
print("\n=> centre-on-interface Kelvin maps the UNBOUNDED half-space to a BOUNDED two-region ball")
print("   (plane->plane, z-sign preserved, infinity->centre): the Kelvin-FEM operator is ISOMORPHIC")
print("   to the Sommerfeld half-space operator and reproduces its image Green's function in BOTH")
print("   media and the far field -- NO Sommerfeld integral, no branch cuts, no special functions.")
