# -*- coding: utf-8 -*-
# DEMO (oo) (verified): the matched HOIBC as a GENUINE 3D SURFACE FE term (Laplace-Beltrami grad_Gamma)
# -- promoting the per-mode eigenvalue of demo_kk/ll/mm to the real surface operator in NGSolve.
#
# demo_mm realised the radiating extended-Kelvin boundary as a radial FE, with the angular
# Laplace-Beltrami Delta_S entering only as its eigenvalue -n(n+1). The genuinely 3D ingredient is the
# SURFACE term itself: the matched HOIBC impedance operator on the inner image sphere is
#       Z_HOIBC = (i kb - 1) I + (i/2kb) Delta_S        (Delta_S = unit-sphere Laplace-Beltrami),
# assembled as the surface bilinear form (weak Delta_S = -grad_Gamma . grad_Gamma):
#       S(u,v) = (i kb - 1) int_Gamma u v ds  -  (i/2kb) rho_b^2 int_Gamma grad_Gamma u . grad_Gamma v ds
# (the rho_b^2 converts the radius-rho_b surface gradient to the unit-sphere Delta_S). This demo
# assembles that term in NGSolve and verifies, two ways, that it IS the matched HOIBC:
#   (1) the surface gradient grad_Gamma reproduces the Laplace-Beltrami SPECTRUM n(n+1) with the
#       correct 2n+1 multiplicities (so grad_Gamma realises Delta_S);
#   (2) the full HOIBC surface form S has generalized spectrum (vs the surface mass) equal to the
#       per-mode matched impedance Lambda_HOIBC,n = i kb - 1 - i n(n+1)/(2 kb), with 2n+1 multiplicities.
# Together with demo_mm (which verified the transformation-optics MEDIUM + the radial coupling to the
# truncation, converging O(h^2) to the closed form) this means the radiating extended-Kelvin boundary
# is now a genuine 3D FE: an isotropic (a/rho)^2 volume medium + this Delta_S surface term. The
# remaining composition (one monolithic volumetric 3D solve) just glues the two verified pieces.
#
# VERIFICATIONS (asserted from computed values; surface-mesh discretization ~1-3% at this maxh):
#  (1) grad_Gamma surface stiffness vs surface mass -> Laplace-Beltrami eigenvalues {0,2,6,12} with
#      multiplicities {1,3,5,7} (i.e. n(n+1), 2n+1).
#  (2) the HOIBC surface form S -> the matched impedance Lambda_HOIBC,n per degree n (2n+1 each).
#
# Needs ngsolve + netgen.occ (a genuine surface FE), numpy/scipy.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
import scipy.sparse as sp
from scipy.linalg import eigh, eig


def _csr(m, ndof):
    r, c, val = m.COO()
    return sp.csr_matrix((np.array(val), (np.array(r), np.array(c))), shape=(ndof, ndof)).toarray()


def surface_forms(radius, order=3, maxh=0.18, complex_fes=False):
    """Assemble the surface-gradient stiffness K = int grad_Gamma u . grad_Gamma v ds and the surface
    mass M = int u v ds on a radius-`radius` sphere (a genuine 2-manifold FE mesh in 3D)."""
    face = Sphere(Pnt(0, 0, 0), radius).faces[0]
    mesh = ng.Mesh(OCCGeometry(face).GenerateMesh(maxh=maxh * radius))
    fes = ng.H1(mesh, order=order, complex=complex_fes)
    u, v = fes.TnT()
    K = ng.BilinearForm(ng.grad(u).Trace() * ng.grad(v).Trace() * ng.ds); K.Assemble()
    M = ng.BilinearForm(u.Trace() * v.Trace() * ng.ds); M.Assemble()
    return _csr(K.mat, fes.ndof), _csr(M.mat, fes.ndof)


print("=" * 80)
print("The matched HOIBC as a GENUINE 3D SURFACE FE term (Laplace-Beltrami grad_Gamma) in NGSolve")
print("=" * 80)

# (1) grad_Gamma realises the Laplace-Beltrami spectrum n(n+1) -------------------------------
print("\n(1) surface-gradient FE spectrum on the unit sphere -> Laplace-Beltrami eigenvalues n(n+1):")
K1, M1 = surface_forms(1.0, order=3, maxh=0.20)
w = np.sort(eigh(K1, M1, eigvals_only=True).real)
print("    lowest eigenvalues:", np.round(w[:9], 3))
groups = [(0, 1, 0.0), (1, 3, 2.0), (2, 5, 6.0), (3, 7, 12.0)]
idx = 0
for n, mult, val in groups:
    grp = w[idx:idx + mult]; idx += mult
    err = np.max(np.abs(grp - val)) / max(val, 1.0)
    print("    n=%d  mult=%d  eig~%.3f (exact %g)  rel.err=%.2e" % (n, mult, grp.mean(), val, err))
    assert err < 5e-2, "grad_Gamma must reproduce the Laplace-Beltrami eigenvalue n(n+1)"
print("    => grad_Gamma on the surface FE IS the unit-sphere Delta_S (eigenvalues n(n+1), 2n+1 each).")

# (2) the HOIBC surface form has spectrum = the matched per-mode impedance Lambda_HOIBC,n ----
a, k, b = 1.0, 4.0, 2.0
kb = k * b
rb = a * a / b
def Lam_hoibc(n):
    return 1j * kb - 1.0 - 1j * n * (n + 1.0) / (2.0 * kb)
print("\n(2) HOIBC SURFACE form S=(ikb-1)M-(i/2kb)rb^2 K on the radius-rb=%.2f sphere (kb=%.0f):" % (rb, kb))
print("    its generalized spectrum vs the surface mass must equal Lambda_HOIBC,n per degree n:")
Kc, Mc = surface_forms(rb, order=3, maxh=0.18, complex_fes=True)
S = (1j * kb - 1.0) * Mc - (1j / (2 * kb)) * rb * rb * Kc       # unit-sphere Delta_S via rb^2 factor
ev = eig(S, Mc, right=False)
print("     n  mult   Lambda_HOIBC,n        FE surface-form eig (mean)   max|err|")
for n in range(0, 4):
    L = Lam_hoibc(n)
    d = np.abs(ev - L); sel = np.argsort(d)[:2 * n + 1]
    err = d[sel].max()
    m = ev[sel].mean()
    print("    %2d   %d    (%7.3f,%7.3f)    (%7.3f,%7.3f)        %.2e"
          % (n, 2 * n + 1, L.real, L.imag, m.real, m.imag, err))
    assert err < 5e-2, "the HOIBC surface form must realise Lambda_HOIBC,n per mode"
print("    => the Delta_S surface form encodes the multipole-dependent matched HOIBC impedance")
print("       (2n+1 multiplicities), to surface-mesh accuracy -- the genuine 3D HOIBC term.")

print("\n" + "=" * 80)
print("RESULT: the matched HOIBC is now a genuine 3D SURFACE FE operator -- the surface gradient")
print("grad_Gamma reproduces Delta_S (eigenvalues n(n+1)) and the assembled HOIBC surface form has the")
print("matched per-mode impedance Lambda_HOIBC,n. With demo_mm's transformation-optics volume medium +")
print("radial coupling (O(h^2) to the closed form), the radiating extended-Kelvin boundary is a")
print("genuine 3D FE: isotropic (a/rho)^2 medium + this Delta_S surface term. Remaining: the single")
print("monolithic volumetric 3D solve gluing the two (each piece now verified).")
print("=" * 80)
