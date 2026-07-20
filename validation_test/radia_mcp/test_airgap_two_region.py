r"""AGE full-machine, step (2) -- TWO FE regions coupled across the un-meshed gap.

Steps 1-3 had ONE FE region + an analytic rotor ring.  A real machine has a rotor FE
region AND a stator FE region, coupled across the air gap.  Here the gap (RI<r<RO) is
NOT meshed; a single H1 space lives on the disconnected rotor annulus (R0<r<RI) and
stator annulus (RO<r<REXT), and the gap's FULL 2x2 DtN couples the two rings:

    [dA/dr|RI ; dA/dr|RO] = M @ [a_RI ; a_RO]     (M = annular_dtn_matrix)

which enters the stiffness as the low-rank block (cr,cs = the cos n*theta boundary
vectors on the rotor / stator rings, norms = pi*R):

    G = [[ -M11/norm_RI, -M12/norm_RO ],
         [  M21/norm_RI,  M22/norm_RO ]]   added as U G U^T,  U = [cr, cs]

Gate: drive the rotor inner ring (R0) and stator outer ring (REXT) with a harmonic;
the whole R0<r<REXT is piecewise-Laplace, so the AGE-coupled FE field in BOTH regions
must equal the 3-region analytic (gap represented exactly, un-meshed).  Single cos
harmonic; real Laplace regions.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx, ds,
                     x, y, atan2, cos, BND, InnerProduct, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_element import annular_dtn_matrix

R0, RI, RO, REXT = 0.05, 0.10, 0.11, 0.20
N = 2
A_ROT, A_STA = 1.0, 0.3


def analytic_3region(n, radii, a_rot, a_sta):
    """6-coeff piecewise-Laplace [rotor | gap | stator] with continuity of A and dA/dr
    at RI and RO; Dirichlet a_rot at R0, a_sta at REXT.  Returns (a1,b1,a2,b2,a3,b3)."""
    R0_, RI_, RO_, REXT_ = radii
    A = lambda r: [r ** n, r ** -n]
    dA = lambda r: [n * r ** (n - 1), -n * r ** (-n - 1)]
    M = np.zeros((6, 6)); rhs = np.zeros(6)
    M[0, 0:2] = A(R0_); rhs[0] = a_rot
    M[1, 0:2] = A(RI_); M[1, 2:4] = [-t for t in A(RI_)]
    M[2, 0:2] = dA(RI_); M[2, 2:4] = [-t for t in dA(RI_)]
    M[3, 2:4] = A(RO_); M[3, 4:6] = [-t for t in A(RO_)]
    M[4, 2:4] = dA(RO_); M[4, 4:6] = [-t for t in dA(RO_)]
    M[5, 4:6] = A(REXT_); rhs[5] = a_sta
    return np.linalg.solve(M, rhs)


def test_two_region_age_matches_3region_analytic():
    geo = SplineGeometry()
    geo.AddCircle((0, 0), REXT, leftdomain=2, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), RO, leftdomain=0, rightdomain=2, bc="stator_ring")
    geo.AddCircle((0, 0), RI, leftdomain=1, rightdomain=0, bc="rotor_ring")
    geo.AddCircle((0, 0), R0, leftdomain=0, rightdomain=1, bc="rotor_inner")
    mesh = Mesh(geo.GenerateMesh(maxh=0.006)); mesh.Curve(6)
    cosn = cos(N * atan2(y, x))
    (M11, M12), (M21, M22) = annular_dtn_matrix(N, RI, RO)
    nri, nro = math.pi * RI, math.pi * RO
    G = np.array([[-M11 / nri, -M12 / nro], [M21 / nri, M22 / nro]])

    fes = H1(mesh, order=6, dirichlet="outer|rotor_inner")
    u, v = fes.TnT()
    a = BilinearForm(fes); a += grad(u) * grad(v) * dx
    with TaskManager():
        a.Assemble()
        Kinv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        gfu = GridFunction(fes)
        gfu.Set(mesh.BoundaryCF({"rotor_inner": A_ROT * cosn, "outer": A_STA * cosn}), BND)
        lr = LinearForm(fes); lr += cosn * v * ds(definedon=mesh.Boundaries("rotor_ring")); lr.Assemble()
        ls = LinearForm(fes); ls += cosn * v * ds(definedon=mesh.Boundaries("stator_ring")); ls.Assemble()
        cr = lr.vec.CreateVector(); cr.data = lr.vec
        cs = ls.vec.CreateVector(); cs.data = ls.vec
        U = [cr, cs]
        rhs = cr.CreateVector(); rhs.data = -(a.mat * gfu.vec)         # no volume source; Dirichlet lift
        # Woodbury: (K + U G U^T) delta = rhs
        Kinv_rhs = rhs.CreateVector(); Kinv_rhs.data = Kinv * rhs
        w = [c.CreateVector() for c in U]
        for i in range(2):
            w[i].data = Kinv * U[i]
        UtKU = np.array([[InnerProduct(U[i], w[j]) for j in range(2)] for i in range(2)])
        small = np.linalg.inv(G) + UtKU
        beta = np.linalg.solve(small, np.array([InnerProduct(U[i], Kinv_rhs) for i in range(2)]))
        delta = Kinv_rhs.CreateVector(); delta.data = Kinv_rhs
        for i in range(2):
            delta.data -= float(beta[i]) * w[i]
        gfu.vec.data += delta

    a1, b1, a2, b2, a3, b3 = analytic_3region(N, (R0, RI, RO, REXT), A_ROT, A_STA)
    print(f"AGE symmetry check: -M12*RI={-M12*RI:.4e} vs M21*RO={M21*RO:.4e} "
          f"(equal => symmetric coupling)")
    diffs, refs = [], []
    for (r_lo, r_hi, ac, bc_) in ((R0, RI, a1, b1), (RO, REXT, a3, b3)):
        for frac in (0.1, 0.5, 0.9):
            rr = r_lo + frac * (r_hi - r_lo)
            for th in (0.3, 1.4, 2.7, -0.9):
                fe = gfu(mesh(rr * math.cos(th), rr * math.sin(th)))
                ex = (ac * rr ** N + bc_ * rr ** -N) * math.cos(N * th)
                diffs.append(abs(fe - ex)); refs.append(abs(ex))
    rel = max(diffs) / max(refs)
    print(f"two-region AGE vs 3-region analytic: max rel err = {rel:.2e}")
    assert rel < 1e-3, f"two-region AGE off analytic by {rel:.2e}"


def main():
    test_two_region_age_matches_3region_analytic()
    print("[OK] AGE step 2: rotor FE region + stator FE region coupled by the gap's full 2x2 "
          "DtN across the UN-MESHED gap, reproducing the 3-region analytic. Next: rotation + torque.")


if __name__ == "__main__":
    main()
