r"""AGE full-machine, steps (3) rotation + (4) torque extraction.

Builds on the two-region AGE coupling (step 2): rotor FE region + stator FE region across
the un-meshed gap.  The rotor excitation is ROTATED by theta_r (no remesh -- a pure phase
of the rotor harmonic, the AGE's reason for being), which requires BOTH cos and sin gap
modes.  Torque is then read straight from the FE gap-ring harmonics via the closed-form
AGE harmonic torque.

  (3) rotation : rotor inner-ring Dirichlet = A_rot cos(n(theta - theta_r)); sweep theta_r,
                 re-solve with the geometry/mesh FIXED.
  (4) torque   : extract the gap-ring phasors (a_cos - i a_sin) from the FE solution,
                 T = airgap_harmonic_torque(n, RI, RO, A_RI, A_RO).

Gate: T_FE(theta_r) equals the analytic torque from the complex 3-region solve at every
angle, and it genuinely VARIES with theta_r (the rotating-field torque, ~ sin(n theta_r)).
"""
import cmath
import math
import os
import sys

import numpy as np
from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx, ds,
                     x, y, atan2, cos, sin, BND, InnerProduct, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_element import (annular_dtn_matrix,
                                                    airgap_harmonic_torque)

R0, RI, RO, REXT = 0.05, 0.10, 0.11, 0.20
N = 2
A_ROT, A_STA, LSTK = 1.0, 0.3, 0.05
NRI, NRO = math.pi * RI, math.pi * RO


def analytic_3region(n, radii, a_rot, a_sta):
    R0_, RI_, RO_, REXT_ = radii
    A = lambda r: [r ** n, r ** -n]
    dA = lambda r: [n * r ** (n - 1), -n * r ** (-n - 1)]
    M = np.zeros((6, 6), complex); rhs = np.zeros(6, complex)
    M[0, 0:2] = A(R0_); rhs[0] = a_rot
    M[1, 0:2] = A(RI_); M[1, 2:4] = [-t for t in A(RI_)]
    M[2, 0:2] = dA(RI_); M[2, 2:4] = [-t for t in dA(RI_)]
    M[3, 2:4] = A(RO_); M[3, 4:6] = [-t for t in A(RO_)]
    M[4, 2:4] = dA(RO_); M[4, 4:6] = [-t for t in dA(RO_)]
    M[5, 4:6] = A(REXT_); rhs[5] = a_sta
    return np.linalg.solve(M, rhs)


def analytic_torque(thr):
    c = analytic_3region(N, (R0, RI, RO, REXT), A_ROT * cmath.exp(-1j * N * thr), A_STA + 0j)
    a2, b2 = c[2], c[3]
    return airgap_harmonic_torque(N, RI, RO, a2 * RI ** N + b2 * RI ** -N,
                                  a2 * RO ** N + b2 * RO ** -N, axial_length=LSTK)


def test_rotation_and_torque():
    geo = SplineGeometry()
    geo.AddCircle((0, 0), REXT, leftdomain=2, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), RO, leftdomain=0, rightdomain=2, bc="stator_ring")
    geo.AddCircle((0, 0), RI, leftdomain=1, rightdomain=0, bc="rotor_ring")
    geo.AddCircle((0, 0), R0, leftdomain=0, rightdomain=1, bc="rotor_inner")
    mesh = Mesh(geo.GenerateMesh(maxh=0.006)); mesh.Curve(6)
    cosn, sinn = cos(N * atan2(y, x)), sin(N * atan2(y, x))
    (M11, M12), (M21, M22) = annular_dtn_matrix(N, RI, RO)
    G = np.array([[-M11 / NRI, -M12 / NRO], [M21 / NRI, M22 / NRO]])
    Gfull = np.zeros((4, 4)); Gfull[0:2, 0:2] = G; Gfull[2:4, 2:4] = G

    fes = H1(mesh, order=6, dirichlet="outer|rotor_inner")
    u, v = fes.TnT()
    a = BilinearForm(fes); a += grad(u) * grad(v) * dx

    def bvec(cf, bc):
        L = LinearForm(fes); L += cf * v * ds(definedon=mesh.Boundaries(bc)); L.Assemble()
        c = L.vec.CreateVector(); c.data = L.vec; return c

    with TaskManager():
        a.Assemble()
        Kinv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        U = [bvec(cosn, "rotor_ring"), bvec(cosn, "stator_ring"),
             bvec(sinn, "rotor_ring"), bvec(sinn, "stator_ring")]   # crc, csc, crs, css
        w = [c.CreateVector() for c in U]
        for i in range(4):
            w[i].data = Kinv * U[i]
        UtKU = np.array([[InnerProduct(U[i], w[j]) for j in range(4)] for i in range(4)])
        smallinv = np.linalg.inv(np.linalg.inv(Gfull) + UtKU)

        def solve(thr):
            gfu = GridFunction(fes)
            gfu.Set(mesh.BoundaryCF({"rotor_inner": A_ROT * cos(N * atan2(y, x) - N * thr),
                                     "outer": A_STA * cosn}), BND)
            rhs = gfu.vec.CreateVector(); rhs.data = -(a.mat * gfu.vec)
            Kr = rhs.CreateVector(); Kr.data = Kinv * rhs
            beta = smallinv @ np.array([InnerProduct(U[i], Kr) for i in range(4)])
            d = Kr.CreateVector(); d.data = Kr
            for i in range(4):
                d.data -= float(beta[i]) * w[i]
            gfu.vec.data += d
            return gfu

        def torque(gfu):
            aRI = InnerProduct(U[0], gfu.vec) / NRI - 1j * InnerProduct(U[2], gfu.vec) / NRI
            aRO = InnerProduct(U[1], gfu.vec) / NRO - 1j * InnerProduct(U[3], gfu.vec) / NRO
            return airgap_harmonic_torque(N, RI, RO, aRI, aRO, axial_length=LSTK)

        thetas = [k * (2 * math.pi / N) / 8 for k in range(8)]      # one torque period
        feT = np.array([torque(solve(thr)) for thr in thetas])

    anT = np.array([analytic_torque(thr) for thr in thetas])
    rel = np.max(np.abs(feT - anT)) / np.max(np.abs(anT))
    print("theta_r [deg]   FE torque [N*m]   analytic [N*m]")
    for thr, f, an in zip(thetas, feT, anT):
        print(f"  {math.degrees(thr):6.1f}      {f:+.5e}     {an:+.5e}")
    assert np.max(np.abs(anT)) > 1e-10, "torque trivially zero"
    assert rel < 2e-3, f"FE torque off analytic by {rel:.2e}"
    assert (feT.max() - feT.min()) / np.max(np.abs(feT)) > 1.0, "torque should swing +/- with rotation"


def main():
    test_rotation_and_torque()
    print("[OK] AGE steps 3+4: the rotor is rotated by a pure PHASE (no remesh) and the torque, "
          "read from the FE gap harmonics, tracks the analytic rotating-field torque T(theta_r). "
          "Two FE regions + mesh-free gap + rotation + torque = the AGE rotating-machine core.")


if __name__ == "__main__":
    main()
