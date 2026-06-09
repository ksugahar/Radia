r"""AGE full-machine, step (1) -- the MULTI-HARMONIC dense gap block.

The single-harmonic AGE (test_airgap_ngsolve_coupling) used a scalar Robin term; a real
machine air-gap carries a full spectrum (slot / pole harmonics), each needing its own DtN
coefficient M22(n).  This wires ALL retained harmonics into ONE FE solve via the spectral
boundary block

    sum_n (M22(n)/norm) c_n c_n^T          (c_n = integral of cos n*theta * v over the inner ring)

added to the FE stiffness as a low-rank correction (solved by Woodbury), with rotor-ring
forcing  sum_n -M21(n) A_in,n c_n.

Gate: drive the (un-meshed) gap with a MULTI-HARMONIC rotor-ring phasor and a
multi-harmonic outer Dirichlet; the AGE-coupled FE field over ro<r<Rext must equal the
SUPERPOSED full-annulus analytic (each harmonic an independent alpha r^n + beta r^-n).
A single scalar M22 would mis-couple the off-harmonics -- this test fails unless each
mode gets its own M22(n).  Real fields, cos modes; single complex-free Laplace FE region.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx, ds,
                     x, y, atan2, cos, InnerProduct, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_element import (airgap_dtn_modes,
                                                    annular_harmonic_coeffs)

RI, RO, REXT = 0.10, 0.11, 0.20
HARM = [1, 2, 3, 5]                                   # retained gap harmonics (cos modes)
A_IN = {1: 0.5, 2: 1.0, 3: 0.3, 5: 0.4}              # rotor-ring cos coefficients
A_EXT = {1: 0.1, 2: 0.3, 3: 0.05, 5: 0.2}            # outer Dirichlet cos coefficients
NORM = math.pi * RO                                  # integral cos^2(n theta) ds over r=ro (n>=1)


def _cf_sum(coeffs):
    cf = coeffs[HARM[0]] * cos(HARM[0] * atan2(y, x))
    for n in HARM[1:]:
        cf = cf + coeffs[n] * cos(n * atan2(y, x))
    return cf


def solve_age_multiharmonic():
    geo = SplineGeometry()
    geo.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), RO, leftdomain=0, rightdomain=1, bc="inner")
    mesh = Mesh(geo.GenerateMesh(maxh=0.006)); mesh.Curve(6)
    fes = H1(mesh, order=6, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes); a += grad(u) * grad(v) * dx
    ds_in = ds(definedon=mesh.Boundaries("inner"))
    modes = airgap_dtn_modes(RI, RO, HARM)
    with TaskManager():
        a.Assemble()
        Kinv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        gfu = GridFunction(fes)
        gfu.Set(_cf_sum(A_EXT), definedon=mesh.Boundaries("outer"))
        # per-mode boundary vectors c_n, projector weights d_n, rotor-ring source
        cvecs, dvals = [], []
        src = gfu.vec.CreateVector(); src[:] = 0.0
        for n in HARM:
            ln = LinearForm(fes); ln += cos(n * atan2(y, x)) * v * ds_in; ln.Assemble()
            c = ln.vec.CreateVector(); c.data = ln.vec; cvecs.append(c)
            m21, m22 = modes[n]
            dvals.append(m22 / NORM)
            src.data += (-m21 * A_IN[n]) * ln.vec               # source += -M21(n) A_in,n c_n
        rhs = src.CreateVector(); rhs.data = src - a.mat * gfu.vec
        # Woodbury for (K + sum d_i c_i c_i^T) delta = rhs  (low-rank dense correction)
        m = len(HARM)
        Kinv_rhs = rhs.CreateVector(); Kinv_rhs.data = Kinv * rhs
        w = [c.CreateVector() for c in cvecs]
        for i in range(m):
            w[i].data = Kinv * cvecs[i]
        M = np.array([[InnerProduct(cvecs[i], w[j]) for j in range(m)] for i in range(m)])
        for i in range(m):
            M[i, i] += 1.0 / dvals[i]
        beta = np.linalg.solve(M, np.array([InnerProduct(cvecs[i], Kinv_rhs) for i in range(m)]))
        delta = Kinv_rhs.CreateVector(); delta.data = Kinv_rhs
        for i in range(m):
            delta.data -= float(beta[i]) * w[i]
        gfu.vec.data += delta
    return mesh, gfu


def test_multiharmonic_age_matches_superposed_analytic():
    mesh, gfu = solve_age_multiharmonic()
    coef = {n: annular_harmonic_coeffs(n, RI, REXT, A_IN[n], A_EXT[n]) for n in HARM}

    def exact(r, th):
        return sum((coef[n][0] * r ** n + coef[n][1] * r ** -n) * math.cos(n * th) for n in HARM)

    diffs, refs = [], []
    for frac in (0.05, 0.4, 0.9):
        rr = RO + frac * (REXT - RO)
        for th in (0.2, 1.0, 2.3, -1.4, 3.0):
            diffs.append(abs(gfu(mesh(rr * math.cos(th), rr * math.sin(th))) - exact(rr, th)))
            refs.append(abs(exact(rr, th)))
    rel = max(diffs) / max(refs)
    print(f"multi-harmonic AGE ({len(HARM)} modes) vs superposed analytic: max rel err = {rel:.2e}")
    assert rel < 1e-4, f"multi-harmonic AGE block off analytic by {rel:.2e}"


def main():
    test_multiharmonic_age_matches_superposed_analytic()
    print("[OK] AGE full-machine step (1): the multi-harmonic spectral gap block carries the "
          "full air-gap spectrum across the un-meshed gap in ONE solve (matches superposed "
          "analytic). Next: two FE regions, rotation, FE torque extraction.")


if __name__ == "__main__":
    main()
