"""Exact eddy-current admittance of a box (slab, square, cube) by heat content.

Scalar model of this lane: (-Lap + s mu sigma) v = -s mu sigma in the box
(0, L)^D, v = 0 on the boundary, Y(s) = Y_DC (1 + <v>).

Let w solve the heat equation w_t = Lap w with w = 0 on the boundary and
w = 1 at tau = 0, and W = int_0^inf exp(-s MS tau) w dtau its Laplace
transform in tau.  Then s MS W - 1 = Lap W, so v = -s MS W and

    <v>(s) = -s MS int_0^inf exp(-s MS tau) m(tau)^D dtau,

because w separates into a product of 1-D solutions whose mean is

    m(tau) = sum_{n odd} 8/(n pi)^2 exp(-(n pi / L)^2 tau)
           = 1 - (4/L) sqrt(tau / pi) + O(erfc(L / (2 sqrt tau))).

For s = j omega the integrand does not decay on the real tau axis; along
the ray tau = rho exp(-j pi/4) both exp(-s MS tau) and m(tau) decay (the
heat trace is analytic for Re tau > 0), so a graded Gauss rule in rho gives
the integral to machine precision.  D = 1 reproduces the slab admittance
(2/(tL)) tanh(tL/2) to 1e-13 (self-test below); D = 3 gives the cube, e.g.
|Y| = 3.431902 at 10 kHz for the 5 mm copper cube, where the Aitken-
accelerated Foster sum of cube3d/05 had 3.431919 and the mixed rank-20
closed-K_ss model 3.443338.  No FEM, no mode truncation.

This replaces the pending NGSolve ground truth as the reference for boxes.
"""
import cmath
import math

import numpy as np

PI = math.pi


def m_heat(tau, L, n_max=201):
    """Mean of the 1-D Dirichlet heat solution with unit initial data at
    complex tau with Re tau > 0 (small |tau|: two-end image formula; else the
    eigen-sum, which converges to 1e-16 with n <= 201 once
    |tau| >= 0.003 L^2 on the pi/4 ray)."""
    tau = np.asarray(tau, dtype=complex)
    out = np.empty_like(tau)
    small = np.abs(tau) < 0.003 * L**2
    out[small] = 1.0 - (4.0 / L) * np.sqrt(tau[small] / PI)
    big = ~small
    if np.any(big):
        n = np.arange(1, n_max + 1, 2, dtype=float)
        lam = (n * PI / L) ** 2
        out[big] = np.sum((8.0 / (n * PI) ** 2)[None, :]
                          * np.exp(-tau[big][:, None] * lam[None, :]), axis=1)
    return out


def ray_rule(rho_min, rho_max, ratio=1.6, nodes=12):
    """Geometric panels on [rho_min, rho_max] with Gauss-Legendre nodes."""
    xg, wg = np.polynomial.legendre.leggauss(nodes)
    edges = [rho_min]
    while edges[-1] < rho_max:
        edges.append(edges[-1] * ratio)
    r, w = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        r.append(0.5 * (b - a) * xg + 0.5 * (b + a))
        w.append(0.5 * (b - a) * wg)
    return np.concatenate(r), np.concatenate(w)


def v_avg_exact(s, D, L, MS):
    """<v>(s) for the D-dimensional box of side L; MS = mu sigma; s = j omega
    (either sign of omega) or 0."""
    if s == 0:
        return 0.0
    sMS = complex(s) * MS
    if sMS.real < 0:
        raise ValueError("s must lie in the closed right half-plane")
    scale = min(1.0 / abs(sMS), L**2 / (3 * PI**2))
    rho, w = ray_rule(1e-13 * scale, 60.0 * max(1.0 / abs(sMS), L**2 / (D * PI**2)))
    phase = cmath.exp(-1j * PI / 4) if sMS.imag > 0 else cmath.exp(1j * PI / 4)
    if abs(sMS.imag) < 1e-12 * abs(sMS):
        phase = 1.0  # real s: the real axis already decays
    tau = rho * phase
    integrand = np.exp(-sMS * tau) * m_heat(tau, L) ** D
    return complex(-sMS * phase * np.sum(w * integrand))


def Y_exact(s, D, L, sigma, mu=4 * PI * 1e-7):
    """Admittance of the box: Y_DC (1 + <v>), Y_DC = sigma L^D."""
    return sigma * L**D * (1.0 + v_avg_exact(s, D, L, mu * sigma))


def v_avg_slab(s, L, MS):
    """Closed form for D = 1: (2/(tL)) tanh(tL/2) - 1, t = sqrt(s MS)."""
    t = cmath.sqrt(complex(s) * MS)
    if abs(t * L) > 200:
        return 2.0 / (t * L) - 1.0
    return (2.0 / (t * L)) * cmath.tanh(t * L / 2) - 1.0


def self_test(L=5e-3, sigma=5.8e7, mu=4 * PI * 1e-7, tol=1e-11):
    """D = 1 against the slab closed form over 1 Hz .. 1 GHz."""
    MS = mu * sigma
    worst = 0.0
    for f in np.logspace(0, 9, 19):
        s = 2j * PI * f
        worst = max(worst, abs(v_avg_exact(s, 1, L, MS) - v_avg_slab(s, L, MS)))
    if worst > tol:
        raise AssertionError(f"heat-content slab check failed: {worst:.2e}")
    return worst


if __name__ == "__main__":
    print(f"slab self-test: max |diff| = {self_test():.1e}")
    for D, name in ((2, "square"), (3, "cube")):
        for f in (1e3, 1e4, 1e5, 1e6):
            y = Y_exact(2j * PI * f, D, 5e-3, 5.8e7)
            print(f"  {name:6s} f = {f:7.0e} Hz: Y = {y.real:.6f} {y.imag:+.6f} j  |Y| = {abs(y):.6f}")
