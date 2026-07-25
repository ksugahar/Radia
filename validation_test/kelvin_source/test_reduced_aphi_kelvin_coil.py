"""Golden: reduced A-Phi with a REAL coil source + Periodic Kelvin at p = 2,
and the curve-order pitfall.

Extends test_aphi_kelvin_eddy.py (uniform background) to the production
configuration: the source is a coaxial circular filament loop (b = 0.6 m,
z0 = 0.5 m, I = 1 kA) inside the inner region -- a genuine decaying
Biot-Savart source, the coil + workpiece pattern.  The solvers are IMPORTED
from the uniform golden (parameterized by A_s), so the two goldens share one
implementation.

Source representation: the reduced weak form needs A_s ONLY inside the
conductor (curl(nu0 curl A_s) = J_coil cancels the coil current exactly, and
the reaction field decays, so nothing is carried into the Kelvin exterior).
A_s is built as the loop's interior multipole expansion -- a polynomial
VectorCF -- and is verified pointwise against the exact elliptic-integral
loop formula to ~1e-10 before use.

Analytic reference: multipole orthogonality means only the n = 1 source
component drives the induced DIPOLE, with the degree-1 reflection
coefficient of the conducting sphere:

    m_z = 2 pi a^3 (B_c / mu0) Gamma_1(x),   B_c = mu0 I b^2 / (2 c0^3),
    Gamma_1 = 3 / (x I_{1/2}(x)/I_{3/2}(x) ... ) - 1
            = -[ 1 - (3/x) coth x + 3/x^2 ]          (verified to 1e-15)

MEASURED (LAB, 2026-07-25), rel = |m_fem - m_ana| / |m_ana|, a/delta = 2:

    lane                                    curve=2      curve=1
    p=1  A* = A-Phi (lanes identical)       2.899%       --
    p=2  A* (plain A-method)                0.475%       2.200%
    p=2  A-Phi                              0.054%       1.870%

Conclusions locked here:

1. Reduced A-Phi with a real coil source + Kelvin at p=2: 0.054% -- the
   uniform-background result (0.053%) carries over unchanged.  The A-method
   floor (0.475% vs 0.473%) is source-independent too.
2. CURVE ORDER MUST MATCH p.  At p=2 on a curve=1 (straight-tet) mesh the
   A-Phi error explodes 35x (0.054% -> 1.870%) and the A*/A-Phi distinction
   is buried under geometry error.  "p=2" always means mesh.Curve(2).

Reference chain: test_aphi_kelvin_eddy.py (formulation, uniform source),
scratchpad loop_sphere_analytic self-checks (Gamma_1 == uniform formula to
1e-15; multipole reconstruction vs elliptic integrals to 1e-10; c_1 = -B_c/2
to 5.8e-15).
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.polynomial import legendre as L
from scipy.special import ellipe, ellipk, iv, lpmv

_HERE = Path(__file__).resolve().parent


def _load_uniform_golden():
    spec = importlib.util.spec_from_file_location(
        "aphi_kelvin_eddy_mod", _HERE / "test_aphi_kelvin_eddy.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


aphi = _load_uniform_golden()

from ngsolve import CoefficientFunction as CF, x, y, z  # noqa: E402

MU_0 = aphi.MU_0
A_C = aphi.A_C
SIGMA = aphi.SIGMA

B_LOOP = 0.6      # loop radius [m]
Z_LOOP = 0.5      # loop height [m]
I_LOOP = 1.0e3    # [A]
C0 = math.hypot(B_LOOP, Z_LOOP)
NMAX = 30

_CACHE: dict = {}


def _a_phi_exact(rho, zz):
    m = 4 * B_LOOP * rho / ((B_LOOP + rho) ** 2 + (zz - Z_LOOP) ** 2)
    return (MU_0 * I_LOOP / math.pi * math.sqrt(B_LOOP / rho) / math.sqrt(m)
            * ((1 - m / 2) * ellipk(m) - ellipe(m)))


def _loop_coeffs():
    """Interior multipole coefficients c_n of the loop A_phi."""
    r0 = 0.35
    th = np.linspace(1e-6, math.pi - 1e-6, 4001)
    u = np.cos(th)
    aphi_vals = np.array([_a_phi_exact(r0 * math.sin(t), r0 * math.cos(t))
                          for t in th])
    c = np.zeros(NMAX + 1)
    for n in range(1, NMAX + 1):
        norm = 2 * n * (n + 1) / (2 * n + 1)
        c[n] = np.trapezoid(aphi_vals * lpmv(1, n, u) * np.sin(th),
                            th) / norm / r0**n
    return c


def _loop_A_s_cf(c):
    """A_s = (A_phi/rho)(-y, x, 0); A_phi/rho = sum_n (-c_n) r^{n-1} P_n'(z/r),
    a polynomial in (z, r^2) -- exact CF, no filament discretization."""
    r2 = x * x + y * y + z * z
    total = None
    for n in range(1, NMAX + 1):
        if abs(c[n]) * A_C**n < 1e-16:
            continue
        dP = L.leg2poly(L.legder([0] * n + [1]))
        term = None
        for k, gk in enumerate(dP):
            if abs(gk) < 1e-14:
                continue
            pw = n - 1 - k
            assert pw % 2 == 0
            mono = gk * (z ** k if k else 1.0)
            for _ in range(pw // 2):
                mono = mono * r2
            term = mono if term is None else term + mono
        contrib = (-c[n]) * term
        total = contrib if total is None else total + contrib
    return CF((-y * total, x * total, 0.0))


def _analytic_moment(omega):
    delta = math.sqrt(2.0 / (omega * MU_0 * SIGMA))
    xk = (1 + 1j) * A_C / delta
    gamma1 = 3.0 / (xk * iv(0.5, xk) / iv(1.5, xk)) - 1.0
    b_c = MU_0 * I_LOOP * B_LOOP**2 / (2 * C0**3)
    return 2 * math.pi * A_C**3 * (b_c / MU_0) * gamma1


def _results():
    if not _CACHE:
        omega, _ = aphi._analytic()          # same a/delta = 2 frequency
        m_ana = _analytic_moment(omega)
        c = _loop_coeffs()
        a_s = _loop_A_s_cf(c)

        mesh2 = aphi._build_mesh(curve=2)
        # pointwise contract: the CF reproduces the exact loop A_phi
        max_rel = 0.0
        for (rr, tt) in ((0.25, 1.2), (0.39, 2.2), (0.2, 0.5)):
            rho, zz = rr * math.sin(tt), rr * math.cos(tt)
            exact = _a_phi_exact(rho, zz)
            got = a_s[1](mesh2(rho, 0.0, zz))    # A_phi at phi=0 is +y comp
            max_rel = max(max_rel, abs(got / exact - 1))
        _CACHE["As_pointwise_rel"] = max_rel

        _CACHE["m_ana"] = m_ana
        _CACHE["astar_p1"] = aphi._solve_astar(mesh2, 1, omega, A_s=a_s)
        _CACHE["aphi_p1"] = aphi._solve_aphi(mesh2, 1, omega, A_s=a_s)
        _CACHE["astar_p2"] = aphi._solve_astar(mesh2, 2, omega, A_s=a_s)
        _CACHE["aphi_p2"] = aphi._solve_aphi(mesh2, 2, omega, A_s=a_s)

        mesh1 = aphi._build_mesh(curve=1)
        _CACHE["aphi_p2_curve1"] = aphi._solve_aphi(mesh1, 2, omega, A_s=a_s)
    return _CACHE


def _rel(m_fem, m_ana):
    return abs(m_fem - m_ana) / abs(m_ana)


@pytest.mark.slow
def test_loop_A_s_cf_matches_exact_elliptic():
    """The multipole CF reproduces the exact loop A_phi in the conductor."""
    res = _results()
    print(f"\n  A_s pointwise max rel = {res['As_pointwise_rel']:.2e}")
    assert res["As_pointwise_rel"] < 1e-6


@pytest.mark.slow
def test_reduced_aphi_p2_matches_analytic_with_coil_source():
    """THE question: reduced A-Phi + Kelvin + real coil source at p=2."""
    res = _results()
    rel = _rel(res["aphi_p2"], res["m_ana"])
    print(f"\n  A-Phi p=2 (coil): m = {res['aphi_p2']:.6e}")
    print(f"  analytic        : m = {res['m_ana']:.6e}")
    print(f"  rel = {rel*100:.3f}%  (measured 0.054%, band 0.2%)")
    assert rel < 0.002, f"reduced A-Phi p=2 off by {rel*100:.3f}%"


@pytest.mark.slow
def test_astar_floor_is_source_independent():
    """The A-method p=2 floor (~0.47%) matches the uniform-source golden."""
    res = _results()
    rel = _rel(res["astar_p2"], res["m_ana"])
    print(f"\n  A* p=2 (coil) rel = {rel*100:.3f}%  "
          f"(uniform golden: 0.473%, band 1.5%)")
    assert rel < 0.015
    rel_aphi = _rel(res["aphi_p2"], res["m_ana"])
    assert rel_aphi * 3.0 < rel, "A-Phi must beat A* by >= 3x with coil source"


@pytest.mark.slow
def test_p1_lanes_agree_with_coil_source():
    """W stays inert at p=1 also for the coil source."""
    res = _results()
    gap = abs(res["astar_p1"] - res["aphi_p1"]) / abs(res["m_ana"])
    rel = _rel(res["aphi_p1"], res["m_ana"])
    print(f"\n  p=1: rel = {rel*100:.3f}%, lane gap = {gap:.2e}")
    assert rel < 0.05 and gap < 1e-4


@pytest.mark.slow
def test_curve_order_must_match_p():
    """p=2 on a curve=1 mesh loses the A-Phi accuracy by an order of magnitude.

    Measured: 0.054% (curve=2) -> 1.870% (curve=1), a 35x degradation --
    geometry error buries the formulation.  'p=2' always means mesh.Curve(2).
    """
    res = _results()
    rel_c2 = _rel(res["aphi_p2"], res["m_ana"])
    rel_c1 = _rel(res["aphi_p2_curve1"], res["m_ana"])
    print(f"\n  A-Phi p=2: curve=2 {rel_c2*100:.3f}% vs curve=1 "
          f"{rel_c1*100:.3f}%  (factor {rel_c1/max(rel_c2,1e-30):.0f})")
    assert rel_c1 > 0.01, (
        f"expected the curve=1 mesh to degrade p=2 accuracy past 1%; "
        f"got {rel_c1*100:.3f}%")
    assert rel_c2 * 10.0 < rel_c1, \
        "curve=2 must beat curve=1 by >= 10x at p=2"


if __name__ == "__main__":
    res = _results()
    m_ana = res["m_ana"]
    print(f"m_analytic = {m_ana:.6e}")
    print(f"A_s pointwise: {res['As_pointwise_rel']:.2e}")
    for key in ("astar_p1", "aphi_p1", "astar_p2", "aphi_p2",
                "aphi_p2_curve1"):
        print(f"  {key:16s} rel = {_rel(res[key], m_ana)*100:.3f}%")
