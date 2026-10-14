"""
Eddy-current magnetic-diffusion weak form vs the analytic SKIN EFFECT.

The time-harmonic eddy-current operator (the core of the A-V and T-Omega
formulations) is the complex magnetic-diffusion equation
    -H'' + j omega mu sigma H = 0.
For a conductor filling x >= 0 with a tangential surface field H0 the exact
solution is the classic skin effect
    H(x) = H0 exp(-(1+j) x / delta),   delta = sqrt(2 / (omega mu sigma)),
i.e. |H| decays as exp(-x/delta) and the phase lags by x/delta radians.

This test solves the complex weak form on a 1-D conductor (depth L = 6 delta,
Neumann at the far end approximating a half-space) and checks the FEM field
against the analytic skin-effect profile -- forging the eddy-current 数式 with
an analytic gate, pure NGSolve (no MCP).
"""
import math
import cmath
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, BilinearForm, GridFunction, grad, dx, BND)
from ngsolve.meshes import Make1DMesh


def _grad0(gfu, mesh):
    g = grad(gfu)(mesh(1e-9))
    return complex(g[0]) if hasattr(g, "__len__") else complex(g)

MU0 = 4e-7 * math.pi
SIGMA = 5.8e7   # copper
F = 1.0e3
OMEGA = 2 * math.pi * F
DELTA = math.sqrt(2.0 / (OMEGA * MU0 * SIGMA))
H0 = 1.0


@pytest.fixture(autouse=True)
def _taskmanager():
    with ng.TaskManager():
        yield


def _solve_skin(order=4, n=800, n_delta=6):
    L = n_delta * DELTA
    mesh = Make1DMesh(n, mapping=lambda t: t * L)
    fes = H1(mesh, order=order, complex=True, dirichlet="left")
    u, v = fes.TnT()
    k2 = 1j * OMEGA * MU0 * SIGMA
    a = BilinearForm(grad(u) * grad(v) * dx + k2 * u * v * dx).Assemble()
    gfu = GridFunction(fes)
    gfu.Set(H0, BND)                      # H0 on the surface (Dirichlet left)
    r = gfu.vec.CreateVector()
    r.data = -a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    return mesh, gfu


def test_skin_depth_decay():
    """|H| decays as exp(-x/delta): at one skin depth |H| = 1/e."""
    mesh, gfu = _solve_skin()
    val = complex(gfu(mesh(DELTA)))
    assert abs(abs(val) - math.exp(-1.0)) < 1e-3, "|H(delta)| = %.5f != 1/e" % abs(val)


def test_skin_phase_lag():
    """Phase lags by x/delta radians: at one skin depth phase = -1 rad."""
    mesh, gfu = _solve_skin()
    ph = cmath.phase(complex(gfu(mesh(DELTA))))
    assert abs(ph - (-1.0)) < 1e-3, "phase(delta) = %.4f != -1 rad" % ph


def test_skin_profile_matches_analytic():
    """Full complex profile matches H0 exp(-(1+j)x/delta) over [0.2, 2.5] delta."""
    mesh, gfu = _solve_skin()
    kc = (1 + 1j) / DELTA
    err = 0.0
    xd = 0.2
    while xd <= 2.5:
        x = xd * DELTA
        err = max(err, abs(complex(gfu(mesh(x))) - H0 * cmath.exp(-kc * x)))
        xd += 0.1
    assert err < 5e-4, "max|H_FEM - H_analytic| = %.3e over [0.2,2.5]delta" % err


def test_surface_impedance():
    """Surface impedance Z_s = E_z(0)/H(0) = -(1/sigma) H'(0)/H(0) equals the
    analytic (1+j)/(sigma*delta) -- the foundation of the SIBC (surface-impedance
    boundary condition). Re(Z_s) = 1/(sigma*delta) gives the eddy power loss/area."""
    mesh, gfu = _solve_skin()
    Zs = -(1.0 / SIGMA) * _grad0(gfu, mesh) / complex(gfu(mesh(0.0)))
    Zs_an = (1 + 1j) / (SIGMA * DELTA)
    assert abs(Zs - Zs_an) / abs(Zs_an) < 1e-3, \
        "Z_s = %.3e+%.3ej != (1+j)/(sigma delta) = %.3e+%.3ej" % (
            Zs.real, Zs.imag, Zs_an.real, Zs_an.imag)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
