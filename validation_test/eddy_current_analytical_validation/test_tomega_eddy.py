"""
T-Omega eddy-current formulation: the SOURCE/REACTION split, gated by the Bessel field.

The T-Omega method writes H = T - grad(Omega): the magnetic scalar potential Omega
carries the curl-free (applied / current-source) part, while the electric vector
potential T carries the eddy reaction (J = curl T) inside the conductor. This is the
multiply-connected eddy core used by COMSOL/JMAG (and the partner of the reduced scalar
potential + cohomology cuts).

We forge that SPLIT on a conducting cylinder (radius a, sigma, mu) in a uniform axial AC
field H0. The applied field is a pure gradient, so it goes entirely into Omega = -H0 z
(=> -grad(Omega) = H0 z_hat everywhere). The reaction T = T_z(r) z_hat then satisfies a
diffusion equation with a HOMOGENEOUS Dirichlet trace T_z(a)=0 (no reaction outside) and a
source driven by the applied field:

    (1/r) d/dr( r dT_z/dr ) = k^2 (T_z + H0),   k=(1+j)/delta,   T_z(a)=0,

i.e. weak form  int(T' v') r dr + k^2 int(T v) r dr = -k^2 H0 int(v) r dr  (rho cancels).
The physical field is then RECONSTRUCTED as H_z = T_z + H0, and must equal the exact
Kelvin/Bessel profile H_z(r) = H0 I0(kr)/I0(ka). So the test proves the T-Omega
decomposition (source in Omega, reaction in T, T=0 on the boundary) reproduces the true
eddy field -- in a 1D radial solve and a genuine 2D unstructured disk solve.

Pure NGSolve + scipy (Bash-robust).
"""
import math
import pytest

ng = pytest.importorskip("ngsolve")
sp = pytest.importorskip("scipy.special")
from scipy.special import iv
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, Integrate,
                     grad, dx, x, CF, BND)
from ngsolve.meshes import Make1DMesh

MU0 = 4e-7 * math.pi
SIGMA = 3.5e7
H0 = 1.0


@pytest.fixture(autouse=True)
def _taskmanager():
    with ng.TaskManager():
        yield


def _k(a, a_over_delta, mu_r=1.0):
    delta = a / a_over_delta
    mu = mu_r * MU0
    omega = 2.0 / (mu * SIGMA * delta ** 2)
    return (1 + 1j) / delta


def _bessel_Hz(r, a, k):
    return H0 * iv(0, k * r) / iv(0, k * a)


# --- 1D radial T-Omega: solve the REACTION T_z, reconstruct H_z = T_z + H0 ---------------
def _solve_T_1d(a, a_over_delta, order=4, n=500):
    k = _k(a, a_over_delta)
    mesh = Make1DMesh(n, mapping=lambda t: t * a)
    fes = H1(mesh, order=order, complex=True, dirichlet="right")  # T_z(a)=0
    u, v = fes.TnT()
    a_f = BilinearForm(grad(u) * grad(v) * x * dx + (k * k) * u * v * x * dx).Assemble()
    f = LinearForm(-(k * k) * H0 * v * x * dx).Assemble()
    T = GridFunction(fes)                 # homogeneous Dirichlet (T_z(a)=0)
    T.vec.data = a_f.mat.Inverse(fes.FreeDofs()) * f.vec
    return mesh, T, k


@pytest.mark.parametrize("a_over_delta", [0.5, 2.0, 5.0])
def test_tomega_1d_recovers_bessel(a_over_delta):
    """Reconstructed H_z = T_z + H0 equals the exact Bessel field at sample radii."""
    a = 0.01
    mesh, T, k = _solve_T_1d(a, a_over_delta)
    for frac in (0.0, 0.25, 0.5, 0.75):
        r = frac * a
        Hz = complex(T(mesh(r))) + H0
        an = _bessel_Hz(r, a, k)
        assert abs(Hz - an) / abs(an) < 1e-5, \
            "a/d=%.1f r/a=%.2f: H_z=%s Bessel=%s" % (a_over_delta, frac, Hz, an)


def test_tomega_split_structure():
    """The split is correct: T_z=0 on the surface (reaction confined), and the centre
    reaction T_z(0) = H0(1/I0(ka) - 1) (= shielded field minus the applied field)."""
    a, a_over_delta = 0.01, 3.0
    mesh, T, k = _solve_T_1d(a, a_over_delta)
    Ta = complex(T(mesh(a)))
    assert abs(Ta) < 1e-6 * H0, "reaction not zero on the surface: T_z(a)=%s" % Ta
    T0 = complex(T(mesh(0.0)))
    T0_an = H0 * (1.0 / iv(0, k * a) - 1.0)
    assert abs(T0 - T0_an) / abs(T0_an) < 1e-4, "centre reaction %s vs %s" % (T0, T0_an)
    # reconstructed centre field is the shielding factor
    assert abs((complex(T(mesh(0.0))) + H0) - H0 / iv(0, k * a)) < 1e-4


# --- 2D unstructured disk T-Omega: same split on a real mesh -----------------------------
def _solve_T_2d(a, a_over_delta, maxh_frac=30.0):
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh
    k = _k(a, a_over_delta)
    geo = SplineGeometry()
    geo.AddCircle((0, 0), a, bc="outer")
    mesh = Mesh(geo.GenerateMesh(maxh=a / maxh_frac))
    fes = H1(mesh, order=4, complex=True, dirichlet="outer")   # T=0 on r=a
    u, v = fes.TnT()
    a_f = BilinearForm(grad(u) * grad(v) * dx + (k * k) * u * v * dx).Assemble()
    f = LinearForm(-(k * k) * H0 * v * dx).Assemble()
    T = GridFunction(fes)
    T.vec.data = a_f.mat.Inverse(fes.FreeDofs()) * f.vec
    return mesh, T, k


def test_tomega_2d_disk_recovers_bessel():
    """The 2D T-Omega split on an unstructured disk reconstructs the Bessel field."""
    a, a_over_delta = 0.01, 3.0
    mesh, T, k = _solve_T_2d(a, a_over_delta)
    Hz0 = complex(T(mesh(0.0, 0.0))) + H0
    an0 = _bessel_Hz(0.0, a, k)
    assert abs(Hz0 - an0) / abs(an0) < 1e-2, "2D centre H_z %s vs Bessel %s" % (Hz0, an0)
    rs = 0.5 * a
    Hzs = complex(T(mesh(rs, 0.0))) + H0
    ans = _bessel_Hz(rs, a, k)
    assert abs(Hzs - ans) / abs(ans) < 1e-2, "2D r/a=0.5 H_z %s vs Bessel %s" % (Hzs, ans)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
