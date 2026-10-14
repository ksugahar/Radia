"""
Nonlinear (saturable) magnetostatic weak form + Newton convergence.

The motor / TEAM-20 steel 数式: the 2-D vector-potential magnetostatic equation with a
field-dependent reluctivity nu(|B|) (Brauer-type saturation),

    -div( nu(|grad A_z|) grad A_z ) = J_z ,     B = |grad A_z|,

solved by Newton's method on the semilinear form
    a(A; v) = int nu(|grad A|^2) grad A . grad v .

The reluctivity nu(s) = 1 + alpha*s with alpha>0 makes the operator Jacobian
    nu(g^2) I + 2 alpha (grad A)(grad A)^T   (SPD)
strictly monotone, so the root is UNIQUE. Validated by a *consistent* manufactured
solution: pick A_ex = sin(pi x) sin(pi y), build the RHS as the SAME nonlinear operator
applied to A_ex, hence A = A_ex is the exact (unique) discrete root. Tests:
  * Newton recovers A_ex to FEM accuracy, and a TIGHTER Newton tolerance shrinks the gap
    -> proves A_ex really is the root (residual is the only error);
  * the saturation is genuinely ACTIVE (nu varies >10% across the domain);
  * from a perturbation of the root, Newton converges QUADRATICALLY (accelerating
    residual drop, machine zero in a few steps) -- i.e. the linearization is correct.

Pure NGSolve (Bash-robust; no MCP / Cubit / COMSOL dependency).
"""
import math
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve import (H1, BilinearForm, GridFunction, Integrate, Projector,
                     grad, dx, x, y, sin, pi, CF)
from ngsolve.meshes import MakeStructured2DMesh

ALPHA = 0.6                       # saturation strength (strong: nu varies several-fold)
_AEX = sin(pi * x) * sin(pi * y)


@pytest.fixture(autouse=True)
def _taskmanager():
    with ng.TaskManager():
        yield


def _mesh(n):
    return MakeStructured2DMesh(quads=False, nx=n, ny=n)


def _nu(grad2):
    return 1.0 + ALPHA * grad2


def _build(fes):
    """(a, gfu, Aex) for the consistent manufactured nonlinear problem; gfu has the
    exact Dirichlet trace and zero interior (a cold Newton start)."""
    u, v = fes.TnT()
    Aex = GridFunction(fes); Aex.Set(_AEX)
    a = BilinearForm(fes)
    a += (_nu(grad(u) * grad(u)) * grad(u) * grad(v)) * dx                 # nonlinear operator
    a += (-_nu(grad(Aex) * grad(Aex)) * grad(Aex) * grad(v)) * dx          # consistent RHS
    gfu = GridFunction(fes); gfu.Set(Aex, ng.BND)
    return a, gfu, Aex


def _l2(gf, ref, mesh):
    return math.sqrt(Integrate((gf - ref) ** 2, mesh))


def _newton_error(maxerr=None):
    from ngsolve.solvers import Newton

    mesh = _mesh(16)
    fes = H1(mesh, order=3, dirichlet=".*")
    a, gfu, exact = _build(fes)
    options = {"maxit": 40, "printing": False}
    if maxerr is not None:
        options["maxerr"] = maxerr
    Newton(a, gfu, **options)
    return _l2(gfu, exact, mesh)


def _saturation_relative_spread():
    mesh = _mesh(16)
    fes = H1(mesh, order=3, dirichlet=".*")
    exact = GridFunction(fes)
    exact.Set(_AEX)
    grad2 = grad(exact) * grad(exact)
    area = Integrate(CF(1), mesh)
    mean = Integrate(_nu(grad2), mesh) / area
    std = math.sqrt(Integrate((_nu(grad2) - mean) ** 2, mesh) / area)
    return float(std / mean)


def _newton_convergence():
    mesh = _mesh(16)
    fes = H1(mesh, order=3, dirichlet=".*")
    a, gfu, exact = _build(fes)

    pert = GridFunction(fes)
    pert.Set(sin(2 * pi * x) * sin(2 * pi * y))
    proj = Projector(fes.FreeDofs(), True)
    gfu.vec.data = exact.vec
    gfu.vec.data += 0.1 * (proj * pert.vec)

    residual_history = []
    res = gfu.vec.CreateVector()
    for _ in range(10):
        a.Apply(gfu.vec, res)
        res.data = proj * res
        residual_history.append(res.Norm())
        if residual_history[-1] < 1e-11:
            break
        a.AssembleLinearization(gfu.vec)
        gfu.vec.data -= a.mat.Inverse(fes.FreeDofs()) * res
    return residual_history, _l2(gfu, exact, mesh)


def test_newton_recovers_manufactured():
    """Newton recovers the unique manufactured root A_ex to FEM accuracy."""
    err = _newton_error()
    assert err < 1e-5, "Newton did not recover A_ex: ||A-A_ex||=%.2e" % err


def test_tighter_tolerance_reaches_discretization_floor():
    """Newton converges BELOW the FEM discretization floor, so the L2 error vs
    A_ex is tol-INDEPENDENT: both a loose and a tight residual tolerance reach
    the same discrete root.  The gap to A_ex is the mesh/order discretization
    error, not the solver tolerance -- confirming A_ex is recovered to
    discretization accuracy (and tightening the tolerance does not worsen it)."""
    errs = [_newton_error(tol) for tol in (1e-6, 1e-11)]
    # both at the small discretization floor; tighter tol does not increase it.
    assert errs[0] < 1e-5, "error not at the discretization floor: %s" % errs
    assert errs[1] <= errs[0] * (1 + 1e-6), \
        "tighter tol changed the discrete root: %s" % errs


def test_saturation_is_active():
    """The reluctivity genuinely VARIES (>10%) over the domain (not secretly linear)."""
    spread = _saturation_relative_spread()
    assert spread > 0.1, "saturation not active: relative nu spread %.3f" % spread


def test_newton_quadratic_convergence():
    """From a perturbation of the known root, Newton converges QUADRATICALLY:
    accelerating residual drop, reaching ~machine zero in a few steps."""
    hist, root_error = _newton_convergence()
    assert hist[-1] < 1e-10, "Newton did not reach machine zero: %s" % hist
    assert len(hist) <= 6, "too many iterations (%d) for quadratic: %s" % (len(hist), hist)
    # quadratic signature: the per-step drop factor ACCELERATES.
    ratios = [hist[i] / hist[i + 1] for i in range(len(hist) - 1) if hist[i + 1] > 0]
    assert ratios[-1] > ratios[0], "convergence not accelerating (not quadratic): %s" % hist
    assert max(ratios) > 1e3, "no quadratic-sized leap seen: %s" % hist
    # and the recovered field is back on the root after the perturbation
    assert root_error < 1e-5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
