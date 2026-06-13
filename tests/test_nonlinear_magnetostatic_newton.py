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


def test_newton_recovers_manufactured():
    """Newton recovers the unique manufactured root A_ex to FEM accuracy."""
    from ngsolve.solvers import Newton
    mesh = _mesh(16)
    fes = H1(mesh, order=3, dirichlet=".*")
    a, gfu, Aex = _build(fes)
    Newton(a, gfu, maxit=40, printing=False)
    err = _l2(gfu, Aex, mesh)
    assert err < 1e-5, "Newton did not recover A_ex: ||A-A_ex||=%.2e" % err


def test_tighter_tolerance_reduces_error():
    """A tighter Newton residual tolerance brings the solution CLOSER to A_ex --
    confirming A_ex is the exact root and the only gap is the solver tolerance."""
    from ngsolve.solvers import Newton
    mesh = _mesh(16)
    fes = H1(mesh, order=3, dirichlet=".*")
    errs = []
    for tol in (1e-6, 1e-11):
        a, gfu, Aex = _build(fes)
        Newton(a, gfu, maxit=40, maxerr=tol, printing=False)
        errs.append(_l2(gfu, Aex, mesh))
    assert errs[1] < errs[0], "tighter tol did not reduce error: %s" % errs


def test_saturation_is_active():
    """The reluctivity genuinely VARIES (>10%) over the domain (not secretly linear)."""
    mesh = _mesh(16)
    fes = H1(mesh, order=3, dirichlet=".*")
    Aex = GridFunction(fes); Aex.Set(_AEX)
    g2 = grad(Aex) * grad(Aex)
    area = Integrate(CF(1), mesh)
    mean = Integrate(_nu(g2), mesh) / area
    std = math.sqrt(Integrate((_nu(g2) - mean) ** 2, mesh) / area)
    assert std / mean > 0.1, "saturation not active: nu spread %.3f (mean %.3f)" % (std, mean)


def test_newton_quadratic_convergence():
    """From a perturbation of the known root, Newton converges QUADRATICALLY:
    accelerating residual drop, reaching ~machine zero in a few steps."""
    mesh = _mesh(16)
    fes = H1(mesh, order=3, dirichlet=".*")
    a, gfu, Aex = _build(fes)

    # start AT the root, then perturb the interior (free) dofs by 10%.
    pert = GridFunction(fes); pert.Set(sin(2 * pi * x) * sin(2 * pi * y))
    proj = Projector(fes.FreeDofs(), True)
    gfu.vec.data = Aex.vec
    gfu.vec.data += 0.1 * (proj * pert.vec)

    res = gfu.vec.CreateVector()
    hist = []
    for _ in range(10):
        a.Apply(gfu.vec, res)
        res.data = proj * res
        hist.append(res.Norm())
        if hist[-1] < 1e-11:
            break
        a.AssembleLinearization(gfu.vec)
        gfu.vec.data -= a.mat.Inverse(fes.FreeDofs()) * res

    assert hist[-1] < 1e-10, "Newton did not reach machine zero: %s" % hist
    assert len(hist) <= 6, "too many iterations (%d) for quadratic: %s" % (len(hist), hist)
    # quadratic signature: the per-step drop factor ACCELERATES.
    ratios = [hist[i] / hist[i + 1] for i in range(len(hist) - 1) if hist[i + 1] > 0]
    assert ratios[-1] > ratios[0], "convergence not accelerating (not quadratic): %s" % hist
    assert max(ratios) > 1e3, "no quadratic-sized leap seen: %s" % hist
    # and the recovered field is back on the root after the perturbation
    assert _l2(gfu, Aex, mesh) < 1e-5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
