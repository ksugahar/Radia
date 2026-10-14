"""The large-problem BDDC+CG route must pass the shared true-residual gate."""
import pytest


ng = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")


def _one_sided_dirichlet_box():
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("air")
    box.faces.Min(occ.X).name = "ground"
    return ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.7))


def test_iterative_cohomology_solve_rejects_bad_true_residual(monkeypatch):
    import ngsolve.krylovspace as krylovspace
    import radia.cohomology_cut as cut

    class ZeroInverse:
        def __init__(self, **_kwargs):
            self.iterations = 1

        def __mul__(self, rhs):
            result = rhs.CreateVector()
            result[:] = 0.0
            return result

    mesh = _one_sided_dirichlet_box()
    solver = cut.CohomologyCutSolver()
    solver._mesh = mesh
    solver._h_basis = [ng.CF((1.0, 0.0, 0.0))]
    monkeypatch.setattr(cut, "_ITERATIVE_DOF_THRESHOLD", 0)
    monkeypatch.setattr(krylovspace, "CGSolver", ZeroInverse)

    with ng.TaskManager(), pytest.raises(
        RuntimeError, match=r"cohomology total-potential solve true relative residual"
    ):
        solver.solve([1.0], order=1, dirichlet="ground")
