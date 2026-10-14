"""Weak source projection and failure contracts for reduced-A Newton."""

import numpy as np
import pytest

def _unit_mesh():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("iron")
    with ng.TaskManager():
        return ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.8))


def test_hdiv_source_projection_is_autodiffdiff_compatible():
    ng = pytest.importorskip("ngsolve")
    from radia.vector_potential_solver import VectorPotentialSolver

    mesh = _unit_mesh()
    solver = VectorPotentialSolver(mesh, iron_domains="iron", order=1)
    solver.set_source_cf(ng.CF((ng.x, 2 * ng.y, 3 * ng.z)))
    with ng.TaskManager():
        source = solver.materialize_source_hdiv(order=1)

    assert source.space.ndof > 0
    assert solver._B_source_cf is source
    assert solver._B_source_gridfunction is source
    assert solver._source_projection_diagnostics["kind"] == "weak-l2"
    assert solver._source_projection_diagnostics["converged"]
    assert solver._source_projection_diagnostics["iterations"] >= 0
    assert solver._source_projection_diagnostics["initial_guess"] == \
        "ngsolve-hdiv-interpolation"
    assert solver._source_projection_diagnostics[
        "initial_algebraic_relative_residual"] >= 0.0
    assert solver._source_projection_diagnostics[
        "algebraic_relative_residual"] < 2.0e-11
    sampled = np.asarray(source(mesh(0.25, 0.25, 0.25)), dtype=float)
    assert sampled == pytest.approx([0.25, 0.5, 0.75], abs=2.0e-12)

    fes = ng.HCurl(mesh, order=1, nograds=True)
    trial = fes.TrialFunction()
    total_b = source + ng.curl(trial)
    energy = ng.BilinearForm(fes, symmetric=True)
    energy += ng.SymbolicEnergy(ng.InnerProduct(total_b, total_b))
    state = ng.GridFunction(fes)
    with ng.TaskManager():
        energy.AssembleLinearization(state.vec)


def test_hdiv_source_projection_satisfies_weak_l2_equations():
    ng = pytest.importorskip("ngsolve")
    from radia.vector_potential_solver import VectorPotentialSolver

    mesh = _unit_mesh()
    source_cf = ng.CF((ng.exp(ng.x), ng.sin(ng.y), ng.cos(ng.z)))
    solver = VectorPotentialSolver(mesh, iron_domains="iron", order=1)
    solver.set_source_cf(source_cf)
    with ng.TaskManager():
        source = solver.materialize_source_hdiv(order=1)
        trial = source.space.TrialFunction()
        test = source.space.TestFunction()
        mass = ng.BilinearForm(source.space, symmetric=True)
        mass += ng.InnerProduct(trial, test) * ng.dx
        rhs = ng.LinearForm(source.space)
        rhs += ng.InnerProduct(source_cf, test) * ng.dx(bonus_intorder=3)
        mass.Assemble()
        rhs.Assemble()

    residual = rhs.vec.CreateVector()
    residual.data = mass.mat * source.vec - rhs.vec
    relative_residual = np.linalg.norm(residual.FV().NumPy()) / max(
        np.linalg.norm(rhs.vec.FV().NumPy()), 1.0e-30)
    assert relative_residual < 2.0e-11


def test_hdiv_source_projection_requires_a_source_and_positive_order():
    from radia.vector_potential_solver import VectorPotentialSolver

    solver = VectorPotentialSolver(_unit_mesh(), iron_domains="iron", order=1)
    with pytest.raises(RuntimeError, match="Set source field first"):
        solver.materialize_source_hdiv()
    import ngsolve as ng

    solver.set_source_cf(ng.CF((0.0, 0.0, 0.0)))
    with pytest.raises(ValueError, match="at least one"):
        solver.materialize_source_hdiv(order=0)
    with pytest.raises(ValueError, match="tolerance"):
        solver.materialize_source_hdiv(tol=0.0)
    with pytest.raises(ValueError, match="maxiter"):
        solver.materialize_source_hdiv(maxiter=0)
    with pytest.raises(ValueError, match="bonus_intorder"):
        solver.materialize_source_hdiv(bonus_intorder=-1)
    solver._kelvin_region = "kelvin"
    with pytest.raises(ValueError, match="physical-only mesh"):
        solver.materialize_source_hdiv()


def test_hdiv_source_projection_persists_failed_cg_diagnostics():
    import ngsolve as ng
    from radia.vector_potential_solver import VectorPotentialSolver

    solver = VectorPotentialSolver(_unit_mesh(), iron_domains="iron", order=1)
    solver.set_source_cf(ng.CF((ng.exp(ng.x), ng.sin(ng.y), ng.cos(ng.z))))
    with ng.TaskManager(), pytest.raises(RuntimeError) as captured:
        solver.materialize_source_hdiv(order=1, tol=1.0e-30, maxiter=1)
    diagnostics = solver._source_projection_diagnostics
    assert diagnostics["converged"] is False
    assert diagnostics["iterations"] >= 1
    assert diagnostics["algebraic_relative_residual"] > 0.0
    assert captured.value.source_projection_diagnostics == diagnostics
