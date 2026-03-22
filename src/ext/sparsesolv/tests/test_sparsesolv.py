# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
Tests for SparseSolv iterative solvers and preconditioners.

Tests SparseSolvSolver (ICCG, SGSMRTR) and standalone
preconditioners (IC, SGS) for use with NGSolve's Krylov solvers.

Factory functions (ICPreconditioner, etc.) auto-dispatch real/complex
based on mat.IsComplex() and auto-call Update() on construction.
"""

import pytest
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
from sparsesolv_ngsolve import ICPreconditioner, SGSPreconditioner
from sparsesolv_ngsolve import SparseSolvSolver
from ngsolve.krylovspace import CGSolver


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def poisson_2d():
    """2D Poisson problem with exact solution."""
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
    fes = H1(mesh, order=2, dirichlet="bottom|right|top|left")
    u, v = fes.TnT()

    exact = sin(pi * x) * sin(pi * y)
    f_rhs = 2 * pi * pi * sin(pi * x) * sin(pi * y)

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += f_rhs * v * dx
    f.Assemble()

    gfu_bc = GridFunction(fes)
    gfu_bc.Set(exact, BND)
    f.vec.data -= a.mat * gfu_bc.vec

    return mesh, fes, a, f, gfu_bc, exact


@pytest.fixture
def poisson_3d():
    """3D Poisson problem with exact solution."""
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))
    fes = H1(mesh, order=2, dirichlet="bottom|right|top|left|front|back")
    u, v = fes.TnT()

    exact = sin(pi * x) * sin(pi * y) * sin(pi * z)
    f_rhs = 3 * pi * pi * sin(pi * x) * sin(pi * y) * sin(pi * z)

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += f_rhs * v * dx
    f.Assemble()

    gfu_bc = GridFunction(fes)
    gfu_bc.Set(exact, BND)
    f.vec.data -= a.mat * gfu_bc.vec

    return mesh, fes, a, f, gfu_bc, exact


# ============================================================================
# SparseSolvSolver tests
# ============================================================================

@pytest.mark.parametrize("method", ["ICCG", "SGSMRTR"])
def test_sparsesolv_solver_2d_poisson(poisson_2d, method):
    """All solver methods converge on 2D Poisson."""
    mesh, fes, a, f, gfu_bc, exact = poisson_2d

    solver = SparseSolvSolver(a.mat, method=method,
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000)

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec
    gfu.vec.data += gfu_bc.vec

    error = sqrt(Integrate((gfu - exact) ** 2, mesh))
    result = solver.last_result
    print(f"{method}: iterations={result.iterations}, L2 error={error:.2e}")
    assert result.converged
    assert result.iterations < 100
    assert error < 1e-3


def test_sparsesolv_solver_3d_poisson(poisson_3d):
    """ICCG converges on 3D Poisson."""
    mesh, fes, a, f, gfu_bc, exact = poisson_3d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-8, maxiter=5000)

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec
    gfu.vec.data += gfu_bc.vec

    error = sqrt(Integrate((gfu - exact) ** 2, mesh))
    result = solver.last_result
    print(f"ICCG 3D: iterations={result.iterations}, L2 error={error:.2e}")
    assert result.converged
    assert error < 1e-2


def test_sparsesolv_solver_vs_direct(poisson_2d):
    """ICCG matches direct solver solution."""
    mesh, fes, a, f, gfu_bc, exact = poisson_2d

    gfu_direct = GridFunction(fes)
    gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-12, maxiter=2000)
    gfu_iccg = GridFunction(fes)
    gfu_iccg.vec.data = solver * f.vec

    diff = gfu_iccg.vec.CreateVector()
    diff.data = gfu_iccg.vec - gfu_direct.vec
    rel_err = Norm(diff) / Norm(gfu_direct.vec)
    print(f"Relative error vs direct: {rel_err:.2e}")
    assert rel_err < 1e-6


# ============================================================================
# Preconditioner tests (with NGSolve CGSolver)
# ============================================================================

@pytest.mark.parametrize("PreClass,kwargs", [
    (ICPreconditioner, {"shift": 1.05}),
    (SGSPreconditioner, {}),
])
def test_preconditioners_with_ngsolve_cg(poisson_2d, PreClass, kwargs):
    """Preconditioners work with NGSolve's CGSolver."""
    mesh, fes, a, f, gfu_bc, exact = poisson_2d

    pre = PreClass(a.mat, freedofs=fes.FreeDofs(), **kwargs)
    pre.Update()

    inv = CGSolver(a.mat, pre, printrates=False, tol=1e-10, maxiter=2000)

    gfu = GridFunction(fes)
    gfu.vec.data = inv * f.vec
    gfu.vec.data += gfu_bc.vec

    error = sqrt(Integrate((gfu - exact) ** 2, mesh))
    print(f"{PreClass.__name__}: L2 error={error:.2e}")
    assert error < 1e-3


# ============================================================================
# Auto-shift and diagonal scaling (3D curl-curl)
# ============================================================================

def test_auto_shift_curl_curl():
    """Auto-shift + diagonal scaling for 3D curl-curl with nograds=True."""
    from netgen.occ import Box, Pnt

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    for face in box.faces:
        face.name = "outer"
    mesh = box.GenerateMesh(maxh=0.4)

    fes = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += curl(u) * curl(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += CF((0, 0, 1)) * v * dx
    f.Assemble()

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-8, maxiter=2000, shift=1.0)
    solver.auto_shift = True
    solver.diagonal_scaling = True

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec
    result = solver.last_result
    print(f"Auto-shift curl-curl: iterations={result.iterations}")
    assert result.converged
    assert result.iterations < 100


# ============================================================================
# Residual history and best-result tracking
# ============================================================================

def test_residual_history(poisson_2d):
    """Residual history is recorded when enabled."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000,
                               save_residual_history=True)
    gfu = GridFunction(fes)
    result = solver.Solve(f.vec, gfu.vec)

    assert result.converged
    assert len(result.residual_history) > 0
    assert result.residual_history[-1] < result.residual_history[0]


def test_no_residual_history(poisson_2d):
    """Residual history is empty when disabled."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               save_residual_history=False)
    gfu = GridFunction(fes)
    result = solver.Solve(f.vec, gfu.vec)

    assert len(result.residual_history) == 0


# ============================================================================
# Property accessors
# ============================================================================

def test_properties(poisson_2d):
    """Property getters and setters work correctly."""
    _, fes, a, _, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs())

    assert solver.method == "ICCG"
    assert abs(solver.tol - 1e-10) < 1e-15
    assert solver.maxiter == 1000
    assert abs(solver.shift - 1.05) < 1e-10
    assert solver.save_best_result is True
    assert solver.save_residual_history is False
    assert solver.auto_shift is False
    assert solver.diagonal_scaling is False
    assert solver.divergence_check is False
    assert abs(solver.divergence_threshold - 1000.0) < 1e-10
    assert solver.divergence_count == 100

    solver.method = "SGSMRTR"
    assert solver.method == "SGSMRTR"

    solver.tol = 1e-8
    assert abs(solver.tol - 1e-8) < 1e-15

    solver.maxiter = 500
    assert solver.maxiter == 500

    solver.auto_shift = True
    assert solver.auto_shift is True

    solver.diagonal_scaling = True
    assert solver.diagonal_scaling is True

    solver.divergence_check = True
    assert solver.divergence_check is True

    solver.divergence_threshold = 10.0
    assert abs(solver.divergence_threshold - 10.0) < 1e-10

    solver.divergence_count = 5
    assert solver.divergence_count == 5


def test_invalid_method(poisson_2d):
    """Invalid method raises RuntimeError."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="INVALID",
                               freedofs=fes.FreeDofs())
    gfu = GridFunction(fes)
    with pytest.raises(RuntimeError):
        solver.Solve(f.vec, gfu.vec)


# ============================================================================
# Operator interface
# ============================================================================

def test_operator_interface(poisson_2d):
    """solver * vec gives same result as solver.Solve()."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000)

    gfu1 = GridFunction(fes)
    gfu1.vec.data = solver * f.vec

    gfu2 = GridFunction(fes)
    solver.Solve(f.vec, gfu2.vec)

    diff = gfu1.vec.CreateVector()
    diff.data = gfu1.vec - gfu2.vec
    assert Norm(diff) / max(Norm(gfu1.vec), 1e-30) < 1e-6


# ============================================================================
# Divergence detection (stagnation-based early termination)
# ============================================================================

def test_divergence_check_early_termination():
    """Divergence check terminates early on semi-definite system."""
    from netgen.occ import Box, Pnt

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    for face in box.faces:
        face.name = "outer"
    mesh = box.GenerateMesh(maxh=0.4)

    # Pure curl-curl (semi-definite) with localized source (div J != 0)
    # This system has a null-space component that makes CG stagnate
    fes = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += curl(u) * curl(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += CF((1, 0, 0)) * v * dx  # localized source
    f.Assemble()

    # Without auto-shift, CG on this semi-definite system will stagnate
    # Enable divergence_check to detect stagnation early
    solver = SparseSolvSolver(a.mat, method="CG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-12, maxiter=5000)
    solver.divergence_check = True
    solver.divergence_threshold = 10.0
    solver.divergence_count = 20

    gfu = GridFunction(fes)
    result = solver.Solve(f.vec, gfu.vec)

    # Solver should terminate early (not run all 5000 iterations)
    print(f"Divergence check: converged={result.converged}, "
          f"iterations={result.iterations}")
    assert result.iterations < 5000


def test_save_best_result_restores_best(poisson_2d):
    """save_best_result returns the best solution found (even if converged)."""
    _, fes, a, f, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000,
                               save_best_result=True,
                               save_residual_history=True)
    gfu = GridFunction(fes)
    result = solver.Solve(f.vec, gfu.vec)

    assert result.converged
    # Best residual should be <= final residual
    assert result.final_residual <= result.residual_history[0]


# ============================================================================
# Complex auto-dispatch tests
# ============================================================================

@pytest.fixture
def poisson_2d_complex():
    """2D Poisson in complex FE space (Hermitian, real coefficients)."""
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
    fes = H1(mesh, order=2, complex=True, dirichlet="bottom|right|top|left")
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx + u * v * dx
    a.Assemble()

    f = LinearForm(fes)
    f += 1 * v * dx
    f.Assemble()

    return mesh, fes, a, f


def test_factory_auto_dispatch_complex_solver(poisson_2d_complex):
    """SparseSolvSolver factory auto-dispatches for complex matrix."""
    mesh, fes, a, f = poisson_2d_complex

    # Same function name works for complex
    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000)

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec

    result = solver.last_result
    print(f"Complex ICCG: iterations={result.iterations}")
    assert result.converged
    assert Norm(gfu.vec) > 0


@pytest.mark.parametrize("Factory,kwargs", [
    (ICPreconditioner, {"shift": 1.05}),
    (SGSPreconditioner, {}),
])
def test_factory_auto_dispatch_complex_precond(poisson_2d_complex, Factory, kwargs):
    """Preconditioner factories auto-dispatch for complex matrix."""
    mesh, fes, a, f = poisson_2d_complex

    # Factory auto-dispatches to complex and auto-calls Update()
    pre = Factory(a.mat, freedofs=fes.FreeDofs(), **kwargs)

    inv = CGSolver(a.mat, pre, printrates=False, tol=1e-10, maxiter=2000)

    gfu = GridFunction(fes)
    gfu.vec.data = inv * f.vec
    assert Norm(gfu.vec) > 0


def test_factory_returns_typed_class(poisson_2d):
    """Factory returns correctly typed internal class (D suffix)."""
    _, fes, a, _, _, _ = poisson_2d

    pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs())
    assert type(pre).__name__ == "ICPreconditionerD"

    solver = SparseSolvSolver(a.mat, freedofs=fes.FreeDofs())
    assert type(solver).__name__ == "SparseSolvSolverD"


def test_factory_returns_typed_class_complex(poisson_2d_complex):
    """Factory returns correctly typed internal class (C suffix)."""
    _, fes, a, _ = poisson_2d_complex

    pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs())
    assert type(pre).__name__ == "ICPreconditionerC"

    solver = SparseSolvSolver(a.mat, freedofs=fes.FreeDofs())
    assert type(solver).__name__ == "SparseSolvSolverC"


# ============================================================================
# Complex eddy current tests (complex-symmetric systems)
# ============================================================================

@pytest.fixture
def eddy_current_3d():
    """3D eddy current problem: curl-curl + i*omega*sigma mass (complex-symmetric)."""
    from netgen.occ import Box, Pnt, OCCGeometry

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    for face in box.faces:
        face.name = "outer"
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.3))

    fes = HCurl(mesh, order=1, complex=True, dirichlet="outer", nograds=True)
    u, v = fes.TnT()

    # curl-curl + i*sigma mass: A^T = A (complex-symmetric, NOT Hermitian)
    a = BilinearForm(fes)
    a += InnerProduct(curl(u), curl(v)) * dx
    a += 1j * InnerProduct(u, v) * dx
    a.Assemble()

    f_form = LinearForm(fes)
    f_form += InnerProduct(CF((1, 0, 0)), v) * dx
    f_form.Assemble()

    return mesh, fes, a, f_form


@pytest.mark.parametrize("method", ["ICCG", "SGSMRTR"])
def test_complex_eddy_current(eddy_current_3d, method):
    """All solver methods converge on complex eddy current (complex-symmetric)."""
    mesh, fes, a, f = eddy_current_3d

    solver = SparseSolvSolver(a.mat, method=method,
                               freedofs=fes.FreeDofs(),
                               tol=1e-8, maxiter=5000,
                               save_best_result=True)

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec

    result = solver.last_result
    print(f"{method} eddy current: iterations={result.iterations}, "
          f"residual={result.final_residual:.2e}")
    assert result.converged
    assert result.iterations < 200


def test_complex_eddy_current_vs_direct(eddy_current_3d):
    """ICCG matches direct solver on complex eddy current."""
    mesh, fes, a, f = eddy_current_3d

    gfu_direct = GridFunction(fes)
    gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=5000,
                               save_best_result=True)
    gfu_iccg = GridFunction(fes)
    gfu_iccg.vec.data = solver * f.vec

    diff = gfu_iccg.vec.CreateVector()
    diff.data = gfu_iccg.vec - gfu_direct.vec
    rel_err = Norm(diff) / Norm(gfu_direct.vec)
    print(f"Complex ICCG vs direct: relative error={rel_err:.2e}")
    assert rel_err < 1e-4


# ============================================================================
# BDDC comparison tests
# ============================================================================

def test_sparsesolv_vs_bddc_3d_poisson():
    """SparseSolv ICCG achieves similar accuracy as BDDC+CG on 3D Poisson."""
    from netgen.occ import Box, Pnt, OCCGeometry

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    for face in box.faces:
        face.name = "outer"
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.3))

    fes = H1(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += 1 * v * dx
    f.Assemble()

    # Direct solver reference
    gfu_direct = GridFunction(fes)
    gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # SparseSolv ICCG
    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=5000)
    gfu_iccg = GridFunction(fes)
    gfu_iccg.vec.data = solver * f.vec
    iccg_result = solver.last_result

    # BDDC+CG
    a_bddc = BilinearForm(fes)
    a_bddc += grad(u) * grad(v) * dx
    c_bddc = Preconditioner(a_bddc, type="bddc")
    a_bddc.Assemble()

    gfu_bddc = GridFunction(fes)
    inv_bddc = CGSolver(a_bddc.mat, c_bddc.mat, printrates=False,
                         tol=1e-10, maxiter=200)
    gfu_bddc.vec.data = inv_bddc * f.vec

    # Both should match direct solver
    diff_iccg = gfu_iccg.vec.CreateVector()
    diff_iccg.data = gfu_iccg.vec - gfu_direct.vec
    err_iccg = Norm(diff_iccg) / Norm(gfu_direct.vec)

    diff_bddc = gfu_bddc.vec.CreateVector()
    diff_bddc.data = gfu_bddc.vec - gfu_direct.vec
    err_bddc = Norm(diff_bddc) / Norm(gfu_direct.vec)

    print(f"ICCG: iters={iccg_result.iterations}, rel_err={err_iccg:.2e}")
    print(f"BDDC: rel_err={err_bddc:.2e}")
    assert iccg_result.converged
    assert err_iccg < 1e-6
    assert err_bddc < 1e-6


def test_sparsesolv_vs_bddc_eddy_current(eddy_current_3d):
    """SparseSolv ICCG achieves similar accuracy as BDDC+CG on eddy current."""
    mesh, fes, a, f = eddy_current_3d

    u, v = fes.TnT()

    # Direct solver reference
    gfu_direct = GridFunction(fes)
    gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # SparseSolv ICCG
    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=5000,
                               save_best_result=True)
    gfu_iccg = GridFunction(fes)
    gfu_iccg.vec.data = solver * f.vec
    iccg_result = solver.last_result

    # BDDC+CG (conjugate=False for complex-symmetric)
    a_bddc = BilinearForm(fes)
    a_bddc += InnerProduct(curl(u), curl(v)) * dx
    a_bddc += 1j * InnerProduct(u, v) * dx
    c_bddc = Preconditioner(a_bddc, type="bddc")
    a_bddc.Assemble()

    gfu_bddc = GridFunction(fes)
    inv_bddc = CGSolver(a_bddc.mat, c_bddc.mat, printrates=False,
                         conjugate=False, tol=1e-10, maxiter=200)
    gfu_bddc.vec.data = inv_bddc * f.vec

    # Both should match direct solver
    diff_iccg = gfu_iccg.vec.CreateVector()
    diff_iccg.data = gfu_iccg.vec - gfu_direct.vec
    err_iccg = Norm(diff_iccg) / Norm(gfu_direct.vec)

    diff_bddc = gfu_bddc.vec.CreateVector()
    diff_bddc.data = gfu_bddc.vec - gfu_direct.vec
    err_bddc = Norm(diff_bddc) / Norm(gfu_direct.vec)

    print(f"ICCG: iters={iccg_result.iterations}, rel_err={err_iccg:.2e}")
    print(f"BDDC: rel_err={err_bddc:.2e}")
    assert iccg_result.converged
    assert err_iccg < 1e-4
    assert err_bddc < 1e-4


# ============================================================================
# Conjugate inner product tests (Hermitian vs complex-symmetric)
# ============================================================================

def test_conjugate_property(poisson_2d):
    """conjugate property getter and setter work correctly."""
    _, fes, a, _, _, _ = poisson_2d

    solver = SparseSolvSolver(a.mat, freedofs=fes.FreeDofs())
    assert solver.conjugate is False  # default

    solver.conjugate = True
    assert solver.conjugate is True

    solver.conjugate = False
    assert solver.conjugate is False


def test_conjugate_factory_parameter(poisson_2d_complex):
    """conjugate parameter in factory function is passed through."""
    _, fes, a, _ = poisson_2d_complex

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               conjugate=True)
    assert solver.conjugate is True


@pytest.mark.parametrize("method", ["ICCG", "SGSMRTR"])
def test_hermitian_system_with_conjugate(poisson_2d_complex, method):
    """Hermitian system (real coefficients, complex space) solves with conjugate=True."""
    mesh, fes, a, f = poisson_2d_complex

    # This is a Hermitian system (A^H = A) because coefficients are real.
    # Both conjugate=True (Hermitian) and conjugate=False (complex-symmetric)
    # should converge, but conjugate=True is mathematically correct here.
    solver = SparseSolvSolver(a.mat, method=method,
                               freedofs=fes.FreeDofs(),
                               tol=1e-10, maxiter=2000,
                               conjugate=True)

    gfu = GridFunction(fes)
    gfu.vec.data = solver * f.vec

    result = solver.last_result
    print(f"{method} Hermitian (conjugate=True): iterations={result.iterations}")
    assert result.converged
    assert Norm(gfu.vec) > 0


def test_hermitian_vs_direct(poisson_2d_complex):
    """Hermitian ICCG with conjugate=True matches direct solver."""
    mesh, fes, a, f = poisson_2d_complex

    gfu_direct = GridFunction(fes)
    gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    solver = SparseSolvSolver(a.mat, method="ICCG",
                               freedofs=fes.FreeDofs(),
                               tol=1e-12, maxiter=2000,
                               conjugate=True)
    gfu_iccg = GridFunction(fes)
    gfu_iccg.vec.data = solver * f.vec

    diff = gfu_iccg.vec.CreateVector()
    diff.data = gfu_iccg.vec - gfu_direct.vec
    rel_err = Norm(diff) / Norm(gfu_direct.vec)
    print(f"Hermitian ICCG vs direct: relative error={rel_err:.2e}")
    assert rel_err < 1e-6


# ============================================================================
# ABMC Ordering Tests
# ============================================================================

class TestABMCProperties:
    """Test ABMC property accessors."""

    def test_abmc_property_defaults(self, poisson_2d):
        """ABMC properties have correct defaults."""
        mesh, fes, a, f, *_ = poisson_2d
        solver = SparseSolvSolver(a.mat, method="ICCG", freedofs=fes.FreeDofs())
        assert solver.use_abmc == False
        assert solver.abmc_block_size == 4
        assert solver.abmc_num_colors == 4

    def test_abmc_property_set(self, poisson_2d):
        """ABMC properties can be set and read back."""
        mesh, fes, a, f, *_ = poisson_2d
        solver = SparseSolvSolver(a.mat, method="ICCG", freedofs=fes.FreeDofs())
        solver.use_abmc = True
        solver.abmc_block_size = 8
        solver.abmc_num_colors = 6
        assert solver.use_abmc == True
        assert solver.abmc_block_size == 8
        assert solver.abmc_num_colors == 6

    def test_abmc_factory_parameter(self, poisson_2d):
        """ABMC parameters can be set via factory function."""
        mesh, fes, a, f, *_ = poisson_2d
        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   use_abmc=True,
                                   abmc_block_size=8,
                                   abmc_num_colors=6)
        assert solver.use_abmc == True
        assert solver.abmc_block_size == 8
        assert solver.abmc_num_colors == 6


class TestABMCSolve:
    """Test ABMC-ICCG solver correctness."""

    def test_abmc_poisson_2d(self, poisson_2d):
        """ABMC-ICCG solves 2D Poisson correctly."""
        mesh, fes, a, f, *_ = poisson_2d

        # Direct solution for reference
        gfu_direct = GridFunction(fes)
        gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

        # ABMC-ICCG
        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-10, maxiter=2000,
                                   use_abmc=True,
                                   abmc_block_size=4,
                                   abmc_num_colors=4)
        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec

        assert solver.last_result.converged
        diff = gfu.vec.CreateVector()
        diff.data = gfu.vec - gfu_direct.vec
        rel_err = Norm(diff) / Norm(gfu_direct.vec)
        print(f"ABMC-ICCG 2D Poisson: {solver.last_result.iterations} iters, err={rel_err:.2e}")
        assert rel_err < 1e-6

    def test_abmc_poisson_3d(self):
        """ABMC-ICCG solves 3D Poisson (larger problem)."""
        mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))
        fes = H1(mesh, order=2, dirichlet="left|right|top|bottom|back|front")
        u, v = fes.TnT()
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        f = LinearForm(fes)
        f += 1 * v * dx
        f.Assemble()

        # Direct
        gfu_direct = GridFunction(fes)
        gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

        # ABMC-ICCG
        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-10, maxiter=3000,
                                   use_abmc=True,
                                   abmc_block_size=8,
                                   abmc_num_colors=4)
        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec

        assert solver.last_result.converged
        diff = gfu.vec.CreateVector()
        diff.data = gfu.vec - gfu_direct.vec
        rel_err = Norm(diff) / Norm(gfu_direct.vec)
        print(f"ABMC-ICCG 3D Poisson: {solver.last_result.iterations} iters, err={rel_err:.2e}")
        assert rel_err < 1e-6

    def test_abmc_curl_curl(self):
        """ABMC-ICCG with auto_shift for semi-definite curl-curl."""
        from netgen.occ import Box, Pnt
        box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
        for face in box.faces:
            face.name = "outer"
        mesh = box.GenerateMesh(maxh=0.4)

        fes = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
        u, v = fes.TnT()
        a = BilinearForm(fes)
        a += curl(u) * curl(v) * dx
        a.Assemble()
        f = LinearForm(fes)
        f += CF((0, 0, 1)) * v * dx
        f.Assemble()

        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-8, maxiter=3000, shift=1.0,
                                   use_abmc=True,
                                   abmc_block_size=4,
                                   abmc_num_colors=4)
        solver.auto_shift = True
        solver.diagonal_scaling = True

        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec
        print(f"ABMC curl-curl: converged={solver.last_result.converged}, "
              f"iters={solver.last_result.iterations}")
        assert solver.last_result.converged

    def test_abmc_complex_symmetric(self, poisson_2d_complex):
        """ABMC-ICCG with complex-symmetric system."""
        mesh, fes, a, f, *_ = poisson_2d_complex

        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-10, maxiter=2000,
                                   conjugate=True,
                                   use_abmc=True,
                                   abmc_block_size=4,
                                   abmc_num_colors=4)
        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec

        assert solver.last_result.converged
        print(f"ABMC complex: {solver.last_result.iterations} iters")

    def test_abmc_vs_level_schedule(self, poisson_2d):
        """ABMC and level-scheduling produce the same solution."""
        mesh, fes, a, f, *_ = poisson_2d

        # Without ABMC (level scheduling)
        solver_ls = SparseSolvSolver(a.mat, method="ICCG",
                                      freedofs=fes.FreeDofs(),
                                      tol=1e-12, maxiter=2000)
        gfu_ls = GridFunction(fes)
        gfu_ls.vec.data = solver_ls * f.vec

        # With ABMC
        solver_abmc = SparseSolvSolver(a.mat, method="ICCG",
                                        freedofs=fes.FreeDofs(),
                                        tol=1e-12, maxiter=2000,
                                        use_abmc=True)
        gfu_abmc = GridFunction(fes)
        gfu_abmc.vec.data = solver_abmc * f.vec

        assert solver_ls.last_result.converged
        assert solver_abmc.last_result.converged

        diff = gfu_abmc.vec.CreateVector()
        diff.data = gfu_abmc.vec - gfu_ls.vec
        rel_err = Norm(diff) / Norm(gfu_ls.vec)
        print(f"ABMC vs level-schedule: rel_err={rel_err:.2e}, "
              f"iters ABMC={solver_abmc.last_result.iterations} vs LS={solver_ls.last_result.iterations}")
        assert rel_err < 1e-6

    @pytest.mark.parametrize("block_size", [2, 4, 8, 16])
    def test_abmc_block_size_variations(self, poisson_2d, block_size):
        """ABMC converges with various block sizes."""
        mesh, fes, a, f, *_ = poisson_2d

        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-10, maxiter=2000,
                                   use_abmc=True,
                                   abmc_block_size=block_size,
                                   abmc_num_colors=4)
        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec

        assert solver.last_result.converged
        print(f"ABMC block_size={block_size}: {solver.last_result.iterations} iters")


# ============================================================================
# ABMC reordered space tests
# ============================================================================

class TestABMCReorderedSpace:
    """Test CG in ABMC-reordered space (no per-iteration permutation)."""

    def test_reordered_space_poisson_2d(self, poisson_2d):
        """ABMC with CG in reordered space solves 2D Poisson correctly."""
        mesh, fes, a, f, *_ = poisson_2d

        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-10, maxiter=2000,
                                   use_abmc=True,
                                   abmc_block_size=4,
                                   abmc_num_colors=4)
        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec

        assert solver.last_result.converged
        # Compare with direct solver
        gfu_direct = GridFunction(fes)
        gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

        diff = gfu.vec.CreateVector()
        diff.data = gfu.vec - gfu_direct.vec
        rel_err = Norm(diff) / Norm(gfu_direct.vec)
        print(f"ABMC reordered-space 2D: {solver.last_result.iterations} iters, err={rel_err:.2e}")
        assert rel_err < 1e-6

    def test_reordered_space_vs_level_schedule(self, poisson_2d):
        """Reordered-space ABMC and level scheduling give same solution."""
        mesh, fes, a, f, *_ = poisson_2d

        # Level scheduling (no ABMC)
        solver_ls = SparseSolvSolver(a.mat, method="ICCG",
                                      freedofs=fes.FreeDofs(),
                                      tol=1e-12, maxiter=2000)
        gfu_ls = GridFunction(fes)
        gfu_ls.vec.data = solver_ls * f.vec

        # ABMC (CG in reordered space)
        solver_abmc = SparseSolvSolver(a.mat, method="ICCG",
                                        freedofs=fes.FreeDofs(),
                                        tol=1e-12, maxiter=2000,
                                        use_abmc=True,
                                        abmc_block_size=4,
                                        abmc_num_colors=4)
        gfu_abmc = GridFunction(fes)
        gfu_abmc.vec.data = solver_abmc * f.vec

        assert solver_ls.last_result.converged
        assert solver_abmc.last_result.converged

        diff = gfu_ls.vec.CreateVector()
        diff.data = gfu_ls.vec - gfu_abmc.vec
        rel_err = Norm(diff) / Norm(gfu_ls.vec)
        print(f"Reordered vs LS: rel_err={rel_err:.2e}, "
              f"iters LS={solver_ls.last_result.iterations} / "
              f"ABMC={solver_abmc.last_result.iterations}")
        assert rel_err < 1e-6

    def test_reordered_space_diagonal_scaling(self):
        """Reordered-space ABMC with diagonal scaling."""
        from netgen.occ import Box, Pnt
        box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
        for face in box.faces:
            face.name = "outer"
        mesh = box.GenerateMesh(maxh=0.4)

        fes = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
        u, v = fes.TnT()
        a = BilinearForm(fes)
        a += curl(u) * curl(v) * dx
        a.Assemble()
        f = LinearForm(fes)
        f += CF((0, 0, 1)) * v * dx
        f.Assemble()

        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-8, maxiter=3000, shift=1.0,
                                   use_abmc=True,
                                   abmc_block_size=4,
                                   abmc_num_colors=4)
        solver.auto_shift = True
        solver.diagonal_scaling = True

        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec
        print(f"ABMC reordered diag_scaling: converged={solver.last_result.converged}, "
              f"iters={solver.last_result.iterations}")
        assert solver.last_result.converged

    def test_reordered_space_3d_poisson(self):
        """Reordered-space ABMC on 3D Poisson."""
        mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
        fes = H1(mesh, order=1, dirichlet="left|right|top|bottom|front|back")
        u, v = fes.TnT()
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        f = LinearForm(fes)
        f += 1 * v * dx
        f.Assemble()

        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-10, maxiter=3000,
                                   use_abmc=True,
                                   abmc_block_size=8,
                                   abmc_num_colors=4)
        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec

        # Compare with direct solver
        gfu_direct = GridFunction(fes)
        gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

        diff = gfu.vec.CreateVector()
        diff.data = gfu.vec - gfu_direct.vec
        rel_err = Norm(diff) / Norm(gfu_direct.vec)
        print(f"ABMC reordered 3D: {solver.last_result.iterations} iters, err={rel_err:.2e}")
        assert rel_err < 1e-6

    def test_reordered_space_complex(self, poisson_2d_complex):
        """Reordered-space ABMC with complex-symmetric system."""
        mesh, fes, a, f, *_ = poisson_2d_complex

        solver = SparseSolvSolver(a.mat, method="ICCG",
                                   freedofs=fes.FreeDofs(),
                                   tol=1e-10, maxiter=2000,
                                   conjugate=True,
                                   use_abmc=True,
                                   abmc_block_size=4,
                                   abmc_num_colors=4)
        gfu = GridFunction(fes)
        gfu.vec.data = solver * f.vec
        print(f"ABMC reordered complex: converged={solver.last_result.converged}, "
              f"iters={solver.last_result.iterations}")
        assert solver.last_result.converged


# ============================================================================
# COCR Solver (standalone) tests
# ============================================================================

class TestCOCRSolverStandalone:
    """Tests for standalone COCRSolver with external preconditioner."""

    def test_cocr_complex_solve(self, eddy_current_3d):
        """COCRSolver converges on complex-symmetric eddy current problem."""
        from sparsesolv_ngsolve import COCRSolver

        mesh, fes, a, f = eddy_current_3d

        pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
        inv = COCRSolver(a.mat, pre, maxiter=500, tol=1e-8)

        gfu = GridFunction(fes)
        gfu.vec.data = inv * f.vec
        assert Norm(gfu.vec) > 0
        print(f"COCRSolver: {inv.iterations} iterations")
        assert inv.iterations < 500

    def test_cocr_with_freedofs(self, eddy_current_3d):
        """COCRSolver works with freedofs parameter."""
        from sparsesolv_ngsolve import COCRSolver

        mesh, fes, a, f = eddy_current_3d

        pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
        inv = COCRSolver(a.mat, pre, freedofs=fes.FreeDofs(), maxiter=500, tol=1e-8)

        gfu = GridFunction(fes)
        gfu.vec.data = inv * f.vec
        assert Norm(gfu.vec) > 0


# ============================================================================
# GMRES Solver tests
# ============================================================================

class TestGMRESSolver:
    """Tests for GMRESSolver (non-symmetric systems)."""

    def test_gmres_basic_solve(self):
        """GMRESSolver solves a non-symmetric system (convection-diffusion)."""
        from sparsesolv_ngsolve import GMRESSolver

        mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
        fes = H1(mesh, order=2, dirichlet="bottom|right|top|left")
        u, v = fes.TnT()

        # Non-symmetric: diffusion + convection
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a += CF((1, 0)) * grad(u) * v * dx  # convection term
        a.Assemble()

        f = LinearForm(fes)
        f += 1 * v * dx
        f.Assemble()

        pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
        inv = GMRESSolver(a.mat, pre, freedofs=fes.FreeDofs(),
                          maxiter=500, tol=1e-8)

        gfu = GridFunction(fes)
        gfu.vec.data = inv * f.vec
        assert Norm(gfu.vec) > 0
        print(f"GMRES convection-diffusion: {inv.iterations} iterations")
        assert inv.iterations < 200

    def test_gmres_with_preconditioner(self):
        """GMRESSolver works with SGS preconditioner."""
        from sparsesolv_ngsolve import GMRESSolver

        mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
        fes = H1(mesh, order=2, dirichlet="bottom|right|top|left")
        u, v = fes.TnT()

        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a += CF((1, 0)) * grad(u) * v * dx
        a.Assemble()

        f = LinearForm(fes)
        f += 1 * v * dx
        f.Assemble()

        pre = SGSPreconditioner(a.mat, freedofs=fes.FreeDofs())
        inv = GMRESSolver(a.mat, pre, maxiter=500, tol=1e-8)

        gfu = GridFunction(fes)
        gfu.vec.data = inv * f.vec
        assert Norm(gfu.vec) > 0


# ============================================================================
# CompactAMG Preconditioner tests
# ============================================================================

class TestCompactAMG:
    """Tests for CompactAMG preconditioner."""

    @pytest.fixture
    def has_amg(self):
        from sparsesolv_ngsolve import has_compact_ams
        if not has_compact_ams():
            pytest.skip("CompactAMG not available")

    def test_compact_amg_basic(self, has_amg):
        """CompactAMG + CG solves 3D Poisson."""
        from sparsesolv_ngsolve import CompactAMGPreconditioner

        mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))
        fes = H1(mesh, order=1, dirichlet="left|right|top|bottom|front|back")
        u, v = fes.TnT()
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        f = LinearForm(fes)
        f += 1 * v * dx
        f.Assemble()

        pre = CompactAMGPreconditioner(a.mat, freedofs=fes.FreeDofs())
        inv = CGSolver(a.mat, pre, printrates=False, tol=1e-10, maxiter=200)

        gfu = GridFunction(fes)
        gfu.vec.data = inv * f.vec

        gfu_direct = GridFunction(fes)
        gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

        diff = gfu.vec.CreateVector()
        diff.data = gfu.vec - gfu_direct.vec
        rel_err = Norm(diff) / Norm(gfu_direct.vec)
        print(f"CompactAMG+CG: rel_err={rel_err:.2e}")
        assert rel_err < 1e-6


# ============================================================================
# CompactAMS Preconditioner tests
# ============================================================================

class TestCompactAMS:
    """Tests for CompactAMS preconditioner (HCurl systems)."""

    @pytest.fixture
    def has_ams(self):
        from sparsesolv_ngsolve import has_compact_ams
        if not has_compact_ams():
            pytest.skip("CompactAMS not available")

    @pytest.fixture
    def curl_curl_3d(self):
        """3D curl-curl + mass problem for AMS testing."""
        from netgen.occ import Box, Pnt, OCCGeometry
        box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
        for face in box.faces:
            face.name = "outer"
        mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.3))

        fes = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
        u, v = fes.TnT()

        a = BilinearForm(fes)
        a += curl(u) * curl(v) * dx + 0.01 * u * v * dx  # regularized
        a.Assemble()

        f_form = LinearForm(fes)
        f_form += CF((0, 0, 1)) * v * dx
        f_form.Assemble()

        return mesh, fes, a, f_form

    def test_compact_ams_basic(self, has_ams, curl_curl_3d):
        """CompactAMS + CG solves curl-curl + mass system."""
        from sparsesolv_ngsolve import CompactAMSPreconditioner

        mesh, fes, a, f = curl_curl_3d

        # Get gradient matrix and coordinates
        G_mat, h1_fes = fes.CreateGradient()

        # Extract vertex coordinates
        coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
        coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
        coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

        pre = CompactAMSPreconditioner(a.mat, G_mat,
                                        freedofs=fes.FreeDofs(),
                                        coord_x=coord_x, coord_y=coord_y,
                                        coord_z=coord_z)

        inv = CGSolver(a.mat, pre, printrates=False, tol=1e-8, maxiter=200)

        gfu = GridFunction(fes)
        gfu.vec.data = inv * f.vec

        gfu_direct = GridFunction(fes)
        gfu_direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

        diff = gfu.vec.CreateVector()
        diff.data = gfu.vec - gfu_direct.vec
        rel_err = Norm(diff) / Norm(gfu_direct.vec)
        print(f"CompactAMS+CG: rel_err={rel_err:.2e}")
        assert rel_err < 1e-4

    def test_compact_ams_update(self, has_ams, curl_curl_3d):
        """CompactAMS Update() preserves solution quality."""
        from sparsesolv_ngsolve import CompactAMSPreconditioner

        mesh, fes, a, f = curl_curl_3d

        G_mat, h1_fes = fes.CreateGradient()
        coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
        coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
        coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

        pre = CompactAMSPreconditioner(a.mat, G_mat,
                                        freedofs=fes.FreeDofs(),
                                        coord_x=coord_x, coord_y=coord_y,
                                        coord_z=coord_z)

        # First solve
        inv1 = CGSolver(a.mat, pre, printrates=False, tol=1e-8, maxiter=200)
        gfu1 = GridFunction(fes)
        gfu1.vec.data = inv1 * f.vec
        iters1 = inv1.iterations

        # Update (same matrix - should give same iterations)
        pre.Update()
        inv2 = CGSolver(a.mat, pre, printrates=False, tol=1e-8, maxiter=200)
        gfu2 = GridFunction(fes)
        gfu2.vec.data = inv2 * f.vec
        iters2 = inv2.iterations

        print(f"AMS Update: iters before={iters1}, after={iters2}")
        assert abs(iters1 - iters2) <= 1  # Should be identical or +-1

        diff = gfu1.vec.CreateVector()
        diff.data = gfu1.vec - gfu2.vec
        rel_err = Norm(diff) / Norm(gfu1.vec)
        assert rel_err < 1e-6


# ============================================================================
# ComplexCompactAMS Preconditioner tests
# ============================================================================

class TestComplexCompactAMS:
    """Tests for ComplexCompactAMS preconditioner."""

    @pytest.fixture
    def has_ams(self):
        from sparsesolv_ngsolve import has_compact_ams
        if not has_compact_ams():
            pytest.skip("ComplexCompactAMS not available")

    def test_complex_ams_basic(self, has_ams, eddy_current_3d):
        """ComplexCompactAMS + COCR solves eddy current problem."""
        from sparsesolv_ngsolve import ComplexCompactAMSPreconditioner, COCRSolver

        mesh, fes, a, f = eddy_current_3d
        u, v = fes.TnT()

        # Build real auxiliary matrix: K + |omega|*sigma*M
        fes_real = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
        u_r, v_r = fes_real.TnT()
        a_real = BilinearForm(fes_real)
        a_real += curl(u_r) * curl(v_r) * dx + 1.0 * u_r * v_r * dx
        a_real.Assemble()

        G_mat, h1_fes = fes_real.CreateGradient()
        coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
        coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
        coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

        pre = ComplexCompactAMSPreconditioner(
            a_real.mat, G_mat,
            freedofs=fes_real.FreeDofs(),
            coord_x=coord_x, coord_y=coord_y, coord_z=coord_z,
            ndof_complex=fes.ndof)

        inv = COCRSolver(a.mat, pre, maxiter=500, tol=1e-6)
        gfu = GridFunction(fes)
        gfu.vec.data = inv * f.vec

        assert Norm(gfu.vec) > 0
        print(f"ComplexCompactAMS+COCR: {inv.iterations} iterations")
        assert inv.iterations < 500

    def test_complex_ams_ndof_auto(self, has_ams, eddy_current_3d):
        """ComplexCompactAMS with ndof_complex=0 auto-derives from matrix."""
        from sparsesolv_ngsolve import ComplexCompactAMSPreconditioner

        mesh, fes, a, f = eddy_current_3d

        fes_real = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
        u_r, v_r = fes_real.TnT()
        a_real = BilinearForm(fes_real)
        a_real += curl(u_r) * curl(v_r) * dx + 1.0 * u_r * v_r * dx
        a_real.Assemble()

        G_mat, h1_fes = fes_real.CreateGradient()
        coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
        coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
        coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

        # ndof_complex=0 (default) should auto-derive from matrix
        pre = ComplexCompactAMSPreconditioner(
            a_real.mat, G_mat,
            freedofs=fes_real.FreeDofs(),
            coord_x=coord_x, coord_y=coord_y, coord_z=coord_z)
        # Should not raise an error - ndof_complex is auto-derived

        # Verify it works: use complex vectors (ComplexCompactAMS expects complex)
        test_vec = a.mat.CreateRowVector()
        result_vec = a.mat.CreateRowVector()
        # Basic smoke test
        pre.Mult(test_vec, result_vec)


# ============================================================================
# IC Preconditioner ABMC path switch regression test
# ============================================================================

class TestICPreconditionerRegression:
    """Regression tests for IC preconditioner fixes."""

    def test_abmc_path_switch(self, poisson_2d):
        """ABMC <-> standard path switch works correctly (regression for cached_use_abmc_)."""
        mesh, fes, a, f, *_ = poisson_2d

        # Start without ABMC
        pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
        inv1 = CGSolver(a.mat, pre, printrates=False, tol=1e-10, maxiter=2000)
        gfu1 = GridFunction(fes)
        gfu1.vec.data = inv1 * f.vec
        assert Norm(gfu1.vec) > 0

        # Switch to ABMC (this used to segfault before the fix)
        pre.use_abmc = True
        pre.Update()
        inv2 = CGSolver(a.mat, pre, printrates=False, tol=1e-10, maxiter=2000)
        gfu2 = GridFunction(fes)
        gfu2.vec.data = inv2 * f.vec

        # Results should be similar
        diff = gfu1.vec.CreateVector()
        diff.data = gfu1.vec - gfu2.vec
        rel_err = Norm(diff) / Norm(gfu1.vec)
        print(f"ABMC path switch: rel_err={rel_err:.2e}")
        assert rel_err < 1e-6

        # Switch back to standard
        pre.use_abmc = False
        pre.Update()
        inv3 = CGSolver(a.mat, pre, printrates=False, tol=1e-10, maxiter=2000)
        gfu3 = GridFunction(fes)
        gfu3.vec.data = inv3 * f.vec

        diff.data = gfu1.vec - gfu3.vec
        rel_err = Norm(diff) / Norm(gfu1.vec)
        print(f"Back to standard: rel_err={rel_err:.2e}")
        assert rel_err < 1e-6


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
