"""Forward B(H) Galerkin Newton without a single-valued H(M) assumption.

NGSolve owns quadrature, sparse assembly and GMRES. The configured HDiv
operator owns the demagnetizing form. Unknown H and projected M share the
HDiv space: G H + N projection(M(H)) = G Hext.
"""
import numpy as np
import ngsolve as ng
from ngsolve.comp import IntegrationRuleSpace
from ngsolve.krylovspace import GMRESSolver

from ._nonlinear import _bh_table_funcs, _MU0


def solve_forward_newton(mesh, fes, table, operator, h_ext, *,
                         tol, maxit, nl_tol, nl_maxit):
    """Return projected M, iteration count and checked residual diagnostics.

    This first route accepts one isotropic material table. Region dictionaries
    fail explicitly. It retains negative differential magnetization; GMRES is
    necessary because the forward Jacobian is not symmetric.
    """
    if isinstance(table, dict):
        raise NotImplementedError("forward-newton currently requires one B(H) table")
    arr = np.asarray(table, dtype=float)
    _, bfun, bder, hmax, mmax = _bh_table_funcs(arr[:, 0], arr[:, 1])
    u, v = fes.TnT()
    space = IntegrationRuleSpace(mesh, order=max(6, 2 * int(fes.globalorder) + 2))
    rules = space.GetIntegrationRules()
    measure = ng.dx(intrules=rules)
    nq = space.ndof
    gh = ng.GridFunction(fes)
    samples = ng.GridFunction(ng.VectorValued(space, dim=3))
    sec = ng.GridFunction(space)
    tangent = ng.GridFunction(ng.VectorValued(space, dim=9))
    mass = ng.BilinearForm(fes)
    mass += ng.InnerProduct(u, v) * measure
    mass.Assemble()
    G = mass.mat
    Gi = G.Inverse(inverse="sparsecholesky")

    def apply(matrix, coefficients):
        x = matrix.CreateRowVector()
        x.FV().NumPy()[:] = coefficients
        y = matrix.CreateColVector()
        matrix.Mult(x, y)
        return y.FV().NumPy().copy()

    def demag(coefficients):
        return np.asarray(operator.apply_configured_demag(
            np.ascontiguousarray(coefficients, dtype=float), True), dtype=float)

    probe = np.cos(np.arange(fes.ndof))
    native_mass = np.asarray(operator.apply_configured_geometry_mass(probe))
    mass_error = float(np.linalg.norm(apply(G, probe) - native_mass)
                       / max(np.linalg.norm(native_mass), 1e-300))
    if not np.isfinite(mass_error) or mass_error > 1e-11:
        raise RuntimeError("forward-newton quadrature and configured geometry masses differ")
    rhs = apply(G, h_ext)
    scale = max(float(np.linalg.norm(rhs)), 1e-300)

    def material(h):
        gh.vec.FV().NumPy()[:] = h
        samples.Interpolate(gh)
        values = samples.vec.FV().NumPy().reshape(3, nq).T
        mag = np.linalg.norm(values, axis=1)
        safe = np.maximum(mag, 1e-30)
        clipped = np.minimum(mag, hmax)
        mm = np.where(mag > hmax, mmax, bfun(clipped) / _MU0 - clipped)
        ss = mm / safe
        ss[mag < 1e-12] = float(bder(0)) / _MU0 - 1
        dd = np.where(mag > hmax, 0., bder(clipped) / _MU0 - 1)
        unit = values / safe[:, None]
        tensors = (ss[:, None, None] * np.eye(3)
                   + (dd - ss)[:, None, None] * unit[:, :, None] * unit[:, None, :])
        if not np.isfinite(tensors).all():
            raise RuntimeError("forward-newton non-finite constitutive tangent")
        sec.vec.FV().NumPy()[:] = ss
        tangent.vec.FV().NumPy().reshape(9, nq)[:] = tensors.reshape(nq, 9).T
        load = ng.LinearForm(fes)
        load += ng.InnerProduct(sec * gh, v) * measure
        load.Assemble()
        return apply(Gi, load.vec.FV().NumPy())

    def state(h):
        m = material(h)
        return apply(G, h) + demag(m) - rhs, m

    h = np.asarray(h_ext, dtype=float).copy()
    residual, m = state(h)
    linear_audit, history = [], []
    for iteration in range(nl_maxit + 1):
        relative = float(np.linalg.norm(residual) / scale)
        history.append(relative)
        if not np.isfinite(relative):
            raise RuntimeError("forward-newton non-finite nonlinear residual")
        print(f"forward-newton {iteration}: true nonlinear residual {relative:.9e}", flush=True)
        if relative <= nl_tol:
            projected = material(h)
            constitutive = float(np.linalg.norm(apply(G, m - projected))
                                 / max(np.linalg.norm(apply(G, m)), 1e-300))
            final = float(np.linalg.norm(apply(G, h) + demag(m) - rhs) / scale)
            if final > nl_tol or constitutive > nl_tol:
                raise RuntimeError("forward-newton final equation verification failed")
            stats = dict(nonlinear_solver="forward-newton",
                         nonlinear_converged_final_stage=True,
                         nonlinear_final_relative_residual=final,
                         nonlinear_newton_iters=iteration,
                         nonlinear_convergence_mode="forward-BH-Galerkin",
                         constitutive_projection_relative_residual=constitutive,
                         quadrature_geometry_mass_relative_error=mass_error,
                         nonlinear_residual_history=history, linear_audit=linear_audit,
                         material_magnetization_cap_applied=False)
            return m, iteration, stats
        if iteration == nl_maxit:
            break
        W = ng.BilinearForm(fes)
        integrator = ng.SymbolicBFI(
            ng.InnerProduct(tangent.Reshape((3, 3)) * u, v), simd_evaluate=False)
        for element_type, rule in rules.items():
            integrator.SetIntegrationRule(element_type, rule)
        W += integrator
        W.Assemble()

        class Jacobian(ng.BaseMatrix):
            def IsComplex(self): return False
            def Height(self): return fes.ndof
            def Width(self): return fes.ndof
            def CreateRowVector(self): return G.CreateRowVector()
            def CreateColVector(self): return G.CreateColVector()
            def CreateVector(self, col):
                return G.CreateColVector() if col else G.CreateRowVector()
            def Mult(self, x, y):
                xx = x.FV().NumPy()
                y.FV().NumPy()[:] = apply(G, xx) + demag(apply(Gi, apply(W.mat, xx)))

        J = Jacobian()
        b = G.CreateRowVector()
        b.FV().NumPy()[:] = -residual
        inner = GMRESSolver(J, pre=Gi, tol=min(float(tol), 1e-8),
                            maxiter=maxit, restart=min(100, maxit), printrates=False)
        delta = G.CreateColVector()
        delta.data = inner * b
        step = delta.FV().NumPy().copy()
        true_linear = float(np.linalg.norm(apply(J, step) + residual)
                            / max(np.linalg.norm(residual), 1e-300))
        linear_audit.append(dict(iterations=int(inner.iterations),
                                 true_relative_residual=true_linear))
        if not np.isfinite(true_linear) or true_linear > 1e-6:
            raise RuntimeError(f"forward-newton true linear residual {true_linear:g} exceeds 1e-6")
        length = 1.
        initial_norm = np.linalg.norm(residual)
        while length >= 2. ** -30:
            trial = h + length * step
            trial_residual, trial_m = state(trial)
            if np.linalg.norm(trial_residual) <= (1. - 1e-4 * length) * initial_norm:
                break
            length *= .5
        else:
            raise RuntimeError("forward-newton residual line search exhausted")
        h, residual, m = trial, trial_residual, trial_m
        # Trial material evaluation owns the tangent for the accepted state.
    raise RuntimeError(f"forward-newton nonlinear iteration limit; residual={history[-1]:g}")
