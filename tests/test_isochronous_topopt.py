"""Adjoint-gradient and reciprocity locks for radia.isochronous_topopt.

Stage-1 gate of docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md: the
per-element density adjoint gradient must match central finite differences to
the 1e-6 class, and the dipole-array reciprocity load must reproduce the
independent C++ analytic charge evaluator.  Bands are set from measured values
with 2-3 decades of margin (promoted research run 2026-07-28:
directional 8.1e-10, per-element
worst 3.9e-7 at the FD noise floor, reciprocity 1.1e-10 at bonus_intorder=10).
"""
from math import pi
from types import SimpleNamespace

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia")

from netgen.occ import OCCGeometry, Pnt, Sphere  # noqa: E402
from ngsolve import HDiv, InnerProduct, Mesh, SetNumThreads, TaskManager  # noqa: E402

from radia.isochronous_topopt import (  # noqa: E402
    CHI_MIN, MU0, DensityAdjointVIM, HelmholtzFilter, HeavisideProjection,
    demag_field_from_solution, density_gradient_from_s_gradient,
    density_discreteness, density_to_s, field_functional_load,
    gradient_pair_points, iron_only_verification_ready,
    iron_only_mesh, optimize_density, orbit_arc_points, uniform_field_load,
    verify_design_iron_only,
)


def test_heaviside_projection_chain_and_grayness_gate():
    projection = HeavisideProjection(beta=4.0, eta=0.45)
    rho = np.array([0.05, 0.3, 0.45, 0.7, 0.95])
    gradient = np.array([1.0, -2.0, 0.5, 0.3, -0.7])
    step = 1e-7
    direction = np.array([0.2, -0.1, 0.3, -0.2, 0.1])
    analytic = projection.chain(rho, gradient)@direction
    fd = gradient@(projection.apply(rho+step*direction)
                   - projection.apply(rho-step*direction))/(2*step)
    np.testing.assert_allclose(analytic, fd, rtol=2e-8, atol=2e-10)
    projected = projection.apply(rho)
    assert np.all(np.diff(projected)>0) and np.all((projected>=0)&(projected<=1))
    ready, metrics = iron_only_verification_ready([0., 0.02, .98, 1.])
    assert ready and metrics["intermediate_fraction"] == 0.0
    ready, metrics = iron_only_verification_ready(np.full(20, .5))
    assert not ready and metrics == density_discreteness(np.full(20, .5))


@pytest.fixture(scope="module")
def problem():
    SetNumThreads(4)
    with TaskManager():
        mesh = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.45))
        fes = HDiv(mesh, order=1)
        prob = DensityAdjointVIM(fes, eps=1e-7)
        rng = np.random.default_rng(1234)
        s0 = np.exp(rng.uniform(np.log(1e-2), 0.0, prob.n_el))
        theta = np.linspace(0.0, 2.0 * pi, 6, endpoint=False)
        orbit = np.stack([1.6 * np.cos(theta), 1.6 * np.sin(theta),
                          0.3 * np.ones_like(theta)], axis=1)
        w_orbit = np.array([1.0, -0.7, 1.3, 0.5, -1.1, 0.9])
        pts, wts = gradient_pair_points(orbit, w_orbit, delta=0.1, axis=0)
        f_state = uniform_field_load(fes, (0.0, 0.0, 1.0))
        f_adj = field_functional_load(fes, pts, wts, axis=2, scale=MU0,
                                      bonus_intorder=10)
        base = prob.objective_and_gradient(s0, f_state, f_adj)
    return SimpleNamespace(mesh=mesh, fes=fes, prob=prob, s0=s0,
                           dipole_points=pts, dipole_weights=wts,
                           f_state=f_state, f_adj=f_adj, base=base, rng=rng)


def _objective_at(problem, s):
    with TaskManager():
        gf, _ = problem.prob.solve(s, problem.f_state)
    return float(InnerProduct(problem.f_adj.vec, gf.vec))


def test_adjoint_gradient_directional_fd(problem):
    """All-element check: grad . ds == central FD along a random direction.

    h = 1e-4: the FD signal |2 h grad.ds| must sit well above the CG solution
    noise (tol * J ~ 1e-19); at h = 1e-5 this direction's derivative is
    cancellation-small and the comparison hits the solver noise floor.
    """
    d = np.random.default_rng(7).standard_normal(problem.prob.n_el)
    ds = problem.s0 * d              # multiplicative direction keeps s > 0
    predicted = float(problem.base.gradient @ ds)
    h = 1e-4
    fd = (_objective_at(problem, problem.s0 + h * ds)
          - _objective_at(problem, problem.s0 - h * ds)) / (2.0 * h)
    rel = abs(fd - predicted) / abs(predicted)
    assert rel < 1e-6, (fd, predicted, rel)


def test_adjoint_gradient_per_element_fd(problem):
    """Per-element FD on the strongest-gradient element plus random picks."""
    g = problem.base.gradient
    picks = {int(np.abs(g).argmax())}
    rng = np.random.default_rng(21)
    while len(picks) < 4:
        picks.add(int(rng.integers(0, problem.prob.n_el)))
    for e in sorted(picks):
        h = 1e-4 * problem.s0[e]
        sp = problem.s0.copy(); sp[e] += h
        sm = problem.s0.copy(); sm[e] -= h
        fd = (_objective_at(problem, sp) - _objective_at(problem, sm)) / (2.0 * h)
        rel = abs(fd - g[e]) / max(abs(g[e]), abs(fd))
        assert rel < 1e-5, (e, g[e], fd, rel)


def test_reciprocity_matches_cpp_field_evaluator(problem):
    """f_adj^T m == MU0 * sum_i w_i H_d[m]_z(x_i) via the independent C++ path."""
    with TaskManager():
        H_d = demag_field_from_solution(problem.prob.demag, problem.base.gfM,
                                        problem.dipole_points)
    independent = float(MU0 * np.dot(problem.dipole_weights, H_d[:, 2]))
    rel = abs(problem.base.objective - independent) / abs(independent)
    assert rel < 1e-7, (problem.base.objective, independent, rel)


def test_warm_start_reproduces_and_saves_iterations(problem):
    """Warm-started re-solve at the same design matches and converges fast."""
    with TaskManager():
        again = problem.prob.objective_and_gradient(
            problem.s0, problem.f_state, problem.f_adj, warm=problem.base)
    rel = abs(again.objective - problem.base.objective) / abs(problem.base.objective)
    assert rel < 1e-9
    assert again.state_iterations <= max(3, problem.base.state_iterations // 5)


def test_native_jacobi_solver_satisfies_true_residual(problem):
    """Native convergence is checked against b-Ax, not its CG recurrence."""
    with TaskManager():
        gf, _ = problem.prob.solve(
            problem.s0, problem.f_state, tol=1e-10,
            solver="native-jacobi")
        _, operator, _ = problem.prob._system(problem.s0)
        residual = problem.f_state.vec.CreateVector()
        residual.data = problem.f_state.vec - operator * gf.vec
        relative = ng.Norm(residual) / ng.Norm(problem.f_state.vec)
    assert relative < 1.1e-10


def test_dipole_point_inside_mesh_raises(problem):
    with pytest.raises(ValueError, match="INSIDE"):
        field_functional_load(problem.fes, [[0.0, 0.0, 0.0]], [1.0])


def test_wrong_design_vector_raises(problem):
    with pytest.raises(ValueError, match="entries"):
        problem.prob.solve(np.ones(3), problem.f_state)
    with pytest.raises(ValueError, match="positive"):
        problem.prob.solve(np.zeros(problem.prob.n_el), problem.f_state)


def test_density_mapping_floor_and_chain_rule():
    chi_iron = 100.0
    s = density_to_s(np.array([0.0, 1.0, 0.5]), chi_iron)
    assert s[0] == pytest.approx(1.0 / CHI_MIN)
    assert s[1] == pytest.approx(1.0 / chi_iron)
    rho = np.array([0.1, 0.4, 0.9])
    grad_s = np.array([2.0, -1.0, 0.5])
    h = 1e-7
    fd = (grad_s * (density_to_s(rho + h, chi_iron)
                    - density_to_s(rho - h, chi_iron)) / (2.0 * h))
    chain = density_gradient_from_s_gradient(rho, grad_s, chi_iron)
    assert np.allclose(chain, fd, rtol=1e-5)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        density_to_s(np.array([1.5]), chi_iron)


def test_gradient_pair_points_contract():
    pts, wts = gradient_pair_points([[1.0, 2.0, 3.0]], [4.0], delta=0.2, axis=1)
    assert pts.shape == (2, 3) and wts.shape == (2,)
    np.testing.assert_allclose(pts[0], [1.0, 2.1, 3.0])
    np.testing.assert_allclose(pts[1], [1.0, 1.9, 3.0])
    np.testing.assert_allclose(wts, [20.0, -20.0])


def test_element_volumes_sum_to_mesh_volume(problem):
    with TaskManager():
        direct = float(ng.Integrate(ng.CoefficientFunction(1.0), problem.mesh))
    total = float(problem.prob.element_volumes.sum())
    assert abs(total - direct) / direct < 1e-12
    deficit = (4.0 * pi / 3.0 - total) / (4.0 * pi / 3.0)
    assert 0.0 < deficit < 0.12    # measured 7.0% straight-facet deficit at maxh=0.45


# ---------------------------------------------------------------- Stage 2
def test_helmholtz_filter_transpose_and_invariants(problem):
    with TaskManager():
        filt = HelmholtzFilter(problem.mesh, radius=0.15)
        n_el = problem.prob.n_el
        rng = np.random.default_rng(3)
        rho = rng.uniform(0.2, 0.8, n_el)
        g = rng.standard_normal(n_el)
        # constants are exact fixed points (pure-Neumann Helmholtz)
        np.testing.assert_allclose(filt.apply(np.ones(n_el)), 1.0, atol=1e-9)
        # chain == transpose of apply: FD of g . apply(rho) per component
        chain = filt.chain(g)
        h = 1e-6
        for e in [0, n_el // 2, n_el - 1]:
            rp = rho.copy(); rp[e] += h
            rm = rho.copy(); rm[e] -= h
            fd = (g @ filt.apply(rp) - g @ filt.apply(rm)) / (2.0 * h)
            assert abs(fd - chain[e]) <= 1e-8 * max(1.0, abs(chain[e]))
        # smoothing: a spike loses amplitude; the P1 Helmholtz realization
        # undershoots slightly (measured -1.0e-3 here) -- the design loop
        # clips the filtered density with the piecewise-exact chain rule.
        spike = np.zeros(n_el); spike[n_el // 3] = 1.0
        smoothed = filt.apply(spike)
        assert smoothed.max() < 0.9 and smoothed.min() > -5e-3


def test_native_linearize_matches_ngsolve_reference(problem):
    with TaskManager():
        cpts, _ = orbit_arc_points(1.4, -0.3, 6)
        f_con = field_functional_load(problem.fes, cpts,
                                      np.full(len(cpts), 1.0 / len(cpts)),
                                      axis=2, scale=MU0, bonus_intorder=10)
        s = problem.s0
        lin = problem.prob.linearize(s, problem.f_state,
                                     [problem.f_adj, f_con])
        single = problem.prob.objective_and_gradient(
            s, problem.f_state, f_con, solver="ngsolve-cg")
    assert lin.values[0] == pytest.approx(problem.base.objective, rel=1e-10)
    np.testing.assert_allclose(lin.jacobians[0], problem.base.gradient,
                               rtol=1e-8, atol=1e-18)
    assert lin.values[1] == pytest.approx(single.objective, rel=1e-10)
    np.testing.assert_allclose(lin.jacobians[1], single.gradient,
                               rtol=1e-8, atol=1e-18)


def test_optimize_density_constrained_monotone(problem):
    with TaskManager():
        filt = HelmholtzFilter(problem.mesh, radius=0.12)
        cpts, _ = orbit_arc_points(1.4, -0.3, 8)
        f_con = field_functional_load(problem.fes, cpts,
                                      np.full(len(cpts), 1.0 / len(cpts)),
                                      axis=2, scale=MU0, bonus_intorder=10)
        rho0 = np.full(problem.prob.n_el, 0.5)
        lin0 = problem.prob.linearize(
            density_to_s(filt.apply(rho0), 100.0), problem.f_state,
            [problem.f_adj, f_con])
        target = float(lin0.values[1])
        result = optimize_density(
            problem.prob, problem.f_state, problem.f_adj, [f_con], [target],
            chi_iron=100.0, volume_fraction=0.5, density_filter=filt,
            move_limit=0.1, max_iterations=6)
    hist = result.history
    assert len(hist) >= 3
    objectives = [h["objective"] for h in hist]
    assert objectives[-1] > float(lin0.values[0])          # real ascent
    assert all(b >= a * (1.0 - 1e-6) for a, b in zip(objectives,
                                                     objectives[1:]))
    volume_max = 0.5 * problem.prob.element_volumes.sum()
    for h in hist:
        assert max(np.array(h["violation"]) / np.array(h["band"])) <= 1.25 + 1e-9
        assert h["volume"] <= volume_max * (1.0 + 1e-9)
    assert result.density.min() >= 0.0 and result.density.max() <= 1.0
    assert result.solves >= len(hist) + 1


def test_optimize_density_unconstrained_runs(problem):
    with TaskManager():
        result = optimize_density(
            problem.prob, problem.f_state, problem.f_adj,
            chi_iron=100.0, volume_fraction=0.5, max_iterations=3)
    objectives = [h["objective"] for h in result.history]
    assert len(objectives) == 3
    assert all(b >= a * (1.0 - 1e-6) for a, b in zip(objectives,
                                                     objectives[1:]))


def test_optimize_density_validation_errors(problem):
    with pytest.raises(ValueError, match="targets"):
        optimize_density(problem.prob, problem.f_state, problem.f_adj,
                         [problem.f_adj], [], chi_iron=100.0,
                         volume_fraction=0.5)
    with pytest.raises(ValueError, match="volume_fraction"):
        optimize_density(problem.prob, problem.f_state, problem.f_adj,
                         chi_iron=100.0, volume_fraction=0.0)
    with pytest.raises(ValueError, match="move_min"):
        optimize_density(problem.prob, problem.f_state, problem.f_adj,
                         chi_iron=100.0, volume_fraction=0.5,
                         move_limit=-0.1)
    with pytest.raises(ValueError, match="volume budget"):
        optimize_density(problem.prob, problem.f_state, problem.f_adj,
                         chi_iron=100.0, volume_fraction=0.5,
                         initial_density=np.ones(problem.prob.n_el))


def test_optimize_density_deep_restoration_prioritizes_feasibility():
    """An infeasible functional start must improve even when J must decrease."""
    chi_iron = 100.0

    class LinearProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            chi = 1.0 / np.asarray(s, dtype=float)
            rho = (chi - CHI_MIN) / (chi_iron - CHI_MIN)
            # Objective and constraint both increase with rho[0].  Reaching
            # the target therefore requires a deliberate objective decrease.
            grad_s = np.zeros(2)
            grad_s[0] = -(chi[0] ** 2) / (chi_iron - CHI_MIN)
            return SimpleNamespace(
                values=np.array([rho[0], rho[0]]),
                jacobians=np.stack([grad_s, grad_s]),
                gfM=None, gfLambdas=(None, None), state_iterations=0,
                adjoint_iterations=(0, 0))

    result = optimize_density(
        LinearProblem(), object(), object(), [object()], [0.2],
        chi_iron=chi_iron, volume_fraction=1.0,
        initial_density=np.array([0.9, 0.1]), band_floor=0.01,
        move_limit=0.2, max_iterations=2)
    assert len(result.history) == 2
    assert all(entry["band_mode"] == "deep" for entry in result.history)
    assert result.history[1]["objective"] < result.history[0]["objective"]
    assert result.history[1]["violation"][0] < result.history[0]["violation"][0]


def test_optimize_density_applies_projection_and_reports_discreteness():
    chi_iron = 100.0

    class LinearProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            chi = 1.0/np.asarray(s, dtype=float)
            material_rho = (chi-CHI_MIN)/(chi_iron-CHI_MIN)
            grad_s = np.array([-(chi[0]**2)/(chi_iron-CHI_MIN), 0.0])
            return SimpleNamespace(values=np.array([material_rho[0]]),
                jacobians=grad_s[None,:],gfM=None,gfLambdas=(None,),
                state_iterations=0,adjoint_iterations=(0,))

    checkpoints=[]
    result = optimize_density(LinearProblem(),object(),object(),
        chi_iron=chi_iron,volume_fraction=.5,
        initial_density=np.array([.4,.6]),density_projection=HeavisideProjection(2.),
        move_limit=.1,max_iterations=1,
        checkpoint_callback=lambda entry,rho:checkpoints.append((entry,rho)))
    assert len(result.history)==1 and result.density[0]>.4
    assert "intermediate_fraction" in result.history[0]
    assert len(checkpoints)==1
    np.testing.assert_allclose(checkpoints[0][1],result.density)


def test_optimize_density_rejects_projected_volume_over_budget():
    class LinearProblem:
        n_el = 2
        element_volumes = np.ones(2)

        def linearize(self, s, state_load, loads, **kwargs):
            raise AssertionError("an infeasible projected design must fail before solve")

    with pytest.raises(ValueError, match="projected iron volume budget"):
        optimize_density(
            LinearProblem(), object(), object(), chi_iron=100.0,
            volume_fraction=0.5, initial_density=np.full(2, 0.5),
            density_projection=HeavisideProjection(beta=4.0, eta=0.3),
            max_iterations=1)


# ---------------------------------------------------------------- Stage 3
def test_iron_only_extraction_and_ersatz_band(problem):
    """Hemisphere extraction: exact volume, sane demag physics (surface
    orientation lock), and the embedded-vs-exact-void ersatz sign.

    Bands from measured values (maxh=0.45 ball, 86 kept tets):
    <Mz> extracted 2.284 (a flipped surface orientation runs away to 18+),
    embedded 0/1 <Mz> 2.089, ersatz -8.6 % (charge clamp under-magnetizes
    the embedded iron; -8.3/-8.8 % on the maxh .35/.22 research meshes).
    """
    with TaskManager():
        vols = problem.prob.element_volumes
        cz = np.asarray(ng.Integrate(ng.z, problem.mesh, element_wise=True),
                        float) / vols
        keep = cz > 0.0
        hemi = iron_only_mesh(problem.mesh, keep)
        V_kept = float(vols[keep].sum())
        V_new = float(ng.Integrate(ng.CoefficientFunction(1.0), hemi))
        assert abs(V_new - V_kept) / V_kept < 1e-12
        assert hemi.ne == int(keep.sum())
        prob_h = DensityAdjointVIM(ng.HDiv(hemi, order=1), eps=1e-7)
        gf, _ = prob_h.solve(np.full(prob_h.n_el, 1.0 / 100.0),
                             uniform_field_load(prob_h.fes, (0, 0, 1.0)))
        Mz_extracted = float(ng.Integrate(gf[2], hemi)) / V_new
        s_bin = np.where(keep, 1.0 / 100.0, 1.0 / CHI_MIN)
        gf_b, _ = problem.prob.solve(s_bin, problem.f_state)
        vals = np.asarray(ng.Integrate(gf_b[2], problem.mesh,
                                       element_wise=True), float)
        Mz_embedded = float(vals[keep].sum() / V_kept)
    assert 2.2 < Mz_extracted < 2.4, Mz_extracted
    ersatz = (Mz_embedded - Mz_extracted) / Mz_extracted
    assert -0.15 < ersatz < -0.03, ersatz


def test_verify_design_iron_only_protocol(problem):
    """Matched-0/1 verification bands from the promoted protocol function.

    Measured on this mesh: band +4.4 % for an external-field functional
    (values -1.77e-7 embedded vs -1.85e-7 iron-only)."""
    with TaskManager():
        vols = problem.prob.element_volumes
        cz = np.asarray(ng.Integrate(ng.z, problem.mesh, element_wise=True),
                        float) / vols
        density = (cz > 0.0).astype(float)
        cpts, _ = orbit_arc_points(1.5, 0.4, 6)

        def state_builder(fes):
            return uniform_field_load(fes, (0.0, 0.0, 1.0))

        def functional_builder(fes):
            return field_functional_load(fes, cpts, np.full(6, 1.0 / 6.0),
                                         axis=2, scale=MU0, bonus_intorder=10)

        ver = verify_design_iron_only(
            problem.prob, density, state_builder, [functional_builder],
            chi_iron=100.0, gram_kwargs=dict(eps=1e-7))
    assert ver.keep.sum() == int((cz > 0.0).sum())
    assert ver.values_embedded.shape == ver.values_iron_only.shape == (1,)
    assert ver.values_embedded[0] < 0.0 and ver.values_iron_only[0] < 0.0
    assert 0.0 < ver.bands[0] < 0.12, ver.bands
    assert ver.embedded_iterations > 0 and ver.iron_only_iterations > 0


def test_iron_only_mesh_guards(problem):
    n_el = problem.prob.n_el
    with pytest.raises(ValueError, match="proper non-empty subset"):
        iron_only_mesh(problem.mesh, np.ones(n_el, bool))
    with pytest.raises(ValueError, match="mask has"):
        iron_only_mesh(problem.mesh, np.ones(3, bool))
    with pytest.raises(ValueError, match="shape"):
        verify_design_iron_only(problem.prob, np.ones(3),
                                lambda fes: None, [], chi_iron=100.0)


def test_density_penalty_mapping():
    chi_iron = 1000.0
    rho = np.array([0.0, 0.5, 1.0])
    s1 = density_to_s(rho, chi_iron, penalty=1.0)
    s3 = density_to_s(rho, chi_iron, penalty=3.0)
    # endpoints penalty-invariant; interior penalized (chi smaller -> s larger)
    assert s3[0] == s1[0] and s3[2] == s1[2]
    assert s3[1] > s1[1]
    # chain rule vs FD at penalty=3
    rho = np.array([0.2, 0.5, 0.9])
    grad_s = np.array([1.0, -2.0, 0.5])
    h = 1e-7
    fd = grad_s * (density_to_s(rho + h, chi_iron, penalty=3.0)
                   - density_to_s(rho - h, chi_iron, penalty=3.0)) / (2 * h)
    chain = density_gradient_from_s_gradient(rho, grad_s, chi_iron,
                                             penalty=3.0)
    np.testing.assert_allclose(chain, fd, rtol=1e-5)
    with pytest.raises(ValueError, match="penalty"):
        density_to_s(rho, chi_iron, penalty=0.5)


def test_orbit_arc_and_direction_pairs():
    pts, radial = orbit_arc_points(2.0, 0.5, 4)
    assert pts.shape == (4, 3) and radial.shape == (4, 3)
    np.testing.assert_allclose(np.linalg.norm(radial, axis=1), 1.0)
    np.testing.assert_allclose(pts[0], [2.0, 0.0, 0.5], atol=1e-14)
    # full circle omits the duplicate endpoint
    assert not np.allclose(pts[-1], pts[0])
    ppts, wts = gradient_pair_points(pts, np.ones(4), delta=0.2,
                                     direction=radial)
    assert ppts.shape == (8, 3)
    np.testing.assert_allclose(ppts[0], [2.1, 0.0, 0.5], atol=1e-14)
    np.testing.assert_allclose(ppts[4], [1.9, 0.0, 0.5], atol=1e-14)
    np.testing.assert_allclose(wts[:4], 5.0)
    np.testing.assert_allclose(wts[4:], -5.0)
    with pytest.raises(ValueError, match="zero direction"):
        gradient_pair_points(pts, np.ones(4), delta=0.2,
                             direction=np.zeros((4, 3)))
