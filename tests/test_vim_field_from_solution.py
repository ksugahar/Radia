"""Golden tests for the C++ radia.vim.FieldFromSolution BDM1 field path.

Locks:
  * C++ materialization directly from the configured charge map and HDiv coefficients;
  * a uniform-M box against Radia's independent analytic cuboid field;
  * the mu_r sphere end-to-end: far field matches the mesh-moment dipole and the
    constant-M-collapse evaluation of the same solution;
  * fail-loud contract (res without gfM / wrong order).
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
rad = pytest.importorskip("radia")

from radia import vim
from radia.vim import _field_batch as fb


def test_custom_tree_options_are_validated_and_cached_by_configuration():
    class FakeGram:
        def __init__(self):
            self.calls = []

        def create_field_evaluator(self, *args):
            self.calls.append(args)
            return object()

    gram = FakeGram()
    result = {
        "gfM": object(),
        "order": 1,
        "_charge_gram": gram,
        "_m_coefficients": np.array([1.0]),
    }
    first = fb._materialize_field_evaluator(
        result, {"theta": 0.2, "leaf_size": 16})
    second = fb._materialize_field_evaluator(
        result, {"leaf_size": 16.0, "theta": 0.2})

    assert first is second
    assert len(gram.calls) == 1
    assert all(isinstance(key, str) for key in result)
    assert isinstance(result["_field_evaluator_custom"], dict)

    with pytest.raises(TypeError, match="must be a dict"):
        fb._materialize_field_evaluator(result, [("theta", 0.2)])
    with pytest.raises(ValueError, match="unknown"):
        fb._materialize_field_evaluator(result, {"opening_angle": 0.2})
    with pytest.raises(ValueError, match="leaf_size"):
        fb._materialize_field_evaluator(result, {"leaf_size": 1.5})
    with pytest.raises(ValueError, match="theta"):
        fb._materialize_field_evaluator(result, {"theta": np.nan})
    with pytest.raises(ValueError, match="probe_count"):
        fb._materialize_field_evaluator(result, {"probe_count": 0})


def test_uniform_box_matches_radia_cpp():
    from netgen.occ import Box, OCCGeometry, Pnt
    rad.UtiDelAll()
    dims = np.array([0.04, 0.05, 0.03])
    center = 0.5 * dims
    M0 = np.array([2.0e5, -1.0e5, 3.0e5])
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(*dims))).GenerateMesh(maxh=0.015))
        gfM = ng.GridFunction(ng.HDiv(mesh, order=1))
        gfM.Set(ng.CoefficientFunction(tuple(M0)))
        _B, gram, _mass = vim.ChargeGram(gfM.space, eps=1e-12)
    box = rad.ObjRecMag(center.tolist(), dims.tolist(), M0.tolist())
    pts = np.array([[0.10, 0.03, 0.02], [-0.05, 0.02, 0.04],
                    [0.02, 0.12, -0.04], [0.08, -0.04, 0.09]])
    H_rad = np.asarray(rad.Fld(box, "h", pts))
    solution = {
        "gfM": gfM,
        "order": 1,
        "curve_order": None,
        "_charge_gram": gram,
        "_m_coefficients": np.ascontiguousarray(gfM.vec.FV().NumPy()),
    }
    H_rt = vim.FieldFromSolution(solution, pts)
    rel = (np.linalg.norm(H_rt - H_rad, axis=1)
           / np.maximum(np.linalg.norm(H_rad, axis=1), 1e-30))
    A_rad = np.asarray(rad.Fld(box, "a", pts))
    with ng.TaskManager():
        A_coefficient = vim.VectorPotentialCoefficientFromSolution(
            solution, construction="quadrature", integration_order=8
        )
        adaptive_targets = np.vstack(
            (pts, np.array([[0.02, 0.025, 0.031]]))
        )
        adaptive_A = vim.VectorPotentialCoefficientFromSolution(
            solution,
            construction="quadrature",
            integration_order=4,
            target_points_m=adaptive_targets,
            maximum_integration_order=8,
        )
        exact_A = vim.VectorPotentialCoefficientFromSolution(solution)
        target = ng.Mesh(
            OCCGeometry(
                Box(Pnt(-0.12, -0.12, -0.12), Pnt(0.22, 0.22, 0.22))
            ).GenerateMesh(maxh=0.08)
        )
        A_hdiv = np.asarray(
            [A_coefficient(target(*point)) for point in pts], dtype=float
        )
        A_adaptive = np.asarray(
            [adaptive_A(target(*point)) for point in pts], dtype=float
        )
        A_exact = np.asarray(
            [exact_A(target(*point)) for point in pts], dtype=float
        )
        H_coefficient = vim.FieldCoefficientFromSolution(solution)
        H_symmetric = vim.FieldCoefficientFromSolution(
            solution, reflection_normal=[0.0, 0.0, 1.0]
        )
        reflected_pts = pts * np.array([1.0, 1.0, -1.0])
        H_raw = np.asarray(
            [H_coefficient(target(*point)) for point in pts], dtype=float
        )
        H_reflected = np.asarray(
            [H_coefficient(target(*point)) for point in reflected_pts],
            dtype=float,
        )
        H_symmetrized = np.asarray(
            [H_symmetric(target(*point)) for point in pts], dtype=float
        )
    A_scale = np.maximum(np.linalg.norm(A_rad, axis=1), 1e-30)
    A_relative_error = np.linalg.norm(A_hdiv - A_rad, axis=1) / A_scale
    adaptive_relative_error = (
        np.linalg.norm(A_adaptive - A_rad, axis=1) / A_scale
    )
    exact_relative_error = np.linalg.norm(A_exact - A_rad, axis=1) / A_scale
    A_stats = solution["vector_potential_coefficient_stats"][8]
    exact_stats = solution["vector_potential_coefficient_stats"][
        ("exact-equivalent-current", None)
    ]
    adaptive_stats = next(
        value
        for value in solution["vector_potential_coefficient_stats"].values()
        if value.get("target_adaptive")
    )
    mixed = np.array([[0.01, 0.02, 0.01], [0.20, 0.20, 0.20],
                      [0.03, 0.01, 0.02]])
    got_m = fb.magnetization_from_solution(
        {"gfM": gfM, "order": 1, "curve_order": None}, mixed)
    expected_m = np.array([M0, [0.0, 0.0, 0.0], M0])
    rad.UtiDelAll()
    assert rel.max() < 2e-8
    assert A_relative_error.max() < 2e-7
    assert adaptive_relative_error.max() < 3e-5
    expected_H_symmetric = 0.5 * (
        H_raw + H_reflected * np.array([-1.0, -1.0, 1.0])
    )
    np.testing.assert_allclose(
        H_symmetrized, expected_H_symmetric, rtol=2e-14, atol=2e-10
    )
    assert H_symmetric.reflection_symmetrized is True
    np.testing.assert_allclose(H_symmetric.reflection_normal, [0.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="reflection_normal must be nonzero"):
        vim.FieldCoefficientFromSolution(
            solution, reflection_normal=[0.0, 0.0, 0.0]
        )
    assert A_stats["source_count"] == A_coefficient.source_count
    assert A_stats["construction"] == (
        "ngsolve-mapped-hdiv-magnetization-volume-integral"
    )
    # The exact equivalent-current construction reproduces the same uniform-M
    # body as Radia's independent analytic cuboid A, to machine precision --
    # the point-dipole cloud above only reaches 1e-7 relative on the same mesh.
    assert exact_relative_error.max() < 1e-12
    assert exact_stats["construction"] == (
        "analytic-equivalent-current-tet-triangle"
    )
    assert exact_stats["requested_construction"] == "auto"
    assert exact_stats["tetrahedron_count"] == mesh.ne
    assert exact_stats["triangle_count"] == 4 * mesh.ne
    assert exact_stats["maximum_relative_curvature_defect"] < 1e-12
    assert exact_stats["maximum_relative_affine_defect"] < 1e-12
    assert exact_A.tetrahedron_count == mesh.ne
    assert exact_A.triangle_count == 4 * mesh.ne
    assert adaptive_stats["target_adaptive"] is True
    assert adaptive_stats["target_point_count"] == len(adaptive_targets)
    assert adaptive_stats["maximum_integration_order"] == 8
    assert sum(
        adaptive_stats["element_integration_order_histogram"].values()
    ) == mesh.ne
    assert max(adaptive_stats["element_integration_order_histogram"]) == 8
    assert np.allclose(got_m, expected_m, rtol=2e-14, atol=2e-10)


def test_exact_vector_potential_equivalent_current():
    """Lock the analytic equivalent-current A against independent oracles.

    A point-dipole quadrature cloud cannot resolve the 1/R^2 kernel when the
    target stands much closer to a source element than that element's own size,
    which is exactly the EarlyTimes beam-tube geometry.  The equivalent-current
    identity closes the integral analytically, so ``curl A`` must reproduce the
    exact charge-kernel H right next to the body, not only in the far field.
    """
    from netgen.occ import Box, OCCGeometry, Pnt
    rad.UtiDelAll()
    mu0 = 4.0e-7 * np.pi
    lower, upper = np.array([0.0, 0.0, 0.0]), np.array([0.04, 0.03, 0.02])
    with ng.TaskManager():
        mesh = ng.Mesh(
            OCCGeometry(Box(Pnt(*lower), Pnt(*upper))).GenerateMesh(maxh=0.012)
        )
        gfM = ng.GridFunction(ng.HDiv(mesh, order=1))
        # Affine M with nonzero curl and nonzero divergence: both the volume
        # equivalent current and the affine face sheet current are exercised.
        affine = ng.CoefficientFunction(
            (2.0e5 + 3.0e6 * ng.x - 1.0e6 * ng.y + 5.0e5 * ng.z,
             -4.0e5 + 8.0e5 * ng.x + 2.0e6 * ng.y - 3.0e6 * ng.z,
             1.0e5 - 2.0e6 * ng.x + 4.0e5 * ng.y + 1.5e6 * ng.z))
        gfM.Set(affine)
        _B, gram, _mass = vim.ChargeGram(gfM.space, eps=1e-12)
    solution = {
        "gfM": gfM,
        "order": 1,
        "curve_order": None,
        "_charge_gram": gram,
        "_m_coefficients": np.ascontiguousarray(gfM.vec.FV().NumPy()),
    }
    with ng.TaskManager():
        exact = vim.VectorPotentialCoefficientFromSolution(solution)
        target = ng.Mesh(
            OCCGeometry(
                Box(Pnt(-0.30, -0.30, -0.30), Pnt(0.30, 0.30, 0.30))
            ).GenerateMesh(maxh=0.10)
        )
        # Far field: an independent mapped-quadrature dipole sum converges here.
        far = np.array([[0.25, 0.20, 0.15], [0.10, -0.12, 0.18]])
        coordinate = ng.CoefficientFunction((ng.x, ng.y, ng.z))
        positions, moments = [], []
        for element in mesh.Elements(ng.VOL):
            rule = ng.IntegrationRule(element.type, 10)
            transformation = mesh.GetTrafo(element)
            mapped = transformation(rule)
            count = len(rule)
            position = fb._mapped_vector_values(coordinate, mapped, count, "points")
            value = fb._mapped_vector_values(gfM, mapped, count, "M")
            weights = np.asarray([float(point.weight)
                                  * float(transformation(point).measure)
                                  for point in rule])
            positions.append(position)
            moments.append(value * weights[:, None])
        positions, moments = np.vstack(positions), np.vstack(moments)
        reference_far = np.array([
            1.0e-7 * np.sum(np.cross(moments, point - positions)
                            / np.linalg.norm(point - positions, axis=1)[:, None] ** 3,
                            axis=0)
            for point in far])
        exact_far = np.asarray([exact(target(*point)) for point in far], float)

        # Near field: 1 mm outside a body meshed at 12 mm, i.e. the regime where
        # the point-dipole construction fails.  curl A must equal mu0 * H.
        near = np.array([[0.020, 0.015, 0.021], [0.041, 0.015, 0.010],
                         [0.020, -0.001, 0.010]])
        step = 5.0e-5
        curl_a = np.empty((len(near), 3))
        for index, point in enumerate(near):
            gradient = np.empty((3, 3))
            for axis in range(3):
                offset = np.zeros(3)
                offset[axis] = step
                plus = np.asarray(exact(target(*(point + offset))), float)
                minus = np.asarray(exact(target(*(point - offset))), float)
                gradient[axis] = (plus - minus) / (2.0 * step)
            curl_a[index] = (gradient[1][2] - gradient[2][1],
                             gradient[2][0] - gradient[0][2],
                             gradient[0][1] - gradient[1][0])
        near_b = mu0 * np.asarray(vim.FieldFromSolution(solution, near), float)

        symmetric = vim.VectorPotentialCoefficientFromSolution(
            solution, reflection_normal=[0.0, 0.0, 1.0])
        mirrored = far * np.array([1.0, 1.0, -1.0])
        symmetrized = np.asarray([symmetric(target(*p)) for p in far], float)
        reflected = np.asarray([exact(target(*p)) for p in mirrored], float)
    rad.UtiDelAll()

    far_scale = np.linalg.norm(reference_far, axis=1)
    assert (np.linalg.norm(exact_far - reference_far, axis=1) / far_scale).max() < 1e-9
    near_scale = np.linalg.norm(near_b, axis=1)
    assert near_scale.min() > 1e-3                      # a real near field, not noise
    # Limited by the central-difference truncation, not by the construction.
    assert (np.linalg.norm(curl_a - near_b, axis=1) / near_scale).max() < 1e-5
    # A is polar where B is axial: A_sym(r) = 0.5 [A(r) + R A(R r)].
    expected = 0.5 * (exact_far + reflected * np.array([1.0, 1.0, -1.0]))
    np.testing.assert_allclose(symmetrized, expected, rtol=1e-9, atol=1e-16)


def test_exact_vector_potential_fails_loud_on_non_affine_space():
    """A quadratic (BDM2) magnetization is not an affine source -- raise."""
    from netgen.occ import Box, OCCGeometry, Pnt
    rad.UtiDelAll()
    with ng.TaskManager():
        mesh = ng.Mesh(
            OCCGeometry(Box(Pnt(0, 0, 0), Pnt(0.04, 0.03, 0.02))).GenerateMesh(
                maxh=0.02)
        )
        gfM = ng.GridFunction(ng.HDiv(mesh, order=2))
        gfM.Set(ng.CoefficientFunction(
            (1.0e6 * ng.x * ng.x, 2.0e6 * ng.y * ng.z, 3.0e5 * ng.z * ng.z)))
    solution = {"gfM": gfM, "order": 2, "curve_order": None}
    with pytest.raises(NotImplementedError, match="affine"):
        vim.VectorPotentialCoefficientFromSolution(solution, construction="exact")
    with pytest.raises(ValueError, match="construction must be"):
        vim.VectorPotentialCoefficientFromSolution(solution, construction="analytic")
    rad.UtiDelAll()


def test_sphere_end_to_end_and_fail_loud():
    from netgen.occ import OCCGeometry, Pnt, Sphere
    rng = np.random.default_rng(3)
    rad.UtiDelAll()
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 0.05)).GenerateMesh(maxh=0.02))
    iron = vim.MeshSoftIron(mesh, mu_r=1000.0)
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, 4e-7 * np.pi * 1000.0])
    top = rad.ObjCnt([iron, bkg])
    res = rad.Solve(top)
    assert "gfM" in res
    assert "_field_evaluator" in res
    evaluator = res["_field_evaluator"]
    assert res["field_evaluator_stats"]["source_kind"] == "analytic-tet-bdm1"
    with ng.TaskManager():
        V_el = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh, ng.VOL,
                                       element_wise=True), float)
    m_dip = float((np.asarray(res["M"]) * V_el[:, None]).sum(axis=0)[2])
    u = rng.normal(size=(20, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    far = 0.15 * u
    r = np.linalg.norm(far, axis=1)
    rh = far / r[:, None]
    mz = np.array([0.0, 0.0, m_dip])
    H_an = (3.0 * (rh @ mz)[:, None] * rh - mz[None, :]) / (4.0 * np.pi * r ** 3)[:, None]
    H_lin = vim.FieldFromSolution(res, far)
    assert res["_field_evaluator"] is evaluator
    H_col = np.asarray(rad.Fld(iron, "h", far))
    scale = np.linalg.norm(H_an, axis=1).mean()
    assert np.linalg.norm(H_lin - H_an, axis=1).max() / scale < 2e-2   # facet multipole tail
    assert np.linalg.norm(H_lin - H_col, axis=1).max() / scale < 2e-3  # same solution
    with pytest.raises(ValueError):
        vim.FieldFromSolution({"order": 1}, far)                       # no gfM
    bad = dict(res)
    bad["order"] = 0
    with pytest.raises(NotImplementedError):
        vim.FieldFromSolution(bad, far)
    rad.UtiDelAll()
