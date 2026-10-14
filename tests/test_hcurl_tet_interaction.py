"""Regression tests for the epsilon-free affine-tet HCurl interaction."""

from __future__ import annotations

import numpy as np
import pytest

import radia._radia_pybind as _rp
import radia.vim as vim
from radia.vim._hcurl_tet_interaction import _affine_polynomial_hmatrix_operator
from radia.vim._vim import (
    _f64_buffer,
    _g01,
    _i32_buffer,
    _monos_vol,
    _outer_tet,
    _outer_tri,
)


def _p2_tet_nodes(vertices):
    vertices = np.asarray(vertices, dtype=float)
    return np.vstack(
        (
            vertices,
            0.5 * (vertices[0] + vertices[1]),
            0.5 * (vertices[1] + vertices[2]),
            0.5 * (vertices[2] + vertices[0]),
            0.5 * (vertices[0] + vertices[3]),
            0.5 * (vertices[1] + vertices[3]),
            0.5 * (vertices[2] + vertices[3]),
        )
    )


def test_reduced_reference_polynomial_gram_matches_charge_gram_oracle():
    vertices = np.array(
        [[[0.0, 0.0, 0.0], [0.7, 0.1, 0.0], [0.1, 0.8, 0.05], [0.05, 0.2, 0.9]]]
    )
    exponents = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
        dtype=np.int32,
    )
    coefficients = np.zeros((2, 1, 4, 3))
    coefficients[0, 0, :, 0] = [1.0, 0.2, -0.1, 0.3]
    coefficients[1, 0, :, 1] = [0.7, -0.4, 0.2, 0.1]
    coefficients[1, 0, :, 2] = [0.2, 0.1, 0.0, -0.2]
    tet_points, tet_weights = _outer_tet(4)
    tri_points, tri_weights = _outer_tri(4)

    reduced = _rp._TetHCurlReducedGram(
        _f64_buffer(vertices),
        _i32_buffer(exponents),
        _f64_buffer(coefficients),
        2,
        _f64_buffer(tet_points),
        _f64_buffer(tet_weights),
    )
    oracle = _rp._ChargeGramHMatrix(
        cell_verts=_f64_buffer(vertices),
        face_verts=np.empty(0),
        n_el=1,
        charge_host=_i32_buffer(np.zeros(4, dtype=int)),
        charge_kind=_i32_buffer(np.zeros(4, dtype=int)),
        charge_expo=_i32_buffer(exponents),
        ref_tet_pts=_f64_buffer(tet_points),
        ref_tet_w=_f64_buffer(tet_weights),
        ref_tri_pts=_f64_buffer(tri_points),
        ref_tri_w=_f64_buffer(tri_weights),
        build=False,
    )
    scalar = np.array([[oracle.entry(i, j) for j in range(4)] for i in range(4)])
    expected = np.zeros((2, 2))
    for component in range(3):
        values = coefficients[:, 0, :, component]
        expected += values @ scalar @ values.T

    assert np.allclose(reduced, expected, rtol=2.0e-14, atol=1.0e-18)


def test_native_component_charge_map_matches_explicit_btgb():
    vertices = np.array(
        [[[0.0, 0.0, 0.0], [0.7, 0.1, 0.0], [0.1, 0.8, 0.05], [0.05, 0.2, 0.9]]]
    )
    exponents = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.int32
    )
    polynomial_coefficients = np.eye(4)
    tet_points, tet_weights = _outer_tet(4)
    gram = _rp._ChargeGramHMatrix.from_local_polynomials(
        cell_verts=_f64_buffer(vertices),
        n_el=1,
        charge_host=_i32_buffer(np.zeros(4, dtype=int)),
        polynomial_coefficients=_f64_buffer(polynomial_coefficients),
        polynomial_exponents=_i32_buffer(exponents),
        ref_tet_pts=_f64_buffer(tet_points),
        ref_tet_w=_f64_buffer(tet_weights),
        eps=1.0e-12,
        leaf=2,
        eta=2.0,
        build=True,
    )
    maps = np.random.default_rng(20260717).normal(size=(3, 4, 3))
    operator = vim.HCurlHMatrixOperator(gram, maps, mu=1.7)
    scalar = np.array([[gram.entry(i, j) for j in range(4)] for i in range(4)])
    expected = 1.7 * sum(charge_map.T @ scalar @ charge_map for charge_map in maps)

    np.testing.assert_allclose(operator.to_dense(), expected, rtol=2.0e-13, atol=1.0e-18)
    np.testing.assert_allclose(
        operator.matvec(np.array([1.0, -0.5, 0.25])),
        expected @ np.array([1.0, -0.5, 0.25]),
        rtol=2.0e-13,
        atol=1.0e-18,
    )
    assert operator.stats()["charge_map_backend"] == "native-cpp-component-csr"


def test_degree18_multitet_operator_matches_dense_reduced_gram():
    vertices = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[1.2, 0.0, 0.0], [2.2, 0.0, 0.0], [1.2, 1.0, 0.0], [1.2, 0.0, 1.0]],
        ]
    )
    exponents = np.asarray(_monos_vol(18), dtype=np.int32)
    coefficients = np.random.default_rng(20260718).normal(
        size=(4, 2, len(exponents), 3)
    )
    coefficients *= np.power(10.0, 0.35 * exponents.sum(axis=1))[None, None, :, None]
    projected = {
        "cell_verts": vertices,
        "cell_count": 2,
        "subtet_count": 2,
        "coefficients": coefficients,
        "exponents": exponents,
    }
    points, weights = _outer_tet(4)
    operator = _affine_polynomial_hmatrix_operator(
        projected,
        outer_points=points,
        outer_weights=weights,
        eps=1.0e-12,
        leafsize=64,
        eta=2.0,
        local_rank_rtol=1.0e-13,
        max_charges=1_000,
        mu=1.0,
    )
    dense = np.asarray(
        _rp._TetHCurlReducedGram(
            _f64_buffer(vertices),
            _i32_buffer(exponents),
            _f64_buffer(coefficients),
            4,
            _f64_buffer(points),
            _f64_buffer(weights),
        )
    )

    np.testing.assert_allclose(operator.to_dense(), dense, rtol=2.0e-12, atol=1.0e-12)
    assert operator.stats()["local_rank_min"] == 12
    assert operator.stats()["local_rank_relative_truncation_error"] == 0.0


def test_curved_gram_reference_density_omits_both_physical_measures():
    vertices = np.asarray(
        ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 4.0))
    )
    nodes = _p2_tet_nodes(vertices)
    tet_points, tet_weights = _outer_tet(4)
    tri_points, tri_weights = _outer_tri(4)
    curve_points, curve_weights = _g01(5)
    common = dict(
        cell_nodes=_f64_buffer(nodes),
        face_nodes=np.empty(0),
        cell_vertices=_i32_buffer([0, 1, 2, 3]),
        face_vertices=np.empty(0, dtype=np.int32),
        n_el=1,
        curve_order=2,
        charge_host=_i32_buffer([0]),
        charge_kind=_i32_buffer([0]),
        charge_expo=_i32_buffer([[0, 0, 0]]),
        ref_tet_pts=_f64_buffer(tet_points),
        ref_tet_w=_f64_buffer(tet_weights),
        ref_tri_pts=_f64_buffer(tri_points),
        ref_tri_w=_f64_buffer(tri_weights),
        curve_gl=_f64_buffer(curve_points),
        curve_gw=_f64_buffer(curve_weights),
        build=False,
    )
    physical = _rp._ChargeGramHMatrix(
        **common, reference_density=False
    ).entry(0, 0)
    reference = _rp._ChargeGramHMatrix(
        **common, reference_density=True
    ).entry(0, 0)

    determinant = 2.0 * 3.0 * 4.0
    assert physical > 0.0
    assert reference > 0.0
    assert physical / reference == pytest.approx(determinant**2, rel=2.0e-14)


def _small_hcurl_problem():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")
    box = occ.Box((0.0, 0.0, 0.0), (1.0, 0.8, 0.6))
    box.mat("cond")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.8))
    fes = ng.HCurl(mesh, order=2, nograds=True)
    vectors = np.random.default_rng(20260717).normal(size=(fes.ndof, 2))
    basis = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        vectors,
        intorder=3,
        materials="cond",
    )
    return ng, mesh, fes, vectors, basis


def test_ngsolve_hcurl_tet_interaction_is_psd_and_epsilon_free():
    ng, mesh, fes, vectors, basis = _small_hcurl_problem()
    with ng.TaskManager():
        interaction = vim.NgsolveHCurlTetVolumeInteraction(
            mesh,
            fes,
            vectors,
            basis,
            materials="cond",
            outer_quad=4,
        )
        dense = vim.NgsolveHCurlTetVolumeInteraction(
            mesh,
            fes,
            vectors,
            basis,
            materials="cond",
            outer_quad=4,
            matrix_free=False,
        )

    diagnostics = interaction.diagnostics()
    assert interaction.matrix.shape == (2, 2)
    assert np.allclose(interaction.matrix, interaction.matrix.T, atol=1.0e-18)
    assert np.linalg.eigvalsh(interaction.matrix)[0] > 0.0
    assert diagnostics["projection_relative_residual"] < 1.0e-12
    assert diagnostics["kernel_epsilon_m"] is None
    assert diagnostics["singular_self_treatment"].startswith("analytic-reference")
    assert interaction(basis, basis) is interaction.matrix
    assert diagnostics["matrix_free"] is True
    assert diagnostics["hmatrix_operator"]["charge_map_backend"] == "native-cpp-component-csr"
    np.testing.assert_allclose(interaction.matrix, dense.matrix, rtol=2.0e-8, atol=1.0e-18)


def test_ngsolve_hcurl_tet_interaction_rejects_underfit_degree():
    ng, mesh, fes, vectors, basis = _small_hcurl_problem()
    with ng.TaskManager(), pytest.raises(ValueError, match="projection residual exceeds"):
        vim.NgsolveHCurlTetVolumeInteraction(
            mesh,
            fes,
            vectors,
            basis,
            degree=0,
            projection_tolerance=1.0e-12,
            materials="cond",
            outer_quad=4,
        )


def test_p2_curved_hcurl_uses_reference_density_and_exact_geometry():
    ng = pytest.importorskip("ngsolve")
    csg = pytest.importorskip("netgen.csg")
    geometry = csg.CSGeometry()
    geometry.Add(csg.Sphere(csg.Pnt(0.0, 0.0, 0.0), 1.0).mat("cond"))
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=2.0))
    mesh.Curve(2)
    fes = ng.HCurl(mesh, order=2, nograds=True)
    vectors = np.random.default_rng(20260717).normal(size=fes.ndof)
    with ng.TaskManager():
        basis = vim.NgsolveHCurlCurlBasis(
            mesh,
            fes,
            vectors,
            intorder=4,
            materials="cond",
        )
        interaction = vim.NgsolveHCurlCellVolumeInteraction(
            mesh,
            fes,
            vectors,
            basis,
            materials="cond",
            projection_quad=4,
            outer_quad=3,
            curve_gauss=4,
            far_quad=2,
            hmatrix_eps=1.0e-7,
            projection_tolerance=1.0e-11,
            geometry_tolerance=1.0e-12,
        )

    diagnostics = interaction.diagnostics()
    assert interaction.matrix[0, 0] > 0.0
    assert diagnostics["backend"] == "curved-p2-reference-density-hcurl-hmatrix-operator"
    assert diagnostics["geometry"] == "curved-p2-reference-density"
    assert diagnostics["polynomial_degree"] == 2
    assert diagnostics["charge_count"] == 10 * mesh.GetNE(ng.VOL)
    assert diagnostics["projection_relative_residual"] < 1.0e-12
    assert diagnostics["geometry_relative_residual"] < 1.0e-12
    assert diagnostics["kernel_epsilon_m"] is None
    assert diagnostics["matrix_free"] is True


def test_public_hcurl_tet_interaction_names_are_exported():
    assert vim.HCurlTetVolumeInteraction.__name__ == "HCurlTetVolumeInteraction"
    assert callable(vim.NgsolveHCurlTetVolumeInteraction)


def test_custom_interaction_does_not_allocate_sampled_epsilon(monkeypatch):
    from radia.vim import _eddy_hybrid

    basis = vim.VolumeCurrentBasis(
        [[0.0, 0.0, 0.0]],
        [1.0],
        [[[1.0, 0.0, 0.0]]],
    )
    monkeypatch.setattr(
        _eddy_hybrid,
        "_default_kernel_epsilon",
        lambda _bases: (_ for _ in ()).throw(AssertionError("sampled epsilon requested")),
    )
    system = vim.AssembleHybridVIM(
        basis,
        sigma=1.0,
        interaction=lambda _left, _right: np.array([[2.0]]),
    )
    assert system.inductance[0, 0] == pytest.approx(2.0)
