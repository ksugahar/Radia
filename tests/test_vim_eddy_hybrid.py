import json

import numpy as np
import pytest

import radia.vim as vim


class _ToySymmetricOperator:
    def __init__(self, matrix):
        self._matrix = np.asarray(matrix, dtype=float)
        self.shape = self._matrix.shape
        self.dtype = self._matrix.dtype
        self.is_hermitian = True

    def matvec(self, vector):
        return self._matrix @ vector

    def matmat(self, matrix):
        return self._matrix @ matrix

    def to_dense(self):
        return self._matrix.copy()

    def stats(self):
        return {"matrix_free": True, "backend": "toy"}


def test_team28_skin_depth_gate_selects_volumetric_hcurl():
    gate = vim.EddySIBCApplicability(
        frequency_hz=50.0,
        sigma=3.4e7,
        characteristic_thickness_m=3.0e-3,
    )

    assert vim.SkinDepth(50.0, 3.4e7) == pytest.approx(12.206e-3, rel=2.0e-3)
    assert gate.skin_depth_m == pytest.approx(vim.SkinDepth(50.0, 3.4e7))
    assert gate.thickness_to_skin_depth == pytest.approx(0.2458, rel=2.0e-3)
    assert gate.sibc_applicable is False
    assert gate.selected_model == "volumetric"
    assert gate.diagnostics()["thickness_separated"] is False

    high_frequency = vim.EddySIBCApplicability(
        frequency_hz=1.0e6,
        sigma=3.4e7,
        characteristic_thickness_m=3.0e-3,
    )
    assert high_frequency.thickness_to_skin_depth > 3.0
    assert high_frequency.sibc_applicable is True
    assert high_frequency.selected_model == "sibc"


def test_empty_surface_basis_keeps_a_volume_only_vim_well_formed():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([[[1.0, 0.0, 0.0]]]),
        names=["bulk"],
    )
    surface = vim.SampledCurrentBasis(
        points=np.zeros((0, 3)),
        weights=np.zeros(0),
        modes=np.zeros((0, 0, 3)),
        kind="surface",
        names=(),
    )

    system = vim.AssembleHybridVIM(
        volume,
        surface,
        sigma=3.4e7,
        kernel_epsilon=0.1,
    )
    rhs = system.block_rhs(volume=np.array([1.0]), surface=np.zeros(0))

    assert system.n_modes == 1
    assert system.blocks["surface"] == (1, 1)
    assert system.block_slice("surface") == slice(1, 1)
    assert system.solve(1j * 2.0 * np.pi * 50.0, rhs).shape == (1,)
    assert system.diagnostics()["passive_blocks"] is True


@pytest.mark.parametrize("solver", ["gmres", "cocr"])
def test_matrix_free_inductance_uses_native_ngsolve_krylov(solver):
    basis = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        weights=np.ones(2),
        current_modes=np.array(
            [
                [[1.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
                [[0.0, 1.0, 0.0], [0.0, 0.25, 0.0]],
            ]
        ),
        names=["bulk0", "bulk1"],
    )
    system = vim.AssembleHybridVIM(
        basis,
        sigma=2.0,
        interaction=vim.HACApKSampledLaplaceInteraction(
            kernel_epsilon=0.2,
            cross_only=False,
            leaf_size=2,
        ),
    )
    s = 1j * 7.0
    rhs = np.array([1.0 + 0.5j, -0.25j])
    inductance = system.inductance.to_dense()
    expected = np.linalg.solve(
        system.resistance + s * inductance,
        rhs,
    )

    actual, solve_info = system.solve(
        s,
        rhs,
        solver=solver,
        return_diagnostics=True,
    )
    np.testing.assert_allclose(actual, expected, rtol=2.0e-10, atol=2.0e-12)
    assert solve_info["backend"] == f"ngsolve-base-matrix-{solver}"
    assert solve_info["native_term_count"] == 1
    assert solve_info["relative_residual_max"] < 1.0e-10
    info = system.diagnostics()
    assert info["inductance_matrix_free"] is True
    assert info["inductance_operator"]["operator_block_count"] == 1
    assert info["passive_blocks"] is True


def test_sampled_planar_log_hacapk_matches_dense_reference():
    basis = vim.VolumeCurrentBasis(
        points=np.array(
            [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0], [0.0, 0.6, 0.0]]
        ),
        weights=np.array([0.2, 0.3, 0.25]),
        current_modes=np.array(
            [
                [[0.0, 0.0, 1.0], [0.0, 0.0, -0.5], [0.0, 0.0, 0.2]],
                [[0.0, 0.0, 0.1], [0.0, 0.0, 0.4], [0.0, 0.0, -0.8]],
            ]
        ),
        names=("jz0", "jz1"),
    )
    backend = vim.HACApKSampledPlanarLogInteraction(
        mu=1.7,
        kernel_epsilon=0.05,
        reference_length=1.0,
        cross_only=False,
        leaf_size=2,
    )
    system = vim.AssembleHybridVIM(basis, sigma=3.0, interaction=backend)

    diff = basis.points[:, None, :2] - basis.points[None, :, :2]
    distance = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff) + 0.05**2)
    kernel = -(1.7 / (2.0 * np.pi)) * (
        basis.weights[:, None] * basis.weights[None, :]
    ) * np.log(distance)
    expected = np.einsum(
        "aik,bjk,ij->ab", basis.modes, basis.modes, kernel
    )
    np.testing.assert_allclose(
        system.inductance.to_dense(), expected, rtol=2.0e-12, atol=2.0e-13
    )


def test_coupled_hdiv_hcurl_keeps_two_independent_native_hmatrices():
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    weights = np.ones(3) / 3.0
    modes = np.array(
        [
            [[1.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-0.2, 0.0, 0.0]],
            [[0.0, 1.0, 0.0], [0.0, -0.3, 0.0], [0.0, 0.4, 0.0]],
        ]
    )
    current = vim.VolumeCurrentBasis(
        points, weights, modes, names=("e0", "e1")
    )
    magnetization = vim.MagnetizationBasis(
        points, weights, modes, names=("m0", "m1")
    )

    def hmatrix_system():
        return vim.AssembleHybridVIM(
            current,
            sigma=2.0,
            interaction=vim.HACApKSampledLaplaceInteraction(
                kernel_epsilon=0.1,
                cross_only=False,
                leaf_size=2,
            ),
        )

    eddy = hmatrix_system()
    magnetic = hmatrix_system().impedance_operator(0.75)
    coupling = np.array([[0.02, 0.01], [-0.01, 0.03]])
    coupled = vim.CoupledHDivHybridVIMSystem(
        magnetization,
        eddy,
        (current,),
        coupling,
        magnetic_operator=magnetic,
    )
    operator = coupled.mixed_operator(s=2j)
    expected = np.array([1.0 + 0.2j, -0.4j, 0.3 + 0.1j, -0.2 + 0.5j])
    rhs = operator.matvec(expected)
    result = coupled.solve(
        s=2j,
        magnetic_rhs=rhs[:2],
        eddy_rhs=rhs[2:],
        tolerance=1.0e-9,
    )

    np.testing.assert_allclose(
        result["solution"][:, 0], expected, rtol=5.0e-8, atol=5.0e-10
    )
    assert result["solver_diagnostics"]["native_term_count"] == 2
    assert result["solver_diagnostics"]["backend"] == "ngsolve-base-matrix-gmres"


def test_coupled_hdiv_hcurl_can_materialize_small_system_for_dense_lu():
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    weights = np.ones(3) / 3.0
    modes = np.array(
        [
            [[1.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-0.2, 0.0, 0.0]],
            [[0.0, 1.0, 0.0], [0.0, -0.3, 0.0], [0.0, 0.4, 0.0]],
        ]
    )
    current = vim.VolumeCurrentBasis(points, weights, modes, names=("e0", "e1"))
    magnetization = vim.MagnetizationBasis(
        points, weights, modes, names=("m0", "m1")
    )

    def hmatrix_system():
        return vim.AssembleHybridVIM(
            current,
            sigma=2.0,
            interaction=vim.HACApKSampledLaplaceInteraction(
                kernel_epsilon=0.1,
                cross_only=False,
                leaf_size=2,
            ),
        )

    eddy = hmatrix_system()
    magnetic = hmatrix_system().impedance_operator(0.75)
    coupled = vim.CoupledHDivHybridVIMSystem(
        magnetization,
        eddy,
        (current,),
        np.array([[0.02, 0.01], [-0.01, 0.03]]),
        magnetic_operator=magnetic,
    )
    operator = coupled.mixed_operator(s=2j)
    expected = np.array([1.0 + 0.2j, -0.4j, 0.3 + 0.1j, -0.2 + 0.5j])
    rhs = operator.matvec(expected)

    result = coupled.solve(
        s=2j,
        magnetic_rhs=rhs[:2],
        eddy_rhs=rhs[2:],
        solver="dense",
    )

    np.testing.assert_allclose(result["solution"][:, 0], expected)
    assert result["solver_diagnostics"] == {
        "backend": "native-dense-reduced-lu",
        "iterations": 0,
        "relative_residual_max": 0.0,
        "materialized_from_matrix_free": True,
    }


def test_ngsolve_bem_laplace_sl_projection_stays_as_base_matrix():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    mesh = ng.Mesh(occ.unit_cube.GenerateMesh(maxh=1.5))
    fes = ng.HDivSurface(mesh, order=0)
    basis = vim.SampledCurrentBasis(
        points=np.array([[0.5, 0.5, 0.0]]),
        weights=np.ones(1),
        modes=np.array(
            [[[1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0]]]
        ),
        kind="surface",
        names=("s0", "s1"),
    )
    projection = np.zeros((fes.ndof, 2))
    projection[0, 0] = 1.0
    projection[min(1, fes.ndof - 1), 1] = 1.0
    with ng.TaskManager():
        interaction = vim.NGSolveProjectedInteraction.from_laplace_sl(
            fes,
            (basis,),
            (projection,),
        )
        operator = interaction.build_operator((basis,))
        value = operator.matvec(np.array([1.0 + 0.2j, -0.3 + 0.1j]))

    assert operator.shape == (2, 2)
    assert np.all(np.isfinite(value))
    assert np.linalg.norm(value) > 0.0
    assert operator.stats()["backend"] == "ngsolve-base-matrix"


def test_hcurl_eddy_cln_model_preserves_vim_response_and_faraday_drive():
    system = vim.HybridVIMSystem(
        resistance=np.array([[2.0]]),
        inductance=np.array([[3.0]]),
        surface_mass=np.array([[0.0]]),
        basis_names=("eddy0",),
        blocks={"volume": (0, 1), "surface": (1, 1)},
    )
    port = np.array([[4.0]])
    model = vim.HCurlEddyCLNFromVIM(system, port)
    s = 1j * 7.0

    np.testing.assert_allclose(model.port_admittance(s), system.port_admittance(s, port))
    expected_current = (-s * 4.0) / (2.0 + 3.0 * s)
    np.testing.assert_allclose(
        model.solve_vector_potential_drive(s, 1.0),
        np.array([expected_current]),
    )
    state = model.derivative_input_state_space()
    np.testing.assert_allclose(state["A"], [[-2.0 / 3.0]])
    np.testing.assert_allclose(state["B"], [[4.0 / 3.0]])
    np.testing.assert_allclose(state["C"], [[4.0]])
    np.testing.assert_allclose(state["D"], [[0.0]])
    assert model.diagnostics()["passive"] is True
    assert model.diagnostics()["finite_rl_state_space"] is True


def test_hcurl_eddy_cln_model_requires_sibc_rationalization_for_state_space():
    model = vim.HCurlEddyCLNModel(
        resistance=np.array([[1.0]]),
        inductance=np.array([[1.0]]),
        surface_mass=np.array([[2.0]]),
        port_rhs=np.array([[1.0]]),
    )

    assert model.has_sibc_termination is True
    with pytest.raises(ValueError, match="rationalized"):
        model.derivative_input_state_space()


def test_hcurl_eddy_cln_matlab_exchange_preserves_row_major_matrices(tmp_path):
    model = vim.HCurlEddyCLNModel(
        resistance=np.array([[2.0, 0.1], [0.1, 1.0]]),
        inductance=np.array([[3.0, 0.2], [0.2, 2.0]]),
        surface_mass=np.zeros((2, 2)),
        port_rhs=np.array([[1.0], [0.5]]),
        basis_names=("eddy0", "eddy1"),
        blocks={"volume": (0, 2)},
    )
    force_operator = np.arange(6.0).reshape(3, 2, 1)
    destination = vim.ExportHCurlEddyCLNJSON(
        model,
        tmp_path / "hcurl_exchange.json",
        force_operator=force_operator,
        metadata={"frequency_hz": 50.0},
    )
    payload = json.loads(destination.read_text(encoding="utf-8"))

    assert payload["schema"] == "radia.hcurl.eddy_cln.exchange.v1"
    assert payload["arrays"]["resistance"]["shape"] == [2, 2]
    np.testing.assert_allclose(
        np.asarray(payload["arrays"]["resistance"]["values"]).reshape(2, 2),
        model.resistance,
    )
    assert payload["arrays"]["force_operator"]["shape"] == [3, 2, 1]
    np.testing.assert_allclose(
        payload["arrays"]["force_operator"]["values"],
        force_operator.ravel(order="C"),
    )
    assert payload["metadata"]["frequency_hz"] == 50.0

    sibc_model = vim.HCurlEddyCLNModel(
        resistance=np.eye(1),
        inductance=np.eye(1),
        surface_mass=np.ones((1, 1)),
        port_rhs=np.ones((1, 1)),
    )
    with pytest.raises(ValueError, match="rationalized"):
        vim.ExportHCurlEddyCLNJSON(sibc_model, tmp_path / "sibc.json")


def test_hcurl_eddy_cln_family_exchange_requires_common_sorted_state_basis(tmp_path):
    def make_model(scale):
        return vim.HCurlEddyCLNModel(
            resistance=scale * np.eye(1),
            inductance=np.eye(1),
            surface_mass=np.zeros((1, 1)),
            port_rhs=np.ones((1, 1)),
        )

    destination = vim.ExportHCurlEddyCLNFamilyJSON(
        [
            {"height_m": 1.0, "model": make_model(2.0)},
            {"height_m": -1.0, "model": make_model(1.0)},
        ],
        tmp_path / "hcurl_family.json",
    )
    payload = json.loads(destination.read_text(encoding="utf-8"))
    assert payload["schema"] == "radia.hcurl.eddy_cln.family.v1"
    assert [item["height_m"] for item in payload["snapshots"]] == [-1.0, 1.0]
    assert payload["shared_state_basis"] is True

    with pytest.raises(ValueError, match="strictly increasing"):
        vim.ExportHCurlEddyCLNFamilyJSON(
            [
                {"height_m": 0.0, "model": make_model(1.0)},
                {"height_m": 0.0, "model": make_model(2.0)},
            ],
            tmp_path / "duplicate.json",
        )


def test_surface_omega_basis_builds_tangential_current():
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    weights = np.array([0.5, 0.5])
    normals = np.array([[0.0, 0.0, 2.0], [0.0, 0.0, 3.0]])
    grad_omega = np.array([
        [[1.0, 0.0, 5.0], [1.0, 0.0, -2.0]],
    ])

    basis = vim.SurfaceOmegaBasis(points, weights, normals, grad_omega, names=["omega_x"])

    assert basis.kind == "surface"
    assert basis.names == ("omega_x",)
    np.testing.assert_allclose(
        basis.modes[0],
        np.array([[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
    )
    np.testing.assert_allclose(basis.mass_matrix(), [[1.0]])


def test_block_krylov_basis_respects_free_dofs_and_metric_orthonormality():
    k = np.diag([2.0, 3.0, 5.0])
    m = np.diag([1.0, 2.0, 4.0])
    ports = np.array([
        [1.0, 0.0],
        [1.0, 1.0],
        [1.0, 2.0],
    ])
    free = np.array([True, False, True])

    basis = vim.BlockKrylovBasis(k, m, ports, steps=2, free_dofs=free)

    assert isinstance(basis, vim.EVRSBasis)
    assert isinstance(basis, vim.ResponseBasis)
    assert np.array_equal(basis.active_dofs, [0, 2])
    assert np.allclose(basis.vectors[1, :], 0.0)
    gram = basis.vectors.conj().T @ m @ basis.vectors
    np.testing.assert_allclose(gram, np.eye(basis.rank), atol=1.0e-12)
    assert basis.diagnostics() == {
        "ndof": 3,
        "active_dofs": 2,
        "rank": basis.rank,
        "port_visible_dofs": basis.rank,
        "eddy_visible_dofs": basis.rank,
        "compression_ratio": basis.rank / 3,
        "inactive_dofs": 1,
        "port_invisible_dofs": 2 - basis.rank,
        "eddy_invisible_dofs": 2 - basis.rank,
        "eliminated_dofs": 3 - basis.rank,
        "port_count": 2,
        "krylov_steps": 2,
        "construction": "block-krylov",
        "parent_space": "HCurl",
    }


def test_block_krylov_relative_rank_test_is_port_scale_invariant():
    stiffness = np.diag([1.0, 1.0 + 1.0e-13])
    mass = np.eye(2)
    port = np.array([1.0, 1.0])

    base = vim.BlockKrylovBasis(
        stiffness,
        mass,
        port,
        steps=2,
        rtol=1.0e-10,
    )
    scaled = vim.BlockKrylovBasis(
        stiffness,
        mass,
        1.0e12 * port,
        steps=2,
        rtol=1.0e-10,
    )

    assert base.rank == 1
    assert scaled.rank == base.rank
    np.testing.assert_allclose(
        np.abs(scaled.vectors),
        np.abs(base.vectors),
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_current_gram_compression_removes_curl_null_response():
    response = vim.EVRSBasis(
        vectors=np.eye(3),
        active_dofs=np.arange(3),
        port_count=1,
        krylov_steps=3,
    )
    current = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        weights=np.ones(2),
        current_modes=np.array(
            [
                [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                [[1.0e-14, 0.0, 0.0], [0.0, 0.0, 0.0]],
            ]
        ),
    )

    compressed_response, compressed_current = (
        vim.CompressHCurlResponseInCurrentGram(
            response,
            current,
            rtol=1.0e-10,
        )
    )

    assert compressed_response.rank == 2
    assert compressed_current.n_modes == 2
    np.testing.assert_allclose(
        compressed_current.mass_matrix(),
        np.eye(2),
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        compressed_current.modes,
        np.einsum(
            "mk,mpc->kpc",
            np.linalg.lstsq(
                response.vectors,
                compressed_response.vectors,
                rcond=None,
            )[0],
            current.modes,
        ),
        atol=1.0e-12,
    )
    diagnostics = compressed_response.diagnostics()
    assert diagnostics["pre_current_gram_rank"] == 3
    assert diagnostics["current_gram_rank"] == 2
    assert diagnostics["construction"].endswith("+current-gram")


def test_evrs_tmethod_algebra_preserves_derham_gauge_and_ports():
    curl_map = np.array([
        [1.0, -1.0, 0.0],
        [-1.0, 1.0, 0.0],
    ])
    div_map = np.array([[1.0, 1.0]])
    grad_map = np.array([[1.0], [1.0], [0.0]])
    evrs_map = np.array([
        [1.0, 0.0],
        [0.0, 1.0],
        [0.0, 1.0],
    ])
    resistance_current = np.array([[2.0, 0.25], [0.25, 3.0]])
    inductance_current = np.array([[5.0, 0.5], [0.5, 7.0]])
    port_current = np.array([[1.0, 0.0], [0.0, 2.0]])

    result = vim.EVRSTMethodAlgebra(
        curl_map,
        div_map,
        grad_map,
        evrs_map,
        resistance_current,
        inductance_current,
        port_current,
        backend="python",
    )

    current_evrs = curl_map @ evrs_map
    np.testing.assert_allclose(result["current_evrs"], current_evrs)
    np.testing.assert_allclose(
        result["resistance_t"],
        curl_map.T @ resistance_current @ curl_map,
    )
    np.testing.assert_allclose(
        result["resistance_evrs"],
        current_evrs.T @ resistance_current @ current_evrs,
    )
    np.testing.assert_allclose(
        result["port_evrs"],
        evrs_map.T @ curl_map.T @ port_current,
    )
    info = result["diagnostics"]
    assert info["div_curl_norm"] == pytest.approx(0.0)
    assert info["div_evrs_norm"] == pytest.approx(0.0)
    assert info["resistance_gauge_norm"] == pytest.approx(0.0)
    assert info["inductance_gauge_norm"] == pytest.approx(0.0)
    assert info["port_gauge_norm"] == pytest.approx(0.0)
    assert info["evrs_resistance_galerkin_residual"] == pytest.approx(0.0)
    assert info["evrs_inductance_galerkin_residual"] == pytest.approx(0.0)

    try:
        cpp_result = vim.EVRSTMethodAlgebra(
            curl_map,
            div_map,
            grad_map,
            evrs_map,
            resistance_current,
            inductance_current,
            port_current,
            backend="cpp",
        )
    except RuntimeError:
        cpp_result = None
    if cpp_result is not None:
        for key in (
            "current_evrs",
            "resistance_t",
            "inductance_t",
            "resistance_evrs",
            "inductance_evrs",
            "port_t",
            "port_evrs",
        ):
            np.testing.assert_allclose(cpp_result[key], result[key])
        for key, value in result["diagnostics"].items():
            assert cpp_result["diagnostics"][key] == pytest.approx(value)


def test_hybrid_vim_mixes_t_volume_and_surface_omega_sibc_blocks():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        weights=np.array([0.25, 0.25]),
        current_modes=np.array([
            [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        ]),
        names=["T_loop"],
    )
    surface = vim.SurfaceOmegaBasis(
        points=np.array([[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]),
        weights=np.array([0.5, 0.5]),
        normals=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        grad_omega_modes=np.array([
            [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        ]),
        names=["Omega_skin"],
    )
    sigma = 5.8e7

    system = vim.AssembleHybridVIM(volume, surface, sigma=sigma, kernel_epsilon=0.2)

    assert system.basis_names == ("T_loop", "Omega_skin")
    assert system.blocks["volume"] == (0, 1)
    assert system.blocks["surface"] == (1, 2)
    np.testing.assert_allclose(system.resistance, system.resistance.conj().T)
    np.testing.assert_allclose(system.inductance, system.inductance.conj().T)
    np.testing.assert_allclose(system.surface_mass, system.surface_mass.conj().T)
    np.testing.assert_allclose(system.resistance[0, 0], 0.5 / sigma)
    np.testing.assert_allclose(system.resistance[1, 1], 0.0)
    np.testing.assert_allclose(system.surface_mass[1, 1], 1.0)
    assert system.inductance[0, 0].real > 0.0
    assert system.inductance[1, 1].real > 0.0
    assert abs(system.inductance[0, 1]) > 0.0

    omega = 2.0 * np.pi * 10_000.0
    zs = vim.SkinImpedance(1j * omega, sigma)
    z = system.impedance(1j * omega, surface_impedance=zs)
    q = np.array([1.0 - 0.3j, 0.4 + 0.2j])
    dissipated = np.vdot(q, z @ q).real
    assert dissipated > 0.0
    assert system.n_modes == 2

    port_rhs = np.eye(2)
    admittance = system.port_admittance(1j * omega, port_rhs, surface_impedance=zs)
    np.testing.assert_allclose(
        admittance,
        vim.ReducedPortAdmittance(system, 1j * omega, port_rhs, surface_impedance=zs),
    )
    port_impedance = vim.ReducedPortImpedance(
        system,
        1j * omega,
        port_rhs,
        surface_impedance=zs,
    )
    np.testing.assert_allclose(port_impedance @ admittance, np.eye(2), atol=1.0e-10)

    rhs_v = vim.ExternalVectorPotentialRHS(volume, [1.0, 0.0, 0.0])
    rhs_s = vim.ExternalVectorPotentialRHS(surface, [1.0, 0.0, 0.0])
    np.testing.assert_allclose(rhs_v, [0.5])
    np.testing.assert_allclose(rhs_s, [-1.0])


def test_local_surface_impedance_gram_preserves_esim_sample_variation():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([[[0.0, 0.0, 1.0]]]),
        names=["bulk"],
    )
    surface = vim.SampledCurrentBasis(
        points=np.array([[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]),
        weights=np.array([1.0, 2.0]),
        modes=np.array([
            [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        ]),
        kind="surface",
        names=("skin_left", "skin_right"),
    )
    system = vim.AssembleHybridVIM(
        volume,
        surface,
        sigma=2.0,
        kernel_epsilon=0.2,
    )
    values = np.array([2.0 + 1.0j, 5.0 + 3.0j])

    gram = vim.AssembleSurfaceImpedanceGram(
        system,
        surface,
        lambda points: values,
        label="esim-per-sample",
    )

    expected = np.diag([0.0, 2.0 + 1.0j, 10.0 + 6.0j])
    np.testing.assert_allclose(gram.matrix, expected)
    np.testing.assert_allclose(
        system.impedance(0.0, surface_impedance=gram),
        system.resistance + expected,
    )
    assert isinstance(gram, vim.SurfaceImpedanceGram)
    info = gram.diagnostics()
    assert info["label"] == "esim-per-sample"
    assert info["samples"] == 2
    assert info["real_min"] == pytest.approx(2.0)
    assert info["real_max"] == pytest.approx(5.0)
    assert info["passive"]

    orthogonalized = system.mixed_galerkin_orthogonalization(
        "surface",
        "volume",
        0.5j,
        surface_impedance=gram,
    )
    np.testing.assert_allclose(
        orthogonalized.reduced_operator,
        system.schur_complement(
            "surface",
            "volume",
            0.5j,
            surface_impedance=gram,
        ),
    )


def test_local_surface_impedance_gram_rejects_implicit_mode_diagonal():
    system = vim.HybridVIMSystem(
        resistance=np.eye(2),
        inductance=np.eye(2),
        surface_mass=np.diag([0.0, 1.0]),
        basis_names=("bulk", "surface"),
        blocks={"volume": (0, 1), "surface": (1, 2)},
    )

    with pytest.raises(TypeError, match="AssembleSurfaceImpedanceGram"):
        system.impedance(1j, surface_impedance=np.array([1.0, 2.0]))


def test_local_surface_impedance_gram_rejects_a_different_surface_mass():
    surface = vim.SampledCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        modes=np.array([[[1.0, 0.0, 0.0]]]),
        kind="surface",
        names=("skin",),
    )
    source = vim.HybridVIMSystem(
        resistance=np.zeros((1, 1)),
        inductance=np.zeros((1, 1)),
        surface_mass=np.ones((1, 1)),
        basis_names=("skin",),
        blocks={"surface": (0, 1)},
    )
    target = vim.HybridVIMSystem(
        resistance=np.zeros((1, 1)),
        inductance=np.zeros((1, 1)),
        surface_mass=np.array([[2.0]]),
        basis_names=("skin",),
        blocks={"surface": (0, 1)},
    )
    gram = vim.AssembleSurfaceImpedanceGram(source, surface, 1.0 + 0.5j)

    with pytest.raises(ValueError, match="different hybrid system"):
        target.impedance(0.0, surface_impedance=gram)
    with pytest.raises(ValueError, match="quadrature/modes"):
        vim.AssembleSurfaceImpedanceGram(target, surface, 1.0 + 0.5j)


def test_local_esim_surface_vim_runs_the_scipy_cell_problem_and_is_consistent():
    surface = vim.SampledCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        modes=np.array([[[1.0, 0.0, 0.0]]]),
        kind="surface",
        names=("skin",),
    )
    system = vim.HybridVIMSystem(
        resistance=np.zeros((1, 1)),
        inductance=np.array([[1.0e-10]]),
        surface_mass=np.ones((1, 1)),
        basis_names=("skin",),
        blocks={"surface": (0, 1)},
    )
    model = vim.LocalESIMSurfaceModel(
        bh_curve=np.array([
            [0.0, 0.0],
            [100.0, 0.1],
            [500.0, 0.45],
            [2_000.0, 1.2],
            [10_000.0, 1.6],
        ]),
        sigma=2.0e6,
        bins=2,
        n_nodes=30,
        cell_tolerance=1.0e-4,
        cell_max_iterations=60,
    )
    rhs = np.array([1.0])
    frequency_hz = 10_000.0

    solved = vim.SolveLocalESIMSurfaceVIM(
        system,
        surface,
        rhs,
        model,
        frequency_hz,
        outer_tolerance=5.0e-3,
        outer_max_iterations=12,
        outer_relaxation=0.7,
    )

    assert solved.converged
    assert solved.iterations >= 1
    assert solved.surface_impedance.diagnostics()["passive"]
    assert solved.history[-1]["cell_model"].endswith("1d-esim")
    assert solved.history[0]["cell_solve_count"] == 2
    assert all(row["cell_solve_count"] == 0 for row in solved.history[1:])
    operator = system.impedance(
        1j * 2.0 * np.pi * frequency_hz,
        surface_impedance=solved.surface_impedance,
    )
    np.testing.assert_allclose(operator @ solved.coefficients, rhs)

    partial = vim.SolveLocalESIMSurfaceVIM(
        system,
        surface,
        rhs,
        model,
        frequency_hz,
        outer_tolerance=1.0e-14,
        outer_max_iterations=1,
        raise_on_nonconvergence=False,
    )
    assert not partial.converged
    partial_operator = system.impedance(
        1j * 2.0 * np.pi * frequency_hz,
        surface_impedance=partial.surface_impedance,
    )
    np.testing.assert_allclose(partial_operator @ partial.coefficients, rhs)

    with pytest.raises(ValueError, match="one nonlinear excitation"):
        vim.SolveLocalESIMSurfaceVIM(
            system,
            surface,
            np.ones((1, 2)),
            model,
            frequency_hz,
        )


def test_local_esim_surface_lut_roundtrip_interpolation_and_model_guard(tmp_path):
    curve = np.array([
        [0.0, 0.0],
        [100.0, 0.1],
        [500.0, 0.45],
        [2_000.0, 1.2],
        [10_000.0, 1.6],
    ])
    cell_model = vim.LocalESIMSurfaceModel(
        bh_curve=curve,
        sigma=2.0e6,
        bins=3,
        n_nodes=30,
        cell_tolerance=1.0e-4,
        cell_max_iterations=60,
        h_floor=1.0,
    )
    frequencies = np.array([10_000.0, 40_000.0])
    fields = np.array([1.0, 100.0, 10_000.0])
    path = tmp_path / "steel-local-esim.npz"

    built = vim.BuildLocalESIMSurfaceLUT(
        cell_model,
        frequencies,
        fields,
        output_path=path,
    )
    loaded = vim.LocalESIMSurfaceLUT.load(path)

    assert path.is_file()
    assert loaded.model_key == built.model_key
    assert loaded.diagnostics()["passive"]
    assert loaded.diagnostics()["interpolation"].endswith("logabs-phase-v1")
    np.testing.assert_allclose(loaded.impedance_ohm, built.impedance_ohm)
    np.testing.assert_allclose(
        loaded.evaluate(fields, frequencies[0]),
        built.impedance_ohm[0],
    )
    quality = vim.ValidateLocalESIMSurfaceLUT(
        cell_model,
        loaded,
        frequencies_hz=[20_000.0],
        h_values_A_per_m=[10.0],
    )
    assert quality["sample_count"] == 1
    assert quality["direct_cells_passive"]
    assert quality["max_relative_error"] < 0.2

    lut_model = cell_model.with_lut(loaded)
    values, diagnostics = lut_model.impedance_samples(
        np.array([10.0, 1_000.0]),
        20_000.0,
    )
    assert np.all(np.isfinite(values))
    assert np.min(values.real) >= -1.0e-12
    assert diagnostics["cell_model"] == "precomputed-2d-local-esim-lut"
    assert diagnostics["cell_solve_count"] == 0
    assert not diagnostics["cell_table_rebuilt"]

    surface = vim.SampledCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        modes=np.array([[[1.0, 0.0, 0.0]]]),
        kind="surface",
        names=("skin",),
    )
    system = vim.HybridVIMSystem(
        resistance=np.zeros((1, 1)),
        inductance=np.array([[1.0e-10]]),
        surface_mass=np.ones((1, 1)),
        basis_names=("skin",),
        blocks={"surface": (0, 1)},
    )
    solved = vim.SolveLocalESIMSurfaceVIM(
        system,
        surface,
        np.array([1.0]),
        lut_model,
        frequencies[0],
        outer_tolerance=5.0e-3,
        outer_max_iterations=12,
        outer_relaxation=0.7,
    )
    assert solved.converged
    assert all(row["cell_solve_count"] == 0 for row in solved.history)
    assert all(row["cell_model"].endswith("local-esim-lut") for row in solved.history)

    with pytest.raises(ValueError, match="outside the local ESIM LUT range"):
        lut_model.impedance_samples([20_000.0], 20_000.0)
    with pytest.raises(ValueError, match="frequency is outside"):
        lut_model.impedance_samples([100.0], 80_000.0)
    with pytest.raises(ValueError, match="model signature does not match"):
        vim.LocalESIMSurfaceModel(
            bh_curve=curve,
            sigma=3.0e6,
            bins=3,
            n_nodes=30,
            cell_tolerance=1.0e-4,
            cell_max_iterations=60,
            h_floor=1.0,
            lut=loaded,
        )


def test_hybrid_vim_accepts_custom_bem_interaction_backend():
    basis = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([
            [[1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
        ]),
        names=["jx", "jy"],
    )

    def fake_bem(left, right):
        assert left is basis
        assert right is basis
        return np.array([[2.0, 0.25], [0.25, 3.0]])

    system = vim.AssembleHybridVIM(basis, sigma=2.0, interaction=fake_bem)

    np.testing.assert_allclose(
        system.inductance,
        np.array([[2.0, 0.25], [0.25, 3.0]]),
    )


def test_explicit_sampled_laplace_interaction_matches_default_backend():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        weights=np.array([0.25, 0.25]),
        current_modes=np.array([
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        ]),
    )
    surface = vim.SurfaceOmegaBasis(
        points=np.array([[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]),
        weights=np.array([0.5, 0.5]),
        normals=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        grad_omega_modes=np.array([
            [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
        ]),
    )

    default = vim.AssembleHybridVIM(volume, surface, sigma=2.0, kernel_epsilon=0.2)
    explicit = vim.AssembleHybridVIM(
        volume,
        surface,
        sigma=2.0,
        interaction=vim.SampledLaplaceInteraction(kernel_epsilon=0.2),
    )

    np.testing.assert_allclose(explicit.inductance, default.inductance)
    np.testing.assert_allclose(explicit.resistance, default.resistance)
    np.testing.assert_allclose(explicit.surface_mass, default.surface_mass)
    assert explicit.interaction_backend == "sampled-laplace"
    info = explicit.diagnostics()
    assert info["interaction_backend"] == "sampled-laplace"
    assert info["passive_blocks"] is True
    assert info["inductance_hermitian_error"] == pytest.approx(0.0)
    assert info["blocks"] == {"volume": [0, 1], "surface": [1, 2]}


def test_hacapk_sampled_interaction_keeps_full_hybrid_operator_matrix_free():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0], [0.8, 0.1, 0.0]]),
        weights=np.array([0.3, 0.4]),
        current_modes=np.array(
            [
                [[1.0, 0.2, 0.0], [0.3, 0.8, 0.1]],
                [[0.0, 0.7, 0.2], [0.5, 0.0, 0.9]],
            ]
        ),
        names=["T0", "T1"],
    )
    surface = vim.SampledCurrentBasis(
        points=np.array([[0.1, 1.0, 0.0], [0.9, 1.1, 0.2]]),
        weights=np.array([0.5, 0.6]),
        modes=np.array([[[0.4, 0.9, 0.0], [0.2, 0.6, 0.3]]]),
        kind="surface",
        names=("Omega0",),
    )
    mu = 1.7
    kernel_epsilon = 0.25
    dense = vim.AssembleHybridVIM(
        volume,
        surface,
        sigma=3.0,
        interaction=vim.SampledLaplaceInteraction(
            mu=mu,
            kernel_epsilon=kernel_epsilon,
        ),
    )
    compressed = vim.AssembleHybridVIM(
        volume,
        surface,
        sigma=3.0,
        interaction=vim.HACApKSampledLaplaceInteraction(
            mu=mu,
            kernel_epsilon=kernel_epsilon,
            aca_eps=1.0e-12,
            leaf_size=64,
            cross_only=False,
        ),
    )

    assert compressed.diagnostics()["inductance_matrix_free"] is True
    stats = compressed.diagnostics()["inductance_operator"]
    assert stats["operator_block_count"] == 1
    assert stats["operator_blocks"][0]["stats"]["cross_only"] is False
    np.testing.assert_allclose(
        compressed.inductance.to_dense(),
        dense.inductance,
        rtol=2.0e-12,
        atol=2.0e-13,
    )


def test_hacapk_cross_operator_preserves_matrix_free_diagonal_backends_and_gmres():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0], [0.7, 0.0, 0.1]]),
        weights=np.array([0.35, 0.45]),
        current_modes=np.array(
            [
                [[1.0, 0.0, 0.2], [0.4, 0.7, 0.0]],
                [[0.1, 0.8, 0.0], [0.0, 0.2, 0.9]],
            ]
        ),
        names=["T0", "T1"],
    )
    surface = vim.SampledCurrentBasis(
        points=np.array([[0.0, 0.9, 0.0], [0.8, 1.0, 0.2]]),
        weights=np.array([0.55, 0.65]),
        modes=np.array([[[0.3, 1.0, 0.0], [0.5, 0.1, 0.4]]]),
        kind="surface",
        names=("Omega0",),
    )
    mu = 2.3
    kernel_epsilon = 0.2
    reference_backend = vim.SampledLaplaceInteraction(
        mu=mu,
        kernel_epsilon=kernel_epsilon,
    )
    reference = vim.AssembleHybridVIM(
        volume,
        surface,
        sigma=4.0,
        interaction=reference_backend,
    )

    exact_volume_operator = vim.HACApKSampledLaplaceInteraction(
        mu=mu,
        kernel_epsilon=kernel_epsilon,
        aca_eps=1.0e-12,
        leaf_size=64,
        cross_only=False,
    ).build_operator((volume,))
    exact_volume = exact_volume_operator.to_dense()

    def exact_diagonal(left, right):
        assert left is volume
        assert right is volume
        return exact_volume_operator

    compressed = vim.AssembleHybridVIM(
        volume,
        surface,
        sigma=4.0,
        interaction=vim.HACApKSampledLaplaceInteraction(
            mu=mu,
            kernel_epsilon=kernel_epsilon,
            aca_eps=1.0e-12,
            leaf_size=64,
            cross_only=True,
            diagonal_interaction=exact_diagonal,
            diagonal_bases=(volume,),
        ),
    )
    expected_l = np.array(reference.inductance, copy=True)
    expected_l[:2, :2] = exact_volume
    actual_l = compressed.inductance.to_dense()
    np.testing.assert_allclose(actual_l, expected_l, rtol=2.0e-12, atol=2.0e-13)
    stats = compressed.inductance.stats()
    assert stats["operator_block_count"] == 2
    assert stats["operator_blocks"][-1]["start"] == 0
    assert stats["operator_blocks"][-1]["stop"] == 3
    assert stats["operator_blocks"][-1]["stats"]["cross_only"] is True

    s = 0.4 + 2.1j
    zs = 0.7 + 0.3j
    rhs = np.array([1.0 + 0.2j, -0.4j, 0.5 - 0.1j])
    expected = np.linalg.solve(
        reference.resistance + s * expected_l + zs * reference.surface_mass,
        rhs,
    )
    np.testing.assert_allclose(
        compressed.solve(s, rhs, surface_impedance=zs),
        expected,
        rtol=3.0e-11,
        atol=3.0e-12,
    )


def test_reduced_interaction_matrix_backend_wires_bem_blocks():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        weights=np.array([0.5, 0.5]),
        current_modes=np.array([
            [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        ]),
        names=["T0", "T1"],
    )
    surface = vim.SurfaceOmegaBasis(
        points=np.array([[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]),
        weights=np.array([0.5, 0.5]),
        normals=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        grad_omega_modes=np.array([
            [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        ]),
        names=["Omega0"],
    )
    reduced_bem = np.array([
        [2.0, 0.1, -0.2],
        [0.1, 3.0, 0.4],
        [-0.2, 0.4, 5.0],
    ])
    backend = vim.ReducedInteractionMatrix(
        (volume, surface),
        reduced_bem,
        name="mock-reduced-bem",
    )

    system = vim.AssembleHybridVIM(volume, surface, sigma=5.0, interaction=backend)

    np.testing.assert_allclose(system.inductance, reduced_bem)
    assert system.interaction_backend == "mock-reduced-bem"
    assert system.diagnostics()["interaction_backend"] == "mock-reduced-bem"
    assert system.blocks["volume"] == (0, 2)
    assert system.blocks["surface"] == (2, 3)
    np.testing.assert_allclose(system.resistance[:2, :2], np.eye(2) / 5.0)
    np.testing.assert_allclose(system.surface_mass[2:, 2:], [[1.0]])


def test_interaction_backend_only_needs_upper_triangular_blocks():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([[[1.0, 0.0, 0.0]]]),
        names=["T0"],
    )
    surface = vim.SurfaceOmegaBasis(
        points=np.array([[0.0, 1.0, 0.0]]),
        weights=np.array([1.0]),
        normals=np.array([[0.0, 0.0, 1.0]]),
        grad_omega_modes=np.array([[[0.0, 1.0, 0.0]]]),
        names=["Omega0"],
    )
    calls = []

    def backend(left, right):
        calls.append((left.names, right.names))
        assert not (left is surface and right is volume)
        if left is volume and right is volume:
            return np.array([[2.0]])
        if left is volume and right is surface:
            return np.array([[0.5]])
        if left is surface and right is surface:
            return np.array([[3.0]])
        raise AssertionError("unexpected block request")

    system = vim.AssembleHybridVIM(volume, surface, sigma=1.0, interaction=backend)

    assert calls == [
        (("T0",), ("T0",)),
        (("T0",), ("Omega0",)),
        (("Omega0",), ("Omega0",)),
    ]
    np.testing.assert_allclose(system.inductance, [[2.0, 0.5], [0.5, 3.0]])


def test_hybrid_vim_impedance_is_positive_real_for_passive_blocks():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        weights=np.array([0.5, 0.5]),
        current_modes=np.array([
            [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        ]),
    )
    surface = vim.SurfaceOmegaBasis(
        points=np.array([[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]),
        weights=np.array([0.5, 0.5]),
        normals=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        grad_omega_modes=np.array([
            [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        ]),
    )
    passive_l = np.array([
        [3.0, 0.2, 0.1],
        [0.2, 2.0, 0.3],
        [0.1, 0.3, 4.0],
    ])
    backend = vim.ReducedInteractionMatrix((volume, surface), passive_l)
    system = vim.AssembleHybridVIM(volume, surface, sigma=5.0, interaction=backend)
    z = system.impedance(0.2 + 1.3j, surface_impedance=0.7 + 0.9j)

    rng = np.random.default_rng(4)
    for _ in range(20):
        q = rng.normal(size=3) + 1j * rng.normal(size=3)
        assert np.vdot(q, z @ q).real > 0.0


def test_magnetization_current_coupling_matches_biot_savart_orientation():
    current = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([2.0]),
        current_modes=np.array([[[1.0, 0.0, 0.0]]]),
        names=["jx"],
    )
    magnetization = vim.MagnetizationBasis(
        points=np.array([[0.0, 1.0, 0.0]]),
        weights=np.array([3.0]),
        magnetization_modes=np.array([[[0.0, 0.0, 1.0]]]),
        names=["mz"],
    )

    fields = vim.CurrentMagneticFluxDensitySamples(
        current,
        magnetization.points,
        kernel_epsilon=0.0,
    )
    expected_bz = vim.MU0 / (4.0 * np.pi) * 2.0
    np.testing.assert_allclose(fields[0, 0], [0.0, 0.0, expected_bz])

    coupling = vim.MagnetizationCurrentCoupling(
        magnetization,
        current,
        kernel_epsilon=0.0,
    )
    np.testing.assert_allclose(coupling, [[3.0 * expected_bz]])


def test_magnetization_current_coupling_accepts_surface_omega_current():
    surface = vim.SurfaceOmegaBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([2.0]),
        normals=np.array([[0.0, 0.0, 1.0]]),
        grad_omega_modes=np.array([[[1.0, 0.0, 0.0]]]),
        names=["omega_x"],
    )
    magnetization = vim.MagnetizationBasis(
        points=np.array([[1.0, 0.0, 0.0]]),
        weights=np.array([1.5]),
        magnetization_modes=np.array([[[0.0, 0.0, -1.0]]]),
        names=["minus_mz"],
    )

    # K = n x grad(Omega) = e_y.  At x=e_x, e_y x e_x = -e_z.
    coupling = vim.MagnetizationCurrentCoupling(
        magnetization,
        surface,
        kernel_epsilon=0.0,
    )

    expected = 1.5 * vim.MU0 / (4.0 * np.pi) * 2.0
    np.testing.assert_allclose(coupling, [[expected]])


def test_coupled_hdiv_evrs_system_names_rectangular_mixed_block():
    eddy = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([2.0]),
        current_modes=np.array([
            [[0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0]],
        ]),
        names=["EVRS0", "EVRS1"],
    )
    magnetization = vim.MagnetizationBasis(
        points=np.array([[1.0, 0.0, 0.0]]),
        weights=np.array([3.0]),
        magnetization_modes=np.array([[[0.0, 0.0, -1.0]]]),
        names=["HDivM0"],
    )
    eddy_system = vim.AssembleHybridVIM(eddy, sigma=4.0, kernel_epsilon=0.5)

    coupled = vim.CoupleHDivMagnetizationToEVRS(
        magnetization,
        eddy,
        eddy_system=eddy_system,
    )

    assert isinstance(coupled, vim.CoupledHDivEVRSSystem)
    assert coupled.n_hdiv_modes == 1
    assert coupled.n_evrs_modes == 2
    np.testing.assert_allclose(coupled.coupling, [[6.0e-7, 0.0]], atol=1.0e-18)
    info = coupled.diagnostics()
    assert info["hdiv_modes"] == 1
    assert info["evrs_modes"] == 2
    assert info["hdiv_mmm_modes"] == 1
    assert info["hcurl_vim_modes"] == 2
    assert info["has_eddy_system"] is True
    assert info["has_shared_material_model"] is False
    assert info["coupling_frobenius_norm"] > 0.0

    material_model = vim.SharedMeshMaterialModel(
        mesh="shared-ngsolve-mesh",
        magnetic_regions="iron",
        conductive_regions=("conductor",),
        mu={"iron": vim.MU0 * 1000.0},
        sigma={"conductor": 5.8e7},
    )
    production_named = vim.CoupleHCurlVIMWithHDivMMM(
        magnetization,
        eddy,
        eddy_system=eddy_system,
        material_model=material_model,
    )

    assert isinstance(production_named, vim.HCurlVIMHDivMMMSystem)
    assert production_named.material_model is material_model
    assert production_named.n_hdiv_mmm_modes == 1
    assert production_named.n_hcurl_vim_modes == 2
    assert production_named.diagnostics()["has_shared_material_model"] is True
    assert production_named.diagnostics()["has_shared_mesh_material_model"] is True
    assert material_model.hdiv_mmm_coefficient() == {"iron": vim.MU0 * 1000.0}
    assert material_model.hcurl_vim_coefficient() == {"conductor": 5.8e7}
    assert material_model.diagnostics()["magnetic_region_count"] == 1
    assert material_model.diagnostics()["conductive_region_count"] == 1
    np.testing.assert_allclose(
        production_named.mixed_energy(np.array([2.0]), np.array([3.0, 5.0])),
        36.0e-7,
        atol=1.0e-18,
    )


def test_eddy_bubble_hcurl_basis_is_vim_and_hdiv_mmm_ready():
    response = vim.EVRSBasis(
        vectors=np.eye(2),
        active_dofs=np.array([0, 1]),
        port_count=1,
        krylov_steps=1,
    )
    current = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([
            [[1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
        ]),
        names=["j0", "j1"],
    )
    free = np.ones(2, dtype=bool)
    policy = vim.EddyDofPolicy(
        free=free,
        sibc_surface=np.zeros(2, dtype=bool),
        surface_candidate=np.zeros(2, dtype=bool),
        loop_bridge=np.zeros(2, dtype=bool),
        local_bubble=np.zeros(2, dtype=bool),
        interface=np.zeros(2, dtype=bool),
        wirebasket=np.zeros(2, dtype=bool),
    )
    bubbling = vim.EddyBubbleReduction(policy, evrs_rank=2, surface_modes=0)
    basis = vim.EddyBubbleHCurlBasis(response, current, bubbling)

    assert basis.n_modes == 2
    assert basis.eddy_basis is current
    system = basis.assemble_vim(sigma=5.8e7, kernel_epsilon=0.2)
    assert isinstance(system, vim.HybridVIMSystem)
    assert system.n_modes == 2

    magnetization = vim.MagnetizationBasis(
        points=np.array([[0.0, 0.0, 1.0]]),
        weights=np.array([1.0]),
        magnetization_modes=np.array([[[0.0, 1.0, 0.0]]]),
        names=["m0"],
    )
    coupled = vim.CoupleEddyBubbleHCurlBasisWithHDivMMM(
        magnetization,
        basis,
        eddy_system=system,
        kernel_epsilon=0.2,
    )
    assert isinstance(coupled, vim.HCurlVIMHDivMMMSystem)
    assert coupled.eddy_basis is current
    info = basis.diagnostics()
    assert info["kind"] == "EddyBubbleHCurlBasis"
    assert info["modes"] == 2
    assert info["eddy_bubbling"]["rule"] == "topology-aware-eddy-bubbling"


def test_hdiv_mmm_couples_to_hybrid_vim_volume_and_surface_blocks():
    volume = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([[[1.0, 0.0, 0.0]]]),
        names=["bulk"],
    )
    surface = vim.SurfaceOmegaBasis(
        points=np.array([[0.0, 1.0, 0.0]]),
        weights=np.array([2.0]),
        normals=np.array([[0.0, 0.0, 1.0]]),
        grad_omega_modes=np.array([[[1.0, 0.0, 0.0]]]),
        names=["skin"],
    )
    magnetization = vim.MagnetizationBasis(
        points=np.array([[1.0, 0.0, 0.0]]),
        weights=np.array([1.5]),
        magnetization_modes=np.array([[[0.0, 0.0, 1.0]]]),
        names=["m0"],
    )
    eddy_system = vim.AssembleHybridVIM(volume, surface, sigma=5.8e7, kernel_epsilon=0.2)

    coupled = vim.CoupleHybridVIMWithHDivMMM(
        magnetization,
        eddy_system,
        (volume, surface),
        kernel_epsilon=0.2,
    )

    assert isinstance(coupled, vim.CoupledHDivHybridVIMSystem)
    assert coupled.n_hdiv_mmm_modes == 1
    assert coupled.n_hcurl_vim_modes == eddy_system.n_modes
    assert coupled.coupling.shape == (1, eddy_system.n_modes)
    op = coupled.mixed_operator(np.array([[2.0]]), 1j * 100.0, surface_impedance=0.3)
    assert op.shape == (1 + eddy_system.n_modes, 1 + eddy_system.n_modes)
    info = coupled.diagnostics()
    assert info["eddy_basis_count"] == 2
    assert info["hybrid_blocks"] == eddy_system.diagnostics()["blocks"]


def test_full_coupled_mixed_galerkin_is_exact_with_eliminated_rhs():
    bulk = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([[[1.0, 0.0, 0.0]]]),
        names=["bulk"],
    )
    bridge = vim.VolumeCurrentBasis(
        points=np.array([[0.5, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([[[0.0, 1.0, 0.0]]]),
        names=["bridge"],
    )
    surface = vim.SurfaceOmegaBasis(
        points=np.array([[0.0, 0.5, 0.0]]),
        weights=np.array([1.0]),
        normals=np.array([[0.0, 0.0, 1.0]]),
        grad_omega_modes=np.array([[[1.0, 0.0, 0.0]]]),
        names=["surface"],
    )
    magnetization = vim.MagnetizationBasis(
        points=np.array([[0.25, 0.25, 0.25]]),
        weights=np.array([1.0]),
        magnetization_modes=np.array([[[0.0, 0.0, 1.0]]]),
        names=["m0"],
    )
    eddy_system = vim.HybridVIMSystem(
        resistance=np.array([
            [4.0, 0.3, -0.2],
            [0.3, 3.0, 0.4],
            [-0.2, 0.4, 2.5],
        ]),
        inductance=np.array([
            [1.0, 0.1, 0.05],
            [0.1, 0.8, -0.1],
            [0.05, -0.1, 0.6],
        ]),
        surface_mass=np.diag([0.0, 0.0, 1.2]),
        basis_names=("bulk", "bridge", "surface"),
        blocks={"bulk": (0, 1), "bridge": (1, 2), "surface": (2, 3)},
    )
    coupled = vim.CoupledHDivHybridVIMSystem(
        magnetization_basis=magnetization,
        eddy_system=eddy_system,
        eddy_bases=(bulk, bridge, surface),
        coupling=np.array([[0.7, -0.25, 0.4]]),
        magnetic_operator=np.array([[5.0]]),
    )
    s = 0.3 + 0.8j
    zs = 0.15 + 0.2j
    full = coupled.mixed_operator(s=s, surface_impedance=zs)
    reduced = coupled.mixed_galerkin_orthogonalization(
        ("bridge", "surface"),
        "bulk",
        s=s,
        surface_impedance=zs,
    )

    keep = np.array([0, 2, 3])
    eliminate = np.array([1])
    expected_schur = (
        full[np.ix_(keep, keep)]
        - full[np.ix_(keep, eliminate)]
        @ np.linalg.solve(
            full[np.ix_(eliminate, eliminate)],
            full[np.ix_(eliminate, keep)],
        )
    )
    assert isinstance(reduced, vim.MixedGalerkinHDivHybridVIMSystem)
    np.testing.assert_allclose(reduced.reduced_operator, expected_schur)
    assert reduced.n_hdiv_modes == 1
    assert reduced.n_hcurl_retained_modes == 2
    assert reduced.n_hcurl_eliminated_modes == 1

    rhs = np.array([1.1, 0.7, -0.3, 0.2])
    solved = reduced.solve(
        magnetic_rhs=rhs[:1],
        eddy_rhs=rhs[1:],
        require_excitation=True,
        return_operator=True,
    )
    expected_reduced_rhs = (
        rhs[keep, np.newaxis]
        - full[np.ix_(keep, eliminate)]
        @ np.linalg.solve(
            full[np.ix_(eliminate, eliminate)],
            rhs[eliminate, np.newaxis],
        )
    )
    np.testing.assert_allclose(solved["reduced_rhs"], expected_reduced_rhs)
    np.testing.assert_allclose(
        solved["full_solution"].ravel(),
        np.linalg.solve(full, rhs),
    )
    assert solved["full_residual_relative_norm"] < 1.0e-14
    assert solved["projected_residual_relative_norm"] < 1.0e-14
    info = reduced.diagnostics()
    assert info["full_coupled_schur"] is True
    assert info["schur_relative_error"] < 1.0e-14


def test_shared_mesh_material_model_validates_positive_scalar_coefficients():
    model = vim.SharedMeshMaterialModel(
        mesh="mesh",
        magnetic_regions=("iron", "air"),
        conductive_regions="coil",
        nu=2.0,
        sigma=5.8e7,
    )

    assert model.has_reluctivity is True
    assert model.has_magnetic_law is True
    assert model.has_conductivity is True
    assert model.hdiv_mmm_coefficient() == 2.0
    assert model.hcurl_vim_coefficient() == 5.8e7

    with pytest.raises(ValueError, match="sigma"):
        vim.SharedMeshMaterialModel(mesh="mesh", sigma=0.0)


def test_hcurl_vim_hdiv_mmm_system_solves_mixed_block_and_schur():
    eddy = vim.VolumeCurrentBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        current_modes=np.array([
            [[1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
        ]),
        names=["j0", "j1"],
    )
    magnetization = vim.MagnetizationBasis(
        points=np.array([[0.0, 0.0, 0.0]]),
        weights=np.array([1.0]),
        magnetization_modes=np.array([[[0.0, 0.0, 1.0]]]),
        names=["m0"],
    )
    eddy_system = vim.HybridVIMSystem(
        resistance=np.diag([2.0, 3.0]),
        inductance=np.diag([0.5, 0.25]),
        surface_mass=np.zeros((2, 2)),
        basis_names=("j0", "j1"),
        blocks={"eddy": (0, 2)},
    )
    coupled = vim.HCurlVIMHDivMMMSystem(
        magnetization_basis=magnetization,
        eddy_basis=eddy,
        coupling=np.array([[0.5, -0.25]]),
        eddy_system=eddy_system,
    )

    magnetic = np.array([[4.0]])
    expected_eddy = np.diag([3.0, 3.5])
    expected = np.array([
        [4.0, 0.5, -0.25],
        [0.5, 3.0, 0.0],
        [-0.25, 0.0, 3.5],
    ])

    np.testing.assert_allclose(coupled.eddy_impedance(2.0), expected_eddy)
    np.testing.assert_allclose(coupled.mixed_operator(magnetic, 2.0), expected)

    rhs = np.array([1.0, 0.2, -0.1])
    solved = coupled.solve(
        magnetic,
        2.0,
        magnetic_rhs=np.array([rhs[0]]),
        eddy_rhs=rhs[1:],
        return_operator=True,
    )
    np.testing.assert_allclose(solved["operator"], expected)
    np.testing.assert_allclose(
        solved["solution"].ravel(),
        np.linalg.solve(expected, rhs),
    )

    k = np.array([[0.5, -0.25]])
    np.testing.assert_allclose(
        coupled.schur_magnetic_operator(magnetic, 2.0),
        magnetic - k @ np.linalg.solve(expected_eddy, k.T),
    )
    np.testing.assert_allclose(
        coupled.schur_eddy_operator(magnetic, 2.0),
        expected_eddy - k.T @ np.linalg.solve(magnetic, k),
    )

    eddy_ports = np.array([[1.0, 0.0], [0.0, 1.0]])
    ports = np.vstack([np.zeros((1, 2)), eddy_ports])
    np.testing.assert_allclose(
        coupled.port_admittance(magnetic, 2.0, eddy_ports),
        ports.T @ np.linalg.solve(expected, ports),
    )


def test_hybrid_vim_named_block_schur_matches_igte_mixed_galerkin_formula():
    system = vim.HybridVIMSystem(
        resistance=np.array([
            [4.0, 0.5],
            [0.5, 3.0],
        ]),
        inductance=np.array([
            [1.0, 0.25],
            [0.25, 0.75],
        ]),
        surface_mass=np.array([
            [0.0, 0.0],
            [0.0, 2.0],
        ]),
        basis_names=("bulk", "surface"),
        blocks={"bulk": (0, 1), "surface": (1, 2)},
    )
    s = 2.0
    zs = 5.0
    z = system.impedance(s, surface_impedance=zs)

    np.testing.assert_allclose(system.block_matrix("bulk", "surface", s, surface_impedance=zs), z[:1, 1:])
    np.testing.assert_allclose(
        system.schur_complement("surface", "bulk", s, surface_impedance=zs),
        z[1:, 1:] - z[1:, :1] @ np.linalg.solve(z[:1, :1], z[:1, 1:]),
    )
    np.testing.assert_allclose(
        system.schur_complement("bulk", "surface", s, surface_impedance=zs),
        z[:1, :1] - z[:1, 1:] @ np.linalg.solve(z[1:, 1:], z[1:, :1]),
    )
    orthogonalized = system.mixed_galerkin_orthogonalization(
        "surface",
        "bulk",
        s,
        surface_impedance=zs,
    )
    assert isinstance(orthogonalized, vim.MixedGalerkinOrthogonalization)
    np.testing.assert_allclose(
        orthogonalized.reduced_operator,
        system.schur_complement("surface", "bulk", s, surface_impedance=zs),
    )
    np.testing.assert_allclose(
        orthogonalized.trial_transform,
        orthogonalized.test_transform,
    )
    rhs = np.array([[1.5, -0.2], [0.7, 1.1]])
    np.testing.assert_allclose(
        orthogonalized.solve(z, rhs),
        np.linalg.solve(z, rhs),
    )
    orthogonal_info = orthogonalized.diagnostics()
    assert orthogonal_info["retained_modes"] == 1
    assert orthogonal_info["eliminated_modes"] == 1
    assert orthogonal_info["trial_orthogonality_relative_defect"] < 1.0e-15
    assert orthogonal_info["test_orthogonality_relative_defect"] < 1.0e-15
    assert orthogonal_info["schur_relative_error"] < 1.0e-15
    np.testing.assert_allclose(
        system.block_rhs(
            bulk=np.array([[1.0, 2.0]]),
            surface=np.array([[3.0, 4.0]]),
        ),
        np.array([[1.0, 2.0], [3.0, 4.0]]),
    )
    np.testing.assert_allclose(
        system.block_rhs(surface=np.array([5.0])),
        np.array([0.0, 5.0]),
    )
    with pytest.raises(KeyError):
        system.block_slice("missing")
    with pytest.raises(ValueError):
        system.block_rhs(bulk=np.ones((2, 1)))
    with pytest.raises(ValueError):
        system.block_rhs(bulk=np.ones((1, 1)), surface=np.ones((1, 2)))


def test_mixed_galerkin_uses_distinct_trial_and_test_for_nonsymmetric_operator():
    system = vim.HybridVIMSystem(
        resistance=np.array([
            [4.0, 1.0],
            [2.0, 3.0],
        ]),
        inductance=np.zeros((2, 2)),
        surface_mass=np.zeros((2, 2)),
        basis_names=("bulk", "surface"),
        blocks={"bulk": (0, 1), "surface": (1, 2)},
    )

    orthogonalized = system.mixed_galerkin_orthogonalization(
        "surface",
        "bulk",
        0.0,
    )

    np.testing.assert_allclose(
        orthogonalized.trial_transform,
        np.array([[-0.25], [1.0]]),
    )
    np.testing.assert_allclose(
        orthogonalized.test_transform,
        np.array([[-0.5], [1.0]]),
    )
    np.testing.assert_allclose(orthogonalized.reduced_operator, np.array([[2.5]]))
    info = orthogonalized.diagnostics()
    assert info["trial_test_relative_difference"] > 0.0
    assert info["trial_orthogonality_relative_defect"] < 1.0e-15
    assert info["test_orthogonality_relative_defect"] < 1.0e-15
    assert info["schur_relative_error"] < 1.0e-15

    with pytest.raises(ValueError, match="operator does not match"):
        orthogonalized.solve(system.impedance(0.0) + np.eye(2), np.ones(2))


def test_hybrid_vim_multi_block_schur_eliminates_evrs_keeps_bridge_and_surface():
    r = np.diag([4.0, 2.0, 3.0, 1.5])
    l = np.array([
        [2.0, 0.1, 0.2, -0.1],
        [0.1, 3.0, 0.25, 0.15],
        [0.2, 0.25, 2.5, -0.05],
        [-0.1, 0.15, -0.05, 1.25],
    ])
    sm = np.diag([0.0, 0.0, 0.0, 2.0])
    system = vim.HybridVIMSystem(
        resistance=r,
        inductance=l,
        surface_mass=sm,
        basis_names=("evrs", "bridge0", "bridge1", "surface"),
        blocks={"volume": (0, 1), "volume1": (1, 3), "surface": (3, 4)},
    )
    s = 0.2 + 1.0j
    zs = 0.7 + 0.3j
    z = system.impedance(s, surface_impedance=zs)
    keep = system.block_indices(("volume1", "surface"))
    elim = system.block_indices("volume")

    np.testing.assert_array_equal(keep, np.array([1, 2, 3]))
    np.testing.assert_allclose(
        system.block_matrix_blocks(("volume1", "surface"), "volume", s, surface_impedance=zs),
        z[np.ix_(keep, elim)],
    )
    np.testing.assert_allclose(
        system.schur_complement_blocks(("volume1", "surface"), "volume", s, surface_impedance=zs),
        z[np.ix_(keep, keep)] - z[np.ix_(keep, elim)] @ np.linalg.solve(
            z[np.ix_(elim, elim)],
            z[np.ix_(elim, keep)],
        ),
    )
    orthogonalized = system.mixed_galerkin_orthogonalization(
        ("volume1", "surface"),
        "volume",
        s,
        surface_impedance=zs,
    )
    np.testing.assert_allclose(
        orthogonalized.reduced_operator,
        system.schur_complement_blocks(
            ("volume1", "surface"),
            "volume",
            s,
            surface_impedance=zs,
        ),
    )
    np.testing.assert_allclose(
        system.schur_complement_blocks("surface", ("volume", "volume1"), s, surface_impedance=zs),
        z[3:, 3:] - z[3:, :3] @ np.linalg.solve(z[:3, :3], z[:3, 3:]),
    )
    with pytest.raises(ValueError):
        system.block_indices(("volume", "volume"))
    with pytest.raises(ValueError):
        system.schur_complement_blocks(("volume", "surface"), "surface", s)


def test_sibc_tail_and_schur_termination_match_digest_asymptote():
    s = 1j * 2.0 * np.pi * 1.0e6
    sigma = 5.8e7
    surface = 2.0 * np.pi * 0.005
    k_sibc = surface * np.sqrt(sigma / vim.MU0)

    np.testing.assert_allclose(
        vim.SIBCAdmittanceTail(s, surface, sigma),
        k_sibc / np.sqrt(s),
    )

    d = 2.0 * np.pi * 1.0e3
    z = vim.SIBCSchurTerminationImpedance(s, k_sibc, d=d)
    y = vim.SIBCSchurTerminationAdmittance(s, k_sibc, d=d)
    np.testing.assert_allclose(y, 1.0 / z)
    np.testing.assert_allclose(
        y / (k_sibc / np.sqrt(s)),
        s / (s + d),
    )


def test_cpp_hybrid_vim_schur_kernel_matches_numpy_when_available():
    radia_cpp = pytest.importorskip("radia._radia_pybind")
    func = getattr(radia_cpp, "_HybridVIMSchurComplement", None)
    if func is None:
        pytest.skip("_HybridVIMSchurComplement is not available in this binary")

    kk = np.array([[4.0 + 0.5j, 0.25], [0.25, 3.0 + 0.1j]], dtype=np.complex128)
    ke = np.array([[0.5 - 0.2j], [-0.25 + 0.1j]], dtype=np.complex128)
    ek = ke.conj().T
    ee = np.array([[2.0 + 1.0j]], dtype=np.complex128)

    np.testing.assert_allclose(
        func(kk, ke, ek, ee),
        kk - ke @ np.linalg.solve(ee, ek),
    )


def test_cpp_hybrid_vim_solve_kernel_matches_numpy_when_available():
    radia_cpp = pytest.importorskip("radia._radia_pybind")
    func = getattr(radia_cpp, "_HybridVIMSolve", None)
    if func is None:
        pytest.skip("_HybridVIMSolve is not available in this binary")

    matrix = np.array(
        [[4.0 + 0.5j, 0.25 - 0.1j], [0.25 + 0.1j, 3.0 + 0.2j]],
        dtype=np.complex128,
    )
    rhs = np.array([[1.0, 0.0 - 0.5j], [0.25j, 2.0]], dtype=np.complex128)
    np.testing.assert_allclose(func(matrix, rhs), np.linalg.solve(matrix, rhs))

    scale = 1.0e-40
    scaled_matrix = scale * np.eye(2, dtype=np.complex128)
    scaled_rhs = scale * np.array([[1.0], [2.0]], dtype=np.complex128)
    np.testing.assert_allclose(func(scaled_matrix, scaled_rhs), [[1.0], [2.0]])


def test_sibc_dc_limits_fail_loud_instead_of_returning_nan():
    sigma = 5.8e7
    assert vim.SkinImpedance(0.0j, sigma) == 0.0j
    with pytest.raises(ValueError, match="undefined at s=0"):
        vim.SIBCAdmittanceTail(0.0j, 1.0, sigma)
    assert vim.SIBCSchurTerminationImpedance(0.0j, 1.0, d=0.0) == 0.0j
    with pytest.raises(ValueError, match="pole at s=0"):
        vim.SIBCSchurTerminationImpedance(0.0j, 1.0, d=1.0)
    assert vim.SIBCSchurTerminationAdmittance(0.0j, 1.0, d=1.0) == 0.0j
    with pytest.raises(ValueError, match="pole at s=0"):
        vim.SIBCSchurTerminationAdmittance(0.0j, 1.0, d=0.0)

    radia_cpp = pytest.importorskip("radia._radia_pybind")
    assert radia_cpp._SkinImpedance(0.0j, sigma, vim.MU0) == 0.0j
    with pytest.raises(RuntimeError, match="undefined at s=0"):
        radia_cpp._SIBCAdmittanceTail(0.0j, 1.0, sigma, vim.MU0)
    assert radia_cpp._SIBCSchurTerminationAdmittance(0.0j, 1.0, 1.0) == 0.0j


@pytest.mark.parametrize("bad", [np.nan, np.inf, complex(1.0, np.nan)])
def test_sibc_helpers_reject_nonfinite_laplace_frequency(bad):
    with pytest.raises(ValueError, match="s must be finite"):
        vim.SkinImpedance(bad, 5.8e7)


def test_hybrid_vim_public_names_are_exported():
    for name in (
        "SkinDepth",
        "EddySIBCApplicability",
        "SampledMagnetizationBasis",
        "VolumeCurrentBasis",
        "MagnetizationBasis",
        "SurfaceOmegaBasis",
        "EddyTracePolynomialDim",
        "EddyParentOrderLedger",
        "EddyFaceTopology",
        "EddyConductorGraphEdge",
        "EddyConductorCycle",
        "EddyConductorGraph",
        "EddyMeshTopology",
        "EddyDofPolicy",
        "EddyReductionPlan",
        "EddyBubbleDecomposition",
        "EddyBubbleHCurlBasis",
        "EddyBubbleReduction",
        "ClassifyNgsolveEddyTopology",
        "NgsolveEddyDofPolicy",
        "NgsolveEddyBubbleReduction",
        "NgsolveEddyBubbleHCurlBasis",
        "NgsolveBridgeCycleCurrentBasis",
        "NgsolveVolumeCurrentBasis",
        "NgsolveMagnetizationBasis",
        "NgsolveHDivMagnetizationBasis",
        "HDivMultipolePortSet",
        "PlanarHarmonicPortSet",
        "NgsolveHDivRegularSolidHarmonicPorts",
        "NgsolvePlanarHarmonicPorts",
        "NgsolveHDivExternalFieldRHS",
        "HDivMMMReducedModel",
        "NgsolveHDivMMMReduction",
        "NgsolveHDivMMMResponseReduction",
        "NgsolveBDMHDivMMMResponseReduction",
        "PlanarHDivMMMReducedSolution",
        "PlanarHDivMMMReducedModel",
        "NgsolvePlanarHDivMMMResponseReduction",
        "NgsolveHCurlCurlBasis",
        "NgsolveSurfaceOmegaBasis",
        "NgsolveMatrixToDense",
        "NgsolveVectorToArray",
        "NgsolveCouplingDofMasks",
        "NgsolveBlockKrylovBasis",
        "NgsolveOperatorBlockKrylovBasis",
        "NgsolveStaticCondensedBlockKrylovBasis",
        "EVRSBasis",
        "CompressHCurlResponseInCurrentGram",
        "BlockKrylovBasis",
        "SampledLaplaceInteraction",
        "HACApKSampledLaplaceInteraction",
        "HACApKSampledPlanarLogInteraction",
        "SampledHACApKOperator",
        "ReducedInteractionMatrix",
        "NGSolveProjectedOperator",
        "NGSolveProjectedInteraction",
        "CoupledReducedOperator",
        "CurrentMagneticFluxDensitySamples",
        "MagnetizationCurrentCoupling",
        "EVRSTMethodAlgebra",
        "ReducedPortAdmittance",
        "ReducedPortImpedance",
        "HCurlEddyCLNModel",
        "HCurlEddyCLNFromVIM",
        "SharedMeshMaterialModel",
        "CoupledHDivEVRSSystem",
        "CoupledHDivHybridVIMSystem",
        "HCurlVIMHDivMMMSolution",
        "HCurlVIMHDivMMMSystem",
        "MixedGalerkinHDivHybridVIMSystem",
        "CoupleHDivMagnetizationToEVRS",
        "CoupleHCurlVIMWithHDivMMM",
        "CoupleHybridVIMWithHDivMMM",
        "CoupleEddyBubbleHCurlBasisWithHDivMMM",
        "AssembleHybridVIM",
        "LocalESIMSurfaceLUT",
        "BuildLocalESIMSurfaceLUT",
        "ValidateLocalESIMSurfaceLUT",
        "TopologyAwareHybridVIM",
        "NgsolveTopologyAwareHybridVIM",
        "NgsolveEddyBubbleHybridVIM",
        "NgsolveHCurlVIMHDivMMM",
        "NgsolveBDMEddyBubbleVIM",
        "SkinImpedance",
        "SIBCAdmittanceTail",
        "SIBCSchurTerminationImpedance",
        "SIBCSchurTerminationAdmittance",
        "ExternalVectorPotentialRHS",
    ):
        assert callable(getattr(vim, name))
        assert name in vim.__all__


def test_eddy_parent_order_ledger_keeps_p_symbolic_not_empirical():
    assert vim.EddyTracePolynomialDim(0) == 1
    assert vim.EddyTracePolynomialDim(2) == 6
    assert vim.EddyTracePolynomialDim(2, face_family="tensor") == 9

    ledger = vim.EddyParentOrderLedger(
        bulk_degree=4,
        bridge_trace_degree=0,
        surface_current_degree=2,
    )
    assert ledger.required_parent_order == 4
    assert ledger.surface_omega_degree == 3
    assert ledger.bridge_trace_dim == 1
    assert not ledger.is_parent_order_admissible(3)
    assert ledger.is_parent_order_admissible(4)
    assert ledger.is_parent_order_admissible(6)

    info = ledger.diagnostics(
        parent_order=6,
        evrs_rank=24,
        cycle_rank=7,
        surface_modes=3,
    )
    assert info["parent_order_admissible"] is True
    assert info["parent_order_excess"] == 2
    assert info["estimated_reduced_modes"] == 24 + 7 + 3

    enriched_bridge = vim.EddyParentOrderLedger(
        bulk_degree=4,
        bridge_trace_degree=2,
        surface_current_degree=2,
    )
    assert enriched_bridge.required_parent_order == 4
    assert enriched_bridge.bridge_modes(cycle_rank=7) == 7 * 6
    assert (
        enriched_bridge.estimated_reduced_modes(
            evrs_rank=24,
            cycle_rank=7,
            surface_modes=3,
            non_sibc_trace_modes=5,
        )
        == 24 + 7 * 6 + 3 + 5
    )

    high_order_bridge = vim.EddyParentOrderLedger(
        bulk_degree=4,
        bridge_trace_degree=6,
        surface_current_degree=2,
    )
    assert high_order_bridge.required_parent_order == 6

    with pytest.raises(ValueError):
        vim.EddyTracePolynomialDim(-1)
    with pytest.raises(ValueError):
        vim.EddyTracePolynomialDim(1, face_family="hexagonal")


def test_ngsolve_sampling_bridge_builds_hybrid_vim_on_unit_box():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))

    volume = vim.NgsolveVolumeCurrentBasis(
        mesh,
        (
            ng.CoefficientFunction((1.0, 0.0, 0.0)),
            ng.CoefficientFunction((0.0, 1.0, 0.0)),
        ),
        intorder=1,
        materials="cond",
        names=["jx", "jy"],
    )
    surface = vim.NgsolveSurfaceOmegaBasis(
        mesh,
        (ng.CoefficientFunction((0.0, 1.0, 0.0)),),
        intorder=1,
        boundaries="skin",
        names=["omega_y"],
    )

    assert volume.n_samples > 0
    assert surface.n_samples > 0
    np.testing.assert_allclose(volume.weights.sum(), 1.0, rtol=1.0e-12)
    np.testing.assert_allclose(volume.mass_matrix(), np.eye(2), atol=1.0e-12)
    assert surface.mass_matrix()[0, 0].real > 0.0

    system = vim.AssembleHybridVIM(volume, surface, sigma=5.8e7, kernel_epsilon=0.1)
    assert system.impedance(1j * 100.0, surface_impedance=vim.SkinImpedance(1j * 100.0, 5.8e7)).shape == (3, 3)


def test_ngsolve_eddy_topology_classifies_skin_faces_and_loop_bridges():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.8))

    topology = vim.ClassifyNgsolveEddyTopology(mesh, conductive_materials="cond")
    info = topology.diagnostics()

    assert isinstance(topology, vim.EddyMeshTopology)
    assert info["conductive_element_count"] > 0
    assert info["conductive_component_count"] == 1
    assert info["surface_face_count"] > 0
    assert info["loop_bridge_face_count"] > 0
    assert info["conductor_exterior_face_count"] == info["surface_face_count"]
    assert info["conductor_conductor_face_count"] == info["loop_bridge_face_count"]
    assert "conductor-exterior" in info["roles"]
    assert "conductor-conductor" in info["roles"]
    assert all(face.requires_surface_basis for face in topology.surface_faces)
    assert all(face.requires_loop_bridge for face in topology.loop_bridge_faces)
    graph = topology.conductor_graph()
    graph_info = graph.diagnostics()
    assert isinstance(graph, vim.EddyConductorGraph)
    assert graph_info["node_count"] == info["conductive_element_count"]
    assert graph_info["edge_count"] == info["loop_bridge_face_count"]
    assert graph_info["cycle_rank"] == info["conductive_graph_cycle_rank"]
    assert graph_info["fundamental_cycle_count"] == graph_info["cycle_rank"]
    assert all(isinstance(cycle, vim.EddyConductorCycle) for cycle in graph.cycle_basis())
    cycle_edges = graph.cycle_edge_matrix()
    assert cycle_edges.shape == (graph_info["cycle_rank"], graph_info["edge_count"])
    assert np.all(np.isin(cycle_edges, (-1.0, 0.0, 1.0)))
    assert np.all(np.count_nonzero(cycle_edges, axis=1) >= 3)


def test_eddy_topology_keeps_insulator_faces_out_of_sibc_set():
    topology = vim.EddyMeshTopology(
        faces=(
            vim.EddyFaceTopology(
                face_nr=1,
                role="conductor-air",
                volume_elements=(1, 2),
                volume_materials=("cond", "air"),
            ),
            vim.EddyFaceTopology(
                face_nr=2,
                role="conductor-exterior",
                volume_elements=(1,),
                volume_materials=("cond",),
                boundary_labels=("skin",),
            ),
            vim.EddyFaceTopology(
                face_nr=3,
                role="conductor-insulator",
                volume_elements=(1, 3),
                volume_materials=("cond", "ceramic"),
            ),
            vim.EddyFaceTopology(
                face_nr=4,
                role="conductor-conductor",
                volume_elements=(1, 4),
                volume_materials=("cond", "cond"),
            ),
        ),
        conductive_materials=("cond",),
        air_materials=("air", "vacuum"),
    )

    assert [face.face_nr for face in topology.surface_faces] == [1, 2, 3]
    assert [face.face_nr for face in topology.sibc_faces] == [1, 2]
    assert [face.face_nr for face in topology.non_sibc_trace_faces] == [3]
    assert [face.face_nr for face in topology.loop_bridge_faces] == [4]
    assert not topology.faces_by_role("conductor-insulator")[0].is_sibc_face
    assert not topology.faces_by_role("conductor-insulator")[0].can_sibc_terminate
    info = topology.diagnostics()
    assert info["surface_face_count"] == 3
    assert info["sibc_face_count"] == 2
    assert info["non_sibc_trace_face_count"] == 1
    assert info["conductor_insulator_face_count"] == 1


def test_eddy_dof_policy_excludes_non_sibc_trace_from_ordinary_evrs():
    free = np.ones(8, dtype=bool)
    sibc = np.zeros(8, dtype=bool)
    surface = np.zeros(8, dtype=bool)
    bridge = np.zeros(8, dtype=bool)
    sibc[[0, 1]] = True
    surface[[0, 1, 2, 3]] = True
    bridge[[4, 5]] = True

    policy = vim.EddyDofPolicy(
        free=free,
        sibc_surface=sibc,
        surface_candidate=surface,
        loop_bridge=bridge,
        local_bubble=np.zeros(8, dtype=bool),
        interface=np.zeros(8, dtype=bool),
        wirebasket=np.zeros(8, dtype=bool),
    )

    np.testing.assert_array_equal(
        policy.non_sibc_trace,
        np.array([False, False, True, True, False, False, False, False]),
    )
    np.testing.assert_array_equal(
        policy.ordinary_evrs_candidate,
        np.array([False, False, False, False, False, False, True, True]),
    )
    info = policy.diagnostics()
    assert info["sibc_surface_only_dofs"] == 2
    assert info["non_sibc_trace_dofs"] == 2
    assert info["loop_bridge_dofs"] == 2
    assert info["ordinary_evrs_candidate_dofs"] == 2
    assert info["partitioned_free_dofs"] == 8

    plan = policy.reduction_plan(evrs_rank=3, surface_modes=2, loop_bridge_modes=1)
    plan_info = plan.diagnostics()
    assert plan_info["non_sibc_trace_dofs"] == 2
    assert plan_info["non_sibc_trace_modes"] == 2
    assert plan_info["ordinary_evrs_candidate_dofs"] == 2
    assert plan_info["estimated_reduced_modes"] == 1 + 2 + 3 + 2

    bubbling = vim.EddyBubbleReduction(
        policy,
        evrs_rank=3,
        surface_modes=2,
        loop_bridge_modes=1,
        bridge_strategy="cycle-basis",
        parent_order=4,
        parent_order_ledger=vim.EddyParentOrderLedger(
            bulk_degree=3,
            bridge_trace_degree=0,
            surface_current_degree=1,
        ),
    )
    bubble_info = bubbling.diagnostics()
    assert isinstance(bubbling, vim.EddyBubbleDecomposition)
    assert bubble_info["rule"] == "topology-aware-eddy-bubbling"
    assert bubble_info["classes"]["sibc_surface"] == 2
    assert bubble_info["classes"]["non_sibc_trace"] == 2
    assert bubble_info["classes"]["loop_bridge"] == 2
    assert bubble_info["classes"]["ordinary_bulk_eddy_bubble"] == 2
    assert bubble_info["parent_order"] == 4
    assert bubble_info["parent_order_ledger"]["required_parent_order"] == 3
    assert bubble_info["parent_order_ledger"]["parent_order_excess"] == 1
    np.testing.assert_array_equal(bubbling.eddy_bubble_candidate, policy.ordinary_evrs_candidate)
    np.testing.assert_array_equal(
        bubbling.structural_keep,
        policy.topology_protected | policy.sibc_surface_only | policy.non_sibc_trace,
    )


def test_ngsolve_eddy_dof_policy_marks_sibc_surface_and_loop_bridge_dofs():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.8))
    fes = ng.HCurl(mesh, order=3, nograds=True)
    topology = vim.ClassifyNgsolveEddyTopology(mesh, conductive_materials="cond")
    policy = vim.NgsolveEddyDofPolicy(mesh, fes, topology)
    info = policy.diagnostics()

    assert isinstance(policy, vim.EddyDofPolicy)
    assert info["free_dofs"] == fes.ndof
    assert info["sibc_surface_dofs"] > 0
    assert info["loop_bridge_dofs"] > 0
    assert info["ordinary_evrs_candidate_dofs"] > 0
    assert info["partitioned_free_dofs"] == info["free_dofs"]
    assert np.count_nonzero(policy.topology_protected & policy.sibc_surface_only) == 0
    assert np.count_nonzero(policy.topology_protected & policy.ordinary_evrs_candidate) == 0
    assert np.count_nonzero(policy.sibc_surface_only & policy.ordinary_evrs_candidate) == 0
    np.testing.assert_array_equal(
        policy.free,
        (
            policy.topology_protected
            | policy.sibc_surface_only
            | policy.non_sibc_trace
            | policy.ordinary_evrs_candidate
        ),
    )

    plan = policy.reduction_plan(evrs_rank=12, surface_modes=3)
    plan_info = plan.diagnostics()
    assert isinstance(plan, vim.EddyReductionPlan)
    assert plan_info["loop_bridge_keep_dofs"] == info["loop_bridge_dofs"]
    assert plan_info["sibc_surface_trace_dofs"] == info["sibc_surface_dofs"]
    assert plan_info["ordinary_evrs_candidate_dofs"] == info["ordinary_evrs_candidate_dofs"]
    assert plan_info["estimated_reduced_modes"] == info["loop_bridge_dofs"] + 12 + 3
    assert plan_info["estimated_reduction_ratio"] < 1.0

    graph = topology.conductor_graph()
    cycle_plan = policy.reduction_plan(
        evrs_rank=12,
        surface_modes=3,
        loop_bridge_modes=graph.cycle_rank,
        bridge_strategy="cycle-basis",
    )
    cycle_info = cycle_plan.diagnostics()
    assert cycle_info["loop_bridge_reduced_modes"] == graph.cycle_rank
    assert cycle_info["loop_bridge_reduction_strategy"] == "cycle-basis"
    assert cycle_info["estimated_reduced_modes"] == graph.cycle_rank + 12 + 3
    assert cycle_info["estimated_reduced_modes"] < plan_info["estimated_reduced_modes"]

    bubbling = vim.NgsolveEddyBubbleReduction(
        mesh,
        fes,
        topology,
        evrs_rank=12,
        surface_modes=3,
        loop_bridge_modes=graph.cycle_rank,
        bridge_strategy="cycle-basis",
        parent_order_ledger=vim.EddyParentOrderLedger(
            bulk_degree=2,
            bridge_trace_degree=0,
            surface_current_degree=1,
        ),
    )
    bubble_info = bubbling.diagnostics()
    assert bubble_info["conductor_graph"]["cycle_rank"] == graph.cycle_rank
    assert bubble_info["classes"]["ordinary_bulk_eddy_bubble"] == info["ordinary_evrs_candidate_dofs"]
    assert bubble_info["plan"]["estimated_reduced_modes"] == graph.cycle_rank + 12 + 3
    assert bubble_info["parent_order"] == 3
    assert bubble_info["parent_order_ledger"]["parent_order_admissible"] is True


def test_ngsolve_bridge_cycle_current_basis_feeds_hybrid_vim():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.8))
    topology = vim.ClassifyNgsolveEddyTopology(mesh, conductive_materials="cond")
    graph = topology.conductor_graph()

    bridge = vim.NgsolveBridgeCycleCurrentBasis(mesh, topology)

    assert bridge.kind == "volume"
    assert bridge.n_modes == graph.cycle_rank
    assert bridge.n_samples == graph.edge_count
    assert np.all(bridge.weights > 0.0)
    assert np.count_nonzero(bridge.modes) > 0
    assert bridge.mass_matrix().shape == (graph.cycle_rank, graph.cycle_rank)
    assert np.min(np.linalg.eigvalsh(bridge.mass_matrix().real)) > 0.0

    system = vim.AssembleHybridVIM(bridge, sigma=5.8e7, kernel_epsilon=0.1)
    assert system.n_modes == graph.cycle_rank
    assert system.diagnostics()["passive_blocks"] is True


def test_ngsolve_curved_sphere_preserves_surface_and_bridge_geometry():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    sphere = occ.Sphere(occ.Pnt(0, 0, 0), 1.0)
    sphere.mat("cond")
    for face in sphere.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(sphere).GenerateMesh(maxh=0.75))
    mesh.Curve(3)

    surface = vim.NgsolveSurfaceOmegaBasis(
        mesh,
        (ng.CoefficientFunction((1.0, 0.0, 0.0)),),
        intorder=8,
        boundaries="skin",
    )
    points, weights, normals = vim.SampleNgsolveVectorCFs(
        mesh,
        (ng.specialcf.normal(3),),
        vb="BND",
        intorder=8,
        boundaries="skin",
    )
    fem_area = ng.Integrate(
        1.0,
        mesh,
        definedon=mesh.Boundaries("skin"),
        order=8,
    )

    np.testing.assert_allclose(surface.points, points, atol=1.0e-14)
    np.testing.assert_allclose(surface.weights, weights, rtol=1.0e-13)
    assert surface.weights.sum() == pytest.approx(fem_area, rel=1.0e-12)
    assert surface.weights.sum() == pytest.approx(4.0 * np.pi, rel=5.0e-4)
    tangential_defect = np.max(
        np.abs(np.einsum("ij,ij->i", surface.modes[0], normals[0]))
    )
    assert tangential_defect < 1.0e-12

    topology = vim.ClassifyNgsolveEddyTopology(mesh, conductive_materials="cond")
    bridge = vim.NgsolveBridgeCycleCurrentBasis(
        mesh,
        topology,
        geometry_intorder=8,
    )
    assert bridge.n_modes == topology.conductor_graph().cycle_rank
    assert np.all(np.isfinite(bridge.points))
    assert np.all(bridge.weights > 0.0)
    assert np.min(np.linalg.eigvalsh(bridge.mass_matrix().real)) > 0.0


def test_ngsolve_topology_aware_hybrid_vim_builder_returns_tri_block_system():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.8))
    fes = ng.HCurl(mesh, order=2, nograds=True)
    vectors = np.zeros((fes.ndof, 2))
    vectors[0, 0] = 1.0
    vectors[1, 1] = 1.0

    built = vim.NgsolveTopologyAwareHybridVIM(
        mesh,
        fes,
        vectors,
        (
            ng.CoefficientFunction((1.0, 0.0, 0.0)),
            ng.CoefficientFunction((0.0, 1.0, 0.0)),
        ),
        sigma=5.8e7,
        conductive_materials="cond",
        surface_boundaries="skin",
        intorder=1,
        geometry_intorder=4,
        kernel_epsilon=0.1,
        port_vector_potentials=(
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
        ),
    )

    assert isinstance(built, vim.TopologyAwareHybridVIM)
    assert built.system.blocks.keys() == {"volume", "volume1", "surface"}
    assert built.volume_basis.n_modes == 2
    assert built.bridge_cycle_basis.n_modes == built.conductor_graph.cycle_rank
    np.testing.assert_allclose(
        np.diag(built.bridge_cycle_basis.mass_matrix()),
        np.ones(built.bridge_cycle_basis.n_modes),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        np.diag(built.surface_basis.mass_matrix()),
        np.ones(built.surface_basis.n_modes),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    assert built.surface_basis.n_modes == 2
    assert built.system.n_modes == 2 + built.conductor_graph.cycle_rank + 2
    assert built.rhs.shape == (built.system.n_modes, 2)
    info = built.diagnostics()
    assert built.system.interaction_backend == "hacapk-sampled-laplace"
    assert info["system"]["inductance_matrix_free"] is True
    assert info["system"]["passive_blocks"] is True
    assert info["reduction_plan"]["loop_bridge_reduction_strategy"] == "cycle-basis"
    assert info["reduction_plan"]["estimated_reduced_modes"] == built.system.n_modes

    fes_hdiv = ng.HDiv(mesh, order=1)
    hdiv_vectors = np.zeros((fes_hdiv.ndof, 1))
    hdiv_vectors[0, 0] = 1.0
    magnetization = vim.NgsolveHDivMagnetizationBasis(
        mesh,
        fes_hdiv,
        hdiv_vectors,
        intorder=1,
        materials="cond",
        names=["M0"],
    )
    coupled = built.couple_hdiv_mmm(magnetization, kernel_epsilon=0.1)
    assert isinstance(coupled, vim.CoupledHDivHybridVIMSystem)
    assert coupled.n_hcurl_vim_modes == built.system.n_modes
    assert coupled.diagnostics()["eddy_basis_count"] == 3


def test_hybrid_vim_dense_solve_is_invariant_to_reduced_basis_scale():
    scale = 1.0e-14
    correlation = 0.1 * np.sqrt(scale)
    resistance = np.array(
        [[1.0, correlation], [correlation, scale]],
        dtype=complex,
    )
    system = vim.HybridVIMSystem(
        resistance=resistance,
        inductance=np.zeros((2, 2), dtype=complex),
        surface_mass=np.zeros((2, 2), dtype=complex),
        basis_names=("magnetic_scale", "small_dual_volume_scale"),
        blocks={"volume": (0, 2)},
        interaction_backend="dense-scale-invariance-test",
    )
    expected = np.array([1.25 - 0.5j, -2.0 + 0.75j])
    rhs = resistance @ expected

    solved = system.solve(0.0, rhs)

    np.testing.assert_allclose(solved, expected, rtol=1.0e-9, atol=1.0e-11)
    relative_residual = np.linalg.norm(resistance @ solved - rhs) / np.linalg.norm(rhs)
    assert relative_residual < 1.0e-14


def test_ngsolve_topology_aware_planar_builder_uses_log_hacapk_by_default():
    ng = pytest.importorskip("ngsolve")
    geom2d = pytest.importorskip("netgen.geom2d")

    geometry = geom2d.SplineGeometry()
    geometry.AddRectangle(
        (0.0, 0.0),
        (1.0, 1.0),
        bcs=("skin", "skin", "skin", "skin"),
        leftdomain=1,
    )
    geometry.SetMaterial(1, "cond")
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=0.7))
    fes = ng.HCurl(mesh, order=2, nograds=True)
    vectors = np.zeros((fes.ndof, 2))
    vectors[0, 0] = 1.0
    vectors[1, 1] = 1.0

    built = vim.NgsolveTopologyAwareHybridVIM(
        mesh,
        fes,
        vectors,
        (ng.CF((1.0, 0.0)),),
        sigma=1.0e6,
        conductive_materials="cond",
        surface_boundaries="skin",
        intorder=2,
        kernel_epsilon=0.05,
    )

    assert built.system.interaction_backend == "hacapk-sampled-planar-log"
    assert built.system.diagnostics()["inductance_matrix_free"] is True
    assert built.bridge_cycle_basis.n_modes == built.conductor_graph.cycle_rank


def test_ngsolve_eddy_bubble_hcurl_basis_builder_feeds_vim_and_hdiv_mmm():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))

    fes = ng.HCurl(mesh, order=2, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    port = ng.LinearForm(fes)
    port += ng.CoefficientFunction((-ng.y, ng.x, 0.0)) * v * ng.dx

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()

    eddy_basis = vim.NgsolveEddyBubbleHCurlBasis(
        mesh,
        fes,
        stiffness,
        mass,
        port,
        steps=2,
        conductive_materials="cond",
        response_backend="dense",
        intorder=1,
        parent_order_ledger=vim.EddyParentOrderLedger(
            bulk_degree=2,
            bridge_trace_degree=0,
            surface_current_degree=1,
        ),
    )

    assert isinstance(eddy_basis, vim.EddyBubbleHCurlBasis)
    assert eddy_basis.current_basis.n_modes == eddy_basis.response_basis.rank
    assert eddy_basis.eddy_bubbling.diagnostics()["plan"]["evrs_rank"] == eddy_basis.rank
    eddy_system = eddy_basis.assemble_vim(sigma=5.8e7, kernel_epsilon=0.1)
    assert eddy_system.n_modes == eddy_basis.n_modes

    fes_hdiv = ng.HDiv(mesh, order=1)
    hdiv_vectors = np.zeros((fes_hdiv.ndof, 1))
    hdiv_vectors[0, 0] = 1.0
    magnetization = vim.NgsolveHDivMagnetizationBasis(
        mesh,
        fes_hdiv,
        hdiv_vectors,
        intorder=1,
        materials="cond",
        names=["M0"],
    )
    coupled = eddy_basis.couple_hdiv_mmm(
        magnetization,
        eddy_system=eddy_system,
        kernel_epsilon=0.1,
    )
    assert isinstance(coupled, vim.HCurlVIMHDivMMMSystem)
    assert coupled.n_hcurl_vim_modes == eddy_basis.n_modes
    assert coupled.n_hdiv_mmm_modes == 1
    assert coupled.diagnostics()["has_eddy_system"] is True


def test_ngsolve_one_call_hcurl_vim_hdiv_mmm_builder_returns_mixed_system():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))

    fes = ng.HCurl(mesh, order=2, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    port = ng.LinearForm(fes)
    port += ng.CoefficientFunction((-ng.y, ng.x, 0.0)) * v * ng.dx

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()

    surface_modes = (ng.CoefficientFunction((1.0, 0.0, 0.0)),)
    hybrid = vim.NgsolveEddyBubbleHybridVIM(
        mesh,
        fes,
        stiffness,
        mass,
        port,
        surface_modes,
        steps=2,
        sigma=5.8e7,
        conductive_materials="cond",
        response_backend="dense",
        intorder=1,
        parent_order_ledger=vim.EddyParentOrderLedger(
            bulk_degree=2,
            bridge_trace_degree=0,
            surface_current_degree=1,
        ),
    )
    assert isinstance(hybrid, vim.TopologyAwareHybridVIM)
    assert hybrid.response_basis.diagnostics()["pre_current_gram_rank"] == 2
    assert hybrid.volume_basis.n_modes == 1
    np.testing.assert_allclose(
        hybrid.volume_basis.mass_matrix(),
        np.eye(1),
        atol=1.0e-12,
    )
    assert hybrid.surface_basis.n_modes == 1
    assert hybrid.system.n_modes == (
        hybrid.volume_basis.n_modes
        + hybrid.bridge_cycle_basis.n_modes
        + hybrid.surface_basis.n_modes
    )

    thin_conductor = vim.EddySIBCApplicability(
        frequency_hz=50.0,
        sigma=5.8e7,
        characteristic_thickness_m=1.0e-3,
    )
    volumetric_hybrid = vim.NgsolveEddyBubbleHybridVIM(
        mesh,
        fes,
        stiffness,
        mass,
        port,
        surface_modes,
        steps=2,
        sigma=5.8e7,
        conductive_materials="cond",
        response_backend="dense",
        intorder=1,
        sibc_applicability=thin_conductor,
    )
    assert volumetric_hybrid.surface_basis.n_modes == 0
    assert volumetric_hybrid.system.blocks["surface"][0] == (
        volumetric_hybrid.system.blocks["surface"][1]
    )
    assert volumetric_hybrid.diagnostics()["surface_model"] == "volumetric"
    assert volumetric_hybrid.diagnostics()["sibc_applicability"][
        "sibc_applicable"
    ] is False

    fes_hdiv = ng.HDiv(mesh, order=1)
    hdiv_vectors = np.zeros((fes_hdiv.ndof, 1))
    hdiv_vectors[0, 0] = 1.0
    with ng.TaskManager():
        hdiv_reduction = vim.NgsolveHDivMMMReduction(
            mesh,
            fes_hdiv,
            hdiv_vectors,
            mu_r=1001.0,
            external_fields=(ng.CoefficientFunction((1.0, 0.0, 0.0)),),
            intorder=1,
            materials="cond",
            names=["M0"],
            demag_eps=1.0e-10,
        )
    assert isinstance(hdiv_reduction, vim.HDivMMMReducedModel)
    np.testing.assert_allclose(
        hdiv_reduction.magnetic_operator,
        hdiv_reduction.mass / 1000.0 + hdiv_reduction.demag,
    )
    assert np.linalg.norm(hdiv_reduction.demag) > 0.0
    assert hdiv_reduction.magnetic_rhs.shape == (1, 1)
    swept_rhs = hdiv_reduction.external_field_rhs(
        ng.CoefficientFunction((0.0, 1.0, 0.0))
    )
    assert swept_rhs.shape == (1, 1)

    def swirl_vector_potential(points):
        points = np.asarray(points)
        return np.column_stack(
            (-points[:, 1], points[:, 0], np.zeros(points.shape[0]))
        )

    mixed = vim.NgsolveBDMEddyBubbleVIM(
        mesh,
        fes,
        stiffness,
        mass,
        port,
        surface_modes,
        hdiv_order=1,
        mu_r=1001.0,
        external_fields=(ng.CoefficientFunction((1.0, 0.0, 0.0)),),
        hdiv_max_modes=1,
        magnetic_materials="cond",
        steps=2,
        sigma=5.8e7,
        conductive_materials="cond",
        response_backend="dense",
        intorder=1,
        port_vector_potentials=(swirl_vector_potential,),
        coupling_kernel_epsilon=0.1,
    )
    assert isinstance(mixed, vim.CoupledHDivHybridVIMSystem)
    assert mixed.n_hdiv_mmm_modes == 1
    assert mixed.n_hcurl_vim_modes == hybrid.system.n_modes
    op = mixed.mixed_operator(None, 1j * 100.0, surface_impedance=0.1)
    assert op.shape == (1 + hybrid.system.n_modes, 1 + hybrid.system.n_modes)
    assert mixed.diagnostics()["has_eddy_rhs"] is True
    assert mixed.diagnostics()["has_response_basis"] is True
    assert mixed.diagnostics()["response_basis"]["current_gram_rank"] == (
        mixed.eddy_bases[0].n_modes
    )
    assert mixed.diagnostics()["has_hdiv_reduction"] is True
    assert mixed.diagnostics()["eddy_block_roles"] == {
        "volume": "bulk",
        "volume1": "bridge",
        "surface": "sibc",
    }
    assert mixed.hdiv_reduction.parent_family == "BDM"
    assert mixed.hdiv_reduction.parent_order == 1
    hdiv_info = mixed.hdiv_reduction.diagnostics()
    assert hdiv_info["demag_hmatrix_active"] is True
    assert "ChargeGramHMatrix" in hdiv_info["demag_hmatrix_backend"]

    solved = mixed.solve_frequency(100.0)
    assert isinstance(solved, vim.HCurlVIMHDivMMMSolution)
    assert solved.eddy_block_roles == ("bulk", "bridge", "sibc")
    assert solved.diagnostics()["eddy_block_roles"] == [
        "bulk",
        "bridge",
        "sibc",
    ]
    assert solved.parent_t_coefficients.shape == (fes.ndof, 1)
    assert solved.parent_magnetization_coefficients.shape == (
        mixed.hdiv_reduction.parent_ndof,
        1,
    )
    assert solved.sampled_magnetization.shape[0] == 1
    assert solved.current_samples("volume").shape[0] == 1
    assert solved.current_samples("surface").shape[0] == 1
    np.testing.assert_allclose(
        solved.current_samples("bulk"),
        solved.current_samples("volume"),
    )
    np.testing.assert_allclose(
        solved.current_samples("bridge"),
        solved.current_samples("volume1"),
    )
    np.testing.assert_allclose(
        solved.current_samples("sibc"),
        solved.current_samples("surface"),
    )
    assert solved.average_joule_loss[0] >= 0.0
    # 1e-8 (not 1e-10): the physical coupling scales (-K/mu0, s K^H) spread the
    # block magnitudes, costing ~1 digit of direct-solve residual (measured
    # 1.6e-10; the backward error stays < 1e-10).
    assert solved.residual_relative_norm < 1.0e-8
    assert mixed.adjacency_class_block_partition() == (
        ("volume1", "surface"),
        ("volume",),
    )
    solved_orthogonalized = mixed.solve_frequency_eddy_bubbled(100.0)
    assert solved_orthogonalized.solver_backend == (
        "native-dense-reduced-lu-mixed-galerkin"
    )
    assert solved_orthogonalized.solver_diagnostics["iterations"] == 0
    with pytest.raises(
        ValueError,
        match="mixed-Galerkin reduction currently requires solver='dense'",
    ):
        mixed.solve_frequency_eddy_bubbled(100.0, solver="gmres")
    np.testing.assert_allclose(
        solved_orthogonalized.magnetization_coefficients,
        solved.magnetization_coefficients,
        rtol=2.0e-5,
        atol=3.0e-8,
    )
    np.testing.assert_allclose(
        solved_orthogonalized.eddy_coefficients,
        solved.eddy_coefficients,
        rtol=2.0e-5,
        atol=3.0e-8,
    )
    np.testing.assert_allclose(
        solved_orthogonalized.parent_t_coefficients,
        solved.parent_t_coefficients,
        rtol=2.0e-5,
        atol=3.0e-8,
    )
    np.testing.assert_allclose(
        solved_orthogonalized.port_response,
        solved.port_response,
        rtol=2.0e-8,
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        solved_orthogonalized.average_joule_loss,
        solved.average_joule_loss,
        rtol=2.0e-8,
        atol=1.0e-10,
    )
    assert solved_orthogonalized.orthogonalized_solution is not None
    assert solved_orthogonalized.mixed_galerkin_diagnostics[
        "full_coupled_schur"
    ] is True
    assert solved_orthogonalized.mixed_galerkin_diagnostics[
        "keep_eddy_blocks"
    ] == ["volume1", "surface"]
    assert solved_orthogonalized.mixed_galerkin_diagnostics[
        "hdiv_demag_hmatrix_backend"
    ] == hdiv_info["demag_hmatrix_backend"]
    assert solved_orthogonalized.residual_relative_norm < 1.0e-10
    assert solved_orthogonalized.residual_backward_error < 1.0e-12
    target = np.array([[2.0, 2.0, 2.0]])
    total_field = solved.eddy_flux_density(target)
    block_field = sum(
        solved.eddy_flux_density(target, block=name)
        for name in solved.eddy_block_names
    )
    assert total_field.shape == (1, 1, 3)
    assert np.all(np.isfinite(total_field))
    np.testing.assert_allclose(total_field, block_field)
    np.testing.assert_allclose(
        solved.eddy_flux_density(target, block="bridge"),
        solved.eddy_flux_density(target, block="volume1"),
    )


def test_ngsolve_hdiv_mmm_response_reduction_uses_multipole_training_ports():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("body")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))
    ports = vim.NgsolveHDivRegularSolidHarmonicPorts(mesh, max_degree=2)

    assert isinstance(ports, vim.HDivMultipolePortSet)
    assert ports.count == 8
    assert ports.diagnostics()["degree_counts"] == {"1": 3, "2": 5}
    with ng.TaskManager():
        reduced = vim.NgsolveBDMHDivMMMResponseReduction(
            mesh,
            order=1,
            mu_r=1001.0,
            external_fields=(
                ng.CoefficientFunction((1.0, 0.0, 0.0)),
                ng.CoefficientFunction((0.0, 1.0, 0.0)),
            ),
            training_fields=ports,
            materials="body",
            max_modes=2,
            solve_tol=1.0e-9,
            demag_eps=1.0e-7,
        )

    generation = reduced.diagnostics()["basis_generation"]
    assert reduced.parent_ndof == reduced.fes.ndof
    assert reduced.diagnostics()["parent_family"] == "BDM"
    assert reduced.diagnostics()["parent_order"] == 1
    assert reduced.n_modes == 2
    assert reduced.magnetic_rhs.shape == (2, 2)
    assert generation["construction"] == "hdiv-mmm-response-energy-pod"
    assert generation["snapshot_backend"] == "radia-cpp-mass-riesz-cg"
    assert generation["snapshot_port_count"] == 10
    assert generation["physical_rhs_columns"] == 2
    assert generation["training_port_count"] == 8
    assert generation["protected_physical_modes"] == 2
    assert generation["training_response_modes"] == 6
    assert generation["available_modes"] == 8
    assert generation["discarded_modes"] == 8
    assert generation["dependent_training_ports"] == ["rsh_l1_x", "rsh_l1_y"]
    assert generation["max_snapshot_relative_residual"] < 2.0e-8
    assert max(generation["response_relative_energy_errors"][:2]) < 1.0e-8
    assert generation["max_response_relative_energy_error"] > 0.5
    assert (
        generation["pod_truncation_curve"][1][
            "max_physical_response_relative_energy_error"
        ]
        < 1.0e-8
    )
    assert generation["energy_orthonormality_error"] < 1.0e-10
    np.testing.assert_allclose(
        reduced.magnetic_operator,
        np.eye(reduced.n_modes),
        atol=1.0e-9,
    )


def test_ngsolve_regular_solid_harmonic_ports_are_divergence_free():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))
    ports = vim.NgsolveHDivRegularSolidHarmonicPorts(mesh, max_degree=3)
    fes = ng.HDiv(mesh, order=2)
    divergence_energies = []
    for field in ports.fields:
        projected = ng.GridFunction(fes)
        projected.Set(field)
        divergence_energies.append(
            float(ng.Integrate(ng.div(projected) ** 2, mesh))
        )

    assert ports.count == 15
    assert ports.diagnostics()["degree_counts"] == {"1": 3, "2": 5, "3": 7}
    assert max(divergence_energies) < 1.0e-20


def test_ngsolve_planar_hdiv_response_reduction_preserves_corner_visible_fields():
    ng = pytest.importorskip("ngsolve")
    geom2d = pytest.importorskip("netgen.geom2d")

    geometry = geom2d.SplineGeometry()
    coordinates = ((0, 0), (1, 0), (1, 0.4), (0.4, 0.4), (0.4, 1), (0, 1))
    points = [geometry.AppendPoint(*point) for point in coordinates]
    for index in range(len(points)):
        geometry.Append(
            ["line", points[index], points[(index + 1) % len(points)]],
            leftdomain=1,
            rightdomain=0,
        )
    geometry.SetMaterial(1, "body")
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=0.35))
    ports = vim.NgsolvePlanarHarmonicPorts(mesh, max_degree=2)

    assert isinstance(ports, vim.PlanarHarmonicPortSet)
    assert ports.count == 4
    assert ports.diagnostics()["degree_counts"] == {"1": 2, "2": 2}
    with ng.TaskManager():
        model = vim.NgsolvePlanarHDivMMMResponseReduction(
            mesh,
            mu_r=1001.0,
            order=1,
            external_fields=(
                ng.CoefficientFunction((1.0, 0.0)),
                ng.CoefficientFunction((0.0, 1.0)),
            ),
            training_fields=ports,
            cg_tol=1.0e-10,
        )
        solution = model.solve()

    generation = model.basis_generation
    assert isinstance(model, vim.PlanarHDivMMMReducedModel)
    assert isinstance(solution, vim.PlanarHDivMMMReducedSolution)
    assert model.body.rt is False
    assert model.diagnostics()["parent_family"] == "BDM"
    assert model.n_modes == 4
    assert generation["protected_physical_modes"] == 2
    assert generation["training_response_modes"] == 2
    assert generation["dependent_training_ports"] == [
        "ph_l1_cos",
        "ph_l1_sin",
    ]
    assert generation["max_snapshot_relative_residual"] < 1.0e-8
    assert generation["max_response_relative_energy_error"] < 1.0e-8
    assert solution.parent_coefficients.shape == (model.parent_ndof, 2)
    assert solution.residual_relative_norm < 1.0e-12
    np.testing.assert_allclose(
        model.magnetic_operator,
        np.eye(model.n_modes),
        atol=1.0e-9,
    )

    degree_three = vim.NgsolvePlanarHarmonicPorts(mesh, max_degree=3)
    with pytest.raises(ValueError, match="order must be at least"):
        vim.NgsolvePlanarHDivMMMResponseReduction(
            mesh,
            mu_r=1001.0,
            body=model.body,
            external_fields=ng.CoefficientFunction((1.0, 0.0)),
            training_fields=degree_three,
        )


def test_ngsolve_hcurl_curl_basis_samples_t_method_current():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))
    fes = ng.HCurl(mesh, order=1, nograds=True)
    vectors = np.zeros((fes.ndof, 1))
    vectors[0, 0] = 1.0

    basis = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        vectors,
        intorder=1,
        materials="cond",
        names=["curl_T0"],
    )

    assert basis.kind == "volume"
    assert basis.n_modes == 1
    assert basis.mass_matrix()[0, 0].real > 0.0


def test_ngsolve_hdiv_magnetization_basis_couples_to_hcurl_current_basis():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("body")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))

    fes_hdiv = ng.HDiv(mesh, order=1)
    hdiv_vectors = np.zeros((fes_hdiv.ndof, 1))
    hdiv_vectors[0, 0] = 1.0
    magnetization = vim.NgsolveHDivMagnetizationBasis(
        mesh,
        fes_hdiv,
        hdiv_vectors,
        intorder=1,
        materials="body",
        names=["M0"],
    )

    fes_hcurl = ng.HCurl(mesh, order=1, nograds=True)
    hcurl_vectors = np.zeros((fes_hcurl.ndof, 1))
    hcurl_vectors[0, 0] = 1.0
    current = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes_hcurl,
        hcurl_vectors,
        intorder=1,
        materials="body",
        names=["curl_T0"],
    )

    assert magnetization.n_modes == 1
    assert magnetization.mass_matrix()[0, 0].real > 0.0
    coupling = vim.MagnetizationCurrentCoupling(
        magnetization,
        current,
        kernel_epsilon=0.1,
    )
    assert coupling.shape == (1, 1)
    assert np.all(np.isfinite(coupling))


def test_ngsolve_coupling_masks_identify_local_static_condensation_dofs():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=3.0))
    fes = ng.HCurl(mesh, order=4, nograds=True)

    masks = vim.NgsolveCouplingDofMasks(fes)

    assert masks["local"].sum() > 0
    assert masks["local_bubble"].sum() == masks["local"].sum()
    assert masks["keep"].sum() == sum(fes.FreeDofs(True))
    assert masks["keep"].sum() + masks["local_bubble"].sum() == sum(fes.FreeDofs(False))


def test_ngsolve_block_krylov_response_feeds_hcurl_t_vim():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))

    fes = ng.HCurl(mesh, order=1, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    port = ng.LinearForm(fes)
    port += ng.CoefficientFunction((-ng.y, ng.x, 0.0)) * v * ng.dx

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()

    response = vim.NgsolveBlockKrylovBasis(
        stiffness,
        mass.mat,
        port,
        steps=2,
        free_dofs=fes.FreeDofs(),
    )
    operator_response = vim.NgsolveOperatorBlockKrylovBasis(
        stiffness,
        mass,
        port,
        steps=2,
        free_dofs=fes.FreeDofs(),
    )
    assert 1 <= response.rank <= 2
    assert operator_response.rank == response.rank

    dense_mass = vim.NgsolveMatrixToDense(mass)
    np.testing.assert_allclose(
        response.vectors.conj().T @ dense_mass @ response.vectors,
        np.eye(response.rank),
        atol=1.0e-11,
    )
    np.testing.assert_allclose(
        operator_response.vectors.conj().T @ dense_mass @ operator_response.vectors,
        np.eye(operator_response.rank),
        atol=1.0e-11,
    )
    assert vim.NgsolveVectorToArray(port).shape == (fes.ndof,)

    volume = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        response.vectors,
        intorder=1,
        materials="cond",
        names=[f"T_resp{i}" for i in range(response.rank)],
    )
    surface = vim.NgsolveSurfaceOmegaBasis(
        mesh,
        (ng.CoefficientFunction((0.0, 1.0, 0.0)),),
        intorder=1,
        boundaries="skin",
        names=["Omega_skin"],
    )

    assert volume.n_modes == response.rank
    assert volume.mass_matrix().trace().real > 0.0

    system = vim.AssembleHybridVIM(volume, surface, sigma=5.8e7, kernel_epsilon=0.1)
    assert system.impedance(
        1j * 100.0,
        surface_impedance=vim.SkinImpedance(1j * 100.0, 5.8e7),
    ).shape == (response.rank + 1, response.rank + 1)

    aext = np.column_stack(
        (
            -volume.points[:, 1],
            volume.points[:, 0],
            np.zeros(volume.n_samples),
        )
    )
    rhs = vim.ExternalVectorPotentialRHS(volume, aext)
    assert rhs.shape == (response.rank,)
    assert np.all(np.isfinite(rhs))


def test_ngsolve_static_condensed_basis_matches_full_parent_response():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=3.0))
    fes = ng.HCurl(mesh, order=3, nograds=True)
    u, v = fes.TnT()
    stiffness_full = ng.BilinearForm(fes)
    stiffness_full += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    stiffness_condensed = ng.BilinearForm(fes, condense=True)
    stiffness_condensed += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    port = ng.LinearForm(fes)
    port += ng.CoefficientFunction((-ng.y, ng.x, 0.0)) * v * ng.dx

    with ng.TaskManager():
        stiffness_full.Assemble()
        stiffness_condensed.Assemble()
        mass.Assemble()
        port.Assemble()

    full = vim.NgsolveBlockKrylovBasis(
        stiffness_full,
        mass,
        port,
        steps=2,
        free_dofs=fes.FreeDofs(False),
    )
    condensed = vim.NgsolveStaticCondensedBlockKrylovBasis(
        stiffness_condensed,
        mass,
        port,
        steps=2,
        free_dofs=fes.FreeDofs(True),
    )
    dense_mass = vim.NgsolveMatrixToDense(mass)

    assert full.rank == condensed.rank
    np.testing.assert_allclose(
        condensed.vectors.conj().T @ dense_mass @ condensed.vectors,
        np.eye(condensed.rank),
        atol=1.0e-11,
    )
    overlap = full.vectors.conj().T @ dense_mass @ condensed.vectors
    np.testing.assert_allclose(np.abs(overlap), np.eye(full.rank), atol=1.0e-10)
    assert condensed.diagnostics()["inactive_dofs"] > 0


@pytest.mark.parametrize("order", [3, 4])
def test_high_order_hcurl_space_reduces_to_low_rank_external_response(order):
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=3.0))

    fes = ng.HCurl(mesh, order=order, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    port = ng.LinearForm(fes)
    port += ng.CoefficientFunction((-ng.y, ng.x, 0.0)) * v * ng.dx

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()

    response = vim.NgsolveBlockKrylovBasis(
        stiffness,
        mass,
        port,
        steps=2,
        free_dofs=fes.FreeDofs(),
    )

    min_active_dofs = 200 if order == 3 else 400
    assert sum(fes.FreeDofs()) >= min_active_dofs
    assert response.rank == 2
    assert response.rank < 0.02 * fes.ndof

    volume = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        response.vectors,
        intorder=1,
        materials="cond",
    )
    assert volume.n_modes == response.rank
    assert volume.mass_matrix().trace().real > 0.0


def test_order6_hcurl_parent_space_compresses_aggressively():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=3.0))

    fes = ng.HCurl(mesh, order=6, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    ports = []
    for cf in (
        ng.CoefficientFunction((-ng.y, ng.x, 0.0)),
        ng.CoefficientFunction((0.0, -ng.z, ng.y)),
    ):
        port = ng.LinearForm(fes)
        port += cf * v * ng.dx
        ports.append(port)

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        for port in ports:
            port.Assemble()

    response = vim.NgsolveBlockKrylovBasis(
        stiffness,
        mass,
        ports,
        steps=3,
        free_dofs=fes.FreeDofs(),
        rtol=1.0e-10,
    )
    info = response.diagnostics()

    assert isinstance(response, vim.EVRSBasis)
    assert info["active_dofs"] >= 1_200
    assert info["rank"] == 6
    assert info["compression_ratio"] < 0.005
    assert info["eliminated_dofs"] > 1_200
    assert info["port_count"] == 2
    assert info["krylov_steps"] == 3
    assert info["construction"] == "ngsolve-dense-block-krylov"


def test_order8_hcurl_static_condensation_then_eliminates_eddy_bubbles():
    ng = pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=3.0))

    fes = ng.HCurl(mesh, order=8, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes, condense=True)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    ports = []
    for cf in (
        ng.CoefficientFunction((-ng.y, ng.x, 0.0)),
        ng.CoefficientFunction((0.0, -ng.z, ng.y)),
    ):
        port = ng.LinearForm(fes)
        port += cf * v * ng.dx
        ports.append(port)

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        for port in ports:
            port.Assemble()

    response = vim.NgsolveStaticCondensedBlockKrylovBasis(
        stiffness,
        mass,
        ports,
        steps=3,
        free_dofs=fes.FreeDofs(True),
        rtol=1.0e-10,
    )
    info = response.diagnostics()

    assert isinstance(response, vim.EVRSBasis)
    assert info["ndof"] >= 2_600
    assert info["active_dofs"] <= 1_100
    assert info["inactive_dofs"] >= 1_500
    assert info["rank"] == 6
    assert info["eddy_visible_dofs"] == 6
    assert info["compression_ratio"] < 0.0025
    assert info["eddy_invisible_dofs"] > 1_000
    assert info["port_count"] == 2
    assert info["krylov_steps"] == 3
    assert info["construction"] == "ngsolve-static-condensed-block-krylov"
