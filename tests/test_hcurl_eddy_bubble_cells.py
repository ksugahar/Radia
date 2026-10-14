import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from netgen.csg import Pnt  # noqa: E402
from netgen.meshing import (  # noqa: E402
    EdgeDescriptor,
    Element1D,
    Element2D,
    Element3D,
    FaceDescriptor,
    Mesh as NetgenMesh,
    MeshPoint,
)

from radia import vim  # noqa: E402
from radia.vim._hcurl_tet_interaction import (  # noqa: E402
    _project_cell_currents_to_subtets,
)


_CELLS_3D = {
    "tet": (
        ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
        ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)),
    ),
    "hex": (
        (
            (0, 0, 0),
            (1, 0, 0),
            (1, 1, 0),
            (0, 1, 0),
            (0, 0, 1),
            (1, 0, 1),
            (1, 1, 1),
            (0, 1, 1),
        ),
        (
            (0, 3, 2, 1),
            (4, 5, 6, 7),
            (0, 1, 5, 4),
            (1, 2, 6, 5),
            (2, 3, 7, 6),
            (3, 0, 4, 7),
        ),
    ),
    "wedge": (
        ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)),
        ((0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
    ),
    "pyramid": (
        ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)),
        ((0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
    ),
}

_CELLS_2D = {
    "trig": (((0, 0, 0), (1, 0, 0), (0, 1, 0)), (0, 1, 2)),
    "quad": (((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)), (0, 1, 2, 3)),
}


def _add_2d_boundary_descriptor(mesh, name="skin"):
    mesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    mesh.SetBCName(0, name)
    edge = EdgeDescriptor()
    edge.edgenr = 1
    edge.surfnr = (1, -1)
    edge.domin = 1
    edge.domout = 0
    edge.name = name
    mesh.Add(edge)


def _single_cell_mesh(family):
    if family in _CELLS_3D:
        vertices, faces = _CELLS_3D[family]
        netgen_mesh = NetgenMesh(dim=3)
        netgen_mesh.SetMaterial(1, "cond")
        boundary = netgen_mesh.Add(
            FaceDescriptor(surfnr=1, domin=1, bc=1)
        )
        netgen_mesh.SetBCName(0, "skin")
        points = [
            netgen_mesh.Add(MeshPoint(Pnt(*point))) for point in vertices
        ]
        netgen_mesh.Add(Element3D(1, points))
        for face in faces:
            netgen_mesh.Add(Element2D(boundary, [points[i] for i in face]))
        return ng.Mesh(netgen_mesh)

    vertices, element = _CELLS_2D[family]
    netgen_mesh = NetgenMesh(dim=2)
    netgen_mesh.SetMaterial(1, "cond")
    _add_2d_boundary_descriptor(netgen_mesh)
    points = [netgen_mesh.Add(MeshPoint(Pnt(*point))) for point in vertices]
    netgen_mesh.Add(Element2D(1, [points[i] for i in element]))
    for first, second in zip(element, element[1:] + element[:1]):
        netgen_mesh.Add(Element1D([points[first], points[second]], index=1))
    return ng.Mesh(netgen_mesh)


def _warped_hex_mesh():
    vertices = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.025, -0.0125, 1.0),
        (1.0375, 0.025, 1.025),
        (0.975, 1.05, 0.9625),
        (-0.025, 0.975, 1.0125),
    )
    netgen_mesh = NetgenMesh(dim=3)
    netgen_mesh.SetMaterial(1, "cond")
    boundary = netgen_mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    netgen_mesh.SetBCName(0, "skin")
    points = [netgen_mesh.Add(MeshPoint(Pnt(*point))) for point in vertices]
    netgen_mesh.Add(Element3D(1, points))
    for face in _CELLS_3D["hex"][1]:
        netgen_mesh.Add(Element2D(boundary, [points[i] for i in face]))
    return ng.Mesh(netgen_mesh)


def _mixed_tri_quad_mesh():
    netgen_mesh = NetgenMesh(dim=2)
    netgen_mesh.SetMaterial(1, "cond")
    _add_2d_boundary_descriptor(netgen_mesh)
    points = [
        netgen_mesh.Add(MeshPoint(Pnt(x, y, 0.0)))
        for x, y in ((0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1))
    ]
    netgen_mesh.Add(Element2D(1, [points[i] for i in (0, 1, 4, 3)]))
    netgen_mesh.Add(Element2D(1, [points[i] for i in (1, 2, 5)]))
    netgen_mesh.Add(Element2D(1, [points[i] for i in (1, 5, 4)]))
    for first, second in ((0, 1), (1, 2), (2, 5), (5, 4), (4, 3), (3, 0)):
        netgen_mesh.Add(Element1D([points[first], points[second]], index=1))
    return ng.Mesh(netgen_mesh)


def _two_triangle_square_mesh():
    netgen_mesh = NetgenMesh(dim=2)
    netgen_mesh.SetMaterial(1, "cond")
    _add_2d_boundary_descriptor(netgen_mesh)
    points = [
        netgen_mesh.Add(MeshPoint(Pnt(*point)))
        for point in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0))
    ]
    netgen_mesh.Add(Element2D(1, [points[i] for i in (0, 1, 2)]))
    netgen_mesh.Add(Element2D(1, [points[i] for i in (0, 2, 3)]))
    for first, second in ((0, 1), (1, 2), (2, 3), (3, 0)):
        netgen_mesh.Add(Element1D([points[first], points[second]], index=1))
    return ng.Mesh(netgen_mesh)


def _eddy_basis(mesh, order=3):
    fes = ng.HCurl(mesh, order=order, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += (ng.curl(u) * ng.curl(v) + 0.1 * u * v) * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    port = ng.LinearForm(fes)
    if int(mesh.dim) == 3:
        port += ng.CF((-ng.y, ng.x, 0.0)) * ng.curl(v) * ng.dx
    else:
        port += ng.curl(v) * ng.dx
    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()
        basis = vim.NgsolveEddyBubbleHCurlBasis(
            mesh,
            fes,
            stiffness,
            mass,
            port,
            steps=1,
            conductive_materials="cond",
            response_backend="dense",
            intorder=6,
            parent_order=order,
        )
    return fes, basis


@pytest.mark.parametrize("family", ("tet", "hex", "wedge", "pyramid", "trig", "quad"))
def test_hcurl_eddy_bubble_builds_for_every_ngsolve_cell_family(family):
    mesh = _single_cell_mesh(family)
    fes, basis = _eddy_basis(mesh, order=6)
    info = basis.diagnostics()

    assert basis.rank == 1
    assert basis.current_basis.n_modes == 1
    assert info["cell_families"]["family_counts"] == {family: 1}
    assert info["cell_families"]["eddy_bubble_supported"] is True
    assert info["eddy_bubbling"]["policy"]["partitioned_free_dofs"] == fes.ndof
    assert np.all(np.isfinite(basis.current_basis.modes))
    if int(mesh.dim) == 2:
        assert np.max(np.abs(basis.current_basis.modes[:, :, :2])) == 0.0
        assert np.max(np.abs(basis.current_basis.modes[:, :, 2])) > 0.0
        with pytest.raises(ValueError, match="planar-log interaction"):
            basis.assemble_vim(sigma=5.8e7)


@pytest.mark.parametrize(
    ("family", "subtet_count"),
    (("tet", 1), ("hex", 6), ("wedge", 3), ("pyramid", 2)),
)
def test_affine_3d_cell_family_uses_epsilon_free_subtet_interaction(
    family,
    subtet_count,
):
    mesh = _single_cell_mesh(family)
    fes = ng.HCurl(mesh, order=2, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x, 0.0)))
    vectors = field.vec.FV().NumPy().copy()
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
        degree=0,
        projection_quad=4,
        outer_quad=4,
        materials="cond",
        projection_tolerance=1.0e-10,
    )
    info = interaction.diagnostics()

    assert interaction.matrix.shape == (1, 1)
    assert interaction.matrix[0, 0] > 0.0
    assert info["family_counts"] == {family: 1}
    assert info["subtet_count"] == subtet_count
    assert info["projection_relative_residual"] < 1.0e-12
    assert info["geometry_relative_residual"] < 1.0e-9
    assert info["kernel_epsilon_m"] is None


@pytest.mark.parametrize("family", ("tet", "hex", "wedge"))
def test_hcurl_cell_cluster_derivative_closes_piola_scaling(family):
    mesh = _single_cell_mesh(family)
    fes = ng.HCurl(mesh, order=2, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x, 0.0)))
    vectors = field.vec.FV().NumPy().copy()
    basis = vim.NgsolveHCurlCurlBasis(
        mesh, fes, vectors, intorder=4, materials="cond"
    )
    interaction = vim.NgsolveHCurlCellVolumeInteraction(
        mesh, fes, vectors, basis, degree=0, projection_quad=4,
        outer_quad=4, materials="cond", projection_tolerance=1.0e-10,
    )
    deformation_space = ng.VectorH1(mesh, order=1)
    scaling = ng.GridFunction(deformation_space)
    scaling.Set(ng.CF((ng.x, ng.y, ng.z)))
    with ng.TaskManager():
        velocities = vim.SampleNgsolveHCurlCellSubtetVelocities(
            mesh, [scaling], interaction, materials="cond"
        )
    left = np.array([0.3+0.2j])
    right = np.array([-0.4+0.1j])
    observed = interaction.matrix.directional_contractions(
        velocities, left, right
    )[0]
    expected = np.vdot(left, interaction.matrix.to_dense() @ right)
    np.testing.assert_allclose(observed, expected, rtol=2.0e-10, atol=2.0e-13)


@pytest.mark.parametrize("family", ("tet", "hex", "wedge"))
def test_hcurl_cell_activation_derivative_groups_subtets_by_parent(family):
    mesh = _single_cell_mesh(family)
    fes = ng.HCurl(mesh, order=2, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x, 0.0)))
    vectors = field.vec.FV().NumPy().copy()
    basis = vim.NgsolveHCurlCurlBasis(
        mesh, fes, vectors, intorder=4, materials="cond"
    )
    interaction = vim.NgsolveHCurlCellVolumeInteraction(
        mesh, fes, vectors, basis, degree=0, projection_quad=4,
        outer_quad=4, materials="cond", projection_tolerance=1.0e-10,
    )
    rho = np.array([0.63]); power = 1.7
    left = np.array([0.3+0.2j]); right = np.array([-0.4+0.1j])
    observed = interaction.matrix.activation_contractions(
        rho, left, right, power=power
    )[0]
    base = np.vdot(left, interaction.matrix.to_dense()@right)
    expected = 2*power*rho[0]**(2*power-1)*base
    np.testing.assert_allclose(observed, expected, rtol=3e-12, atol=2e-14)
    step = 2e-6
    plus = np.vdot(left, interaction.matrix.activation_to_dense(
        rho+step, power=power)@right)
    minus = np.vdot(left, interaction.matrix.activation_to_dense(
        rho-step, power=power)@right)
    np.testing.assert_allclose(observed, (plus-minus)/(2*step),
        rtol=2e-10, atol=2e-12)


def test_warped_hex_uses_uniform_h_refinement_until_both_residuals_pass():
    mesh = _warped_hex_mesh()
    fes = ng.HCurl(mesh, order=2, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x, ng.z * ng.x)))
    vectors = field.vec.FV().NumPy().copy()
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
        degree=6,
        projection_quad=5,
        outer_quad=3,
        projection_tolerance=1.0e-5,
        geometry_tolerance=7.0e-3,
        hex_geometry_backend="affine-subtet",
        max_subdivision_levels=1,
        materials="cond",
    )
    info = interaction.diagnostics()

    assert info["subdivision_strategy"] == "uniform"
    assert info["subdivision_level"] == 1
    assert info["subtet_count"] == 48
    assert info["projection_residual_history"][1] < info["projection_residual_history"][0]
    assert info["geometry_residual_history"][1] < info["geometry_residual_history"][0]
    assert info["projection_relative_residual"] < 1.0e-5
    assert info["geometry_relative_residual"] < 7.0e-3
    assert interaction.matrix[0, 0] > 0.0


def test_warped_order1_hex_uses_direct_q2_reference_density_by_default():
    mesh = _warped_hex_mesh()
    fes = ng.HCurl(mesh, order=1, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x, ng.z * ng.x)))
    vectors = field.vec.FV().NumPy().copy()
    basis = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        vectors,
        intorder=5,
        materials="cond",
    )

    interaction = vim.NgsolveHCurlCellVolumeInteraction(
        mesh,
        fes,
        vectors,
        basis,
        materials="cond",
        matrix_free=True,
    )
    info = interaction.diagnostics()

    assert info["geometry_backend"] == "direct-q2-hex-reference-density"
    assert info["subdivision_strategy"] == "direct-q2-hex-reference-density"
    assert info["polynomial_degree"] == 2
    assert info["required_polynomial_degree"] == 2
    assert info["subdivision_level"] == 0
    assert info["subtet_count"] == 6
    assert info["charge_count"] == 27
    assert info["projection_relative_residual"] < 1.0e-12
    assert info["geometry_relative_residual"] < 1.0e-12
    assert info["hmatrix_operator"]["scalar_gram_backend"] == (
        "direct-q2-hex-reference-density"
    )
    assert info["hmatrix_operator"]["tensor_degree"] == 2
    assert interaction.matrix[0, 0] > 0.0


def test_affine_order1_hex_uses_direct_backend_and_matches_legacy_subtet():
    mesh = _single_cell_mesh("hex")
    fes = ng.HCurl(mesh, order=1, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x, 0.0)))
    vectors = field.vec.FV().NumPy().copy()
    basis = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        vectors,
        intorder=4,
        materials="cond",
    )

    direct = vim.NgsolveHCurlCellVolumeInteraction(
        mesh,
        fes,
        vectors,
        basis,
        degree=0,
        projection_quad=4,
        materials="cond",
        hex_geometry_backend="direct",
        matrix_free=False,
        projection_tolerance=1.0e-12,
        geometry_tolerance=1.0e-12,
    )
    subtet = vim.NgsolveHCurlCellVolumeInteraction(
        mesh,
        fes,
        vectors,
        basis,
        degree=0,
        projection_quad=4,
        outer_quad=4,
        materials="cond",
        hex_geometry_backend="affine-subtet",
        matrix_free=False,
        projection_tolerance=1.0e-12,
        geometry_tolerance=1.0e-12,
    )
    automatic = vim.NgsolveHCurlCellVolumeInteraction(
        mesh,
        fes,
        vectors,
        basis,
        materials="cond",
        matrix_free=False,
        projection_tolerance=1.0e-12,
        geometry_tolerance=1.0e-12,
    )

    assert direct.diagnostics()["geometry_backend"] == (
        "direct-q2-hex-reference-density"
    )
    assert subtet.diagnostics()["geometry_backend"] == "piecewise-affine-subtet"
    assert automatic.diagnostics()["geometry_backend"] == (
        "direct-q2-hex-reference-density"
    )
    np.testing.assert_allclose(
        direct.matrix,
        subtet.matrix,
        rtol=2.0e-4,
        atol=1.0e-16,
    )


def test_hex_geometry_backend_rejects_unknown_policy():
    mesh = _single_cell_mesh("hex")
    fes = ng.HCurl(mesh, order=1, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x, 0.0)))
    vectors = field.vec.FV().NumPy().copy()
    basis = vim.NgsolveHCurlCurlBasis(mesh, fes, vectors, intorder=4)

    with pytest.raises(ValueError, match="hex_geometry_backend"):
        vim.NgsolveHCurlCellVolumeInteraction(
            mesh,
            fes,
            vectors,
            basis,
            hex_geometry_backend="silently-refine-forever",
        )


@pytest.mark.parametrize(
    ("family", "degree", "residual_limit"),
    (
        ("hex", 18, 1.0e-8),
        ("wedge", 12, 1.0e-8),
        ("pyramid", 18, 1.0e-4),
    ),
)
def test_p6_non_tet_uses_high_order_analytic_moments(
    family,
    degree,
    residual_limit,
):
    mesh = _single_cell_mesh(family)
    fes, basis = _eddy_basis(mesh, order=6)
    interaction = basis.cell_volume_interaction(mesh, fes, materials="cond")
    info = interaction.diagnostics()

    assert info["polynomial_degree"] == degree
    assert info["projection_relative_residual"] < residual_limit
    assert info["subdivision_level"] == 0
    assert info["kernel_epsilon_m"] is None
    assert info["matrix_free"] is True
    assert info["charge_count"] <= 3 * info["subtet_count"]
    assert (
        info["hmatrix_operator"]["uncompressed_charge_count"]
        > info["charge_count"]
    )
    assert interaction.matrix[0, 0] > 0.0


def test_p6_pyramid_apex_refinement_reaches_strict_projection_gate():
    mesh = _single_cell_mesh("pyramid")
    fes, basis = _eddy_basis(mesh, order=6)
    projected = _project_cell_currents_to_subtets(
        mesh,
        fes,
        basis.response_basis.vectors,
        degree=18,
        projection_quad=10,
        subdivision_level=8,
        materials="cond",
    )

    assert projected["subtet_count"] == 114
    assert projected["relative_residual"] < 1.0e-8
    assert projected["relative_geometry_error"] < 1.0e-12


def test_p7_hex_caps_analytic_degree_and_accepts_measured_residual():
    mesh = _single_cell_mesh("hex")
    fes, basis = _eddy_basis(mesh, order=7)
    interaction = basis.cell_volume_interaction(mesh, fes, materials="cond")
    info = interaction.diagnostics()

    assert info["required_polynomial_degree"] == 21
    assert info["polynomial_degree"] == 18
    assert info["degree_capped"] is True
    assert info["projection_relative_residual"] < 1.0e-4
    assert info["subdivision_level"] == 0
    assert info["dense_moment_pairs"] < 20_000_000
    assert interaction.matrix[0, 0] > 0.0


def test_generic_tet_cell_interaction_matches_native_tet_path():
    mesh = _single_cell_mesh("tet")
    fes = ng.HCurl(mesh, order=2, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x, 0.0)))
    vectors = field.vec.FV().NumPy().copy()
    basis = vim.NgsolveHCurlCurlBasis(mesh, fes, vectors, intorder=4)

    generic = vim.NgsolveHCurlCellVolumeInteraction(
        mesh, fes, vectors, basis, degree=0, projection_quad=4, outer_quad=4
    )
    native = vim.NgsolveHCurlTetVolumeInteraction(
        mesh, fes, vectors, basis, degree=0, projection_quad=4, outer_quad=4
    )

    np.testing.assert_allclose(generic.matrix, native.matrix, rtol=2.0e-14)


def test_mixed_trig_quad_eddy_bubble_preserves_internal_edge_adjacency():
    mesh = _mixed_tri_quad_mesh()
    fes, basis = _eddy_basis(mesh)
    topology = vim.ClassifyNgsolveEddyTopology(mesh, "cond")
    policy = vim.NgsolveEddyDofPolicy(mesh, fes, topology)
    info = basis.diagnostics()

    assert info["cell_families"]["family_counts"] == {"quad": 1, "trig": 2}
    assert info["cell_families"]["mixed_cell_families"] is True
    assert len(topology.loop_bridge_faces) == 2
    assert topology.conductor_graph().node_count == 3
    assert topology.conductor_graph().edge_count == 2
    assert np.count_nonzero(policy.loop_bridge) > 0
    assert policy.diagnostics()["partitioned_free_dofs"] == fes.ndof
    assert np.max(np.abs(basis.current_basis.modes[:, :, :2])) == 0.0


def _constant_planar_current_interaction(mesh):
    fes = ng.HCurl(mesh, order=2, nograds=True)
    field = ng.GridFunction(fes)
    field.Set(ng.CF((-ng.y, ng.x)))
    vectors = field.vec.FV().NumPy().copy()
    basis = vim.NgsolveHCurlCurlBasis(
        mesh,
        fes,
        vectors,
        intorder=6,
        materials="cond",
    )
    return vim.NgsolveHCurlPlanarVolumeInteraction(
        mesh,
        fes,
        vectors,
        basis,
        degree=0,
        projection_order=6,
        projection_tolerance=1.0e-12,
        materials="cond",
    )


def test_planar_log_interaction_agrees_for_one_quad_and_two_trigs():
    quad = _constant_planar_current_interaction(_single_cell_mesh("quad"))
    trig = _constant_planar_current_interaction(_two_triangle_square_mesh())

    assert quad.matrix[0, 0] > 0.0
    assert trig.matrix[0, 0] > 0.0
    np.testing.assert_allclose(quad.matrix, trig.matrix, rtol=1.0e-4)
    assert quad.diagnostics()["kernel_epsilon_m"] is None
    assert quad.diagnostics()["matrix_free"] is True
    assert trig.diagnostics()["family_counts"] == {"trig": 2}


@pytest.mark.parametrize("family", ("trig", "quad"))
def test_planar_eddy_bubble_feeds_log_vim_system(family):
    mesh = _single_cell_mesh(family)
    fes, basis = _eddy_basis(mesh)
    interaction = basis.planar_volume_interaction(
        mesh,
        fes,
        degree=3,
        projection_order=10,
        projection_tolerance=1.0e-9,
        materials="cond",
    )
    system = basis.assemble_vim(sigma=5.8e7, interaction=interaction)

    assert system.n_modes == basis.rank
    assert system.inductance[0, 0] > 0.0
    assert system.diagnostics()["passive_blocks"] is True
    assert interaction.diagnostics()["projection_relative_residual"] < 1.0e-10
