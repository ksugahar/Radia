"""Fast contract tests for the independent HDiv-to-H1 Hodge operator."""

import ngsolve as ng
import netgen.occ as occ
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from ngsolve.meshes import MakeStructured3DMesh

from radia import vim


_OUTER = "back|left|front|right|bottom|top"


def _warped_hex_mesh():
    return MakeStructured3DMesh(
        hexes=True,
        nx=2,
        ny=2,
        nz=2,
        mapping=lambda x, y, z: (
            x + 0.08 * y * z,
            y + 0.05 * x * z,
            z + 0.06 * x * y,
        ),
    )


def test_h1_hodge_operator_is_a_bounded_mapped_hex_bdm2_projection():
    mesh = _warped_hex_mesh()
    hdiv = ng.HDiv(mesh, order=2)
    h1 = ng.H1(mesh, order=2, dirichlet=_OUTER)
    with ng.TaskManager():
        operator = vim.H1HodgeDemagOperator(
            hdiv, h1, definedon=mesh.Materials("default"))

    active = np.flatnonzero(np.asarray(hdiv.FreeDofs(), dtype=bool))
    dense = np.asarray(operator.mat.ToDense())[np.ix_(active, active)]
    rows, columns, values = operator.mass.COO()
    mass = sp.csr_matrix(
        (values, (rows, columns)), shape=(hdiv.ndof, hdiv.ndof))
    mass = mass[active, :][:, active].toarray()
    spectrum = sla.eigh(0.5 * (dense + dense.T), mass, eigvals_only=True)

    assert spectrum[0] >= -2.0e-11
    assert spectrum[-1] <= 1.0 + 2.0e-10
    assert np.max(np.abs(dense - dense.T)) <= 2.0e-11
    diagnostics = operator.Diagnostics()
    assert diagnostics["operator"] == "C.T @ K^-1 @ C"
    assert diagnostics["hdiv_active_dofs"] == len(active)
    assert diagnostics["unit_stiffness"] is True
    assert diagnostics["contraction_contract"] == "standard-unit-H1 metric"
    assert "finite-domain reference" in diagnostics["claim_boundary"]


def test_h1_hodge_operator_exposes_potential_and_rayleigh_quotient():
    mesh = _warped_hex_mesh()
    hdiv = ng.HDiv(mesh, order=2)
    h1 = ng.H1(mesh, order=2, dirichlet=_OUTER)
    with ng.TaskManager():
        operator = vim.H1HodgeDemagOperator(
            hdiv, h1, definedon="default")
        source = ng.GridFunction(hdiv)
        source.Set(ng.CF((ng.x, 0.2 * ng.y, -0.1 * ng.z)))
        potential = operator.Potential(source)
        factor = operator.DemagFactor(source)

    assert potential.space is h1
    assert np.all(np.isfinite(potential.vec.FV().NumPy()))
    assert 0.0 <= factor <= 1.0 + 2.0e-10


def test_h1_hodge_operator_rejects_an_ungrounded_potential_space():
    mesh = _warped_hex_mesh()
    hdiv = ng.HDiv(mesh, order=1)
    h1 = ng.H1(mesh, order=1)
    with np.testing.assert_raises_regex(ValueError, "grounded/constrained"):
        vim.H1HodgeDemagOperator(
            hdiv, h1, definedon=mesh.Materials("default"))


def test_h1_hodge_operator_rejects_a_nonpositive_scalar_metric():
    mesh = _warped_hex_mesh()
    hdiv = ng.HDiv(mesh, order=1)
    h1 = ng.H1(mesh, order=1, dirichlet=_OUTER)
    with np.testing.assert_raises_regex(ValueError, "positive and finite"):
        vim.H1HodgeDemagOperator(
            hdiv,
            h1,
            definedon="default",
            stiffness_coefficient=0.0,
        )


def test_h1_hodge_bdm2_reduction_feeds_hcurl_eddy_bubble_on_mapped_hexes():
    mesh = _warped_hex_mesh()
    h1 = ng.H1(mesh, order=2, dirichlet=_OUTER)
    hcurl = ng.HCurl(mesh, order=2, nograds=True)
    trial, test = hcurl.TnT()
    stiffness = ng.BilinearForm(hcurl)
    stiffness += (
        ng.curl(trial) * ng.curl(test) + 0.05 * trial * test
    ) * ng.dx
    mass = ng.BilinearForm(hcurl)
    mass += trial * test * ng.dx
    port = ng.LinearForm(hcurl)
    port += ng.CF((-ng.y, ng.x, 0.0)) * test * ng.dx

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()
        mixed = vim.NgsolveBDMEddyBubbleVIM(
            mesh,
            hcurl,
            stiffness,
            mass,
            port,
            (),
            hdiv_order=2,
            mu_r=5.0,
            external_fields=(ng.CF((ng.x, 0.0, 0.0)),),
            external_names=("linear_Hx",),
            hdiv_max_modes=1,
            hdiv_intorder=3,
            hdiv_solve_tol=1.0e-10,
            magnetic_materials="default",
            hdiv_definedon="default",
            demag_operator_factory=lambda hdiv: vim.H1HodgeDemagOperator(
                hdiv,
                h1,
                definedon="default",
                boundary_contract="finite-dirichlet-regression",
            ),
            steps=2,
            sigma=1.0e6,
            conductive_materials="default",
            response_backend="dense",
            intorder=3,
            kernel_epsilon=0.1,
            coupling_kernel_epsilon=0.1,
            interaction=vim.HACApKSampledLaplaceInteraction(
                mu=4.0e-7 * np.pi,
                kernel_epsilon=0.1,
                cross_only=False,
            ),
        )
        solution = mixed.solve_frequency(100.0)

    reduction = mixed.hdiv_reduction
    generation = reduction.basis_generation
    assert generation["snapshot_backend"] == "ngsolve-mass-preconditioned-cg"
    assert generation["max_snapshot_relative_residual"] < 1.0e-8
    assert reduction.diagnostics()["demag_backend"] == "H1HodgeDemagOperator"
    assert reduction.parent_family == "BDM"
    assert reduction.parent_order == 2
    assert mixed.hdiv_reduction is reduction
    assert mixed.eddy_system.interaction_backend == "hacapk-sampled-laplace"
    assert solution.residual_relative_norm < 1.0e-10
    assert solution.average_joule_loss[0] > 0.0


def test_restricted_bdm_response_requires_an_exact_space_operator_factory():
    mesh = _warped_hex_mesh()
    field = (ng.CF((ng.x, 0.0, 0.0)),)

    with np.testing.assert_raises_regex(ValueError, "requires demag_operator_factory"):
        vim.NgsolveBDMHDivMMMResponseReduction(
            mesh,
            order=2,
            mu_r=5.0,
            external_fields=field,
            materials="default",
            hdiv_definedon="default",
        )
    with np.testing.assert_raises_regex(ValueError, "mutually exclusive"):
        vim.NgsolveBDMHDivMMMResponseReduction(
            mesh,
            order=2,
            mu_r=5.0,
            external_fields=field,
            materials="default",
            demag_operator=object(),
            demag_operator_factory=lambda _hdiv: object(),
        )


def test_h1_hodge_response_cg_respects_partial_region_free_dofs():
    body = occ.Box(occ.Pnt(-1.0, -1.0, -1.0), occ.Pnt(0.0, 1.0, 1.0))
    air = occ.Box(occ.Pnt(0.0, -1.0, -1.0), occ.Pnt(1.0, 1.0, 1.0))
    body.mat("body")
    air.mat("air")
    for face in (*body.faces, *air.faces):
        face.name = "outer"
    mesh = ng.Mesh(
        occ.OCCGeometry(occ.Glue([body, air])).GenerateMesh(maxh=1.5)
    )
    hdiv = ng.HDiv(mesh, order=1, definedon=mesh.Materials("body"))
    h1 = ng.H1(mesh, order=2, dirichlet="outer")

    with ng.TaskManager():
        demag = vim.H1HodgeDemagOperator(hdiv, h1, definedon="body")
        reduction = vim.NgsolveHDivMMMResponseReduction(
            mesh,
            hdiv,
            mu_r=5.0,
            external_fields=(ng.CF((ng.x, 0.0, 0.0)),),
            max_modes=1,
            materials="body",
            demag_operator=demag,
            solve_tol=1.0e-10,
        )

    assert sum(hdiv.FreeDofs()) < hdiv.ndof
    assert reduction.basis_generation["max_snapshot_relative_residual"] < 1.0e-8
    assert reduction.diagnostics()["min_demag_eigenvalue"] >= -1.0e-10
    solution = reduction.solve()
    assert isinstance(solution, vim.HDivMMMReducedSolution)
    assert solution.reduced_coefficients.shape == (1, 1)
    assert solution.parent_coefficients.shape == (hdiv.ndof, 1)
    assert solution.sampled_magnetization.shape == (
        1,
        reduction.magnetization_basis.n_samples,
        3,
    )
    assert solution.average_magnetization.shape == (1, 3)
    assert solution.residual_relative_norm < 1.0e-12
    manual_average = np.einsum(
        "ik,i->k",
        solution.sampled_magnetization[0],
        reduction.magnetization_basis.weights,
    ) / np.sum(reduction.magnetization_basis.weights)
    np.testing.assert_allclose(solution.average_magnetization[0], manual_average)
