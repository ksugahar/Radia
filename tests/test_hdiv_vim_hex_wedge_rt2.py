"""Fast public-pipeline regression for flat HEX/WEDGE BDM2 HDiv-VIM.

The filename is a legacy label from before the NGSolve BDM/RT distinction was
made explicit in Radia documentation.
"""

import numpy as np
import pytest
import scipy.linalg as sla
import scipy.sparse as sp

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402
import radia._radia_pybind as _rp  # noqa: E402

from radia.vim import (  # noqa: E402
    ChargeGram, FieldFromSolution, MagnetizationSource, Solve,
)


def _cube(kind):
    mapping = lambda x, y, z: (  # noqa: E731
        0.02 * (x - 0.5), 0.02 * (y - 0.5), 0.02 * (z - 0.5))
    kwargs = {"hexes": True} if kind == "hex" else {"prism": True}
    try:
        mesh = MakeStructured3DMesh(
            nx=1, ny=1, nz=1, mapping=mapping, **kwargs)
    except TypeError:
        pytest.skip("this NGSolve build cannot generate a structured prism mesh")
    expected_vertices = 8 if kind == "hex" else 6
    if {len(el.vertices) for el in mesh.Elements(ng.VOL)} != {expected_vertices}:
        pytest.skip(f"structured {kind} generator returned a different topology")
    return mesh


def _dense_demag(B, gram):
    B = sp.csr_matrix(B)
    operator = np.empty((B.shape[1], B.shape[1]))
    for column in range(B.shape[1]):
        basis = np.zeros(B.shape[1])
        basis[column] = 1.0
        charge = np.ascontiguousarray(B @ basis, dtype=np.float64)
        operator[:, column] = B.T @ np.asarray(gram.matvec(charge))
    return 0.5 * (operator + operator.T)


@pytest.mark.parametrize(
    ("kind", "expected_fes_dofs", "expected_charge_dofs"),
    [("hex", 108, 81), ("wedge", 141, 96)],
)
def test_rt2_charge_solve_and_field_pipeline(kind, expected_fes_dofs,
                                             expected_charge_dofs):
    mesh = _cube(kind)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=2)
        B, gram, mass = ChargeGram(fes, eps=1e-14, leafsize=256)
        result = Solve(
            mesh, mu_r=100.0, H_ext=ng.CF((0, 0, 1000.0)),
            order=2, tol=1e-9)

    assert fes.ndof == expected_fes_dofs
    assert B.shape == (expected_charge_dofs, expected_fes_dofs)
    assert mass.shape == (expected_fes_dofs, expected_fes_dofs)
    assert gram.ndof() == expected_charge_dofs
    assert gram.hex_state_check()["ctor"] == gram.hex_state_check()["now"]
    assert result["order"] == 2
    assert abs(result["demag"] - 1.0 / 3.0) < 5e-3
    assert result["iters"] < 50

    field = np.asarray(FieldFromSolution(
        result, np.asarray([[0.04, 0.0, 0.0], [0.0, 0.0, 0.04]]),
        algorithm="direct"))
    assert field.shape == (2, 3)
    assert np.isfinite(field).all()
    assert np.max(np.abs(field[:, 2])) > 1.0


@pytest.mark.parametrize("kind", ["hex", "wedge"])
def test_curved_rt2_nonlinear_material_uses_energy_newton(kind):
    """Curve(2) BDM2 polyhedra use the production C++ nonlinear path."""
    mesh = _cube(kind)
    mesh.Curve(2)
    mu0 = 4.0e-7*np.pi
    linear_bh = np.asarray([
        [0.0, 0.0],
        [1.0e3, mu0*100.0*1.0e3],
        [1.0e5, mu0*100.0*1.0e5],
    ])
    with ng.TaskManager():
        linear = Solve(
            mesh, mu_r=100.0, H_ext=ng.CF((0, 0, 1000.0)),
            order=2, curve_order=2, gram_eps=1e-7, tol=1e-9)
        nonlinear = Solve(
            mesh, bh_table=linear_bh, H_ext=ng.CF((0, 0, 1000.0)),
            order=2, curve_order=2, gram_eps=1e-7, tol=1e-9,
            nl_maxit=30)

    assert nonlinear["order"] == 2
    assert nonlinear["curve_order"] == 2
    assert nonlinear["nonlinear"] is True
    assert nonlinear["linear_solver"] == "energy-newton-cpp"
    assert nonlinear["iters"] < 30
    linear_m = np.asarray(linear["M_avg"], dtype=float)
    nonlinear_m = np.asarray(nonlinear["M_avg"], dtype=float)
    relative = np.linalg.norm(nonlinear_m-linear_m)/np.linalg.norm(linear_m)
    assert relative < 1e-6, relative


@pytest.mark.parametrize("kind", ["hex", "wedge"])
def test_rt2_prescribed_magnetization_uses_the_native_field_path(kind):
    mesh = _cube(kind)
    with ng.TaskManager():
        source = MagnetizationSource(mesh, (0.0, 0.0, 1.0e5), order=2)

    assert source.stats["order"] == 2
    assert source.stats["hmatrix_built"] is False
    assert source.stats["projection_relative_residual"] < 1e-12
    field = np.asarray(source.Field(
        [[0.04, 0.0, 0.0], [0.0, 0.0, 0.04]], algorithm="direct"))
    assert field.shape == (2, 3)
    assert np.isfinite(field).all()
    assert np.max(np.abs(field[:, 2])) > 1.0


@pytest.mark.parametrize("kind", ["hex", "wedge"])
def test_rt2_demag_spectrum_is_physical(kind):
    mesh = _cube(kind)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=2)
        B, gram, mass = ChargeGram(fes, eps=1e-14, leafsize=256)
        operator = _dense_demag(B, gram)
    eigenvalues = sla.eigh(
        operator, sp.csr_matrix(mass).toarray(), eigvals_only=True)
    assert eigenvalues.min() > -1e-10
    assert eigenvalues.max() <= 1.0 + 2e-5


def test_mapped_hex_rt2_material_solve_fails_loud_before_wrong_physics():
    mesh = MakeStructured3DMesh(
        hexes=True, nx=2, ny=1, nz=1,
        mapping=lambda x, y, z: (
            0.02*(x + 0.35*y*z),
            0.02*(y + 0.20*x*z),
            0.02*(z + 0.25*x*y),
        ),
    )
    with ng.TaskManager(), pytest.raises(
            NotImplementedError, match="mapped/non-affine HEX BDM2"):
        Solve(
            mesh, mu_r=100.0, H_ext=ng.CF((0, 0, 1000.0)),
            order=2, tol=1e-9)


def test_mapped_hex_bdm1_material_solve_remains_the_supported_lane():
    mesh = MakeStructured3DMesh(
        hexes=True, nx=2, ny=1, nz=1,
        mapping=lambda x, y, z: (
            0.02*(x + 0.35*y*z),
            0.02*(y + 0.20*x*z),
            0.02*(z + 0.25*x*y),
        ),
    )
    with ng.TaskManager():
        result = Solve(
            mesh, mu_r=10.0, H_ext=ng.CF((0, 0, 1000.0)),
            order=1, tol=1e-8)
    assert result["order"] == 1
    assert result["iters"] < 100
    assert np.isfinite(result["M_avg"]).all()


def test_rt2_hex_far_block_uses_accurate_complete_tensor_rule():
    """A separated affine pair exercises the fast far block, not the exact near recurrence."""
    x5, w5 = np.polynomial.legendre.leggauss(5)
    x5, w5 = 0.5*(x5 + 1.0), 0.5*w5

    def nodes(offset):
        return np.asarray([
            (offset + ix/2, iy/2, iz/2)
            for iz in range(3) for iy in range(3) for ix in range(3)
        ], dtype=float)

    gram = _rp._ChargeGramHMatrix(
        hex_cell_nodes=np.concatenate([nodes(0.0), nodes(4.0)]).ravel(),
        quad_face_nodes=np.empty(0), n_el=2, n_bf=0,
        charge_host=np.asarray([0, 1], dtype=np.int32),
        charge_kind=np.asarray([0, 0], dtype=np.int32),
        charge_expo=np.zeros(6, dtype=np.int32),
        sym_tet_pts=np.asarray([0.25, 0.25, 0.25]),
        sym_tet_w=np.asarray([1.0/6.0]),
        sym_tri_pts=np.asarray([1.0/3.0, 1.0/3.0]),
        sym_tri_w=np.asarray([0.5]),
        gl_out=x5, gw_out=w5, gl_in=x5, gw_in=w5,
        far_tet_pts=np.asarray([0.25, 0.25, 0.25]),
        far_tet_w=np.asarray([1.0/6.0]),
        far_tri_pts=np.asarray([1.0/3.0, 1.0/3.0]),
        far_tri_w=np.asarray([0.5]),
        near_grade=0.5, far_inner_factor=1.0, build=False,
    )

    x10, w10 = np.polynomial.legendre.leggauss(10)
    x10, w10 = 0.5*(x10 + 1.0), 0.5*w10
    q = np.asarray([(x, y, z) for z in x10 for y in x10 for x in x10])
    qw = np.asarray([wx*wy*wz for wz in w10 for wy in w10 for wx in w10])
    reference = 0.0
    source = q + np.asarray([4.0, 0.0, 0.0])
    for target, weight in zip(q, qw):
        reference += weight*np.sum(qw/np.linalg.norm(target - source, axis=1))
    reference /= 4.0*np.pi

    assert gram.stats()["hex_affine_exact_near_factor"] == 1.0
    assert abs(gram.entry(0, 1) - reference) <= 2e-12*abs(reference)
