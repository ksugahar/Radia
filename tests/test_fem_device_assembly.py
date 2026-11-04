"""Device-assembled H1 and HCurl operators against NGSolve assembly on identical inputs.

Guards: on straight tets, trilinear hexes, prisms, warped pyramids, a tet/hex/pyramid mesh and curved tets, with volume and boundary terms at
p=1,2,3, the device sparsity pattern and values act on random vectors like the NGSolve matrix of
the same form. A wrong metric, Piola map, orientation class, quadrature rule or reference tensor
moves the action far beyond the detector. The ledger names only permitted payload families, so a
matrix crossing fails here.
Assumes: NGSolve assembly is the reference; this is parity of one discretization, not accuracy.
ACTION_REL is a PROPOSED detector, not an adopted criterion. Optional GPU lane: skips visibly
unless RADIA_FEM_NATIVE_TEST_DIR is set; a configured lane fails if the binding cannot load.
"""
import os

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

ACTION_REL = 1e-11
BONUS = (2, 2, 2)
COEFFICIENTS = dict(stiff=2.0, mass=3.0, mass_im=5.0, trace={"back": 7.0, "top": 2.0, "exterior": 4.0},
                    trace_im={"back": 1.5, "exterior": 0.5}, stiff_aux=2.0, mass_aux=4.0, trace_aux={"back": 8.0, "exterior": 3.0})


@pytest.fixture(scope="module")
def route():
    if not os.environ.get("RADIA_FEM_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_FEM_NATIVE_TEST_DIR")
    import radia.fem_gpu as module
    return module


def _tets():
    from ngsolve.meshes import MakeStructured3DMesh
    return MakeStructured3DMesh(hexes=False, nx=3, ny=2, nz=2,
                                mapping=lambda x, y, z: (2 * x + 0.1 * y, y + 0.2 * z, 0.5 * z + 0.1 * x))


def _hexes():
    from ngsolve.meshes import MakeStructured3DMesh
    return MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2,
                                mapping=lambda x, y, z: (x * (1 + 0.2 * y), y + 0.1 * x * z, z + 0.05 * x * y))


def _curved():
    from netgen.csg import CSGeometry, Sphere, Pnt
    geometry = CSGeometry()
    geometry.Add(Sphere(Pnt(0, 0, 0), 1).bc("back"))
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=0.6))
    mesh.Curve(2)
    return mesh


def _prisms():
    from ngsolve.meshes import MakeStructured3DMesh
    return MakeStructured3DMesh(hexes=False, prism=True, nx=2, ny=2, nz=2,
                                mapping=lambda x, y, z: (x * (1 + 0.2 * y), y + 0.1 * x * z, z + 0.05 * x * y))


def _pyramids():
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import pyramid_cube
    return pyramid_cube(0.3)


def _mixed():
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import mixed_grid
    return mixed_grid((2, 2, 2))


MESHES = {"tets": _tets, "hexes": _hexes, "curved": _curved, "prisms": _prisms, "pyramids": _pyramids, "mixed": _mixed}


def _coefficients(route, mesh):
    known = set(mesh.GetBoundaries())
    values = dict(COEFFICIENTS)
    for field in ("trace", "trace_im", "trace_aux"):
        values[field] = {k: v for k, v in COEFFICIENTS[field].items() if k in known}
    return route.Coefficients(**values)


def _relative(actual, expected):
    return float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))


def _check_space(route, fes, real_fes, coefficients):
    system = route.DeviceSystem(fes, coefficients, bonus=BONUS)
    assert system.native.nonzeros() > 0
    reference = route.cpu_form(fes, coefficients, bonus=BONUS)
    reference.Assemble()
    auxiliary = route.cpu_form(real_fes, coefficients, bonus=BONUS, part="aux")
    auxiliary.Assemble()
    generator = np.random.default_rng(7)
    for _ in range(2):
        x = generator.standard_normal(fes.ndof)
        vector = reference.mat.CreateRowVector()
        vector.FV().NumPy()[:] = x
        out = vector.CreateVector()
        out.data = reference.mat * vector
        expected = np.array(out.FV().NumPy(), copy=True)
        real_part = _relative(system.action("re", x), expected.real)
        assert real_part <= ACTION_REL, f"real action off by {real_part}, limit {ACTION_REL}"
        imaginary = _relative(system.action("im", x), expected.imag)
        assert imaginary <= ACTION_REL, f"imaginary action off by {imaginary}, limit {ACTION_REL}"
        vector = auxiliary.mat.CreateRowVector()
        vector.FV().NumPy()[:] = x
        out = vector.CreateVector()
        out.data = auxiliary.mat * vector
        aux = _relative(system.action("aux", x), out.FV().NumPy())
        assert aux <= ACTION_REL, f"auxiliary action off by {aux}, limit {ACTION_REL}"
    payloads = {row["payload"] for row in system.ledger()}
    assert not any("matrix" in name or "operator" in name for name in payloads), payloads
    assert any(name.startswith("mesh.dofs") for name in payloads), payloads
    assert not any(row["direction"] == "D2H" and row["payload"].startswith(("mesh.", "table.", "coef.")) for row in system.ledger())


@pytest.mark.parametrize("mesh_name", sorted(MESHES))
@pytest.mark.parametrize("order", [1, 2, 3])
def test_h1_action_matches_ngsolve(route, mesh_name, order):
    if mesh_name == "curved" and order == 3:
        pytest.skip("the curved case is covered at p=1 and p=2 to stay within the runtime target")
    mesh = MESHES[mesh_name]()
    _check_space(route, ng.H1(mesh, order=order, complex=True), ng.H1(mesh, order=order), _coefficients(route, mesh))


@pytest.mark.parametrize("mesh_name", sorted(MESHES))
@pytest.mark.parametrize("order", [1, 2, 3])
def test_hcurl_action_matches_ngsolve(route, mesh_name, order):
    if mesh_name == "curved" and order == 3:
        pytest.skip("the curved case is covered at p=1 and p=2 to stay within the runtime target")
    mesh = MESHES[mesh_name]()
    nograds = order == 1
    _check_space(route, ng.HCurl(mesh, order=order, complex=True, nograds=nograds), ng.HCurl(mesh, order=order, nograds=nograds),
                 _coefficients(route, mesh))


def test_unknown_region_and_backend_raise(route):
    mesh = _tets()
    with pytest.raises(ValueError, match="regions the mesh does not have"):
        route.DeviceSystem(ng.H1(mesh, order=1), route.Coefficients(stiff={"nowhere": 1.0}))
    with pytest.raises(ValueError, match="Backend must be"):
        route.solve_h1(ng.H1(mesh, order=1), route.Coefficients(stiff=1.0), np.zeros(mesh.nv), Backend="auto")
