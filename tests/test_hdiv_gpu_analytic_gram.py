"""Device-filled analytic constant-density Gram against the CPU charge Gram.

Guards: ``_ResidentHDiv.from_analytic_constant_block`` fills, scales and
normalizes the block on the device from geometry alone: its action equals
``matvec_sym`` of the independent CPU ``_ChargeGramHMatrix`` built with the same
all-analytic near route, on tetrahedra alone and on tetrahedra with triangles.
Refuses a CPU backend, more than 128 charges, degenerate geometry, mismatched
vectors and material solves. Records the geometry transfer and no matrix.
Assumes: the CPU Gram is the reference on identical geometry; RELATIVE is the
limit of the retired resident apply check. Parity, not independent accuracy.
Optional GPU lane: skips visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set; a
configured lane loads the binding through the provenance-checked loader.
"""

import os

import numpy as np
import pytest

RELATIVE = 1e-12
TETRA = np.array([[0., 0., 0.], [.4, 0., 0.], [0., .4, 0.], [0., 0., .4]])
TETRA2 = .5*TETRA + np.array([.6, .1, .1])
CELLS = np.concatenate([TETRA.ravel(), TETRA2.ravel()])
FACES = np.concatenate([TETRA[[0, 1, 2]].ravel(), TETRA2[[0, 2, 3]].ravel()])


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


@pytest.mark.parametrize("faces", [np.empty(0), FACES], ids=["tetrahedra", "tetrahedra-triangles"])
def test_analytic_gram_action_matches_cpu_charge_gram(request, cuda, faces):
    import ngsolve as ng
    import radia._radia_pybind as rp
    with ng.TaskManager():
        cpu = rp._ChargeGramHMatrix(CELLS, faces, 2, eps=1e-12, leaf=64, eta=2., near_factor=1e30, build=True, far_quad=0)
    n = CELLS.size//12 + faces.size//9
    gram = cuda._ResidentHDiv.from_analytic_constant_block(CELLS, faces, backend="gpu")
    try:
        rng = np.random.default_rng(n)
        worst = 0.
        for _ in range(4):
            x = rng.standard_normal(n)
            with ng.TaskManager():
                expected = np.asarray(cpu.matvec_sym(x))
            worst = max(worst, np.linalg.norm(gram.apply(x) - expected)/np.linalg.norm(expected))
        request.node.user_properties.append(("relative", float(worst)))
        assert worst <= RELATIVE, f"relative={worst:.3e} limit={RELATIVE:.0e}"
        record = gram.assembly()
        assert record["route"] == "device-analytic-constant-block" and record["charges"] == n, record
        assert {row["payload"] for row in record["transfer_ledger"]} == {"resident_plan", "analytic_geometry"}, record
        assert gram.stats()["applications"] == 4
        with pytest.raises(ValueError, match="dimension"):
            gram.apply(np.ones(n + 1))
    finally:
        gram.close()
    with pytest.raises(RuntimeError, match="closed"):
        gram.apply(np.ones(n))


def test_analytic_gram_has_no_material_solve(tmp_path, cuda):
    from netgen.csg import unit_cube
    from radia import hdiv_gpu
    path = tmp_path/"cube.vol"
    unit_cube.GenerateMesh(maxh=1.).Save(str(path))
    setup = hdiv_gpu._prepare_device_components(path, 0, cuda, aca=None, compression="dense")
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=setup["dofs"], backend="gpu")
    gram = cuda._ResidentHDiv.from_analytic_constant_block(CELLS, np.empty(0), backend="gpu")
    try:
        with pytest.raises(ValueError, match="no charge map or material solve"):
            gram.solve_linear(mass, np.ones(2), 1., 1e-10, 10, "mass-diagonal")
    finally:
        gram.close()
        mass.close()


@pytest.mark.parametrize("cells,faces,backend,message", [
    (CELLS, np.empty(0), "cpu", "backend"),
    (np.tile(CELLS, 65), np.empty(0), "gpu", "budget"),
    (CELLS[:11], np.empty(0), "gpu", "budget"),
    (np.empty(0), np.empty(0), "gpu", "budget"),
    (np.zeros(12), np.empty(0), "gpu", "Degenerate"),
], ids=["cpu", "too-many", "partial", "empty", "degenerate"])
def test_analytic_gram_refuses_unsupported_input(cuda, cells, faces, backend, message):
    with pytest.raises(ValueError, match=message):
        cuda._ResidentHDiv.from_analytic_constant_block(cells, faces, backend=backend)
