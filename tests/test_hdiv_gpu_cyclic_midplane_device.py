"""A connected HEX sector with the midplane mirror on the device reproduces the ring and the CPU route.

Guards: `image="-z"` with `image_cyclic` and the azimuthal periodic trace on the device, whose image terms T = M o R
combine a mirror with a rotation, on the dense route at BDM1 and BDM2.
Assumes: an annular HEX yoke, a field along +z, a mesh that is the mirror of its upper half. The ring is the same device
route with every image removed and the CPU route shares the discretization, so this is parity, not independent accuracy.
The limits are PROPOSED detectors, not adopted criteria. The checked-mesh receipt gate is not under test and is bypassed
for scratch meshes. Optional GPU lane: skips visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set; a configured lane loads the
binding through the provenance-checked loader.
"""

import contextlib
import importlib.util
import os
from pathlib import Path

import numpy as np
import pytest

RING_REL = 1e-6
CPU_REL = 1e-6
INV_CHI = 1.0/9.0
FOLD = 8
_RECEIPT = dict(label_contract="unused.json", mesh_check_receipt="unused.receipt", mesh_check_receipt_sha256="0"*64,
                mesh_check_provenance={})


def _cpu_test():
    spec = importlib.util.spec_from_file_location("hdiv_cyclic_midplane", Path(__file__).with_name("test_hdiv_cyclic_midplane.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextlib.contextmanager
def _unchecked_snapshot(path, *_):
    yield Path(path), "0"*64, "0"*64


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


@pytest.fixture
def prepared(monkeypatch):
    from radia import hdiv_gpu
    monkeypatch.setattr(hdiv_gpu, "_checked_mesh_snapshot", _unchecked_snapshot)
    return hdiv_gpu


def _mean_z(device, x):
    """Mean of M_z over the meshed part, which is the full-domain mean on the images of this field."""
    import ngsolve as ng
    magnetization = ng.GridFunction(device.fespace)
    magnetization.vec.FV().NumPy()[:] = x
    return float(ng.Integrate(magnetization[2], device.fespace.mesh))/float(ng.Integrate(ng.CF(1.0), device.fespace.mesh))


def _device_mean(prepared, path, order, **images):
    import ngsolve as ng
    source = ng.CF((0.0, 0.0, 1e5))
    with prepared.PreparedHDivLinear(path, backend="gpu", order=order, gram_compression="dense", **_RECEIPT, **images) as device:
        result = device.solve(source, INV_CHI, tolerance=1e-10, max_iterations=4000)
        return _mean_z(device, result["x"]), device.dofs


@pytest.mark.parametrize("order", (1, 2))
def test_the_half_sector_reproduces_the_ring_and_the_cpu_route(tmp_path, request, prepared, cuda, order):
    cpu = _cpu_test()
    seam = cpu.SEAM
    ring_path, half_path = tmp_path/"ring.vol", tmp_path/"half.vol"
    half_mesh = cpu._yoke(half_path, 2, FOLD, (0.0, 0.004))
    cpu._yoke(ring_path, 2, FOLD, (-0.004, 0.0, 0.004), ring=True)
    images = dict(image="-z", image_cyclic=FOLD, cyclic_periodic_boundaries=seam)
    ring, ring_dofs = _device_mean(prepared, ring_path, order)
    half, half_dofs = _device_mean(prepared, half_path, order, **images)
    lone, _ = _device_mean(prepared, half_path, order)
    reference = cpu._solve(half_mesh, order, **images)
    cpu_half = float(reference["M_avg"][2])
    ring_error, cpu_error, lone_error = abs(half - ring)/abs(ring), abs(half - cpu_half)/abs(cpu_half), abs(lone - ring)/abs(ring)
    request.node.user_properties += [("half_dofs", int(half_dofs)), ("ring_dofs", int(ring_dofs)), ("ring_error", float(ring_error)),
                                     ("cpu_error", float(cpu_error)), ("lone_error", float(lone_error))]
    assert half_dofs == int(reference["ndof"]), (half_dofs, reference["ndof"])
    assert ring_error <= RING_REL, f"order={order} half={half!r} ring={ring!r} error={ring_error:.3e} limit={RING_REL:.0e}"
    assert cpu_error <= CPU_REL, f"order={order} device={half!r} cpu={cpu_half!r} error={cpu_error:.3e} limit={CPU_REL:.0e}"
    assert lone_error > 1e3*RING_REL, f"a half sector without images passes: {lone!r} ring {ring!r}"


def test_a_plus_midplane_with_rotations_refuses(prepared):
    with pytest.raises(NotImplementedError, match="midplane"):
        prepared._device_image_terms("+z", FOLD)
    terms = prepared._device_image_terms("-z", FOLD)
    assert len(terms["image_masks"]) == 2*FOLD - 1 and np.count_nonzero(terms["image_masks"] == 4) == FOLD, terms
    assert terms["image_rotations"][0] != 0.0, terms
