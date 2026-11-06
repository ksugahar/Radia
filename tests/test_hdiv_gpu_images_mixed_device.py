"""Mirror and cyclic images on the connected tet, pyramid and hex device route.

Guards: the CPU route refuses images on a mixed mesh, so the reference is the CPU typed route on the explicitly meshed
copies. On a transition cube set off from its rotated copies and on a polar sector of hexes, pyramids and tets whose rotated
copies touch it at the seam, BDM1 and BDM2, the device entries with images equal the sum of the CPU entries over the copies,
alternating signs included, and the device operator and linear solve equal the dense CPU equation folded from those entries.
On those meshes and on a transition cube set off from or cut by its mirror, the dense and compressed linear solves and the
nonlinear solve reproduce the element means of the CPU solve of the explicit mesh, averaged over its copies, as closely as
those copies agree. A classification that ignores the image or a typed inner integral on the unmapped point moves them far
beyond the detectors. A space with a compressed trace, a connected periodic sector and a field request refuse.
Assumes: the CPU typed route on a mesh of explicit copies is the reference; a solve in the mirror-symmetric or cyclic
subspace is the same discretization only when the copies carry the same elements, so the copies are built from one cell
list. That route integrates a pair of two hosts from one side only, chosen by cost and then by number, so the explicit matrix
is not invariant under the group: the copies of one solve differ (recorded as reference_spread), and the subspace equation is
one row of that matrix, not its mean. Its difference from the explicit solve is therefore first order in that one-sided rule,
measured up to about the spread; only the folded equation is the same equation at the seam. The REL values and the factor on
the spread are PROPOSED detectors, not adopted criteria. The checked-mesh receipt gate is not under test and is
bypassed for scratch meshes. Optional GPU lane: skips visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set; a configured lane
loads the binding through the provenance-checked loader.
"""

import contextlib
import math
import os
from pathlib import Path

import numpy as np
import pytest

ENTRY_REL = 1e-9
OPERATOR_REL = 1e-9
SOLVE_REL = 1e-6
ACA_SOLVE_REL = 1e-5
SPREAD_FACTOR = 2.0
INV_CHI = 1.0/499.0
PAIR_SAMPLE = 1500
_BH = [[0.0, 0.0], [200.0, 0.8], [1000.0, 1.4], [5000.0, 1.7], [50000.0, 1.95], [500000.0, 2.6]]
_GEOMETRY = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents",
             "tet_rule", "tri_rule", "tet_low_rule", "tri_low_rule", "far_factor", "family",
             "hex_rules", "hex_bands", "host_types")
_FOLD = dict(disjoint=6, sector=8)


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


def _parts(mesh):
    ngmesh = mesh.ngmesh
    points = [tuple(float(value) for value in point.p) for point in ngmesh.Points()]
    return points, [[vertex.nr - 1 for vertex in element.vertices] for element in ngmesh.Elements3D()]


def _rotated_copies(points, cells, fold):
    """The copies of one cell list, disjoint in numbering: a rotated copy touches its neighbour without sharing a vertex."""
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import build_mesh
    all_points, all_cells = [], []
    for k in range(fold):
        c, s = math.cos(2.0*math.pi*k/fold), math.sin(2.0*math.pi*k/fold)
        offset = len(all_points)
        all_points += [(c*x - s*y, s*x + c*y, z) for x, y, z in points]
        all_cells += [[vertex + offset for vertex in cell] for cell in cells]
    return build_mesh(all_points, all_cells)


def _mirrored(points, cells):
    """The cell list and its x -> -x copy; vertices on the plane are shared, so the normal trace is continuous there."""
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import build_mesh, oriented
    all_points, index = list(points), {}
    for number, (x, y, z) in enumerate(points):
        if abs(x) < 1e-12:
            index[number] = number
        else:
            index[number] = len(all_points)
            all_points.append((-x, y, z))
    copies = []
    for cell in cells:
        image = [index[vertex] for vertex in cell]
        copies.append([image[i] for i in (1, 0, 3, 2, 5, 4, 7, 6)] if len(image) == 8 else oriented(all_points, image))
    return build_mesh(all_points, cells + copies)


def _save(mesh, path):
    mesh.ngmesh.Save(str(path))
    return path


@pytest.fixture(scope="module")
def meshes(tmp_path_factory):
    from validation_test.hmatrix_cuda.hdiv_mixed_meshes import build_mesh, mixed_grid, transition_cube
    root = tmp_path_factory.mktemp("mixed-images")
    points, cells = _parts(transition_cube())
    box = lambda x0: [(x0 + 0.02*x, 0.02*(y - 0.5), 0.02*(z - 0.5)) for x, y, z in points]
    alpha = 2.0*math.pi/_FOLD["sector"]
    polar = lambda radius, theta, z: (radius*math.cos(theta), radius*math.sin(theta), z)
    sector = mixed_grid((2, 1, 1), lengths=(0.02, alpha, 0.008), origin=(0.03, -0.5*alpha, -0.004), mapping=polar)
    target = dict(disjoint=_save(build_mesh(box(0.04), cells), root/"disjoint.vol"),
                  disjoint_copies=_save(_rotated_copies(box(0.04), cells, _FOLD["disjoint"]), root/"disjoint_copies.vol"),
                  sector=_save(sector, root/"sector.vol"),
                  gap=_save(build_mesh(box(0.004), cells), root/"gap.vol"),
                  gap_copies=_save(_mirrored(box(0.004), cells), root/"gap_copies.vol"),
                  cut=_save(build_mesh(box(0.0), cells), root/"cut.vol"),
                  cut_copies=_save(_mirrored(box(0.0), cells), root/"cut_copies.vol"))
    target["sector_copies"] = _save(_rotated_copies(*_parts(sector), _FOLD["sector"]), root/"sector_copies.vol")
    return target


@contextlib.contextmanager
def _unchecked_snapshot(path, *_):
    yield Path(path), "0"*64, "0"*64


_RECEIPT = dict(label_contract="unused.json", mesh_check_receipt="unused.receipt", mesh_check_receipt_sha256="0"*64,
                mesh_check_provenance={})


@pytest.fixture
def prepared(monkeypatch):
    from radia import hdiv_gpu
    monkeypatch.setattr(hdiv_gpu, "_checked_mesh_snapshot", _unchecked_snapshot)
    return hdiv_gpu


def _hosts(setup):
    table = setup["charge_gram"]
    return (table["cell_vertices"].size//81, table["face_vertices"].size//27,
            {(int(kind), int(host), tuple(int(e) for e in expo)): index for index, (kind, host, expo) in
             enumerate(zip(table["charge_kind"], table["charge_host"], table["charge_exponents"].reshape(-1, 3)))})


def _entries(cuda, arguments, rows, cols):
    rows, cols = np.ascontiguousarray(rows, dtype=np.uint64), np.ascontiguousarray(cols, dtype=np.uint64)
    return np.concatenate([cuda._diagnostic_flat_charge_entries(**arguments, rows=rows[s:s + 65536], cols=cols[s:s + 65536],
                                                                backend="gpu") for s in range(0, len(rows), 65536)])


@pytest.mark.parametrize("name,order,alternating", [("disjoint", 1, False), ("disjoint", 2, False), ("sector", 1, False),
                                                    ("sector", 2, False), ("sector", 1, True), ("sector", 2, True)])
def test_image_entries_equal_the_sum_over_the_explicit_copies(meshes, request, cuda, name, order, alternating):
    import ngsolve as ng
    from radia import hdiv_gpu
    from radia.vim._vim import _build_charge_gram_mixed
    fold = _FOLD[name]
    with ng.TaskManager():
        base = hdiv_gpu._prepare_device_components(meshes[name], order, cuda, aca=None, compression="dense")
        copies = hdiv_gpu._prepare_device_components(meshes[name + "_copies"], order, cuda, aca=None, compression="dense")
        _, gram, _, _ = _build_charge_gram_mixed(ng.HDiv(ng.Mesh(str(meshes[name + "_copies"])), order=order), build_hmatrix=False)
    assert base["family"] == "mixed" and copies["family"] == "mixed"
    table = base["charge_gram"]
    cells, faces, _ = _hosts(base)
    _, _, ring = _hosts(copies)
    # Copy k of a host follows the cell list: its lattice is the rotated base lattice, which keeps the reference monomials.
    for kind, count, key in ((0, cells, "cell_vertices"), (1, faces, "face_vertices")):
        lattice = table[key].reshape(count, -1, 3)
        for k in range(fold):
            c, s = math.cos(2.0*math.pi*k/fold), math.sin(2.0*math.pi*k/fold)
            moved = np.stack((c*lattice[..., 0] - s*lattice[..., 1], s*lattice[..., 0] + c*lattice[..., 1], lattice[..., 2]), axis=-1)
            other = copies["charge_gram"][key].reshape(fold*count, -1, 3)[k*count:(k + 1)*count]
            assert np.abs(moved - other).max() <= 1e-12, f"{name} copy {k} kind {kind}: {np.abs(moved - other).max():.3e}"
    n = table["charge_host"].size
    rng = np.random.default_rng(20261105)
    upper = np.c_[np.triu_indices(n)]
    pairs = upper if len(upper) <= PAIR_SAMPLE else np.unique(np.r_[np.c_[np.arange(n), np.arange(n)],
                                                                    upper[rng.integers(0, len(upper), PAIR_SAMPLE)]], axis=0)
    signs = [(-1.0)**k if alternating else 1.0 for k in range(fold)]
    expected = np.zeros(len(pairs))
    charge = lambda index: (int(table["charge_kind"][index]), int(table["charge_host"][index]), tuple(int(e) for e in table["charge_exponents"][3*index:3*index + 3]))
    copies_sum = lambda target, source: sum(signs[k]*gram.entry(ring[target], ring[source[0], source[1] + k*(faces if source[0] else cells), source[2]]) for k in range(fold))
    for index, (a, b) in enumerate(pairs):
        expected[index] = copies_sum(charge(a), charge(b))
        # A host against its own images has no lower host to integrate from: the device folds both directed blocks.
        if charge(a)[:2] == charge(b)[:2]:
            expected[index] = 0.5*(expected[index] + copies_sum(charge(b), charge(a)))
    device = {key: table[key] for key in _GEOMETRY}
    device.update(hdiv_gpu._device_image_terms(None, fold, alternating))
    actual = _entries(cuda, device, pairs[:, 0], pairs[:, 1])
    plain = _entries(cuda, {key: table[key] for key in _GEOMETRY}, pairs[:, 0], pairs[:, 1])
    direct = np.abs(_entries(cuda, {key: table[key] for key in _GEOMETRY}, np.arange(n), np.arange(n)))
    scale = np.sqrt(direct[pairs[:, 0]]*direct[pairs[:, 1]])
    relative = np.abs(actual - expected)/scale
    worst = int(np.argmax(relative))
    moved = np.abs(plain - expected)/scale
    request.node.user_properties += [("charges", int(n)), ("pairs", len(pairs)), ("max_relative", float(relative.max())), ("moved", float(moved.max()))]
    assert relative.max() <= ENTRY_REL, (f"{name} order={order} alternating={alternating} pair={pairs[worst].tolist()} {charge(pairs[worst, 0])} {charge(pairs[worst, 1])} actual={actual[worst]!r} "
                                         f"expected={expected[worst]!r} relative={relative[worst]:.3e} limit={ENTRY_REL:.0e}")
    assert moved.max() > 1e3*max(ENTRY_REL, float(relative.max())), f"{name} images are a no-op: moved={moved.max():.3e} relative={relative.max():.3e}"


def _element_means(field, mesh):
    import ngsolve as ng
    return np.column_stack([np.asarray(ng.Integrate(field[k], mesh, order=8, element_wise=True)) for k in range(3)])


def _copy_frames(means, cells, images):
    """The element means of each copy, mapped back onto the cells of copy 0."""
    if "image_cyclic" in images:
        frames = []
        for k in range(images["image_cyclic"]):
            c, s = math.cos(2.0*math.pi*k/images["image_cyclic"]), math.sin(2.0*math.pi*k/images["image_cyclic"])
            block = means[k*cells:(k + 1)*cells]
            frames.append(np.column_stack((c*block[:, 0] + s*block[:, 1], -s*block[:, 0] + c*block[:, 1], block[:, 2])))
        return np.array(frames)
    assert images == dict(image="-x"), images
    return np.array([means[:cells], means[cells:2*cells]*np.array([1.0, -1.0, -1.0])])


def _cpu_means(path, cells, order, applied, images, **kw):
    """The mean over the copies of one CPU solve of the explicit mesh, and the largest distance of a copy from that mean."""
    import ngsolve as ng
    import radia as rad
    from radia.vim import _solve
    rad.UtiDelAll()
    with ng.TaskManager():
        mesh = ng.Mesh(str(path))
        result = _solve.hdiv_demag_solve(mesh, H_ext=ng.CF(applied), order=order, gram_eps=1e-10, tol=1e-11, **kw)
        frames = _copy_frames(_element_means(result["gfM"], mesh), cells, images)
    mean = frames.mean(axis=0)
    return mean, max(float(np.linalg.norm(frame - mean)/np.linalg.norm(mean)) for frame in frames)


def _device_means(device, state):
    import ngsolve as ng
    field = ng.GridFunction(device.fespace)
    field.vec.FV().NumPy()[:] = state
    with ng.TaskManager():
        return _element_means(field, device.fespace.mesh)


_SOLVE_CASES = [("disjoint", 1, "dense", dict(image_cyclic=6), (0.0, 0.0, 1e5)), ("disjoint", 2, "aca", dict(image_cyclic=6), (0.0, 0.0, 1e5)),
                ("sector", 1, "dense", dict(image_cyclic=8), (0.0, 0.0, 1e5)), ("sector", 2, "dense", dict(image_cyclic=8), (0.0, 0.0, 1e5)),
                ("sector", 2, "aca", dict(image_cyclic=8), (0.0, 0.0, 1e5)), ("gap", 2, "dense", dict(image="-x"), (1e5, 0.0, 0.0)),
                ("cut", 1, "dense", dict(image="-x"), (1e5, 0.0, 0.0)), ("cut", 2, "dense", dict(image="-x"), (1e5, 0.0, 0.0))]


@pytest.mark.parametrize("name,order,compression,images,applied", _SOLVE_CASES)
def test_linear_solve_matches_the_explicit_copies(meshes, request, prepared, cuda, name, order, compression, images, applied):
    with prepared.PreparedHDivLinear(meshes[name], backend="gpu", order=order, gram_compression=compression, **_RECEIPT, **images) as device:
        solved = device.solve(applied, INV_CHI, tolerance=1e-10, max_iterations=4000)
        means = _device_means(device, solved["x"])
        rows = {(row["payload"], row["direction"]): row["bytes"] for row in device.stats()["assembly"]["transfer_ledger"]}
    reference, spread = _cpu_means(meshes[name + "_copies"], len(means), order, applied, images, mu_r=500.0)
    relative = np.linalg.norm(means - reference)/np.linalg.norm(reference)
    request.node.user_properties += [("solve_relative", float(relative)), ("reference_spread", spread), ("iterations", int(solved["iterations"]))]
    assert ("image_transforms", "H2D") in rows
    limit = (SOLVE_REL if compression == "dense" else ACA_SOLVE_REL) + SPREAD_FACTOR*spread
    assert relative <= limit, f"{name} order={order} {compression} {images} relative={relative:.3e} spread={spread:.3e} limit={limit:.3e}"


def _folded_gram(base, copies, gram, fold, signs):
    """The Gram of the sector in its subspace, from the exact dense CPU Gram of the explicit copies, on the rule the entries test holds."""
    table = base["charge_gram"]
    cells, faces, _ = _hosts(base)
    _, _, ring = _hosts(copies)
    kinds, hosts = np.asarray(table["charge_kind"], dtype=np.int64), np.asarray(table["charge_host"], dtype=np.int64)
    keys = [(int(kind), int(host), tuple(int(e) for e in table["charge_exponents"][3*index:3*index + 3])) for index, (kind, host) in enumerate(zip(kinds, hosts))]
    first = np.array([ring[key] for key in keys])
    rows = np.zeros((len(keys), len(keys)))
    for column, (kind, host, exponents) in enumerate(keys):
        lifted = np.zeros(fold*len(keys))
        for k in range(fold):
            lifted[ring[kind, host + k*(faces if kind else cells), exponents]] = signs[k]
        rows[:, column] = np.asarray(gram.matvec_sym(lifted))[first]
    folded = np.triu(rows) + np.triu(rows, 1).T
    same = (kinds[:, None] == kinds[None, :]) & (hosts[:, None] == hosts[None, :])
    folded[same] = 0.5*(rows + rows.T)[same]
    return folded


@pytest.mark.parametrize("order,alternating", [(1, False), (1, True), (2, False)])
def test_operator_and_solve_equal_the_folded_cpu_equation(meshes, request, prepared, cuda, order, alternating):
    import ngsolve as ng
    from scipy.sparse import csr_matrix
    from radia.vim._vim import _build_charge_gram_mixed
    fold = _FOLD["sector"]
    applied = (0.0, 0.0, 1e5)
    with ng.TaskManager():
        base = prepared._prepare_device_components(meshes["sector"], order, cuda, aca=None, compression="dense",
                                                   images=prepared._DeviceImages(None, fold, alternating))
        copies = prepared._prepare_device_components(meshes["sector_copies"], order, cuda, aca=None, compression="dense")
        B, _, M, _ = _build_charge_gram_mixed(ng.HDiv(ng.Mesh(str(meshes["sector"])), order=order), build_hmatrix=False)
        _, gram, _, _ = _build_charge_gram_mixed(ng.HDiv(ng.Mesh(str(meshes["sector_copies"])), order=order), build_hmatrix=False)
        gram.build_exact_dense_normalized_gram(maximum_memory_mb=2048)
    folded = _folded_gram(base, copies, gram, fold, [(-1.0)**k if alternating else 1.0 for k in range(fold)])
    B = csr_matrix(B)
    nd = base["dofs"]
    assert nd == B.shape[1] and folded.shape[0] == B.shape[0]
    operator = cuda._ResidentHDiv.from_flat_charges(**base["charge_gram"], dofs=nd, backend="gpu")
    try:
        rng = np.random.default_rng(order)
        worst = 0.0
        for _ in range(4):
            x = rng.standard_normal(nd)
            expected = B.T @ (folded @ (B @ x))
            worst = max(worst, float(np.linalg.norm(operator.apply(np.ascontiguousarray(x)) - expected)/np.linalg.norm(expected)))
    finally:
        operator.close()
    reference = np.linalg.solve((B.T @ (folded @ B.toarray())) + INV_CHI*csr_matrix(M).toarray(), base["rhs_operator"] @ np.array(applied))
    with prepared.PreparedHDivLinear(meshes["sector"], backend="gpu", order=order, gram_compression="dense", image_cyclic=fold,
                                     image_cyclic_alternating=alternating, **_RECEIPT) as device:
        solved = device.solve(applied, INV_CHI, tolerance=1e-10, max_iterations=4000)
    relative = float(np.linalg.norm(solved["x"] - reference)/np.linalg.norm(reference))
    request.node.user_properties += [("dofs", nd), ("operator_relative", worst), ("solve_relative", relative), ("iterations", int(solved["iterations"]))]
    assert worst <= OPERATOR_REL, f"order={order} alternating={alternating} operator relative={worst:.3e} limit={OPERATOR_REL:.0e}"
    assert relative <= SOLVE_REL, f"order={order} alternating={alternating} solve relative={relative:.3e} limit={SOLVE_REL:.0e}"


@pytest.mark.parametrize("name,order,images,applied", [("sector", 1, dict(image_cyclic=8), (0.0, 0.0, 1e5)), ("cut", 2, dict(image="-x"), (1e5, 0.0, 0.0))])
def test_nonlinear_solve_matches_the_explicit_copies(meshes, request, prepared, cuda, name, order, images, applied):
    with prepared.PreparedHDivNonlinear.from_bh_table(meshes[name], _BH, backend="gpu", order=order, **_RECEIPT, **images) as device:
        solved = device.solve(applied, tolerance=1e-9, linear_tolerance=1e-10)
        means = _device_means(device, solved["x"])
    reference, spread = _cpu_means(meshes[name + "_copies"], len(means), order, applied, images, bh_table=_BH, nl_tol=1e-11)
    relative = np.linalg.norm(means - reference)/np.linalg.norm(reference)
    request.node.user_properties += [("solve_relative", float(relative)), ("reference_spread", spread), ("newton", int(solved["newton_iterations"]))]
    limit = SOLVE_REL + SPREAD_FACTOR*spread
    assert relative <= limit, f"{name} order={order} {images} relative={relative:.3e} spread={spread:.3e} limit={limit:.3e}"


def test_the_kernel_expansion_refuses_images(meshes, cuda):
    """Guards: the expansion reads the direct term only; images fail in Python and in native code and are never dropped.

    Assumes: the route's own parameters reach the refusal; its accuracy is not read here.
    """
    from radia import hdiv_gpu
    aca = dict(hdiv_gpu.GRAM_ACA_PRESET, expansion_order=6)
    with pytest.raises(NotImplementedError, match="no image route"):
        hdiv_gpu._prepare_device_components(meshes["sector"], 1, cuda, aca=aca, compression="hca",
                                            images=hdiv_gpu._DeviceImages(image_cyclic=8))
    setup = hdiv_gpu._prepare_device_components(meshes["sector"], 1, cuda, aca=aca, compression="hca")
    with pytest.raises(ValueError, match="no image route"):
        cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], **hdiv_gpu._device_image_terms(image_cyclic=8),
                                             dofs=setup["dofs"], backend="gpu")


def test_a_compressed_trace_a_periodic_sector_and_a_field_request_refuse(meshes, prepared, cuda):
    with pytest.raises(NotImplementedError, match="compressed normal trace"):
        prepared.PreparedHDivLinear(meshes["cut"], backend="gpu", order=1, gram_compression="dense", **_RECEIPT, image="+x")
    with pytest.raises(NotImplementedError, match="pure HEX"):
        prepared._DeviceImages(None, 6, False, ("periodic_min", "periodic_max")).checked_terms(None, 1, prepared._MIXED)
    with prepared.PreparedHDivLinear(meshes["gap"], backend="gpu", order=1, gram_compression="dense", **_RECEIPT, image="-x") as device:
        with pytest.raises(ValueError, match="no route"):
            device.field(np.zeros(device.dofs), np.array([[0.05, 0.0, 0.0]]))
