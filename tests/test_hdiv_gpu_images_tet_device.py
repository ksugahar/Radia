"""Mirror and cyclic images on the affine and curved TET device routes.

Guards: on tets touching a mirror plane or set off from it, in the corner of two planes, and on a sector whose rotated
copies touch it at the seam, RT0, BDM1 and BDM2, affine and curved P2, the device entries, the energy operator on the
dense and the compressed routes and the linear and nonlinear solves with images agree with the CPU route on identical
inputs, and so does the direct field of the same state, with and without images, at probes off, near and inside the
sources. Dropping the image sum, a sign, the graded rule of a touching P1 image, the Duffy rule of a touching curved image
or the inverse point map moves them far beyond the detectors.
Assumes: the CPU TET route is the reference on identical NGSolve inputs; this is parity of one discretization, not
independent accuracy. The REL values are PROPOSED detectors, not adopted criteria. The checked-mesh receipt gate is not
under test and is bypassed for scratch meshes. Optional GPU lane: skips visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set; a
configured lane loads the binding through the provenance-checked loader.
"""

import contextlib
import math
import os
from pathlib import Path

import numpy as np
import pytest

ENTRY_REL = 1e-10
OPERATOR_REL = 1e-9
ACA_OPERATOR_REL = 1e-5
SOLVE_REL = 1e-7
ACA_SOLVE_REL = 1e-5
FIELD_REL = 1e-9
INV_CHI = 1.0/499.0
# Off, near (0.5 mm above a top face of the boxes, 0.6 mm above the sector) and inside the sources; no image lands on a face.
_PROBES = np.array([[0.03, 0.015, 0.012], [0.02, 0.031, 0.0173], [0.0137, 0.0031, 0.0105], [-0.04, 0.01, 0.02], [0.0404, 0.0011, 0.0046]])
PAIR_SAMPLE = 20000
_BH = [[0.0, 0.0], [200.0, 0.8], [1000.0, 1.4], [5000.0, 1.7], [50000.0, 1.95], [500000.0, 2.6]]
_KEYS = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents", "tet_rule", "tri_rule",
         "tet_low_rule", "tri_low_rule", "far_factor", "family")
_CURVED_KEYS = ("cell_numbers", "face_numbers", "curve_gl", "curve_gw")
# Kuhn tets of a hex around its diagonal 0-6, in the order of a generated mesh, and the triangles of its faces on the same
# diagonals: bottom, top, outer, inner, then the two azimuthal faces.
_KUHN = ((0, 2, 1, 6), (0, 3, 2, 6), (0, 7, 3, 6), (0, 4, 7, 6), (0, 5, 4, 6), (0, 1, 5, 6))
_FACES = (((0, 3, 2), (0, 2, 1)), ((4, 5, 6), (4, 6, 7)), ((1, 2, 6), (1, 6, 5)), ((0, 4, 7), (0, 7, 3)),
          ((0, 1, 5), (0, 5, 4)), ((3, 7, 6), (3, 6, 2)))


def _occ(path, shape, maxh):
    from netgen.occ import OCCGeometry
    OCCGeometry(shape).GenerateMesh(maxh=maxh).Save(str(path))
    return path


def _tet_sector(path, cells, fold, inner=0.03, outer=0.05, half=0.004):
    """Kuhn tets of a polar sector over [-alpha/2, alpha/2]; rotated copies touch at the seam vertices. Every face is skin."""
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint, Pnt
    alpha = 2.0*math.pi/fold
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "yoke")
    skin = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    mesh.SetBCName(0, "skin")
    point = {}
    for it, ir, iz in np.ndindex(cells + 1, 2, 2):
        theta, radius = -0.5*alpha + it*alpha/cells, (inner, outer)[ir]
        point[ir, it, iz] = mesh.Add(MeshPoint(Pnt(radius*math.cos(theta), radius*math.sin(theta), (-half, half)[iz])))
    for it in range(cells):
        corners = [point[0, it, 0], point[1, it, 0], point[1, it + 1, 0], point[0, it + 1, 0],
                   point[0, it, 1], point[1, it, 1], point[1, it + 1, 1], point[0, it + 1, 1]]
        for tet in _KUHN:
            mesh.Add(Element3D(1, [corners[k] for k in tet]))
        for index, triangles in enumerate(_FACES):
            if index < 4 or (index == 4 and it == 0) or (index == 5 and it == cells - 1):
                for triangle in triangles:
                    mesh.Add(Element2D(skin, [corners[k] for k in triangle]))
    mesh.Save(str(path))
    return path


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


@pytest.fixture(scope="module")
def meshes(tmp_path_factory):
    from netgen.occ import Box, Pnt, Sphere
    root = tmp_path_factory.mktemp("tet-images")
    half_ball = Sphere(Pnt(0.0, 0.0, 0.0), 0.02)*Box(Pnt(0.0, -0.05, -0.05), Pnt(0.05, 0.05, 0.05))
    return dict(cut=_occ(root/"cut.vol", Box(Pnt(0.0, -0.01, -0.01), Pnt(0.02, 0.01, 0.01)), 0.012),
                corner=_occ(root/"corner.vol", Box(Pnt(0.0, 0.0, -0.01), Pnt(0.02, 0.02, 0.01)), 0.014),
                gap=_occ(root/"gap.vol", Box(Pnt(0.004, -0.01, -0.01), Pnt(0.024, 0.01, 0.01)), 0.012),
                far=_occ(root/"far.vol", Box(Pnt(0.2, -0.01, -0.01), Pnt(0.22, 0.01, 0.01)), 0.012),
                half_ball=_occ(root/"half_ball.vol", half_ball, 0.016),
                sector=_tet_sector(root/"sector.vol", 1, 8),
                sector2=_tet_sector(root/"sector2.vol", 2, 8))


def _terms(image=None, cyclic=None):
    from radia import hdiv_gpu
    return hdiv_gpu._device_image_terms(image, cyclic)


def _cpu_gram(fes, terms, curve_order):
    from radia.vim._vim import build_charge_gram
    kw = dict(image_masks=terms["image_masks"].tolist(), image_signs=terms["image_signs"].tolist())
    if "image_rotations" in terms:
        kw["image_rot_angle"] = terms["image_rotations"].tolist()
    return build_charge_gram(fes, eps=1e-10, curve_order=curve_order, _materialize_mass=False, _build_hmatrix=False, **kw)


def _case(path, order, curve_order, terms, cuda):
    """CPU charge map and Gram, and the device arguments of the same mesh; the images are added to the geometry."""
    import ngsolve as ng
    from radia import hdiv_gpu
    from radia.vim._vim import _curve_mesh
    with ng.TaskManager():
        setup = hdiv_gpu._prepare_device_components(path, order, cuda, aca=None, compression="dense", curve_order=curve_order)
        mesh = ng.Mesh(str(path))
        if curve_order == 2:
            _curve_mesh(mesh, 2)
        fes = ng.HDiv(mesh, order=order)
        B, gram, _ = _cpu_gram(fes, terms, 2 if curve_order == 2 else None)
    keys = _KEYS + (_CURVED_KEYS if curve_order == 2 else ())
    arguments = {key: setup["charge_gram"][key] for key in keys}
    arguments.update(terms)
    return arguments, B, gram, setup


def _entries(cuda, arguments, rows, cols):
    rows, cols = np.ascontiguousarray(rows, dtype=np.uint64), np.ascontiguousarray(cols, dtype=np.uint64)
    return np.concatenate([cuda._diagnostic_flat_charge_entries(**arguments, rows=rows[s:s + 65536], cols=cols[s:s + 65536],
                                                                backend="gpu") for s in range(0, len(rows), 65536)])


_ENTRY_CASES = [("cut", 0, 1, "-x", None), ("cut", 1, 1, "-x", None), ("cut", 2, 1, "-x", None), ("cut", 1, 1, "+x", None),
                ("cut", 2, 1, "+x", None), ("corner", 1, 1, "-x-y", None), ("corner", 2, 1, "+x-y", None),
                ("gap", 2, 1, "-x", None), ("far", 2, 1, "-x", None), ("sector", 0, 1, None, 8), ("sector", 1, 1, None, 8),
                ("sector", 2, 1, None, 8), ("sector2", 1, 1, None, 8), ("sector2", 2, 1, None, 8),
                ("half_ball", 1, 2, "-x", None), ("half_ball", 2, 2, "-x", None), ("half_ball", 2, 2, "+x", None),
                ("sector2", 1, 2, None, 8), ("sector2", 2, 2, None, 8)]


@pytest.mark.parametrize("name,order,curve_order,image,cyclic", _ENTRY_CASES)
def test_image_entries_match_cpu_charge_gram(meshes, request, cuda, name, order, curve_order, image, cyclic):
    arguments, _, gram, _ = _case(meshes[name], order, curve_order, _terms(image, cyclic), cuda)
    n = arguments["charge_host"].size
    rng = np.random.default_rng(20261105)
    upper = np.c_[np.triu_indices(n)]
    pairs = upper if len(upper) <= PAIR_SAMPLE else np.unique(np.r_[np.c_[np.arange(n), np.arange(n)],
                                                                    upper[rng.integers(0, len(upper), PAIR_SAMPLE)]], axis=0)
    rows, cols = pairs[:, 0], pairs[:, 1]
    actual = _entries(cuda, arguments, rows, cols)
    expected = np.array([gram.entry(int(i), int(j)) for i, j in pairs])
    diagonal = np.array([gram.entry(i, i) for i in range(n)])
    bare = {key: value for key, value in arguments.items() if not key.startswith("image_")}
    plain = _entries(cuda, bare, rows, cols)
    # A charge an image of opposite sign cancels has a roundoff diagonal: it is measured at the largest scale, and a
    # cancelled entry keeps the roundoff of its terms, the bound the CPU charge scale itself uses.
    scale = np.sqrt(np.maximum(np.abs(diagonal), 1e-6*np.abs(diagonal).max()))
    direct = np.abs(_entries(cuda, bare, np.arange(n), np.arange(n)))
    roundoff = 64*np.finfo(float).eps*(1 + arguments["image_masks"].size)*np.sqrt(direct[rows]*direct[cols])
    relative = np.maximum(np.abs(actual - expected) - roundoff, 0.0)/(scale[rows]*scale[cols])
    worst = int(np.argmax(relative))
    request.node.user_properties += [("charges", int(n)), ("pairs", len(pairs)), ("max_relative", float(relative.max()))]
    assert relative.max() <= ENTRY_REL, (
        f"{name} {image or cyclic} order={order} curve={curve_order} pair=({rows[worst]},{cols[worst]}) actual={actual[worst]!r} "
        f"expected={expected[worst]!r} relative={relative[worst]:.3e} limit={ENTRY_REL:.0e}")
    # Without the images the same pairs differ: a no-op image route cannot pass.
    moved = np.abs(plain - expected)/(scale[rows]*scale[cols])
    assert moved.max() > 1e4*max(ENTRY_REL, float(relative.max())), (
        f"{name} images are a no-op: moved={moved.max():.3e} relative={relative.max():.3e}")


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


@pytest.mark.parametrize("name,order,curve_order,compression,images", [
    ("cut", 2, 1, "dense", dict(image="-x")), ("cut", 1, 1, "aca", dict(image="-x")), ("corner", 2, 1, "dense", dict(image="-x-y")),
    ("gap", 2, 1, "aca", dict(image="-x")), ("sector", 2, 1, "dense", dict(image_cyclic=8)), ("sector2", 1, 1, "aca", dict(image_cyclic=8)),
    ("half_ball", 2, 2, "dense", dict(image="-x")), ("half_ball", 1, 2, "aca", dict(image="-x"))])
def test_image_operator_matches_cpu(meshes, request, cuda, name, order, curve_order, compression, images):
    from radia import hdiv_gpu
    terms = _terms(images.get("image"), images.get("image_cyclic"))
    import ngsolve as ng
    with ng.TaskManager():
        setup = hdiv_gpu._prepare_device_components(meshes[name], order, cuda, aca=None, compression=compression, curve_order=curve_order,
                                                    images=hdiv_gpu._DeviceImages(images.get("image"), images.get("image_cyclic")))
    _, B, gram, _ = _case(meshes[name], order, curve_order, terms, cuda)
    gram.build_exact_dense_normalized_gram(maximum_memory_mb=2048)
    nd = setup["dofs"]
    assert nd == B.shape[1] and len(setup["charge_gram"]["charge_host"]) == B.shape[0]
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    try:
        assembly = hdiv_gpu._checked_device_assembly(operator.assembly(), setup["charge_gram"])
        rng = np.random.default_rng(order)
        worst = 0.0
        for _ in range(4):
            x = rng.standard_normal(nd)
            expected = B.T @ np.asarray(gram.matvec_sym(np.ascontiguousarray(B @ x)))
            actual = operator.apply(np.ascontiguousarray(x))
            worst = max(worst, float(np.linalg.norm(actual - expected)/np.linalg.norm(expected)))
    finally:
        operator.close()
    limit = OPERATOR_REL if compression == "dense" else ACA_OPERATOR_REL
    rows = {(row["payload"], row["direction"]): row["bytes"] for row in assembly["transfer_ledger"]}
    request.node.user_properties += [("dofs", nd), ("max_relative", worst)]
    assert rows[("image_transforms", "H2D")] == 72*terms["image_masks"].size
    assert worst <= limit, f"{name} order={order} {compression} {images} relative={worst:.3e} limit={limit:.0e}"


def _cpu_solve(path, order, curve_order, **kw):
    import ngsolve as ng
    import radia as rad
    from radia import vim
    rad.UtiDelAll()
    with ng.TaskManager():
        mesh = ng.Mesh(str(path))
        extra = dict(curve_order=2) if curve_order == 2 else {}
        if curve_order == 2:
            from radia.vim._vim import _curve_mesh
            _curve_mesh(mesh, 2)
        result = vim.Solve(mesh, H_ext=ng.CF((0.0, 0.0, 1e5)), order=order, gram_eps=1e-10, tol=1e-11, **extra, **kw)
    return np.asarray(result["_m_coefficients"], dtype=np.float64), result


@pytest.mark.parametrize("name,order,curve_order,compression,images", [
    ("gap", 1, 1, "dense", dict(image="-x")), ("cut", 2, 1, "dense", dict(image="-x")), ("cut", 2, 1, "aca", dict(image="-x")),
    ("sector2", 2, 1, "dense", dict(image_cyclic=8)), ("sector2", 1, 1, "aca", dict(image_cyclic=8)),
    ("half_ball", 2, 2, "dense", dict(image="-x"))])
def test_linear_solve_matches_cpu(meshes, request, prepared, cuda, name, order, curve_order, compression, images):
    x_cpu, _ = _cpu_solve(meshes[name], order, curve_order, mu_r=500.0, **images)
    with prepared.PreparedHDivLinear(meshes[name], backend="gpu", order=order, gram_compression=compression, curve_order=curve_order,
                                     **_RECEIPT, **images) as device:
        solved = device.solve((0.0, 0.0, 1e5), INV_CHI, tolerance=1e-10, max_iterations=4000)
    relative = np.linalg.norm(solved["x"] - x_cpu)/np.linalg.norm(x_cpu)
    request.node.user_properties += [("solve_relative", float(relative)), ("iterations", int(solved["iterations"]))]
    limit = SOLVE_REL if compression == "dense" else ACA_SOLVE_REL
    assert relative <= limit, f"{name} order={order} curve={curve_order} {compression} {images} relative={relative:.3e} limit={limit:.0e}"


@pytest.mark.parametrize("name,order,images", [("cut", 2, dict(image="-x")), ("sector2", 1, dict(image_cyclic=8))])
def test_nonlinear_solve_matches_cpu(meshes, request, prepared, cuda, name, order, images):
    x_cpu, _ = _cpu_solve(meshes[name], order, 1, bh_table=_BH, nl_tol=1e-9, **images)
    with prepared.PreparedHDivNonlinear.from_bh_table(meshes[name], _BH, backend="gpu", order=order, **_RECEIPT, **images) as device:
        solved = device.solve((0.0, 0.0, 1e5), tolerance=1e-9, linear_tolerance=1e-10)
    relative = np.linalg.norm(solved["x"] - x_cpu)/np.linalg.norm(x_cpu)
    request.node.user_properties += [("solve_relative", float(relative)), ("newton", int(solved["newton_iterations"]))]
    assert relative <= SOLVE_REL, f"{name} order={order} {images} relative={relative:.3e} limit={SOLVE_REL:.0e}"


@pytest.mark.parametrize("name,order,curve_order,images", [
    ("gap", 1, 1, dict(image="-x")), ("gap", 2, 1, {}), ("cut", 2, 1, dict(image="-x")), ("corner", 1, 1, dict(image="-x-y")),
    ("sector2", 1, 1, dict(image_cyclic=8)), ("sector2", 2, 1, dict(image_cyclic=8)), ("half_ball", 1, 2, dict(image="-x")),
    ("half_ball", 2, 2, dict(image="-x")), ("half_ball", 2, 2, {})])
def test_field_of_the_cpu_state_matches_cpu(meshes, request, prepared, cuda, name, order, curve_order, images):
    from radia import vim
    x_cpu, cpu = _cpu_solve(meshes[name], order, curve_order, mu_r=500.0, **images)
    field_cpu = np.asarray(vim.FieldFromSolution(cpu, _PROBES, algorithm="direct"), float)
    with prepared.PreparedHDivLinear(meshes[name], backend="gpu", order=order, gram_compression="dense", curve_order=curve_order,
                                     **_RECEIPT, **images) as device:
        same_state = device.field(x_cpu, _PROBES)
    relative = np.linalg.norm(same_state["field"] - field_cpu)/np.linalg.norm(field_cpu)
    request.node.user_properties += [("field_relative", float(relative)), ("field_route", same_state["route"])]
    assert relative <= FIELD_REL, (f"{name} order={order} curve={curve_order} {images} field relative={relative:.3e} "
                                   f"actual={same_state['field'].tolist()} expected={field_cpu.tolist()}")
