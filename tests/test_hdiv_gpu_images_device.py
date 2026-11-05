"""Mirror and cyclic images on the pure-HEX device route.

Guards: on boxes with one or two mirror planes (touching or with a gap, affine or distorted), and on a polar sector
whose rotated copies touch it at the seam (linear and Q2), BDM1 and BDM2, the device entries, the energy operator on the
dense and the compressed routes, the linear and nonlinear solves and the direct field with images agree with the CPU
route on identical inputs. Dropping the image sum, a sign, the symmetrized image block or the inverse point map moves them
far beyond the detectors. The self branch of an on-plane host under its image is not separated, affine or trapezoid: a
point on the plane is its own image, so both branches read one point. A sector of several cells is compared with the
direct route of the explicitly meshed ring, since the CPU image dispatch reads those touching pairs through the inverse
rotation of the source and selects a far rule; that CPU entry error stays visible in the recorded properties. Every
list the CPU route refuses is refused here. A probe whose image lands on a source face raises on the closed-form route.
Assumes: the CPU hex route is the reference on identical NGSolve inputs, and the ring route on its rotated copies; this
is parity of one discretization, not independent accuracy. The REL values are PROPOSED detectors, not adopted criteria.
The checked-mesh receipt gate is not under test and is bypassed for scratch meshes. Optional GPU lane: skips visibly
unless RADIA_HDIV_NATIVE_TEST_DIR is set; a configured lane loads the binding through the provenance-checked loader.
"""

import contextlib
import math
import os
from pathlib import Path

import numpy as np
import pytest

ENTRY_REL = 1e-10
# The CPU dispatch reads a rotated pair through the inverse rotation of its source, so a pair of a one-cell sector near
# a rule band may take another rule there; the ring test holds the device to the direct route on the same pairs.
SECTOR_ENTRY_REL = 1e-8
OPERATOR_REL = 1e-9
ACA_OPERATOR_REL = 1e-5
RING_REL = 1e-5
FIELD_REL = 1e-9
SOLVE_REL = 1e-7
# The CPU Newton stops on its residual: at 1e-9 it stands 1.6e-7 from its own converged state on the HEX -x cut, and two
# runs on one input can take different paths there. At 1e-11 it is converged, so the reference is that state.
REFERENCE_NL_TOL = 1e-11
ACA_SOLVE_REL = 1e-5
INV_CHI = 1.0/499.0
SEAM = ("periodic_min", "periodic_max")
# Off every source plane of the cases: the closed forms are singular on a face (see the refusal test).
_PROBES = np.array([[0.05, 0.0, 0.0], [0.02, 0.03, 0.012], [0.03, -0.02, 0.02], [0.0, 0.0, 0.03], [0.028, 0.0, 0.0]])
_BH = [[0.0, 0.0], [200.0, 0.8], [1000.0, 1.4], [5000.0, 1.7], [50000.0, 1.95], [500000.0, 2.6]]
_KEYS = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents", "tet_rule", "tri_rule",
         "tet_low_rule", "tri_low_rule", "far_factor", "family", "hex_rules", "hex_bands")


def _hex_faces(corners, boundary):
    return [[corners[k] for k in quad] for quad, keep in zip(
        ((3, 2, 1, 0), (4, 5, 6, 7), (0, 1, 5, 4), (2, 3, 7, 6), (3, 0, 4, 7), (1, 2, 6, 5)), boundary) if keep]


def _box(path, lo, hi, cells, mapping=None):
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint, Pnt
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "iron")
    face = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    mesh.SetBCName(0, "skin")
    axes = [np.linspace(lo[k], hi[k], cells[k] + 1) for k in range(3)]
    points = {}
    for i, j, k in np.ndindex(*(n + 1 for n in cells)):
        xyz = (float(axes[0][i]), float(axes[1][j]), float(axes[2][k]))
        points[i, j, k] = mesh.Add(MeshPoint(Pnt(*(mapping(*xyz) if mapping else xyz))))
    for i, j, k in np.ndindex(*cells):
        corners = [points[i + a, j + b, k + c] for c in (0, 1) for a, b in ((0, 0), (1, 0), (1, 1), (0, 1))]
        mesh.Add(Element3D(1, corners))
        for quad in _hex_faces(corners, (k == 0, k == cells[2] - 1, j == 0, j == cells[1] - 1, i == 0, i == cells[0] - 1)):
            mesh.Add(Element2D(face, quad))
    mesh.Save(str(path))
    return path


def _sector(path, cells, fold, copies=1, inner=0.03, outer=0.05, half=0.004):
    """Polar HEX sector over [-alpha/2, alpha/2], or the ring of `fold` rotated copies sharing their seam nodes."""
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint, Pnt
    alpha, total, closed = 2.0*math.pi/fold, cells*copies, copies == fold
    count = total if closed else total + 1
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "yoke")
    skin = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    mesh.SetBCName(0, "skin")
    if not closed:
        low, high = mesh.Add(FaceDescriptor(surfnr=2, domin=1, domout=0, bc=2)), mesh.Add(FaceDescriptor(surfnr=3, domin=1, domout=0, bc=3))
        mesh.SetBCName(1, SEAM[0])
        mesh.SetBCName(2, SEAM[1])
    point = {}
    for it, ir, iz in np.ndindex(count, 2, 2):
        theta, radius = -0.5*alpha + it*alpha/cells, (inner, outer)[ir]
        point[ir, it, iz] = mesh.Add(MeshPoint(Pnt(radius*math.cos(theta), radius*math.sin(theta), (-half, half)[iz])))
    for it in range(total):
        next_it = (it + 1) % count
        corners = [point[0, it, 0], point[1, it, 0], point[1, next_it, 0], point[0, next_it, 0],
                   point[0, it, 1], point[1, it, 1], point[1, next_it, 1], point[0, next_it, 1]]
        mesh.Add(Element3D(1, corners))
        for quad in ((3, 2, 1, 0), (4, 5, 6, 7), (1, 2, 6, 5), (3, 0, 4, 7)):
            mesh.Add(Element2D(skin, [corners[k] for k in quad]))
        if not closed and it == 0:
            mesh.Add(Element2D(low, [corners[k] for k in (0, 1, 5, 4)]))
        if not closed and it == total - 1:
            mesh.Add(Element2D(high, [corners[k] for k in (2, 3, 7, 6)]))
    mesh.Save(str(path))
    return path


def _shear(x, y, z):
    return (x*(1.0 + 8.0*y), y + 0.2*x, z*(1.0 + 6.0*x))


def _trapezoid(x, y, z):
    """Keeps the plane x = 0 and makes the host on it a trapezoid there, so its self rule is not the affine one."""
    return (x*(1.0 + 8.0*y), y + 0.2*x, z*(1.0 + 30.0*y))


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


@pytest.fixture(scope="module")
def meshes(tmp_path_factory):
    root = tmp_path_factory.mktemp("images")
    return dict(
        cut=_box(root/"cut.vol", (0.0, -0.01, -0.01), (0.02, 0.01, 0.01), (2, 1, 1)),
        corner=_box(root/"corner.vol", (0.0, 0.0, -0.01), (0.02, 0.02, 0.01), (2, 1, 1)),
        gap=_box(root/"gap.vol", (0.004, -0.01, -0.01), (0.024, 0.01, 0.01), (2, 1, 1)),
        distorted=_box(root/"distorted.vol", (0.0, -0.01, -0.01), (0.02, 0.01, 0.01), (2, 1, 1), _shear),
        trapezoid=_box(root/"trapezoid.vol", (0.0, -0.01, -0.01), (0.02, 0.01, 0.01), (2, 1, 1), _trapezoid),
        sector=_sector(root/"sector.vol", 1, 8),
        sector2=_sector(root/"sector2.vol", 2, 8),
        ring1=_sector(root/"ring1.vol", 1, 8, copies=8),
        ring2=_sector(root/"ring2.vol", 2, 8, copies=8),
        poles=_box(root/"poles.vol", (0.032, -0.006, 0.004), (0.048, 0.006, 0.010), (2, 1, 1)))


def _terms(image=None, cyclic=None, alternating=False):
    from radia import hdiv_gpu
    return hdiv_gpu._device_image_terms(image, cyclic, alternating)


def _cpu_kw(terms):
    kw = dict(image_masks=terms["image_masks"].tolist(), image_signs=terms["image_signs"].tolist())
    if "image_rotations" in terms:
        kw["image_rot_angle"] = terms["image_rotations"].tolist()
    return kw


def _hex_rule_rows(fes):
    from radia.vim import _vim
    glout, glin = _vim._hex_gram_rule_counts(fes)
    rules = [_vim._g01(glout), _vim._g01(glin), _vim._g01(_vim._HEX_GLPAIR_N), _vim._g01(_vim._HEX_GLPAIR_N)]
    if fes.globalorder == 1:
        rules += [_vim._g01(_vim._HEX_GLNEAR_P1_N), _vim._g01(_vim._HEX_GLIN_SELF_P1_N)]
    return np.ascontiguousarray(np.concatenate([np.c_[np.full(len(x), float(i)), np.asarray(x), np.asarray(w)]
                                                for i, (x, w) in enumerate(rules)]))


def _gram_case(path, order, terms, *, seam=(), deform=None):
    """CPU Gram and the device arguments of the same charge basis and lattices."""
    import ngsolve as ng
    from radia.vim._vim import _build_charge_gram_hex, _charge_basis_hex, _curve_mesh
    mesh = ng.Mesh(str(path))
    deformation = None
    if deform is not None:
        _curve_mesh(mesh, 2)
        deformation = deform(mesh)
    fes = ng.HDiv(mesh, order=order)
    with ng.TaskManager():
        B, gram, _, _ = _build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False, materialize_mass=False,
                                               excluded_boundaries=seam, **_cpu_kw(terms))
        basis = _charge_basis_hex(fes, cob_quad=max(3, order + 1), materialize_mass=False, excluded_boundaries=seam)
    empty = np.zeros((0, 4))
    arguments = dict(cell_vertices=np.ascontiguousarray(basis["cell_nodes"], dtype=np.float64).reshape(-1),
                     face_vertices=np.ascontiguousarray(basis["face_nodes"], dtype=np.float64).reshape(-1),
                     charge_host=np.ascontiguousarray(basis["host"], dtype=np.int32),
                     charge_kind=np.ascontiguousarray(basis["kind"], dtype=np.int32),
                     charge_exponents=np.ascontiguousarray(basis["expo"], dtype=np.int32).reshape(-1),
                     tet_rule=empty, tri_rule=empty, tet_low_rule=empty, tri_low_rule=empty, far_factor=1e30,
                     family="hex", hex_rules=_hex_rule_rows(fes), hex_bands=np.array([1.0, 1.0]), **terms)
    return arguments, gram, deformation


def _entries(cuda, arguments, rows, cols):
    rows, cols = np.ascontiguousarray(rows, dtype=np.uint64), np.ascontiguousarray(cols, dtype=np.uint64)
    return np.concatenate([cuda._diagnostic_flat_charge_entries(**arguments, rows=rows[s:s + 65536], cols=cols[s:s + 65536],
                                                                backend="gpu") for s in range(0, len(rows), 65536)])


def _radial_q2(mesh):
    import ngsolve as ng
    deformation = ng.GridFunction(ng.VectorH1(mesh, order=2))
    square = ng.x*ng.x + ng.y*ng.y
    deformation.Interpolate(ng.CF((60.0*ng.x*square, 60.0*ng.y*square, 8.0*ng.z*square)))
    mesh.SetDeformation(deformation)
    return deformation


_ENTRY_CASES = [("cut", "-x", None, False, (), None), ("cut", "+x", None, False, (), None), ("corner", "-x-y", None, False, (), None),
                ("corner", "+x-y", None, False, (), None), ("gap", "-x", None, False, (), None),
                ("distorted", "-x", None, False, (), None), ("trapezoid", "-x", None, False, (), None), ("poles", None, 6, True, (), None),
                ("sector", None, 8, False, SEAM, None), ("sector", None, 8, False, SEAM, _radial_q2)]


@pytest.mark.parametrize("order", (1, 2))
@pytest.mark.parametrize("name,image,cyclic,alternating,seam,deform", _ENTRY_CASES)
def test_image_entries_match_cpu_charge_gram(meshes, request, cuda, name, image, cyclic, alternating, seam, deform, order):
    arguments, gram, _deformation = _gram_case(meshes[name], order, _terms(image, cyclic, alternating), seam=seam, deform=deform)
    n = arguments["charge_host"].size
    rows, cols = np.triu_indices(n)
    actual = _entries(cuda, arguments, rows, cols)
    expected = np.array([gram.entry(int(i), int(j)) for i, j in zip(rows, cols)])
    diagonal = np.array([gram.entry(i, i) for i in range(n)])
    bare = {key: arguments[key] for key in _KEYS}
    plain = _entries(cuda, bare, rows, cols)
    # A charge an image of opposite sign cancels has a roundoff diagonal: it is measured at the largest scale, and a
    # cancelled entry keeps the roundoff of its terms, the bound the CPU charge scale itself uses.
    scale = np.sqrt(np.maximum(np.abs(diagonal), 1e-6*np.abs(diagonal).max()))
    direct = np.abs(plain[rows == cols])
    roundoff = 64*np.finfo(float).eps*(1 + arguments["image_masks"].size)*np.sqrt(direct[rows]*direct[cols])
    relative = np.maximum(np.abs(actual - expected) - roundoff, 0.0)/(scale[rows]*scale[cols])
    worst = int(np.argmax(relative))
    limit = SECTOR_ENTRY_REL if cyclic else ENTRY_REL
    request.node.user_properties += [("charges", int(n)), ("max_relative", float(relative.max()))]
    assert relative.max() <= limit, (
        f"{name} {image or cyclic} order={order} pair=({rows[worst]},{cols[worst]}) actual={actual[worst]!r} "
        f"expected={expected[worst]!r} relative={relative[worst]:.3e} limit={limit:.0e}")
    # Without the images the same pairs differ: a no-op image route cannot pass.
    moved = np.abs(plain - expected)/(scale[rows]*scale[cols])
    assert moved.max() > 1e4*max(ENTRY_REL, float(relative.max())), (
        f"{name} images are a no-op: moved={moved.max():.3e} relative={relative.max():.3e}")


@pytest.mark.parametrize("order", (1, 2))
@pytest.mark.parametrize("sector,ring_name", (("sector", "ring1"), ("sector2", "ring2")))
def test_sector_cells_match_the_explicit_ring(meshes, request, cuda, sector, ring_name, order):
    """Guards: touching pairs of a one- and a two-cell sector take the rule of the explicit ring.

    Assumes: the ring's direct route on exactly rotated copies is the reference; both are this discretization.
    """
    fold = 8
    arguments, cpu, _ = _gram_case(meshes[sector], order, _terms(None, fold), seam=SEAM)
    import ngsolve as ng
    from radia.vim._vim import _build_charge_gram_hex, _charge_basis_hex
    ring = ng.HDiv(ng.Mesh(str(meshes[ring_name])), order=order)
    with ng.TaskManager():
        _, direct, _, _ = _build_charge_gram_hex(ring, eps=1e-10, build_hmatrix=False, materialize_mass=False)
        basis = _charge_basis_hex(ring, cob_quad=max(3, order + 1), materialize_mass=False)
    host, kind = arguments["charge_host"], arguments["charge_kind"]
    cells, faces = int(kind.tolist().count(0)), int(kind.tolist().count(1))
    per_cell, per_face = int(np.sum((kind == 0) & (host == 0))), int(np.sum((kind == 1) & (host == 0)))
    cell_nodes = arguments["cell_vertices"].reshape(-1, 27, 3)
    face_nodes = arguments["face_vertices"].reshape(-1, 9, 3)
    ring_cells = np.asarray(basis["cell_nodes"]).reshape(-1, 27, 3)
    ring_faces = np.asarray(basis["face_nodes"]).reshape(-1, 9, 3)
    ring_cell_hosts = len(ring_cells)

    def rotated(nodes, targets, k):
        c, s = math.cos(k*2*math.pi/fold), math.sin(k*2*math.pi/fold)
        image = nodes @ np.array(((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0))).T
        found = [int(np.argmin(np.abs(targets - row).max(axis=(1, 2)))) for row in image]
        assert max(np.abs(targets[j] - row).max() for j, row in zip(found, image)) < 1e-12
        return found
    maps = [(rotated(cell_nodes, ring_cells, k), rotated(face_nodes, ring_faces, k)) for k in range(fold)]

    def ring_charge(a, k):
        h = int(host[a])
        if kind[a] == 0:
            return maps[k][0][h]*per_cell + (a - h*per_cell)
        return ring_cell_hosts*per_cell + maps[k][1][h]*per_face + (a - cells - h*per_face)
    n = host.size
    rows, cols = np.triu_indices(n)
    expected = np.array([sum(direct.entry(ring_charge(int(a), 0), ring_charge(int(b), k)) for k in range(fold))
                         for a, b in zip(rows, cols)])
    actual = _entries(cuda, arguments, rows, cols)
    cpu_values = np.array([cpu.entry(int(a), int(b)) for a, b in zip(rows, cols)])
    diagonal = expected[rows == cols]
    scale = np.sqrt(np.abs(diagonal))
    relative = np.abs(actual - expected)/(scale[rows]*scale[cols])
    cpu_relative = np.abs(cpu_values - expected)/(scale[rows]*scale[cols])
    worst = int(np.argmax(relative))
    request.node.user_properties += [("charges", int(n)), ("device_vs_ring", float(relative.max())),
                                     ("cpu_vs_ring", float(cpu_relative.max()))]
    assert relative.max() <= RING_REL, (f"{sector} order={order} pair=({rows[worst]},{cols[worst]}) actual={actual[worst]!r} "
                                        f"expected={expected[worst]!r} relative={relative[worst]:.3e} limit={RING_REL:.0e}")


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


@pytest.mark.parametrize("name,order,compression,images", [
    ("cut", 2, "dense", dict(image="-x")), ("cut", 1, "aca", dict(image="-x")), ("corner", 2, "dense", dict(image="-x-y")),
    ("gap", 2, "aca", dict(image="-x")), ("poles", 2, "dense", dict(image_cyclic=6)), ("poles", 2, "aca", dict(image_cyclic=6)),
    ("poles", 1, "dense", dict(image_cyclic=6))])
def test_image_operator_matches_cpu(meshes, request, cuda, name, order, compression, images):
    import ngsolve as ng
    from radia import hdiv_gpu
    from radia.vim._vim import _build_charge_gram_hex
    with ng.TaskManager():
        setup = hdiv_gpu._prepare_device_components(meshes[name], order, cuda, aca=None, compression=compression,
                                                    images=hdiv_gpu._DeviceImages(images.get("image"), images.get("image_cyclic")))
        fes = ng.HDiv(ng.Mesh(str(meshes[name])), order=order)
        B, gram, _, _ = _build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False, materialize_mass=False,
                                               **_cpu_kw(_terms(images.get("image"), images.get("image_cyclic"))))
    gram.build_exact_dense_normalized_gram(maximum_memory_mb=512)
    nd = setup["dofs"]
    assert nd == B.shape[1] and len(setup["charge_gram"]["charge_host"]) == B.shape[0]
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    try:
        assembly = hdiv_gpu._checked_device_assembly(operator.assembly(), setup["charge_gram"])
        assert assembly["route"] == "device-hex"
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
    request.node.user_properties += [("dofs", nd), ("max_relative", worst)]
    assert worst <= limit, f"{name} order={order} {compression} {images} relative={worst:.3e} limit={limit:.0e}"


_SOLVE_CASES = [("gap", 1, "dense", dict(image="-x")), ("cut", 2, "dense", dict(image="-x")),
                ("poles", 1, "dense", dict(image_cyclic=6)), ("poles", 2, "dense", dict(image_cyclic=6)),
                ("poles", 2, "aca", dict(image_cyclic=6)), ("cut", 2, "aca", dict(image="-x")),
                ("distorted", 1, "dense", dict(image="-x")),
                ("poles", 2, "dense", dict(image_cyclic=6, image_cyclic_alternating=True))]


def _cpu_equation(path, order, terms, rhs):
    """The dense CPU equation of one image list on the device's mass rule E^T W E, IntegrationRule(HEX, 2*max(6, 2p+2))."""
    import ngsolve as ng
    from scipy.sparse import coo_matrix, csr_matrix
    from radia.vim._vim import _build_charge_gram_hex
    with ng.TaskManager():
        fes = ng.HDiv(ng.Mesh(str(path)), order=order)
        B, gram, _, _ = _build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False, materialize_mass=False, **_cpu_kw(terms))
        u, v = fes.TnT()
        mass = ng.BilinearForm(fes)
        mass += u*v*ng.dx(intrules={ng.HEX: ng.IntegrationRule(ng.HEX, 2*max(6, 2*order + 2))})
        mass.Assemble()
    gram.build_exact_dense_normalized_gram(maximum_memory_mb=256)
    B = csr_matrix(B)
    columns = [B.T @ np.asarray(gram.matvec_sym(np.ascontiguousarray(B[:, k].toarray().ravel()))) for k in range(fes.ndof)]
    rows, cols, values = mass.mat.COO()
    M = coo_matrix((values, (rows, cols)), shape=(fes.ndof, fes.ndof)).toarray()
    return np.linalg.solve(np.column_stack(columns) + INV_CHI*M, rhs)


def _cpu_solve(path, order, **kw):
    import ngsolve as ng
    import radia as rad
    from radia import vim
    rad.UtiDelAll()
    with ng.TaskManager():
        mesh = ng.Mesh(str(path))
        result = vim.Solve(mesh, H_ext=ng.CF((0.0, 0.0, 1e5)), order=order, gram_eps=1e-10, tol=1e-11, **kw)
    return result, np.asarray(result["_m_coefficients"], dtype=np.float64)


@pytest.mark.parametrize("name,order,compression,images", _SOLVE_CASES)
def test_linear_solve_and_field_match_cpu(meshes, request, prepared, cuda, name, order, compression, images):
    import ngsolve as ng
    from radia import vim
    cpu, x_cpu = _cpu_solve(meshes[name], order, mu_r=500.0, **images)
    field_cpu = np.asarray(vim.FieldFromSolution(cpu, _PROBES, algorithm="direct"), float)
    with prepared.PreparedHDivLinear(meshes[name], backend="gpu", order=order, gram_compression=compression,
                                     **_RECEIPT, **images) as device:
        solved = device.solve((0.0, 0.0, 1e5), INV_CHI, tolerance=1e-10, max_iterations=4000)
        same_state = device.field(x_cpu, _PROBES)
        rows = {(row["payload"], row["direction"]): row["bytes"] for row in device.stats()["assembly"]["transfer_ledger"]}
    terms = prepared._device_image_terms(images.get("image"), images.get("image_cyclic"))
    reference = x_cpu
    if name == "distorted":
        # On a non-affine host the two routes integrate the mass on different rules: the CPU equation of the device's.
        with ng.TaskManager():
            setup = prepared._prepare_device_components(meshes[name], order, cuda, aca=None, compression=compression,
                                                        images=prepared._DeviceImages(images.get("image"), images.get("image_cyclic")))
        reference = _cpu_equation(meshes[name], order, terms, setup["rhs_operator"] @ np.array([0.0, 0.0, 1e5]))
        request.node.user_properties += [("cpu_solve_vs_equation", float(np.linalg.norm(x_cpu - reference)/np.linalg.norm(reference)))]
    relative = np.linalg.norm(solved["x"] - reference)/np.linalg.norm(reference)
    field_relative = np.linalg.norm(same_state["field"] - field_cpu)/np.linalg.norm(field_cpu)
    request.node.user_properties += [("solve_relative", float(relative)), ("field_relative", float(field_relative)),
                                     ("iterations", int(solved["iterations"])), ("field_route", same_state["route"])]
    assert rows[("image_transforms", "H2D")] == 72*terms["image_masks"].size
    limit = SOLVE_REL if compression == "dense" else ACA_SOLVE_REL
    assert relative <= limit, f"{name} order={order} {compression} {images} relative={relative:.3e} limit={limit:.0e}"
    assert field_relative <= FIELD_REL, (f"{name} order={order} {images} field relative={field_relative:.3e} "
                                         f"actual={same_state['field'].tolist()} expected={field_cpu.tolist()}")


@pytest.mark.parametrize("name,order,images", [("cut", 2, dict(image="-x")), ("poles", 2, dict(image_cyclic=6))])
def test_nonlinear_solve_matches_cpu(meshes, request, prepared, cuda, name, order, images):
    _, x_cpu = _cpu_solve(meshes[name], order, bh_table=_BH, nl_tol=REFERENCE_NL_TOL, **images)
    with prepared.PreparedHDivNonlinear.from_bh_table(meshes[name], _BH, backend="gpu", order=order, **_RECEIPT, **images) as device:
        solved = device.solve((0.0, 0.0, 1e5), tolerance=1e-9, linear_tolerance=1e-10)
    relative = np.linalg.norm(solved["x"] - x_cpu)/np.linalg.norm(x_cpu)
    request.node.user_properties += [("solve_relative", float(relative)), ("newton", int(solved["newton_iterations"]))]
    assert relative <= SOLVE_REL, f"{name} order={order} {images} relative={relative:.3e} limit={SOLVE_REL:.0e}"


@pytest.mark.parametrize("masks,signs,angles,text", [
    ([0], [1.0], [], "identity"),
    ([0], [1.0], [0.0], "identity"),
    ([0], [1.0], [0.5], "inversion"),
    ([0, 0], [1.0, -1.0], [0.5, -0.5], "inversion"),
    ([1], [1.0, 1.0], [], "mismatch"),
    ([8], [1.0], [], "0..7"),
    ([0]*64, [1.0]*64, [], "63"),
])
def test_invalid_image_lists_fail_loud(meshes, cuda, masks, signs, angles, text):
    arguments, _, _ = _gram_case(meshes["gap"], 1, _terms("-x"))
    arguments.update(image_masks=np.array(masks, dtype=np.int32), image_signs=np.array(signs, dtype=np.float64),
                     image_rotations=np.array(angles, dtype=np.float64))
    with pytest.raises(ValueError, match=text):
        _entries(cuda, arguments, np.array([0]), np.array([0]))


def test_a_field_probe_on_a_source_face_refuses(meshes, prepared):
    """Guards: the order-1 closed form is singular on a source face; the device raises where the CPU route returns NaN."""
    with prepared.PreparedHDivLinear(meshes["poles"], backend="gpu", order=1, gram_compression="dense", **_RECEIPT, image_cyclic=6) as device:
        state = np.ones(device.fespace.ndof)
        with pytest.raises(RuntimeError, match="nonfinite"):
            device.field(state, np.array([[0.02, 0.03, 0.01]]))


def test_a_plus_plane_solve_and_a_labeled_seam_refuse(meshes, prepared):
    with pytest.raises(NotImplementedError, match="compressed normal trace"):
        prepared.PreparedHDivLinear(meshes["cut"], backend="gpu", order=1, gram_compression="dense", **_RECEIPT, image="+x")
    with pytest.raises(NotImplementedError, match="azimuthal"):
        prepared.PreparedHDivLinear(meshes["sector"], backend="gpu", order=1, gram_compression="dense", **_RECEIPT, image_cyclic=8)


def test_mirror_and_cyclic_together_refuse():
    from radia import hdiv_gpu
    with pytest.raises(NotImplementedError, match="cannot be combined"):
        hdiv_gpu._device_image_terms("-x", 4)
    with pytest.raises(ValueError, match="even pole count"):
        hdiv_gpu._device_image_terms(None, 5, True)
