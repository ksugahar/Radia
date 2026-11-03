"""Pure-HEX device route against the CPU hex charge Gram.

Guards: on an affine 2x2x2 cube and a bilinearly sheared 2x2x1 lattice
(every cell a non-parallelepiped trilinear hex), both left-handed as NGSolve
structures them (det J < 0), and the sheared lattice mirrored in z
(det J > 0), BDM1 and BDM2, the staged
device input reproduces the CPU hex charge order and Q2 lattices, and the
device-built hex Gram entries (Duffy, product, far, affine analytic and
composite near pairs), the energy operator B^T G B (device B from the
NGSolve moment tables), the reference mass and the linear material solve
agree with the CPU route ``_build_charge_gram_hex`` (its default rules) and
the NGSolve mass on the material rule. Dropping the symmetrized pair block,
the pair-domain Duffy rule, the distorted composite near rule or the face
normal orientation, or the det J sign of the cell charge rows, moves these
far beyond the detectors. Thin (aspect above 12) affine and trilinear lattices
and a non-conforming box contact take the CPU graded near family at BDM1, whose
device rule is the graded tensor outer with the analytic or cone/fan inner.
Assumes: the CPU hex route is the reference on identical NGSolve inputs;
this is parity of one discretization, not independent accuracy. The REL
values are PROPOSED detectors, not adopted criteria. Optional GPU lane: skips
visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set; a configured lane loads the
binding through the provenance-checked loader and fails if it cannot.
"""

import os
from pathlib import Path

import numpy as np
import pytest

ENTRY_REL = 1e-10
OPERATOR_REL = 1e-9
MASS_REL = 1e-12
LINEAR_REL = 1e-8
INV_CHI = 0.1
PAIR_SAMPLE = 20000
_GEOMETRY = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents",
             "tet_rule", "tri_rule", "tet_low_rule", "tri_low_rule", "far_factor", "family",
             "hex_rules", "hex_bands")


def _shear(x, y, z):
    return (x + 0.2*x*y, y + 0.2*y*z, z + 0.2*z*x)


def _mirrored_shear(x, y, z):
    return _shear(x, y, -z)


_MESHES = {"cube": dict(nx=2, ny=2, nz=2), "sheared": dict(nx=2, ny=2, nz=1, mapping=_shear),
           "mirrored": dict(nx=2, ny=2, nz=1, mapping=_mirrored_shear)}
# Structured numbering keeps one vertex order per reference face; these share an order across reference faces.
_NAMES = sorted(_MESHES) + ["pair", "renumbered"]
# BDM1 pairs of these meshes take the graded near family (CPU stats counters below).
_GRADED = ("contact", "thin", "thin-sheared")
_CASES = [(name, order) for name in _NAMES for order in (1, 2)] + [(name, 1) for name in _GRADED]


def _thin(x, y, z):
    return (x, y, 0.02*z)


def _thin_shear(x, y, z):
    return _thin(*_shear(x, y, z))


def _box_mesh(boxes):
    """Boxes as separate HEX cells: a contact between them shares lattice nodes but no whole entity."""
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint, Pnt
    mesh = Mesh(dim=3)
    face = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    for lo, hi in boxes:
        corners = [mesh.Add(MeshPoint(Pnt(*[(hi if bit else lo)[axis] for axis, bit in enumerate(bits)])))
                   for bits in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))]
        mesh.Add(Element3D(1, corners))
        for quad in ((3, 2, 1, 0), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)):
            mesh.Add(Element2D(face, [corners[k] for k in quad]))
    return mesh


def _renumbered_mesh(options):
    """The structured mesh with a fixed random point numbering, as a general mesher would number it."""
    import numpy as np
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint, Pnt
    from ngsolve.meshes import MakeStructured3DMesh
    source = MakeStructured3DMesh(hexes=True, **options).ngmesh
    points = [point.p for point in source.Points()]
    target = Mesh(dim=3)
    face = target.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    new = np.random.default_rng(7).permutation(len(points))
    handles = [None]*len(points)
    for index in np.argsort(new):
        handles[index] = target.Add(MeshPoint(Pnt(*points[index])))
    for element in source.Elements3D():
        target.Add(Element3D(1, [handles[vertex.nr - 1] for vertex in element.vertices]))
    for element in source.Elements2D():
        target.Add(Element2D(face, [handles[vertex.nr - 1] for vertex in element.vertices]))
    return target


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


@pytest.fixture(scope="module")
def meshes(tmp_path_factory):
    from ngsolve.meshes import MakeStructured3DMesh
    root = tmp_path_factory.mktemp("hex")
    paths = {}
    for name, options in _MESHES.items():
        paths[name] = root / (name + ".vol")
        MakeStructured3DMesh(hexes=True, **options).ngmesh.Save(str(paths[name]))
    paths["pair"] = Path(__file__).resolve().parent / "fixtures" / "hdiv" / "hex_pair.vol"
    paths["renumbered"] = root / "renumbered.vol"
    _renumbered_mesh(_MESHES["sheared"]).Save(str(paths["renumbered"]))
    for name, mapping in (("thin", _thin), ("thin-sheared", _thin_shear)):
        paths[name] = root / (name + ".vol")
        MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=1, mapping=mapping).ngmesh.Save(str(paths[name]))
    paths["contact"] = root / "contact.vol"
    _box_mesh((((0., 0., 0.), (2., 2., 2.)), ((2., 0., 0.), (3., 1., 1.)))).Save(str(paths["contact"]))
    return paths


def _cpu(path, order):
    import ngsolve as ng
    from radia.vim._vim import _build_charge_gram_hex
    with ng.TaskManager():
        mesh = ng.Mesh(str(path))
        fes = ng.HDiv(mesh, order=order)
        B, gram, _, _ = _build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False, materialize_mass=False)
        # The device mass is E^T W E on the material rule IntegrationRule(HEX, 2*max(6, 2p+2)).
        u, v = fes.TnT()
        mass = ng.BilinearForm(fes)
        mass += u*v*ng.dx(intrules={ng.HEX: ng.IntegrationRule(ng.HEX, 2*max(6, 2*order + 2))})
        mass.Assemble()
    gram.build_exact_dense_normalized_gram(maximum_memory_mb=256)
    rows, cols, values = mass.mat.COO()
    from scipy.sparse import coo_matrix, csr_matrix
    return csr_matrix(B), gram, coo_matrix((values, (rows, cols)), shape=(fes.ndof, fes.ndof)).tocsr()


def _case(path, order, cuda, compression="dense"):
    import ngsolve as ng
    from radia import hdiv_gpu
    with ng.TaskManager():
        setup = hdiv_gpu._prepare_device_components(path, order, cuda, aca=None, compression=compression)
    B, gram, M = _cpu(path, order)
    assert setup["family"] == "hex" and setup["dofs"] == B.shape[1] == M.shape[0]
    assert len(setup["charge_gram"]["charge_host"]) == B.shape[0]
    return setup, B, gram, M


@pytest.mark.parametrize("name,order", _CASES)
def test_hex_device_input_reproduces_cpu_charges_and_lattices(meshes, name, order):
    import ngsolve as ng
    from radia.vim._vim import _charge_basis_hex, _prepare_gpu_hex_charge_inputs
    with ng.TaskManager():
        fes = ng.HDiv(ng.Mesh(str(meshes[name])), order=order)
        staged = _prepare_gpu_hex_charge_inputs(fes, backend="gpu")
        cpu = _charge_basis_hex(fes, cob_quad=max(3, order + 1), materialize_mass=False)
    for key, cpu_key in (("charge_host", "host"), ("charge_kind", "kind"), ("charge_expo", "expo"),
                         ("cell_nodes", "cell_nodes"), ("face_nodes", "face_nodes")):
        assert np.array_equal(np.asarray(staged[key]).ravel(), np.asarray(cpu[cpu_key]).ravel()), key


@pytest.mark.parametrize("name,order", _CASES)
def test_hex_device_entries_match_cpu_charge_gram(meshes, request, cuda, name, order):
    setup, _, gram, _ = _case(meshes[name], order, cuda)
    if name in _GRADED:
        stats = gram.stats()
        routed = stats["hex_pair_anisotropy_graded"] + stats["hex_pair_nonconforming"]
        assert routed > 0, f"{name} took no graded near pair: {stats}"
    device = {key: setup["charge_gram"][key] for key in _GEOMETRY}
    n = len(device["charge_host"])
    rng = np.random.default_rng(20261103)
    pairs = np.unique(np.r_[np.c_[np.arange(n), np.arange(n)], rng.integers(0, n, size=(PAIR_SAMPLE, 2))], axis=0)
    rows, cols = (np.ascontiguousarray(pairs[:, k], dtype=np.uint64) for k in (0, 1))
    actual = np.concatenate([cuda._diagnostic_flat_charge_entries(**device, rows=rows[s:s + 65536],
                                                                  cols=cols[s:s + 65536], backend="gpu")
                             for s in range(0, len(rows), 65536)])
    expected = np.array([gram.entry(int(i), int(j)) for i, j in pairs])
    diagonal = np.array([gram.entry(i, i) for i in range(n)])
    relative = np.abs(actual - expected)/np.sqrt(np.abs(diagonal[pairs[:, 0]]*diagonal[pairs[:, 1]]))
    worst = int(np.argmax(relative))
    request.node.user_properties += [("charges", n), ("pairs", len(pairs)), ("max_relative", float(relative.max()))]
    assert relative.max() <= ENTRY_REL, (
        f"{name} order={order} pair={pairs[worst].tolist()} actual={actual[worst]!r} "
        f"expected={expected[worst]!r} relative={relative[worst]:.3e} limit={ENTRY_REL:.0e}")


@pytest.mark.parametrize("name,order", _CASES)
def test_hex_device_operator_and_mass_match_cpu(meshes, request, cuda, name, order):
    setup, B, gram, M = _case(meshes[name], order, cuda)
    nd = setup["dofs"]
    assert setup["evaluation_tables"]["vertices"].shape[1] == 81
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    try:
        from radia import hdiv_gpu
        assert hdiv_gpu._checked_device_assembly(operator.assembly(), setup["charge_gram"])["route"] == "device-hex"
        assert hdiv_gpu._checked_material_assembly(mass.assembly(), setup["evaluation_tables"],
                                                   mass=True)["route"] == "device-hex-mass"
        rng = np.random.default_rng(order)
        worst = dict(operator=0., mass=0.)
        for _ in range(4):
            x = rng.standard_normal(nd)
            expected = B.T @ np.asarray(gram.matvec_sym(np.ascontiguousarray(B @ x)))
            actual = operator.apply(np.ascontiguousarray(x))
            worst["operator"] = max(worst["operator"], np.linalg.norm(actual - expected)/np.linalg.norm(expected))
            y = mass.reference_mass_apply(np.ascontiguousarray(x))["y"]
            worst["mass"] = max(worst["mass"], np.linalg.norm(y - M @ x)/np.linalg.norm(M @ x))
        request.node.user_properties += [("dofs", nd)] + [("max_relative_" + k, float(v)) for k, v in worst.items()]
        assert worst["operator"] <= OPERATOR_REL and worst["mass"] <= MASS_REL, (
            f"{name} order={order} {worst} limits operator={OPERATOR_REL:.0e} mass={MASS_REL:.0e}")
    finally:
        mass.close()
        operator.close()


@pytest.mark.parametrize("name,order,compression", [("cube", 1, "aca"), ("sheared", 1, "dense"),
                                                    ("sheared", 2, "aca"), ("mirrored", 2, "dense"),
                                                    ("thin", 1, "aca"), ("thin-sheared", 1, "dense")])
def test_hex_device_linear_solve_matches_cpu_equation(meshes, request, cuda, name, order, compression):
    setup, B, gram, M = _case(meshes[name], order, cuda, compression)
    nd = setup["dofs"]
    columns = [B.T @ np.asarray(gram.matvec_sym(np.ascontiguousarray(B[:, k].toarray().ravel()))) for k in range(nd)]
    rhs = setup["rhs_operator"] @ np.array([0., 0., 1e5])
    expected = np.linalg.solve(np.column_stack(columns) + INV_CHI*M.toarray(), rhs)
    operator = cuda._ResidentHDiv.from_flat_charges(**setup["charge_gram"], dofs=nd, backend="gpu")
    mass = cuda._ResidentMass.from_flat_tables(**setup["evaluation_tables"], dofs=nd, backend="gpu")
    try:
        result = operator.solve_linear(mass, np.ascontiguousarray(rhs), INV_CHI, 1e-11, 4000, "mass-diagonal")
        relative = np.linalg.norm(result["x"] - expected)/np.linalg.norm(expected)
        request.node.user_properties += [("dofs", nd), ("relative", float(relative)), ("iterations", result["iterations"])]
        assert result["relative_residual"] <= 1e-11 and operator.assembly()["compression"] == compression
        assert relative <= LINEAR_REL, (f"{name} order={order} {compression} relative={relative:.3e} "
                                        f"limit={LINEAR_REL:.0e} iterations={result['iterations']}")
    finally:
        mass.close()
        operator.close()
