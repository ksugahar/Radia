"""Guards: directed persisted-ID verification, namespace and immutable load boundary.
Assumes: SYNTHETIC protocol objects below; actual Netgen loading is a separate
native qualification, never simulated by production code.
"""
from copy import deepcopy
from types import SimpleNamespace, ModuleType
import sys

import pytest

# Source-only protocol contract: load the actual module without importing the
# package's compiled solver initialization. Native load is separately qualified.
import importlib.util
from pathlib import Path
_spec = importlib.util.spec_from_file_location(
    "kelvin_contract_under_test",
    Path(__file__).resolve().parents[1] / "src/radia/kelvin_identify_ngsolve.py")
_contract = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_contract)
verify_kelvin_identification = _contract.verify_kelvin_identification
load_verified_kelvin_mesh = _contract.load_verified_kelvin_mesh
_native_rows = _contract._native_rows
_verify_kelvin_structure = _contract._verify_kelvin_structure


def _pid(n):
    return SimpleNamespace(nr=n)


@pytest.fixture
def synthetic(monkeypatch):
    meshing = ModuleType("netgen.meshing")
    meshing.IdentificationType = SimpleNamespace(PERIODIC=2)
    monkeypatch.setitem(sys.modules, "netgen.meshing", meshing)
    ng = SimpleNamespace()
    ng.descriptors = [SimpleNamespace(bc=2, bcname="master"),
                      SimpleNamespace(bc=1, bcname="slave")]
    ng.surfaces = [SimpleNamespace(index=1, vertices=list(map(_pid, [4, 5, 6]))),
                   SimpleNamespace(index=2, vertices=list(map(_pid, [1, 2, 3])))]
    ng.rows = [(_pid(4), _pid(1), 1), (_pid(5), _pid(2), 1), (_pid(6), _pid(3), 1)]
    ng.GetIdentifications = lambda: ng.rows
    ng.GetNrIdentifications = lambda: 1
    ng.GetIdentificationType = lambda nr: 2
    ng.FaceDescriptors = lambda: ng.descriptors
    ng.Elements2D = lambda: ng.surfaces
    mesh = SimpleNamespace(ngmesh=ng, GetBoundaries=lambda: ("slave", "master"))
    return mesh


def _verify(mesh):
    return verify_kelvin_identification(mesh, master="master", slave="slave", group=1)


def test_fd_and_bc_namespaces_and_pointids(synthetic):
    before = deepcopy([(a.nr, b.nr, nr) for a, b, nr in synthetic.ngmesh.rows])
    info = _verify(synthetic)
    assert info["n_pairs"] == 3
    assert info["fd_bc_table"][0] == {"fd": 1, "bc": 2, "label": "master"}
    assert info["name_status"] == "NOT_VERIFIABLE_FROM_THIS_API"
    assert [(a.nr, b.nr, nr) for a, b, nr in synthetic.ngmesh.rows] == before


@pytest.mark.parametrize("case", ["partial", "duplicate", "reversed", "foreign", "self", "wrong_group", "type", "alias", "fd_label"])
def test_invalid_structure_rejected(synthetic, case):
    ng = synthetic.ngmesh
    if case == "partial": ng.rows.pop()
    elif case == "duplicate": ng.rows.append(ng.rows[0])
    elif case == "reversed": ng.rows = [(b, a, nr) for a, b, nr in ng.rows]
    elif case == "foreign": ng.rows[0] = (_pid(99), _pid(1), 1)
    elif case == "self": ng.rows[0] = (_pid(4), _pid(4), 1)
    elif case == "wrong_group": ng.rows = [(a, b, 2) for a, b, nr in ng.rows]
    elif case == "type": ng.GetIdentificationType = lambda nr: 1
    elif case == "alias":
        ng.descriptors.append(SimpleNamespace(bc=3, bcname="MASTER"))
        synthetic.GetBoundaries = lambda: ("slave", "master", "MASTER")
    elif case == "fd_label": ng.descriptors[0].bcname = "unknown"
    with pytest.raises(ValueError): _verify(synthetic)


def test_no_old_pair_shape_or_integer_conversion_fallback(synthetic):
    synthetic.ngmesh.rows = [(_pid(4), _pid(1))]
    with pytest.raises(ValueError, match="triple"): _native_rows(synthetic.ngmesh)


RAW = b"SYNTHETIC protocol fixture, not a native mesh\nidentifications\n3\n4 1 1\n5 2 1\n6 3 1\nidentificationtypes\n1\n2\nidentificationnames\n1\nkt\n"


def test_named_loader_binds_actual_file_and_native_rows(tmp_path, monkeypatch, synthetic):
    path = tmp_path / "synthetic.vol"
    path.write_bytes(RAW)
    seen = []
    module = ModuleType("ngsolve")
    module.Mesh = lambda actual: (seen.append(actual) or synthetic)
    monkeypatch.setitem(sys.modules, "ngsolve", module)
    mesh, evidence = load_verified_kelvin_mesh(path, "master", "slave", 1, "kt")
    assert mesh is synthetic
    assert seen == [str(path.resolve())]
    assert evidence["name_status"] == "VERIFIED_PERSISTED_FILE_AND_NATIVE_ROWS"
    assert path.read_bytes() == RAW


@pytest.mark.parametrize("change", ["name", "type", "duplicate_section", "native_rows", "input_changed"])
def test_named_load_boundary_rejects_drift(tmp_path, monkeypatch, synthetic, change):
    path = tmp_path / "synthetic.vol"
    raw = RAW
    if change == "name": raw = raw.replace(b"kt\n", b"other\n")
    elif change == "type": raw = raw.replace(b"types\n1\n2", b"types\n1\n1")
    elif change == "duplicate_section": raw += b"identificationnames\n1\nkt\n"
    path.write_bytes(raw)
    module = ModuleType("ngsolve")
    def load(actual):
        if change == "native_rows": synthetic.ngmesh.rows.pop()
        if change == "input_changed": path.write_bytes(raw+b"changed")
        return synthetic
    module.Mesh = load
    monkeypatch.setitem(sys.modules, "ngsolve", module)
    with pytest.raises(ValueError):
        load_verified_kelvin_mesh(path, "master", "slave", 1, "kt")


def test_in_memory_direction_observed_not_named_proof(synthetic):
    forward = _verify_kelvin_structure(synthetic, boundaries=("master", "slave"))
    assert forward["master"] == "master"
    synthetic.ngmesh.rows = [(b, a, nr) for a, b, nr in synthetic.ngmesh.rows]
    reverse = _verify_kelvin_structure(synthetic, boundaries=("master", "slave"))
    assert reverse["master"] == "slave"
    assert reverse["direction_status"] == "OBSERVED_NATIVE_NOT_PRODUCER_ORIENTATION_PROOF"
    assert reverse["name_status"] == "NOT_VERIFIABLE_FROM_THIS_API"
    with pytest.raises(ValueError): _verify(synthetic)


def test_exact_label_union_across_split_bc_patches(synthetic):
    synthetic.ngmesh.descriptors.append(SimpleNamespace(bc=3, bcname="master"))
    synthetic.GetBoundaries = lambda: ("slave", "master", "master")
    synthetic.ngmesh.surfaces.append(SimpleNamespace(index=3, vertices=list(map(_pid, [4, 5, 6]))))
    info = _verify(synthetic)
    assert info["master_points"] == 3
    assert info["n_pairs"] == 3


def test_boundary_names_are_not_assumed_indexed_by_bc_or_fd(synthetic):
    synthetic.GetBoundaries = lambda: ("master", "slave")
    info = _verify(synthetic)
    assert info["fd_bc_table"][0] == {"fd": 1, "bc": 2, "label": "master"}
