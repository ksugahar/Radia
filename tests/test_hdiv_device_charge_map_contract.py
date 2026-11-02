"""Guards: the device charge map B/B^T is built from NGSolve tables, never uploaded as a CSR.

Assumes: no GPU, CUDA context or native artifact. Static source scans and
host-side consistency checks; device-vs-CPU B agreement is a separate GPU run.
"""
import importlib.util
from pathlib import Path
import re

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
HDIV = ROOT / "src/core/gpu/hdiv"
_spec = importlib.util.spec_from_file_location("hdiv_device_charge_map_contract", ROOT / "src/radia/hdiv_gpu.py")
hdiv_gpu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hdiv_gpu)

CELL_EXPO = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]], np.int32)
FACE_EXPO = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [2, 0]], np.int32)


def staged(cells=2, faces=1):
    """BDM2-shaped staged inputs in the CPU charge order and consistent tables."""
    host = np.r_[np.repeat(np.arange(cells), 4), np.repeat(np.arange(faces), 6)]
    kind = np.r_[np.zeros(4*cells, int), np.ones(6*faces, int)]
    expo = np.r_[np.tile(CELL_EXPO.ravel(), cells), np.tile(np.c_[FACE_EXPO, np.zeros(6, int)].ravel(), faces)]
    inputs = dict(cell_verts=np.arange(12.*cells), face_verts=np.arange(9.*faces), n_el=cells,
                  charge_host=host, charge_kind=kind, charge_expo=expo,
                  ref_tet_pts=np.full(3, .25), ref_tet_w=[1/6], ref_tri_pts=np.full(2, 1/3), ref_tri_w=[.5])
    tables = dict(cell_class=np.zeros(cells, np.int32), cell_dofs=np.arange(30*cells, dtype=np.int32).reshape(cells, 30),
                  cell_tables=np.ones((2, 4, 30)), cell_exponents=CELL_EXPO,
                  face_class=np.zeros(faces, np.int32), face_dofs=np.arange(6*faces, dtype=np.int32).reshape(faces, 6),
                  face_tables=np.ones((1, 6, 6)), face_exponents=FACE_EXPO)
    return inputs, tables, (inputs["cell_verts"].copy(), inputs["face_verts"].copy())


def test_table_arguments_join_the_gram_binding_payload():
    inputs, tables, geometry = staged()
    arguments = hdiv_gpu._device_gram_arguments(tables, inputs, geometry)
    assert set(hdiv_gpu._CHARGE_TABLE_KEYS) <= set(arguments)
    expected = hdiv_gpu._expected_device_ledger(arguments)
    assert expected[("charge_cell_tables", "H2D")] == tables["cell_tables"].nbytes
    assert expected[("charge_face_dofs", "H2D")] == tables["face_dofs"].nbytes
    assert expected[("charge_map_flag", "D2H")] == 4
    assert not {p for p, _ in expected if "exponents" in p and p != "charge_exponents"}, sorted(expected)
    assert ("charge_map_csr_host_built", "H2D") not in expected


@pytest.mark.parametrize("damage,match", [
    (lambda i, t, g: g[0].__setitem__(0, -1.), "geometry"),
    (lambda i, t, g: g[1].__setitem__(8, -1.), "geometry"),
    (lambda i, t, g: t.update(cell_exponents=CELL_EXPO[::-1].copy()), "charge order"),
    (lambda i, t, g: t.update(face_class=np.zeros(2, np.int32)), "charge order"),
    (lambda i, t, g: i.update(charge_host=np.r_[i["charge_host"][:-1], 1]), "charge order")])
def test_staged_tables_must_reproduce_the_gram_geometry_and_order(damage, match):
    inputs, tables, geometry = staged()
    damage(inputs, tables, geometry)
    with pytest.raises(ValueError, match=match):
        hdiv_gpu._device_gram_arguments(tables, inputs, geometry)


def test_device_route_reports_charge_map_as_device_built():
    source = (ROOT / "src/radia/hdiv_gpu.py").read_text(encoding="utf-8")
    prepare = source[source.index("def _prepare_device_components"):source.index("def _prepare_host_components")]
    assert "_device_charge_tables" in prepare and "charge_map=" not in prepare, prepare
    initialize = source[source.index("    def _initialize"):source.index("    def solve")]
    assert "indptr" not in initialize and "B.data" not in initialize, initialize


def test_flat_binding_takes_tables_not_charge_map_values():
    source = (HDIV / "hdiv_resident_python.cpp").read_text(encoding="utf-8")
    start = source.index('def_static("from_flat_charges"')
    names = re.findall(r'py::arg\("([a-z_]+)"\)', source[start:source.index('.def("assembly"', start)])
    assert names[10:] == ["dofs", "cell_class", "cell_dofs", "cell_tables", "cell_exponents", "face_class",
                          "face_dofs", "face_tables", "face_exponents", "backend"], names
    assert not {"offsets", "indices", "values"} & set(names), names
    flat = source[source.index("std::shared_ptr<DeviceOperator> constructFlat("):start]
    assert not re.search(r"\bMapInput\b", flat[:flat.index("PYBIND11_MODULE")]), "flat binding still builds a host CSR"


def test_device_charge_builder_downloads_only_its_validity_flag():
    source = (HDIV / "hdiv_flat_charge_map.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 1 and "sizeof(int)" in downloads[0], downloads
    assert "cudaMemcpyHostToDevice" not in source and "SortPairs" in source
    operator = (HDIV / "hdiv_resident_python_cuda.cu").read_text(encoding="utf-8")
    flat = operator[operator.index("const FlatInput& flat,const ChargeMapInput& charges)"):operator.index("void planTransfer")]
    assert "buildChargeMap" in flat and "SparseMap" not in flat and "charge_map_csr_host_built" not in flat, flat


def test_device_only_charge_headers_stay_out_of_shared_native_closures():
    for header, owner in (("hdiv_flat_charge_map.cuh", "hdiv_resident_python_cuda.cu"),
                          ("hdiv_charge_staging.hpp", "hdiv_material_python.cpp")):
        users = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / "src").rglob("*")
                       if p.suffix in {".cu", ".cuh", ".cpp", ".hpp", ".h"} and f'#include "{header}"' in p.read_text(encoding="utf-8", errors="replace"))
        users += sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / "packages").rglob("*")
                        if p.suffix in {".cu", ".cuh", ".cpp", ".hpp", ".h"} and header in p.read_text(encoding="utf-8", errors="replace"))
        assert users == ["src/core/gpu/hdiv/" + owner], (header, users)


def test_staging_checks_every_cell_and_face_against_its_class():
    source = (HDIV / "hdiv_charge_staging.hpp").read_text(encoding="utf-8")
    staging = source[source.index("inline FlatChargeTables stageFlatChargeTables"):]
    for guard in ("not determined by the element vertex order", "not determined by the face vertex order",
                  "TET map differs from the device affine convention", "TET Jacobian differs from the device affine convention",
                  "normal trace is not carried by the boundary-element DOFs", "non-conforming or corrupt",
                  "ET_TET", "ET_TRIG", "GetCurveOrder()>1"):
        assert guard in staging, guard
    validate = (HDIV / "hdiv_flat_charge_map.cuh").read_text(encoding="utf-8")
    assert "Charge map rows disagree with the Gram charge order" in validate
