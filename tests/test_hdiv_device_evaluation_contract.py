"""Guards: the device material route builds E from reference tables and reports what it uploads.

Assumes: no GPU, CUDA context or native artifact. Static source scans and
host-side ledger checks; device-vs-NGSolve E agreement is a separate GPU run.
"""
import importlib.util
from pathlib import Path
import re

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
HDIV = ROOT / "src/core/gpu/hdiv"
_spec = importlib.util.spec_from_file_location("hdiv_device_evaluation_contract", ROOT / "src/radia/hdiv_gpu.py")
hdiv_gpu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hdiv_gpu)


def tables(elements=3, rule=4, local=6, classes=2):
    return dict(rule=np.zeros((rule, 4)), reference_shapes=np.zeros((classes, rule, local, 3)),
                vertices=np.zeros((elements, 12)), element_class=np.zeros(elements, np.int32),
                element_dofs=np.zeros((elements, local), np.int32),
                element_points=np.zeros((elements, rule), np.int32))


def ledger(t, *extra, drop=(), law=24):
    sizes = {"element_vertices_m": t["vertices"].nbytes, "reference_shapes": t["reference_shapes"].nbytes,
             "material_rule": t["rule"].nbytes, "element_class": t["element_class"].nbytes,
             "element_dofs": t["element_dofs"].nbytes, "element_points": t["element_points"].nbytes,
             "evaluation_invalid_flag": 4}
    rows = [dict(payload=p, direction="H2D", bytes=b) for p, b in sizes.items() if p not in drop]
    rows += [dict(payload="evaluation_flag", direction="D2H", bytes=4),
             dict(payload="inverse_bh_law", direction="H2D", bytes=law)]
    rows += [dict(payload=p, direction=d, bytes=b) for p, d, b in extra]
    return dict(route="device-flat-tet-evaluation", transfer_ledger=rows)


def test_material_ledger_accepts_only_the_exact_table_rows():
    t = tables()
    assert hdiv_gpu._checked_material_assembly(ledger(t), t)["route"] == "device-flat-tet-evaluation"
    resized = ledger(t)
    resized["transfer_ledger"][1]["bytes"] += 8
    forged = [ledger(t, ("evaluation_csr_host_built", "H2D", 64)), ledger(t, ("material_weights_host", "H2D", 8)),
              ledger(t, ("reference_shapes", "D2H", 8)), ledger(t, ("element_dofs", "H2D", 72)),
              ledger(t, drop=("element_points",)), ledger(t, law=0), resized,
              dict(ledger(t), route="host-evaluation")]
    for case in forged:
        with pytest.raises(RuntimeError):
            hdiv_gpu._checked_material_assembly(case, t)


def test_device_route_reports_evaluation_as_device_built():
    assert hdiv_gpu.PENDING_HOST_OPERATORS == ("charge_map", "mass")
    source = (ROOT / "src/radia/hdiv_gpu.py").read_text(encoding="utf-8")
    prepare = source[source.index("def _prepare_device_components"):source.index("def _prepare_host_components")]
    assert "_device_evaluation_tables" in prepare and "_fe_components" not in prepare, prepare
    initialize = source[source.index("    def _initialize"):source.index("    def solve")]
    assert "from_flat_tables" in initialize and "E.indptr" not in initialize, initialize


def test_device_evaluation_builder_downloads_only_its_validity_flag():
    source = (HDIV / "hdiv_flat_evaluation.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 1 and "sizeof(int)" in downloads[0], downloads
    assert "cudaMemcpyHostToDevice" not in source and ".upload(" in source


def test_flat_material_binding_takes_tables_not_evaluation_values():
    source = (HDIV / "hdiv_resident_python.cpp").read_text(encoding="utf-8")
    block = source[source.index('def_static("from_flat_tables"'):source.index('.def("assembly",[](radia_hdiv_bh::MaterialOperator')]
    names = re.findall(r'py::arg\("([a-z_]+)"\)', block)
    assert names == ["rule", "reference_shapes", "vertices", "element_class", "element_dofs",
                     "element_points", "dofs", "law", "backend"], names


def test_staging_checks_every_element_against_its_class_and_the_affine_map():
    source = (HDIV / "hdiv_material_staging.hpp").read_text(encoding="utf-8")
    staging = source[source.index("inline FlatEvaluationTables stageFlatEvaluationTables"):]
    for guard in ("not determined by the element vertex order", "TET map differs from the device affine convention",
                  "TET Jacobian differs from the device affine convention", "Material global sample coverage mismatch",
                  "ET_TET", "GetCurveOrder()>1"):
        assert guard in staging, guard
