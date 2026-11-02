"""Guards: the device material route builds E and the reference mass from tables and reports what it moves.

Assumes: no GPU, CUDA context or native artifact. Static source scans and
host-side ledger checks; device-vs-NGSolve E and M agreement is a separate GPU run.
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


def ledger(t, *extra, drop=(), law=24, persistent=4096, transient=2048):
    sizes = {"element_vertices_m": t["vertices"].nbytes, "reference_shapes": t["reference_shapes"].nbytes,
             "material_rule": t["rule"].nbytes, "element_class": t["element_class"].nbytes,
             "element_dofs": t["element_dofs"].nbytes, "element_points": t["element_points"].nbytes,
             "evaluation_invalid_flag": 4}
    rows = [dict(payload=p, direction="H2D", bytes=b) for p, b in sizes.items() if p not in drop]
    rows += [dict(payload=p, direction="D2H", bytes=b) for p, b in
             (("evaluation_flag", 4), ("mass_nonzeros", 4), ("mass_flag", 4), ("mass_spectral_bound", 8)) if p not in drop]
    rows += [dict(payload="inverse_bh_law", direction="H2D", bytes=law)]
    rows += [dict(payload=p, direction=d, bytes=b) for p, d, b in extra]
    return dict(route="device-flat-tet-evaluation", transfer_ledger=rows,
                persistent_device_bytes=persistent, setup_transient_device_bytes=transient)


def test_material_ledger_accepts_only_the_exact_table_rows():
    t = tables()
    assert hdiv_gpu._checked_material_assembly(ledger(t), t)["route"] == "device-flat-tet-evaluation"
    resized = ledger(t)
    resized["transfer_ledger"][1]["bytes"] += 8
    forged = [ledger(t, ("evaluation_csr_host_built", "H2D", 64)), ledger(t, ("material_weights_host", "H2D", 8)),
              ledger(t, ("reference_shapes", "D2H", 8)), ledger(t, ("element_dofs", "H2D", 72)),
              ledger(t, ("mass_csr_host_built", "H2D", 96)), ledger(t, ("mass_values", "D2H", 8)),
              ledger(t, drop=("element_points",)), ledger(t, drop=("mass_spectral_bound",)),
              ledger(t, drop=("mass_nonzeros",)), ledger(t, law=0), resized,
              ledger(t, persistent=0), ledger(t, transient=True), ledger(t, persistent=4096.),
              dict(ledger(t), route="host-evaluation")]
    for case in forged:
        with pytest.raises(RuntimeError):
            hdiv_gpu._checked_material_assembly(case, t)


def test_device_route_reports_evaluation_and_mass_as_device_built():
    assert hdiv_gpu.PENDING_HOST_OPERATORS == ("charge_map",)
    source = (ROOT / "src/radia/hdiv_gpu.py").read_text(encoding="utf-8")
    prepare = source[source.index("def _prepare_device_components"):source.index("def _prepare_host_components")]
    assert "_device_evaluation_tables" in prepare and "_fe_components" not in prepare, prepare
    assert "_fe_sources(" in prepare and "_fe_forms(" not in prepare, prepare
    sources = source[source.index("def _fe_sources"):source.index("def _fe_forms")]
    assert "BilinearForm" not in sources, sources
    initialize = source[source.index("    def _initialize"):source.index("    def solve")]
    assert "from_flat_tables" in initialize and "E.indptr" not in initialize and "_mass" not in initialize, initialize
    solve = source[source.index("    def solve"):source.index("    def stats")]
    assert "solve_nonlinear_material(" in solve and "indptr" not in solve and "_host_reference" not in solve, solve


def test_device_material_builder_downloads_only_scalar_flags_and_counts():
    source = (HDIV / "hdiv_flat_evaluation.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 4, downloads
    assert all("sizeof(int)" in d or "sizeof(double)" in d for d in downloads), downloads
    assert "cudaMemcpyHostToDevice" not in source and ".upload(" in source
    for name in ("mass_nonzeros", "mass_flag", "mass_spectral_bound"):
        assert '"%s","D2H"' % name in source, name


def test_device_mass_is_symmetric_by_construction_and_checked_on_device():
    source = (HDIV / "hdiv_flat_evaluation.cuh").read_text(encoding="utf-8")
    entries = source[source.index("__global__ void flatMassEntries"):source.index("__global__ void flatMassReduce")]
    # w*(E_r*E_s): the inner product commutes exactly, so M_rs and M_sr are bitwise equal.
    assert "sum.addProduct(w,values[base+r]*values[base+s])" in entries, entries
    diagnostics = source[source.index("__global__ void flatMassDiagnostics"):source.index("struct FlatEvaluationTransfer")]
    assert "1.e-13*fmax(1.,fmax(fabs(v),fabs(values[t])))" in diagnostics and "values[at]>0" in diagnostics
    build = source[source.index("void buildMass"):]
    for guard in ("not symmetric", "nonpositive diagonal", "invalid sparsity pattern", "Invalid device mass spectral bound"):
        assert guard in build, guard
    assert "Flat mass entry budget exceeded" in source


def test_newton_reuses_one_prepared_mass_operand():
    solver = (HDIV / "hmatrix_material_solver.cuh").read_text(encoding="utf-8")
    operand = solver[solver.index("inline MaterialSolution solveMaterialOperand"):solver.index("inline MaterialSolution solveLinearMaterial")]
    assert ".upload(mass" not in operand and "HostMass" not in operand, "PCG core must not upload a host mass"
    newton = (HDIV / "hmatrix_nonlinear_material.cuh").read_text(encoding="utf-8")
    loop = newton[newton.index("solveNonlinearMaterialOperand("):newton.index("// Checked host reference mass")]
    assert "solveMaterialOperand(op,referenceMass" in loop and "HostMass" not in loop and "solveLinearMaterial(" not in loop
    assert "material.massOperand()" in newton


def test_flat_material_binding_takes_tables_not_evaluation_values():
    source = (HDIV / "hdiv_resident_python.cpp").read_text(encoding="utf-8")
    block = source[source.index('def_static("from_flat_tables"'):source.index('.def("assembly",[](radia_hdiv_bh::MaterialOperator')]
    names = re.findall(r'py::arg\("([a-z_]+)"\)', block)
    assert names == ["rule", "reference_shapes", "vertices", "element_class", "element_dofs",
                     "element_points", "dofs", "law", "backend"], names


def test_newton_binding_takes_no_mass_and_host_reference_is_separate():
    source = (HDIV / "hdiv_resident_python.cpp").read_text(encoding="utf-8")
    device = source[source.index('.def("solve_nonlinear_material",'):source.index('.def("solve_nonlinear_material_host_reference"')]
    names = re.findall(r'py::arg\("([a-z_]+)"\)', device)
    assert names == ["material", "rhs", "initial", "nonlinear_tolerance", "linear_tolerance", "maximum_iterations",
                     "maximum_linear_iterations", "preconditioner", "degree", "lower_fraction"], names
    owner = (HDIV / "hdiv_resident_python_cuda.cu").read_text(encoding="utf-8")
    host = owner[owner.index("solveHostReference(radia_hmatrix_cuda::Resident& op"):]
    assert "owns its reference mass; host reference rejected" in host.split("SparseMap reference")[0]


def test_staging_checks_every_element_against_its_class_and_the_affine_map():
    source = (HDIV / "hdiv_material_staging.hpp").read_text(encoding="utf-8")
    staging = source[source.index("inline FlatEvaluationTables stageFlatEvaluationTables"):]
    for guard in ("not determined by the element vertex order", "TET map differs from the device affine convention",
                  "TET Jacobian differs from the device affine convention", "Material global sample coverage mismatch",
                  "ET_TET", "GetCurveOrder()>1"):
        assert guard in staging, guard
