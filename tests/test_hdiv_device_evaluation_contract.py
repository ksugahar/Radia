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
    return dict(route="device-flat-tet-matrix-free", transfer_ledger=rows,
                persistent_device_bytes=persistent, setup_transient_device_bytes=transient)


def test_material_ledger_accepts_only_the_exact_table_rows():
    t = tables()
    assert hdiv_gpu._checked_material_assembly(ledger(t), t)["route"] == "device-flat-tet-matrix-free"
    resized = ledger(t)
    resized["transfer_ledger"][1]["bytes"] += 8
    forged = [ledger(t, ("evaluation_csr_host_built", "H2D", 64)), ledger(t, ("material_weights_host", "H2D", 8)),
              ledger(t, ("reference_shapes", "D2H", 8)), ledger(t, ("element_dofs", "H2D", 72)),
              ledger(t, ("mass_csr_host_built", "H2D", 96)), ledger(t, ("mass_values", "D2H", 8)),
              ledger(t, drop=("element_points",)), ledger(t, drop=("mass_spectral_bound",)),
              ledger(t, drop=("mass_nonzeros",)), ledger(t, law=0), resized,
              ledger(t, persistent=0), ledger(t, transient=True), ledger(t, persistent=4096.),
              dict(ledger(t), route="host-evaluation"), dict(ledger(t), route="device-flat-tet-explicit-oracle")]
    for case in forged:
        with pytest.raises(RuntimeError):
            hdiv_gpu._checked_material_assembly(case, t)


def test_device_route_reports_evaluation_and_mass_as_device_built():
    assert not hasattr(hdiv_gpu, "PENDING_HOST_OPERATORS"), "no host-built operator remains on this route"
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
    # w*(phi_i.phi_j): the same Piola values in commuted order, so M_rs and M_sr are bitwise equal.
    assert "sum.addProduct(w,piolaValue(J,det,Sk+3*i,c)*piolaValue(J,det,Sk+3*j,c))" in entries, entries
    diagnostics = source[source.index("__global__ void flatMassDiagnostics"):source.index("struct FlatEvaluationTransfer")]
    assert "1.e-13*fmax(1.,fmax(fabs(v),fabs(values[t])))" in diagnostics and "values[at]>0" in diagnostics
    build = source[source.index("void buildMass"):]
    for guard in ("not symmetric", "nonpositive diagonal", "invalid sparsity pattern", "Invalid device mass spectral bound"):
        assert guard in build, guard
    assert "Flat mass triples exceed the checked int CUB item count" in source


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
    block = source[source.index('def_static("from_flat_tables"'):source.index('.def_static("validate_flat_tables"')]
    names = re.findall(r'py::arg\("([a-z_]+)"\)', block)
    assert names == ["rule", "reference_shapes", "vertices", "element_class", "element_dofs",
                     "element_points", "dofs", "law", "backend"], names
    assert "&constructFlatMaterial," in block, "production binding must build the matrix-free context"
    # The explicit-E diagnostic oracle is retired; coverage moved to independent NGSolve references.
    assert "_explicit_evaluation_oracle" not in source and "createFlatMaterialOracle" not in source
    scope = source[source.index('.def_static("validate_flat_tables"'):source.index('.def("assembly",[](radia_hdiv_bh::MaterialOperator')]
    assert re.findall(r'py::arg\("([a-z_]+)"\)', scope) == names[:-2], "scope check takes the same tables, no law"
    helper = source[source.index("py::dict validateFlatTables("):source.index("std::shared_ptr<DeviceOperator> constructFlat(")]
    assert "validateFlatEvaluationScope(flatEvaluationInput(" in helper and "create" not in helper, helper


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


def test_reference_mass_diagnostic_moves_vectors_only():
    source = (HDIV / "hmatrix_bh_material.cuh").read_text(encoding="utf-8")
    action = source[source.index("MassAction massAction("):source.index("// Diagnostic y=E x")]
    assert action.count("cudaMemcpy(") == 2 and "nd*sizeof(double)" in action, action
    evaluation = source[source.index("// Diagnostic y=E x"):source.index("const std::string& route()")]
    assert evaluation.count("cudaMemcpy(") == 2 and "in*sizeof(double)" in evaluation and "out*sizeof(double)" in evaluation, evaluation
    assert "frozen=false" in evaluation, "a diagnostic action overwrites the point field, so it must drop the linearization"
    assert "has no device-built reference mass" in action and "M.apply(" in action


def _block(source, start, end):
    return source[source.index(start):source.index(end)]


def test_flat_material_scope_has_structural_bounds_not_fixed_size_caps():
    evaluation = (HDIV / "hdiv_flat_evaluation.cuh").read_text(encoding="utf-8")
    material = (HDIV / "hmatrix_bh_material.cuh").read_text(encoding="utf-8")
    staging = (HDIV / "hdiv_material_staging.hpp").read_text(encoding="utf-8")
    blocks = dict(scope=_block(evaluation, "inline void validateFlatEvaluation(", "// CUB scratch of FlatEvaluationBuilder"),
                  flat_context=_block(material, "MappedMaterial(const FlatEvaluationInput& in", "radia_hmatrix_cuda::MassOperand massOperand()"),
                  staging_context=_block(staging, "MaterialStagingContext(const std::string& path", "inline MaterialQuadraturePayload stageMaterialQuadrature"))
    for name, block in blocks.items():
        for cap in ("16384", "65536", "4000000", "2000000", "512)<<20", "ne>4096"):
            assert cap not in block, (name, cap)
    for name in ("scope", "staging_context"):
        assert "std::numeric_limits<int>::max()" in blocks[name], name
    assert "requireFlatEvaluationDevice(in," in blocks["flat_context"], "flat context must run the live budget check"


def test_mass_keys_are_64_bit_through_sort_reduce_and_row_search():
    source = (HDIV / "hdiv_flat_evaluation.cuh").read_text(encoding="utf-8")
    assert "using MassKey=unsigned long long;" in source
    for kernel, nxt in (("flatMassEntries", "flatMassReduce"), ("flatMassReduce", "flatMassRowOffsets"),
                        ("flatMassRowOffsets", "flatMassFind")):
        body = _block(source, "__global__ void " + kernel, nxt)
        assert "const MassKey* keys" in body or "MassKey* keys" in body, kernel
        assert "unsigned* keys" not in body and "unsigned(" not in body, kernel
    mass = source[source.index("void buildMass"):]
    assert "Buffer<MassKey> keys,sortedKeys,uniqueKeys;" in mass and "Buffer<unsigned> keys" not in mass


def test_every_cub_item_count_is_a_dimension_the_scope_check_bounds():
    source = (HDIV / "hdiv_flat_evaluation.cuh").read_text(encoding="utf-8")
    builder = _block(source, "struct FlatEvaluationBuilder", "// Host bounds before any upload")
    calls = re.findall(r"cub::(\w+::\w+)\(([^;]*)\)\);", builder)
    assert len(calls) == 10, calls
    counts = {call[1].split(",")[6 if call[0].endswith("SortPairs") else -1].strip() for call in calls}
    # slots <= triples; nnz <= triples.
    assert counts == {"int(slots)", "int(triples)", "nnz", "int(nd)"}, counts
    scope = _block(source, "inline void validateFlatEvaluation(", "// CUB scratch of FlatEvaluationBuilder")
    for guard in ("Flat mass triples exceed the checked int CUB item count",
                  "exceed int32 element, DOF or point numbering"):
        assert guard in scope, guard
    assert "explicitOracle" not in source and "flatMassEntriesExplicit" not in source, "explicit-E oracle path survives"
    scratch = _block(source, "inline std::size_t flatEvaluationSortScratch(", "inline std::size_t requireFlatEvaluationDevice(")
    queried = sorted(re.findall(r"cub::(\w+::\w+)\(nullptr", scratch))
    assert queried == sorted(c[0] for c in calls if c[1].startswith("nullptr")), "budget scratch must query every builder CUB call"


def test_matrix_free_evaluation_has_one_piola_formula_and_no_float_atomics():
    source = (HDIV / "hdiv_flat_evaluation.cuh").read_text(encoding="utf-8")
    assert source.count("(J[c][0]*S[0]+J[c][1]*S[1]+J[c][2]*S[2])/det") == 1, "one Piola formula for every E consumer"
    for kernel, nxt in (("flatEvaluationEntries", "flatEvaluationApply"), ("flatEvaluationApply", "flatEvaluationTransposeApply"),
                        ("flatEvaluationTransposeApply", "slotKeys"), ("flatMassEntries", "flatMassReduce")):
        body = _block(source, "__global__ void " + kernel + "(", "__global__ void " + nxt + "(")
        assert "piolaValue(J,det," in body and "affineTet(vertices+12*e,J)" in body, kernel
    transpose = _block(source, "__global__ void flatEvaluationTransposeApply(", "__global__ void slotKeys(")
    assert "atomic" not in transpose and "incidence[at]" in transpose, "E^T must be a deterministic gather"


def test_production_material_stores_no_global_evaluation_matrix():
    source = (HDIV / "hmatrix_bh_material.cuh").read_text(encoding="utf-8")
    flat = _block(source, "MappedMaterial(const FlatEvaluationInput& in", "radia_hmatrix_cuda::MassOperand massOperand()")
    assert "builder.build(in,flat,weights,M,massDiagonal);" in flat and "&E" not in flat and "&ET" not in flat
    assert 'matrixFree=true;routeName="device-flat-tet-matrix-free";' in flat
    for method in ("MaterialEvaluation evaluate(", "void applyFrozenDevice("):
        body = source[source.index(method):]
        body = body[:body.index("\n    }\n")]
        assert "E.apply(" not in body and "ET.apply(" not in body and "applyE(" in body and "applyET(" in body, method


def test_python_budget_terms_match_the_native_publication():
    header = (HDIV / "hdiv_inverse_bh.hpp").read_text(encoding="utf-8")
    names = re.search(r"names\[count\]=\{([^}]*)\}", header).group(1)
    native = tuple(re.findall(r'"([a-z_]+)"', names))
    assert native == hdiv_gpu._MATERIAL_BUDGET_TERMS, (native, hdiv_gpu._MATERIAL_BUDGET_TERMS)
    assert int(re.search(r"static constexpr std::size_t count=(\d+);", header).group(1)) == len(native)


def test_explicit_evaluation_oracle_is_retired_everywhere():
    """Guards: the explicit-E diagnostic oracle does not return as a second global E path.
    Assumes: independent NGSolve mass/gradient references and the E/E^T adjoint replace it."""
    retired = ("_explicit_evaluation_oracle", "createFlatMaterialOracle", "explicitOracle",
               "flatMassEntriesExplicit", "device-flat-tet-explicit-oracle")
    sources = [*HDIV.glob("*.cuh"), *HDIV.glob("*.hpp"), *HDIV.glob("*.cu"), *HDIV.glob("*.cpp"),
               ROOT / "src/radia/hdiv_gpu.py"]
    found = [(path.name, name) for path in sources for name in retired
             if name in path.read_text(encoding="utf-8")]
    assert len(sources) > 5 and not found, found
