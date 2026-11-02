"""Guards: the device charge-Gram route takes primitives only and reports what it uploads.

Assumes: no GPU, CUDA context or native artifact. Static source scans and
host-side mapping checks; CPU/GPU Gram agreement is a separate validation run.
"""
import importlib.util
from pathlib import Path
import re
import sys
import types

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
HDIV = ROOT / "src/core/gpu/hdiv"
_spec = importlib.util.spec_from_file_location("hdiv_device_gram_contract", ROOT / "src/radia/hdiv_gpu.py")
hdiv_gpu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hdiv_gpu)


def staged(far=True, n_el=1):
    inputs = dict(cell_verts=np.arange(12.), face_verts=np.arange(9.), n_el=n_el,
                  charge_host=[0, 0], charge_kind=[0, 1], charge_expo=[0, 0, 0, 1, 0, 0],
                  ref_tet_pts=np.full(6, .25), ref_tet_w=[.1, .2],
                  ref_tri_pts=np.full(2, 1/3), ref_tri_w=[.5],
                  eps=1e-10, leaf=64, eta=2.)
    if far:
        inputs.update(ref_tet_pts_lo=np.full(3, .25), ref_tet_w_lo=[1/6],
                      ref_tri_pts_lo=np.full(2, 1/3), ref_tri_w_lo=[.5], ho_far_factor=2.)
    return inputs


def tables(inputs):
    """Staged charge tables consistent with ``inputs`` (one cell, one face monomial)."""
    expo = np.asarray(inputs["charge_expo"], dtype=np.int32)
    return (dict(cell_class=np.zeros(1, np.int32), cell_dofs=np.arange(6, dtype=np.int32)[None],
                 cell_tables=np.ones((1, 1, 6)), cell_exponents=expo[None, :3],
                 face_class=np.zeros(1, np.int32), face_dofs=np.arange(3, dtype=np.int32)[None],
                 face_tables=np.ones((1, 1, 3)), face_exponents=expo[None, 3:5]),
            (np.asarray(inputs["cell_verts"], dtype=np.float64), np.asarray(inputs["face_verts"], dtype=np.float64)))


def gram(inputs, aca=None):
    charge, geometry = tables(inputs)
    return hdiv_gpu._device_gram_arguments(charge, inputs, geometry, aca)


def test_staged_inputs_map_to_primitive_rows_with_explicit_aca_parameters():
    arguments = gram(staged())
    preset = hdiv_gpu.GRAM_ACA_PRESET
    assert arguments["compression"] == "aca" and arguments["aca"] == preset, arguments["aca"]
    assert arguments["aca"] is not preset, "callers must not share the module constant"
    tuned = dict(preset, leaf_size=8, eta=1.5)
    assert gram(staged(), tuned)["aca"] == tuned, "caller parameters must reach the binding unchanged"
    assert arguments["tet_rule"].shape == (2, 4) and arguments["tri_rule"].shape == (1, 4)
    assert arguments["tri_rule"][0, 2] == 0 and arguments["tet_low_rule"][0, 3] == pytest.approx(1/6)
    assert arguments["far_factor"] == 2.
    assert {k: str(v.dtype) for k, v in arguments.items() if k.startswith("charge_")} == {
        "charge_host": "int32", "charge_kind": "int32", "charge_exponents": "int32"}
    assert not {"eps", "leaf", "eta", "leaves", "coefficients"} & set(arguments), sorted(arguments)


def test_missing_low_rules_disable_the_far_split():
    arguments = gram(staged(far=False))
    assert arguments["far_factor"] == 1e30 and arguments["tet_low_rule"].shape == (0, 4)


@pytest.mark.parametrize("change,error", [
    (dict(n_el=2), ValueError), (dict(image_masks=[1]), NotImplementedError)])
def test_inconsistent_or_image_payloads_fail_before_the_device(change, error):
    inputs = staged()
    inputs.update(change)
    with pytest.raises(error):
        gram(inputs)


def charge_budget(compression, **change):
    terms = dict(charge_prepare=900, charge_map=300, charge_map_scratch=200, resident_vectors=100,
                 gram_dense=0, aca_batch_scratch=0)
    terms["gram_dense" if compression == "dense" else "aca_batch_scratch"] = 400
    terms.update(change)
    return dict(schema="radia.hdiv-flat-charge-budget.v1", terms=terms, total_bytes=sum(terms.values()),
                limiting_term=max(terms, key=terms.get), excludes="ACA factor storage")


def ledger(arguments, *extra, drop=()):
    rows = [dict(payload=p, direction=d, bytes=8 if b is None else b)
            for (p, d), b in hdiv_gpu._expected_device_ledger(arguments).items() if p not in drop]
    return dict(route="device-flat-tet", charges=int(arguments["charge_host"].size),
                compression=arguments["compression"], hmatrix={k: 0 for k in hdiv_gpu._ACA_COUNTERS},
                charge_budget=charge_budget(arguments["compression"]),
                transfer_ledger=rows + [dict(payload=p, direction=d, bytes=b) for p, d, b in extra])


def test_charge_budget_requires_the_named_terms_of_the_selected_storage():
    arguments = gram(staged())
    accepted = hdiv_gpu._checked_device_assembly(ledger(arguments), arguments)["charge_budget"]
    assert tuple(accepted["terms"]) == hdiv_gpu._CHARGE_BUDGET_TERMS
    assert accepted["limiting_term"] == "charge_prepare" and accepted["total_bytes"] == 1900
    assert hdiv_gpu._checked_charge_budget(charge_budget("dense"), "dense")["terms"]["gram_dense"] == 400
    wrong_total = charge_budget("aca")
    wrong_total["total_bytes"] += 1
    wrong_limit = dict(charge_budget("aca"), limiting_term="charge_map")
    renamed = charge_budget("aca")
    renamed["terms"] = {("charge_scratch" if k == "charge_map_scratch" else k): v for k, v in renamed["terms"].items()}
    for bad in (None, dict(charge_budget("aca"), schema="v0"), wrong_total, wrong_limit, renamed,
                charge_budget("aca", gram_dense=7), charge_budget("aca", aca_batch_scratch=0),
                charge_budget("aca", charge_map=0), charge_budget("aca", resident_vectors=True),
                charge_budget("aca", charge_prepare=9.0), charge_budget("dense")):
        with pytest.raises(RuntimeError, match="charge budget"):
            hdiv_gpu._checked_charge_budget(bad, "aca")


def test_assembly_ledger_requires_the_exact_primitive_rows():
    arguments = gram(staged())
    assert hdiv_gpu._checked_device_assembly(ledger(arguments), arguments)["charges"] == 2
    flags_only = dict(route="device-flat-tet", charges=2, compression="aca",
                      hmatrix={k: 0 for k in hdiv_gpu._ACA_COUNTERS}, transfer_ledger=[
        dict(payload="prepare_flag", direction="D2H", bytes=4),
        dict(payload="energy_flag", direction="D2H", bytes=4)])
    duplicate = ledger(arguments, ("tet_rule", "H2D", 64))
    resized = ledger(arguments)
    resized["transfer_ledger"][0]["bytes"] += 8
    forged = [flags_only, duplicate, resized, ledger(arguments, drop=("resident_plan",)),
              ledger(arguments, ("host_hmatrix_coefficients", "H2D", 8)),
              ledger(arguments, ("cell_vertices_m", "D2H", 96)), dict(ledger(arguments), route="host-leaves"),
              ledger(arguments, ("charge_map_csr_host_built", "H2D", 64)), ledger(arguments, drop=("charge_cell_tables",)),
              dict(ledger(arguments), compression="dense"), dict(ledger(arguments), hmatrix={}),
              dict(ledger(arguments), hmatrix=dict({k: 0 for k in hdiv_gpu._ACA_COUNTERS}, rank_sum=-1)),
              ledger(arguments, drop=("aca_permutation",)), ledger(arguments, ("aca_lowrank_factors", "D2H", 64)),
              ledger(arguments, ("aca_leaf_status", "D2H", 0)),
              ledger(arguments, ("aca_leaf_status", "D2H", 8), ("aca_leaf_status", "D2H", 8))]
    for bad in ("8", -8, True, 0):
        row = ledger(arguments)
        next(r for r in row["transfer_ledger"] if r["payload"] == "resident_plan")["bytes"] = bad
        forged.append(row)
    for case in forged:
        with pytest.raises(RuntimeError):
            hdiv_gpu._checked_device_assembly(case, arguments)


@pytest.mark.parametrize("change", [
    dict(eps=0.), dict(eps=1.), dict(eps=float("nan")), dict(eta=0.), dict(eta=float("inf")),
    dict(leaf_size=0), dict(leaf_size=4097), dict(leaf_size=8.), dict(leaf_size=True),
    dict(max_rank=0), dict(max_rank=1025), dict(rank_cap="drop"), dict(eps="1e-10"), dict(extra=1)])
def test_gram_aca_parameters_are_explicit_and_range_checked(change):
    with pytest.raises(ValueError, match="gram_aca"):
        hdiv_gpu._checked_gram_aca(dict(hdiv_gpu.GRAM_ACA_PRESET, **change))
    missing = dict(hdiv_gpu.GRAM_ACA_PRESET)
    del missing["leaf_size"]
    with pytest.raises(ValueError, match="gram_aca"):
        hdiv_gpu._checked_gram_aca(missing)


def test_preset_is_copied_and_explicit_values_are_owned():
    first, second = hdiv_gpu._checked_gram_aca(None), hdiv_gpu._checked_gram_aca(None)
    assert first == hdiv_gpu.GRAM_ACA_PRESET and first is not second is not hdiv_gpu.GRAM_ACA_PRESET
    caller = dict(hdiv_gpu.GRAM_ACA_PRESET, leaf_size=np.int64(32))
    checked = hdiv_gpu._checked_gram_aca(caller)
    caller["leaf_size"] = 1
    assert checked["leaf_size"] == 32 and type(checked["leaf_size"]) is int, checked


def _fake_device_staging(monkeypatch, calls, validate):
    """Stand-ins for NGSolve, CPU charge staging and the native module (no FE work)."""
    fes = types.SimpleNamespace(ndof=10)
    ng = types.ModuleType("ngsolve")
    ng.Mesh = lambda path: calls.append("mesh") or "mesh"
    ng.HDiv = lambda mesh, order: fes
    vim = types.ModuleType("radia.vim._vim")
    vim._prepare_gpu_charge_inputs = lambda *a, **k: calls.append("charge-staging")
    vim._prepare_gpu_curved_charge_inputs = lambda *a, **k: calls.append("curved-charge-staging")
    vim._curve_mesh = lambda mesh, order: calls.append(("curve", order))
    vim._volume_vertex_counts = lambda mesh: frozenset({4})
    vim._prepare_gpu_hex_charge_inputs = lambda *a, **k: calls.append("hex-charge-staging")
    for name, module in (("ngsolve", ng), ("radia", types.ModuleType("radia")),
                         ("radia.vim", types.ModuleType("radia.vim")), ("radia.vim._vim", vim)):
        monkeypatch.setitem(sys.modules, name, module)
    monkeypatch.setattr(hdiv_gpu, "_device_evaluation_tables",
                        lambda *a: calls.append("evaluation-staging") or dict(rule="tables"))
    monkeypatch.setattr(hdiv_gpu, "_device_charge_tables", lambda *a: calls.append("charge-tables"))
    return types.SimpleNamespace(_ResidentBHMaterial=types.SimpleNamespace(validate_flat_tables=validate))


def test_material_scope_check_precedes_charge_staging_and_gram(monkeypatch):
    calls = []

    def validate(**arguments):
        calls.append(("material-scope", arguments))
        raise ValueError("Flat evaluation entry budget exceeded")
    native = _fake_device_staging(monkeypatch, calls, validate)
    with pytest.raises(ValueError, match="entry budget"):
        hdiv_gpu._prepare_device_components("mesh.vol", 2, native, aca=None)
    assert calls == ["mesh", "evaluation-staging", ("material-scope", dict(rule="tables", dofs=10))], calls


def test_invalid_aca_parameters_fail_before_any_mesh_work(monkeypatch):
    calls = []
    native = _fake_device_staging(monkeypatch, calls, lambda **a: calls.append("material-scope"))
    with pytest.raises(ValueError, match="gram_aca leaf_size"):
        hdiv_gpu._prepare_device_components("mesh.vol", 2, native,
                                            aca=dict(hdiv_gpu.GRAM_ACA_PRESET, leaf_size=0))
    assert calls == [], calls


def _budget(**change):
    terms = dict(reference_tables=10, incidence=400, incidence_sort=50,
                 mass_triples=300, mass_csr=40, material_state=110)
    budget = dict(schema="radia.hdiv-flat-material-budget.v1", terms=terms, total_bytes=sum(terms.values()),
                  limiting_term="incidence", excludes="CUB scratch and the inverse-BH law")
    budget.update(change)
    return budget


def test_scope_check_budget_is_returned_as_an_owned_checked_copy(monkeypatch):
    calls, published = [], _budget()
    native = _fake_device_staging(monkeypatch, calls, lambda **a: calls.append("material-scope") or published)
    sys.modules["radia.vim._vim"]._prepare_gpu_charge_inputs = lambda *a, **k: calls.append("charge-staging") or (None, None)
    monkeypatch.setattr(hdiv_gpu, "_device_charge_tables", lambda *a: calls.append("charge-tables") or (None, None))
    monkeypatch.setattr(hdiv_gpu, "_device_gram_arguments", lambda *a: "gram")
    monkeypatch.setattr(hdiv_gpu, "_fe_sources", lambda mesh, fes: dict(dofs=fes.ndof))
    setup = hdiv_gpu._prepare_device_components("mesh.vol", 2, native, aca=None)
    assert calls == ["mesh", "evaluation-staging", "material-scope", "charge-staging", "charge-tables"], calls
    assert setup["material_budget"]["limiting_term"] == "incidence" and setup["material_budget"]["total_bytes"] == 910
    published["terms"]["incidence"] = 1
    assert setup["material_budget"]["terms"]["incidence"] == 400, "budget must be an owned copy"


@pytest.mark.parametrize("change", [
    dict(schema="radia.hdiv-flat-material-budget.v0"), dict(total_bytes=909), dict(limiting_term="mass_triples"),
    dict(terms=dict(reference_tables=10)), dict(terms=dict(_budget()["terms"], material_state=0)),
    dict(terms=dict(_budget()["terms"], material_state=True)), dict(terms=dict(_budget()["terms"], material_state=2.)),
    dict(terms=dict(reversed(list(_budget()["terms"].items()))))])
def test_material_budget_publication_is_rejected_when_inconsistent(change):
    budget = _budget(**change)
    if "terms" in change and "total_bytes" not in change and all(type(v) is int for v in budget["terms"].values()):
        budget["total_bytes"] = sum(budget["terms"].values())
    with pytest.raises(RuntimeError, match="HDiv material budget"):
        hdiv_gpu._checked_material_budget(budget)


class _Owner:
    def __init__(self, calls, name, fail=None, stats=None):
        calls.append(name)
        if fail:
            raise fail
        self.calls, self.name, self._stats = calls, name, stats or {}

    def assembly(self):
        return {}

    def stats(self):
        return self._stats

    def close(self):
        self.calls.append(self.name + "-close")


def _fake_workflow(monkeypatch, calls, material_fail=None, used=900, bound=910):
    native = types.SimpleNamespace(
        _ResidentBHMaterial=types.SimpleNamespace(from_flat_tables=lambda **a: _Owner(
            calls, "material", material_fail, dict(device_payload_bound_bytes=bound))),
        _ResidentHDiv=types.SimpleNamespace(from_flat_charges=lambda **a: _Owner(calls, "operator")))
    package = types.ModuleType("hdivfake")
    gpu = types.ModuleType("hdivfake.gpu")
    gpu._native = lambda name: native
    monkeypatch.setitem(sys.modules, "hdivfake", package)
    monkeypatch.setitem(sys.modules, "hdivfake.gpu", gpu)
    monkeypatch.setattr(hdiv_gpu, "__package__", "hdivfake")
    monkeypatch.setattr(hdiv_gpu, "_prepare_device_components", lambda *a, **k: dict(
        dofs=4, family="affine-tet", rhs_operator=None, material_budget=_budget(), charge_gram={},
        evaluation_tables={}))
    monkeypatch.setattr(hdiv_gpu, "_checked_material_assembly", lambda assembly, tables: dict(
        persistent_device_bytes=used - 100, setup_transient_device_bytes=100))
    monkeypatch.setattr(hdiv_gpu, "_checked_device_assembly", lambda assembly, gram: {})
    prepared = hdiv_gpu.PreparedHDivNonlinear.__new__(hdiv_gpu.PreparedHDivNonlinear)
    law = dict(knots=[0., 1.], coefficients=[0.] * 4, integral_base=[0., 0.])
    return prepared, law


def test_material_context_and_its_live_budget_precede_the_gram(monkeypatch):
    calls = []
    prepared, law = _fake_workflow(monkeypatch, calls)
    prepared._initialize("mesh.vol", law, 1, dict(hdiv_gpu.GRAM_ACA_PRESET), 1)
    assert calls == ["material", "operator"], calls
    assert prepared.stats()["material_budget"]["limiting_term"] == "incidence"
    calls.clear()
    budget_error = RuntimeError("Flat material device budget exceeded: limiting term incidence")
    prepared, law = _fake_workflow(monkeypatch, calls, material_fail=budget_error)
    with pytest.raises(RuntimeError, match="limiting term incidence"):
        prepared._initialize("mesh.vol", law, 1, dict(hdiv_gpu.GRAM_ACA_PRESET), 1)
    assert calls == ["material"], "no Gram may be built after a material budget failure"


def test_allocations_beyond_the_preflight_bound_fail_and_release(monkeypatch):
    calls = []
    prepared, law = _fake_workflow(monkeypatch, calls, used=911, bound=910)
    with pytest.raises(RuntimeError, match="exceed the preflight budget"):
        prepared._initialize("mesh.vol", law, 1, dict(hdiv_gpu.GRAM_ACA_PRESET), 1)
    assert calls == ["material", "material-close"], calls


def test_aca_metadata_rows_are_optional_but_positive():
    arguments = gram(staged())
    assert hdiv_gpu._checked_device_assembly(ledger(arguments, *[(p, d, 8) for p, d in hdiv_gpu._ACA_OPTIONAL_ROWS]),
                                             arguments)["compression"] == "aca"
    assert ("aca_permutation", "H2D") in hdiv_gpu._expected_device_ledger(arguments)
    assert hdiv_gpu._expected_device_ledger(arguments)[("aca_permutation", "H2D")] == 8*2


def test_disabled_far_split_has_no_low_rule_rows():
    arguments = gram(staged(far=False))
    payloads = {p for p, _ in hdiv_gpu._expected_device_ledger(arguments)}
    assert not {"tet_low_rule", "tri_low_rule"} & payloads, sorted(payloads)
    hdiv_gpu._checked_device_assembly(ledger(arguments), arguments)
    with pytest.raises(RuntimeError, match="rejected"):
        hdiv_gpu._checked_device_assembly(ledger(arguments, ("tet_low_rule", "H2D", 0)), arguments)


@pytest.mark.parametrize("expo", [
    [2**31 - 1, 2**31 - 1, 0], [2**31 - 1, 0, 0], [-1, 0, 0], [0, -2**31, 0], [1, 1, 0], [0, 0, 2]])
def test_cell_exponents_are_bounded_before_any_sum(expo):
    inputs = staged()
    inputs["charge_expo"] = expo + [0, 0, 0]
    with pytest.raises(ValueError, match="cell degree <= 1"):
        gram(inputs)


@pytest.mark.parametrize("expo,ok", [
    ([2, 0, 0], True), ([1, 1, 0], True), ([0, 0, 1], False), ([2, 1, 0], False), ([3, 0, 0], False)])
def test_face_exponents_accept_the_degree_two_border_only(expo, ok):
    inputs = staged()
    inputs["charge_expo"] = [1, 0, 0] + expo
    if ok:
        assert gram(inputs)["charge_exponents"][3:].tolist() == expo
    else:
        with pytest.raises(ValueError, match="face degree <= 2"):
            gram(inputs)


def test_native_guards_bound_components_before_the_degree_sum():
    source = (HDIV / "hdiv_flat_gram.cuh").read_text(encoding="utf-8")
    kernel = source[source.index("__global__ void prepareFlatCharges"):source.index("__global__ void assembleFlatGram")]
    assert kernel.index("q.expo[2]>limit") < kernel.index("q.degree=q.expo[0]+q.expo[1]+q.expo[2]"), "device sums before bounding"
    validate = source[source.index("void validate()const"):source.index("__device__ inline double flatPow")]
    assert "expo[3*a+c]>limit" in validate, "host FlatBlock::validate lacks per-component exponent bounds"


def test_device_gram_source_downloads_only_validity_flags():
    source = (HDIV / "hdiv_flat_gram.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 1 and "sizeof(int)" in downloads[0], downloads
    assert "cudaMemcpyToSymbol" not in source and "cudaMemcpyHostToDevice" not in source


def test_diagnostic_entries_stay_off_the_construction_path():
    source = (HDIV / "hdiv_flat_diagnostic.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 1 and "count*sizeof(double)" in downloads[0], downloads
    assert "cudaMemcpyHostToDevice" not in source and "Plan" not in source and "materializeDeferred" not in source
    owner = (HDIV / "hdiv_resident_python_cuda.cu").read_text(encoding="utf-8")
    create = owner[owner.index("std::shared_ptr<DeviceOperator> createFlat("):]
    assert "flatDiagnosticEntries" not in create, "construction must not evaluate diagnostic entries"
    users = sorted(p.name for p in (ROOT / "src").rglob("*") if p.is_file() and p.suffix in {".cu", ".cuh", ".cpp", ".hpp", ".h"}
                   and '"hdiv_flat_diagnostic.cuh"' in p.read_text(encoding="utf-8", errors="replace"))
    assert users == ["hdiv_resident_python_cuda.cu"], users


def test_device_aca_source_downloads_only_integer_metadata():
    source = (HDIV / "hdiv_flat_aca.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert len(downloads) == 2 and all("sizeof(int)" in d for d in downloads), downloads
    assert {d.split(",")[1] for d in downloads} == {"active.data", "report.data"}, downloads
    assert re.search(r"Buffer<int> active;", source) and re.search(r"Buffer<int> report;", source)
    report = source[source.index("__global__ void acaReport"):source.index("struct AcaDenseSpan")]
    assert "frobenius2" not in report and "maxAbs" not in report, "ACA report leaks norms"
    assert "cudaMemcpyHostToDevice" not in source and "cudaMemcpyToSymbol" not in source


def test_flat_binding_accepts_no_matrix_payload():
    source = (HDIV / "hdiv_resident_python.cpp").read_text(encoding="utf-8")
    start = source.index('def_static("from_flat_charges"')
    block = source[start:source.index('.def("assembly"', start)]
    names = re.findall(r'py::arg\("([a-z_]+)"\)', block)
    assert names[:10] == ["cell_vertices", "face_vertices", "charge_host", "charge_kind",
                          "charge_exponents", "tet_rule", "tri_rule", "tet_low_rule",
                          "tri_low_rule", "far_factor"], names
    assert not {"leaves", "permutation", "scale", "a1", "a2"} & set(names), names


def test_python_counter_schema_matches_the_native_publication():
    source = (HDIV / "hdiv_resident_python_cuda.cu").read_text(encoding="utf-8")
    block = source[source.index("result->record.counters={"):]
    block = block[:block.index("};")]
    native = re.findall(r'\{"([a-z_]+)",s\.', block)
    assert native == list(hdiv_gpu._ACA_COUNTERS), (native, hdiv_gpu._ACA_COUNTERS)


def test_python_charge_budget_terms_match_the_native_publication():
    source = (HDIV / "hdiv_resident_python_cuda.cu").read_text(encoding="utf-8")
    budget = source[source.index("struct ChargeBudget"):source.index("class Operator")]
    native = tuple(re.findall(r'"([a-z_]+)"', re.search(r"names\[count\]=\{([^}]*)\}", budget).group(1)))
    assert native == hdiv_gpu._CHARGE_BUDGET_TERMS, (native, hdiv_gpu._CHARGE_BUDGET_TERMS)
    assert int(re.search(r"static constexpr std::size_t count=(\d+);", budget).group(1)) == len(native)


def test_flat_route_checks_the_budget_before_allocation_without_fixed_caps():
    source = (HDIV / "hdiv_resident_python_cuda.cu").read_text(encoding="utf-8")
    create = source[source.index("std::shared_ptr<DeviceOperator> createFlat("):]
    create = create[:create.index("\n}\n")]
    assert "4608" not in create and "16384" not in create, "fixed charge or DOF cap on the flat route"
    assert create.index("chargeBudget(") < min(create.index("FlatGramOracle oracle"), create.index("FlatAcaOracle oracle")), \
        "charge budget runs after oracle construction"
    flat = source[source.index("Operator(const radia_hmatrix_cuda::Plan& plan,Oracle& oracle"):]
    flat = flat[:flat.index("\n    }\n")]
    assert "plan,nullptr,0)" in flat, "flat route keeps the combined 2 GiB setup cap"
    assert "builder.nonzeros(charges),radia_hdiv_bh::FlatChargeMapBuilder::scratchBound(flat,charges)," in flat \
        and "flatReserve(" not in flat, "charge-map scratch must be reserved together with the map, not separately"
    resident = (HDIV / "hmatrix_resident.cuh").read_text(encoding="utf-8")
    build_map = resident[resident.index("void buildChargeMap("):resident.index("void materializeDeferred(")]
    # B, B^T, the DOF vector and the builder scratch are live together: one checked sum, one reserve test.
    assert "fitsAboveReserve(free,total,checkedSum({indices,values,scratch}))" in build_map, build_map
    assert build_map.index("fitsAboveReserve(") < build_map.index("build(chargeMap,chargeTranspose)"), build_map
    assert "add_test(NAME radia.hdiv_device_reserve.cpp" in (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    for name, cap in (("hdiv_flat_gram.cuh", "n>4608"), ("hdiv_flat_charge_map.cuh", "cells>4096"),
                      ("hdiv_flat_charge_map.cuh", "dofs>16384"), ("hdiv_charge_staging.hpp", ">4608")):
        assert cap not in (HDIV / name).read_text(encoding="utf-8"), (name, cap)
    aca = (HDIV / "hdiv_flat_aca.cuh").read_text(encoding="utf-8")
    compress = aca[aca.index("void compress(){"):aca.index("perm.upload(")]
    assert "flatReserve(live" in compress and "2*sizeof(double)" in compress, "ACA permutation uploaded outside the reserve"


def hex_staged(cell=(2, 2, 2), face=(2, 2, 0)):
    """One HEX cell and one exterior quad with the CPU hex rules (lattices agree exactly)."""
    rule = (np.array([.25, .75]), np.array([.5, .5]))
    inputs = dict(family="hex", n_el=1, cell_nodes=np.arange(81.), face_nodes=np.arange(27.),
                  charge_host=[0, 0], charge_kind=[0, 1], charge_expo=list(cell) + list(face),
                  rules=dict(out=rule, inner=rule, pair=rule, pair_affine=rule), near_grade=1., far_inner=1.)
    charge = dict(cell_class=np.zeros(1, np.int32), cell_dofs=np.arange(6, dtype=np.int32)[None],
                  cell_tables=np.ones((1, 1, 6)), cell_exponents=np.array([cell], np.int32),
                  face_class=np.zeros(1, np.int32), face_dofs=np.arange(3, dtype=np.int32)[None],
                  face_tables=np.ones((1, 1, 3)), face_exponents=np.array([face[:2]], np.int32))
    return inputs, charge, (np.arange(81.), np.arange(27.))


def test_hex_rules_cross_as_indexed_rows_with_bands_and_no_tet_rule():
    inputs, charge, geometry = hex_staged()
    arguments = hdiv_gpu._device_gram_arguments(charge, inputs, geometry, None)
    assert arguments["family"] == "hex" and arguments["far_factor"] == 1e30
    assert arguments["hex_rules"].shape == (8, 3) and arguments["hex_rules"][:, 0].tolist() == [0, 0, 1, 1, 2, 2, 3, 3]
    assert arguments["hex_bands"].tolist() == [1., 1.]
    assert all(arguments[name].shape == (0, 4) for name in ("tet_rule", "tri_rule", "tet_low_rule", "tri_low_rule"))
    ledger = hdiv_gpu._expected_device_ledger(dict(arguments, compression="dense"))
    assert ("cell_q2_nodes_m", "H2D") in ledger and ("hex_invalid_flag", "H2D") in ledger
    assert not {p for p, _ in ledger} & {"cell_vertices_m", "tet_rule", "charge_kind", "charge_map_cell_vertices_m"}


@pytest.mark.parametrize("cell,face,ok", [
    ((2, 2, 2), (2, 2, 0), True), ((3, 0, 0), (0, 0, 0), False), ((0, 0, 0), (0, 0, 1), False),
    ((-1, 0, 0), (0, 0, 0), False), ((0, 0, 0), (3, 0, 0), False)])
def test_hex_exponents_are_bounded_per_axis_not_by_total_degree(cell, face, ok):
    inputs, charge, geometry = hex_staged(cell, face)
    if ok:
        assert hdiv_gpu._device_gram_arguments(charge, inputs, geometry, None)["charge_exponents"].tolist() == \
            list(cell) + list(face)
    else:
        with pytest.raises(ValueError, match="HEX axis degree <= 2"):
            hdiv_gpu._device_gram_arguments(charge, inputs, geometry, None)


def test_hex_lattice_disagreement_fails_before_the_device():
    inputs, charge, (cells, faces) = hex_staged()
    with pytest.raises(ValueError, match="hex geometry and the CPU route disagree"):
        hdiv_gpu._device_gram_arguments(charge, inputs, (cells + 1e-6, faces), None)


@pytest.mark.parametrize("counts,order,error", [
    ({4, 8}, 1, "pure TET or pure HEX"), ({6}, 1, "pure TET or pure HEX"), ({8}, 0, "RT0 is affine TET only")])
def test_mixed_meshes_and_hex_rt0_fail_loudly(monkeypatch, counts, order, error):
    vim = types.ModuleType("radia.vim._vim")
    vim._volume_vertex_counts = lambda mesh: frozenset(counts)
    for name, module in (("radia", types.ModuleType("radia")), ("radia.vim", types.ModuleType("radia.vim")),
                         ("radia.vim._vim", vim)):
        monkeypatch.setitem(sys.modules, name, module)
    with pytest.raises(NotImplementedError, match=error):
        hdiv_gpu._mesh_family(None, 1, order)


def test_hex_gram_downloads_only_flags_and_table_sizes():
    source = (HDIV / "hdiv_hex_gram.cuh").read_text(encoding="utf-8")
    downloads = re.findall(r"cudaMemcpy\w*\(([^;]*)DeviceToHost", source)
    assert sorted(d.split(",")[2].strip() for d in downloads) == ["sizeof(int)", "sizeof(x)", "sizeof(y)"], downloads
    assert "cudaMemcpyToSymbol" not in source and "cudaMemcpyHostToDevice" not in source
