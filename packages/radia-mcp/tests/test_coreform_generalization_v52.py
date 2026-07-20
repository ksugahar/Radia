from copy import deepcopy

from radia_mcp.cubit.imprint_import_identity_v52 import ACIS, APREPRO, SCHEME, SHEET, validate_public_identity, validate_source_identity


CASE_IDS = {
    "v52_public_sheet_imprint_merge_tolerance_curvesplit_topology_owner_mismatch",
    "v52_public_volume_scheme_assignment_autosmooth_curveinterval_seed_owner_mismatch",
    "v52_source_tool_acis_import_tolerance_healing_bodyname_layer_owner_mismatch",
    "v52_source_tool_aprepro_variable_scope_include_order_expression_unit_owner_mismatch",
}


def _generation(prefix: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": prefix, **{name: prefix for name in names}}


def _records():
    split = {"curve:1": ["curve:2", "curve:3"]}
    sheet = {**_generation("sheet-v52", ("tolerance_generation", "split_generation", "topology_generation", "owner_generation", "result_generation")), "merge_tolerance": 1e-7, "result_merge_tolerance": 1e-7, "curve_split_map": split, "result_curve_split_map": split, "topology_counts": {"sheets": 2, "surfaces": 2, "curves": 4, "vertices": 4}, "result_topology_counts": {"sheets": 2, "surfaces": 2, "curves": 4, "vertices": 4}, "body_owner": "headless:sheet-v52", "result_body_owner": "headless:sheet-v52", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64}
    schemes = {"volume:1": "map"}; intervals = {"curve:1": 8}; seeds = {"volume:1": 52}
    scheme = {**_generation("scheme-v52", ("scheme_generation", "smooth_generation", "interval_generation", "seed_generation", "owner_generation", "result_generation")), "volume_schemes": schemes, "result_volume_schemes": schemes, "autosmooth": True, "result_autosmooth": True, "curve_intervals": intervals, "result_curve_intervals": intervals, "volume_seeds": seeds, "result_volume_seeds": seeds, "volume_owner": "headless:scheme-v52", "result_volume_owner": "headless:scheme-v52", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64}
    healing = ["stitch", "remove_sliver"]; names = ["rotor", "airgap"]; layers = {"rotor": "moving", "airgap": "interface"}
    acis = {**_generation("acis-v52", ("tolerance_generation", "healing_generation", "body_generation", "layer_generation", "owner_generation", "result_generation")), "import_tolerance": 1e-6, "replayed_import_tolerance": 1e-6, "healing_operations": healing, "replayed_healing_operations": healing, "body_names": names, "replayed_body_names": names, "body_layers": layers, "replayed_body_layers": layers, "geometry_owner": "headless:acis-v52", "replayed_geometry_owner": "headless:acis-v52", "result_sha256": "c" * 64, "accepted_result_sha256": "c" * 64}
    scopes = {"global": {"gap": {"value": 0.8, "unit": "mm"}}}; includes = ["units.apr", "mesh.apr"]; expressions = {"airgap_total": {"expression": "2*gap", "unit": "mm", "resolved": 1.6}}
    aprepro = {**_generation("apr-v52", ("scope_generation", "include_generation", "expression_generation", "unit_generation", "owner_generation", "result_generation")), "variable_scopes": scopes, "replayed_variable_scopes": scopes, "include_order": includes, "replayed_include_order": includes, "expressions": expressions, "replayed_expressions": expressions, "journal_owner": "headless:apr-v52", "replayed_journal_owner": "headless:apr-v52", "result_sha256": "d" * 64, "accepted_result_sha256": "d" * 64}
    return {SHEET: sheet, SCHEME: scheme, ACIS: acis, APREPRO: aprepro}


def test_v52_positive_public_and_source_replays_are_accepted():
    records = _records()
    assert validate_public_identity(records)["status"] == "ok"
    assert validate_source_identity(records)["status"] == "ok"


def test_v52_public_mutations_are_rejected():
    records = deepcopy(_records()); records[SHEET]["result_merge_tolerance"] = 1e-3
    assert validate_public_identity(records)["status"] == "needs_attention"


def test_v52_source_mutations_are_rejected():
    records = deepcopy(_records()); records[APREPRO]["replayed_include_order"] = list(reversed(records[APREPRO]["include_order"]))
    assert validate_source_identity(records)["status"] == "needs_attention"


def test_v52_invalid_canonical_records_are_rejected():
    records = deepcopy(_records()); records[SCHEME]["volume_schemes"] = {"volume:1": "tetmesh"}; records[SCHEME]["result_volume_schemes"] = {"volume:1": "tetmesh"}
    assert validate_public_identity(records)["status"] == "needs_attention"
