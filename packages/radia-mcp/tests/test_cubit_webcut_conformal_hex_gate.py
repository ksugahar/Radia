import copy, json
from radia_mcp.cubit.server import cubit_webcut_conformal_hex_gate, cubit_webcut_journal_execution_gate


def good():
    return {"source_kind":"source_native_local_journal","source_sha256":"a"*64,"execution_mode":"python_api_headless","headless_flags":["-nographics","-batch"],"gui_daemon_enabled":False,"commands":["reset","cylinder","webcut xplane","webcut yplane","imprint all","merge all","mesh volume all","block add face in surface with area > 1.56 and z_min < -0.4 and z_max > 0.4"],"volume_ids":[1,2,3,4],"element_counts":{"hex":24,"tet":0,"wedge":0,"pyramid":0},"webcut_volume_relative_drift":5.3e-6,"quarter_volume_relative_spread":0.0,"interfaces":[{"adjacent_volumes":[1,2],"face_count":4,"area":1.0} for _ in range(4)],"boundary_block_face_count":16,"quality":{"scaled_jacobian":{"min":.85},"shape":{"min":.83}},"process_exit_code":3,"startup_diagnostics":["ERROR: Could not open file: C:/x/plugins","ERROR: Could not open file: -commandplugindir","ERROR: Could not open file: -nojournal"],"script_error_lines":[],"result_artifact_fresh":True,"owned_processes_remaining":0,"public_gate_status":"ok"}


def test_accepts_live_webcut_geometry_and_source_lanes():
    row=good(); assert json.loads(cubit_webcut_conformal_hex_gate(row))["status"]=="ok"; assert json.loads(cubit_webcut_journal_execution_gate(row))["status"]=="ok"


def test_rejects_boolean_drift_and_unshared_interface():
    row=copy.deepcopy(good()); row["webcut_volume_relative_drift"]=2e-5; row["interfaces"][0]["adjacent_volumes"]=[1]; result=json.loads(cubit_webcut_conformal_hex_gate(row)); assert result["status"]=="needs_attention"; assert result["checks"]["webcut_volume_drift_bounded"] is False


def test_rejects_wrong_journal_order_and_script_error():
    row=copy.deepcopy(good()); row["commands"]=["mesh volume all","webcut xplane","webcut yplane","imprint all","merge all"]; row["script_error_lines"]=["SyntaxError"]; result=json.loads(cubit_webcut_journal_execution_gate(row)); assert result["status"]=="needs_attention"; assert result["checks"]["imprint_merge_before_mesh"] is False; assert result["checks"]["process_exit_acceptable"] is False
