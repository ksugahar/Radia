import copy,json
from radia_mcp.build123d.server import build123d_cross_kernel_mass_topology_diagnosis_gate, build123d_upstream_source_external_cad_contract_gate

def good():
    return {"source_kind":"upstream_native_example","source_sha256":"a"*64,"step_sha256":"b"*64,"native":{"volume":1015.939,"area":3671.793,"face_count":48,"edge_count":114},"self_roundtrip":{"volume":1015.939,"area":3671.793},"external_imports":[{"mode":m,"volume":804.370,"area":3812.209,"volume_count":1,"surface_count":48,"curve_count":108} for m in ("noheal","heal")]}

def test_classifies_invariant_external_kernel_loss_as_complete_evidence():
    result=json.loads(build123d_cross_kernel_mass_topology_diagnosis_gate(json.dumps(good()))); assert result["status"]=="ok"; assert result["portable"] is False; assert result["diagnosis"]=="external_kernel_mass_topology_loss"; assert result["healing_not_root_cause"] is True

def test_rejects_missing_heal_mode_evidence():
    row=copy.deepcopy(good()); row["external_imports"]=row["external_imports"][:1]; result=json.loads(build123d_cross_kernel_mass_topology_diagnosis_gate(json.dumps(row))); assert result["status"]=="needs_attention"

def test_accepts_upstream_source_contract_with_explicit_rejection():
    row={"source_kind":"upstream_native_example","source_sha256":"a"*64,"source_preserved":True,"official_assertion_reproduced":True,"step_sha256":"b"*64,"self_roundtrip_valid":True,"external_headless":True,"diagnosis_gate_status":"ok","diagnosis":"external_kernel_mass_topology_loss"}; assert json.loads(build123d_upstream_source_external_cad_contract_gate(json.dumps(row)))["status"]=="ok"
