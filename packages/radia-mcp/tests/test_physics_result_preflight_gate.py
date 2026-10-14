import json
from radia_mcp.radia_ngsolve.physics_result_preflight_gate import physics_result_preflight_gate
from radia_mcp.radia_ngsolve.server import physics_result_preflight_gate as mcp_gate
def good(): return {"physics_tag":"wave","physics_type":"transient_wave","selection_entities":[1],"expressions":["wave.normE","wave.normB"],"study_tag":"std1","solution_tag":"sol1","failed_expression_probes":[{"expression":"other.normE","status":"rejected"}],"required_license_features_available":True,"solve_completed":True,"result_evaluation_completed":True,"result_metrics":{"wave.normE":2.0,"wave.normB":1.0},"owned_scope_cleaned":True}
def test_accepts(): assert physics_result_preflight_gate(json.dumps(good()))["status"]=="ok"; assert json.loads(mcp_gate(json.dumps(good())))["status"]=="ok"
def test_rejects_namespace():
 r=good();r["expressions"]=["other.normE"];assert physics_result_preflight_gate(json.dumps(r))["status"]=="needs_attention"
def test_rejects_license():
 r=good();r["required_license_features_available"]=False;assert physics_result_preflight_gate(json.dumps(r))["status"]=="needs_attention"
