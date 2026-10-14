"""Solver-neutral preflight for physics-scoped result expressions."""
from __future__ import annotations
import json

def physics_result_preflight_gate(summary_json: str) -> dict:
    row=json.loads(summary_json)
    if not isinstance(row,dict): raise ValueError("summary_json must decode to an object")
    tag=str(row.get("physics_tag") or ""); expr=[str(x) for x in row.get("expressions") or []]
    entities=row.get("selection_entities") or []; failed=row.get("failed_expression_probes") or []
    metrics=row.get("result_metrics") or {}
    checks={"physics_tag_recorded":bool(tag),"physics_type_recorded":bool(row.get("physics_type")),"selection_nonempty":bool(entities) and all(isinstance(x,int) and x>0 for x in entities),"expressions_nonempty":bool(expr),"expressions_use_exact_physics_namespace":bool(tag) and all(x.startswith(tag+".") for x in expr),"study_and_solution_recorded":bool(row.get("study_tag")) and bool(row.get("solution_tag")),"failed_probes_distinguished_from_results":all(isinstance(x,dict) and x.get("status")=="rejected" for x in failed),"required_license_features_available":row.get("required_license_features_available") is True,"solve_and_evaluation_completed":row.get("solve_completed") is True and row.get("result_evaluation_completed") is True,"result_metric_keys_match_expressions":set(metrics)==set(expr),"result_metrics_finite_nonnegative":bool(metrics) and all(isinstance(v,(int,float)) and v>=0 and v<float("inf") for v in metrics.values()),"owned_result_scope_cleaned":row.get("owned_scope_cleaned") is True}
    return {"policy":"physics_result_preflight_gate_v1","status":"ok" if all(checks.values()) else "needs_attention","checks":checks,"issues":[k for k,v in checks.items() if not v]}
