"""Cross-kernel STEP mass/topology diagnosis with evidence-status separation."""
from __future__ import annotations
import math


def cross_kernel_mass_topology_diagnosis_gate(summary, *, self_rtol=1e-10, external_rtol=1e-6, mode_rtol=1e-12):
    if not isinstance(summary,dict): raise ValueError("summary must be a mapping")
    native=summary.get("native") or {}; self_row=summary.get("self_roundtrip") or {}; imports=summary.get("external_imports") or []
    def num(row,key):
        value=float(row[key])
        if not math.isfinite(value) or value<=0: raise ValueError(f"{key} must be positive and finite")
        return value
    nv,na=num(native,"volume"),num(native,"area"); sv,sa=num(self_row,"volume"),num(self_row,"area")
    rows=[{"mode":str(row.get("mode") or "").lower(),"volume":num(row,"volume"),"area":num(row,"area"),"volume_count":int(row.get("volume_count",-1)),"surface_count":int(row.get("surface_count",-1)),"curve_count":int(row.get("curve_count",-1))} for row in imports]
    if not rows: raise ValueError("external_imports must not be empty")
    self_v=abs(sv-nv)/nv; self_a=abs(sa-na)/na; ext_v=max(abs(x["volume"]-nv)/nv for x in rows); ext_a=max(abs(x["area"]-na)/na for x in rows); vspread=(max(x["volume"] for x in rows)-min(x["volume"] for x in rows))/nv; aspread=(max(x["area"] for x in rows)-min(x["area"] for x in rows))/na
    modes={x["mode"] for x in rows}; face=int(native["face_count"]); edge=int(native["edge_count"]); topology=all(x["volume_count"]==1 and x["surface_count"]==face and x["curve_count"]==edge for x in rows)
    self_ok=max(self_v,self_a)<=self_rtol; modes_ok=max(vspread,aspread)<=mode_rtol; external_ok=max(ext_v,ext_a)<=external_rtol and topology
    diagnosis="step_self_roundtrip_loss" if not self_ok else "external_import_mode_inconsistency" if not modes_ok else "external_kernel_mass_topology_loss" if not external_ok else "portable"
    evidence={"source_identity_bound":summary.get("source_kind")=="upstream_native_example" and len(str(summary.get("source_sha256","")))==64,"step_digest_bound":len(str(summary.get("step_sha256","")))==64,"self_roundtrip_closes":self_ok,"heal_noheal_present":{"heal","noheal"}.issubset(modes),"external_modes_invariant":modes_ok,"single_external_volume":all(x["volume_count"]==1 for x in rows)}
    return {"policy":"build123d_cross_kernel_mass_topology_diagnosis_gate_v1","status":"ok" if all(evidence.values()) else "needs_attention","portable":diagnosis=="portable","diagnosis":diagnosis,"healing_not_root_cause":diagnosis=="external_kernel_mass_topology_loss" and modes_ok,"evidence_checks":evidence,"issues":[k for k,v in evidence.items() if not v],"portability_checks":{"external_mass_matches":max(ext_v,ext_a)<=external_rtol,"external_topology_matches":topology},"metrics":{"self_volume_relative_error":self_v,"self_area_relative_error":self_a,"external_volume_relative_error":ext_v,"external_area_relative_error":ext_a,"external_volume_mode_spread":vspread,"external_area_mode_spread":aspread,"native_face_count":face,"native_edge_count":edge,"external_surface_counts":[x["surface_count"] for x in rows],"external_curve_counts":[x["curve_count"] for x in rows]}}


def upstream_source_external_cad_contract_gate(summary):
    if not isinstance(summary,dict): raise ValueError("summary must be a mapping")
    checks={"upstream_source_bound":summary.get("source_kind")=="upstream_native_example" and len(str(summary.get("source_sha256","")))==64,"source_preserved":summary.get("source_preserved") is True,"official_assertion_reproduced":summary.get("official_assertion_reproduced") is True,"step_digest_bound":len(str(summary.get("step_sha256","")))==64,"self_roundtrip_valid":summary.get("self_roundtrip_valid") is True,"external_headless_import_recorded":summary.get("external_headless") is True,"diagnosis_gate_passed":summary.get("diagnosis_gate_status")=="ok","portability_decision_recorded":summary.get("diagnosis") in {"portable","external_kernel_mass_topology_loss"}}
    return {"policy":"build123d_upstream_source_external_cad_contract_gate_v1","status":"ok" if all(checks.values()) else "needs_attention","checks":checks,"issues":[k for k,v in checks.items() if not v]}
