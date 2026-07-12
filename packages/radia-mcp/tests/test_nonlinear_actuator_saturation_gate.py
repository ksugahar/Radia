import copy, json
from radia_mcp.radia_ngsolve.server import nonlinear_actuator_saturation_knee_gate


def rows():
    currents=[0.1,0.2,0.4,0.8,1.2,1.6]
    inductance=[6.31,7.12,6.93,5.01,3.72,2.96]
    force_i2=[94.95,120.47,114.22,60.83,34.10,21.87]
    return [{"current_a":i,"flux_linkage_wb_turn":i*l,"secant_inductance_h":l,"magnetic_energy_j":0.1*(k+1),"magnetic_coenergy_j":0.1*(k+1) if k<5 else 0.8,"weighted_stress_radial_force_n":0.0,"weighted_stress_axial_force_n":i*i*f,"axial_force_per_i2_n_per_a2":f,"fixed_iron_b_t":0.2*(k+1),"plunger_b_t":0.27*(k+1),"node_count":11097,"element_count":21941} for k,(i,l,f) in enumerate(zip(currents,inductance,force_i2))]


def test_accepts_initial_permeability_rise_then_shared_saturation_knee():
    result=json.loads(nonlinear_actuator_saturation_knee_gate(json.dumps({"rows":rows()})))
    assert result["status"]=="ok"
    assert result["metrics"]["knee_current_a"]==0.2


def test_rejects_linear_i_squared_response_without_saturation():
    data=rows()
    for row in data:
        row["secant_inductance_h"]=7.0; row["flux_linkage_wb_turn"]=7.0*row["current_a"]; row["axial_force_per_i2_n_per_a2"]=100.0; row["weighted_stress_axial_force_n"]=100.0*row["current_a"]**2
    result=json.loads(nonlinear_actuator_saturation_knee_gate(json.dumps({"rows":data})))
    assert result["status"]=="needs_attention"
    assert result["checks"]["shared_nonendpoint_knee"] is False


def test_rejects_different_inductance_and_force_knees():
    data=copy.deepcopy(rows()); data[2]["axial_force_per_i2_n_per_a2"]=130.0
    result=json.loads(nonlinear_actuator_saturation_knee_gate(json.dumps({"rows":data})))
    assert result["status"]=="needs_attention"
    assert result["checks"]["shared_nonendpoint_knee"] is False
