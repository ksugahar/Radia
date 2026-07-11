import copy
import json

from radia_mcp.radia_ngsolve.capacitance_matrix_gate import two_conductor_capacitance_matrix_gate
from radia_mcp.radia_ngsolve.server import two_conductor_capacitance_matrix_gate as mcp_gate


def good():
    return {"capacitance_unit":"pF","maxwell_matrix":[[13.7725494371,-8.3702899988],[-8.3702899988,38.6629718563]],"mutual_matrix":[[5.4022594383,8.3702899988],[8.3702899988,30.2926818575]]}


def test_accepts_reciprocal_matrix_representations():
    result=two_conductor_capacitance_matrix_gate(good()); assert result["status"]=="ok"; assert result["checks"]["maxwell_positive_definite"] is True; assert json.loads(mcp_gate(json.dumps(good())))["status"]=="ok"


def test_rejects_stale_mutual_representation():
    row=copy.deepcopy(good()); row["mutual_matrix"][0][0]*=1.2; result=two_conductor_capacitance_matrix_gate(row); assert result["status"]=="needs_attention"; assert result["checks"]["representations_agree"] is False


def test_rejects_nonreciprocal_maxwell_matrix():
    row=copy.deepcopy(good()); row["maxwell_matrix"][1][0]*=0.8; result=two_conductor_capacitance_matrix_gate(row); assert result["status"]=="needs_attention"; assert result["checks"]["maxwell_reciprocal"] is False
