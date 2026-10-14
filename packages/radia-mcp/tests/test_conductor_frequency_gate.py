import copy
import json

from radia_mcp.radia_ngsolve.conductor_frequency_gate import twin_conductor_skin_effect_frequency_gate
from radia_mcp.radia_ngsolve.server import twin_conductor_skin_effect_frequency_gate as mcp_gate


FREQUENCIES = [20, 1000, 5000, 10000, 20000, 30000, 40000, 50000, 100000]
RESISTANCE = [
    [0.005098729053313812, 0.005098724447575972],
    [0.005103618873241277, 0.00510361136846435],
    [0.005218289652215128, 0.005218210725374226],
    [0.005546964797643319, 0.005546676493814849],
    [0.0065681829874930965, 0.0065672058089034325],
    [0.007749092785118857, 0.007747263565854886],
    [0.008919714548224987, 0.00891705340555678],
    [0.010026305878767499, 0.01002296268141749],
    [0.014511053490392088, 0.01450764825143959],
]
INDUCTANCE = [
    [1.2954448216181467e-07, 1.2953157659179423e-07],
    [1.295223959050115e-07, 1.2950888701941803e-07],
    [1.2900847227983817e-07, 1.2899526551267324e-07],
    [1.275754041250983e-07, 1.2756301706769884e-07],
    [1.2349784015980382e-07, 1.2348799611564974e-07],
    [1.1942660827160411e-07, 1.194197599234571e-07],
    [1.1594099552731967e-07, 1.159371690607628e-07],
    [1.130452713210417e-07, 1.1304423394497461e-07],
    [1.0425141748742023e-07, 1.04258354268295e-07],
]


def test_twin_conductor_live_shape_passes_and_dispatches():
    result = twin_conductor_skin_effect_frequency_gate(FREQUENCIES, RESISTANCE, INDUCTANCE)
    assert result["status"] == "ok"
    assert result["metrics"]["frequency_count"] == 9
    assert json.loads(mcp_gate(FREQUENCIES, RESISTANCE, INDUCTANCE))["status"] == "ok"


def test_twin_conductor_gate_rejects_active_asymmetric_and_nonmonotone_rows():
    bad_r = copy.deepcopy(RESISTANCE)
    bad_l = copy.deepcopy(INDUCTANCE)
    bad_r[4][1] = -bad_r[4][1]
    bad_r[5][0] = bad_r[4][0] * 0.5
    bad_l[6][1] = bad_l[5][1] * 1.2
    result = twin_conductor_skin_effect_frequency_gate(FREQUENCIES, bad_r, bad_l)
    assert result["status"] == "needs_attention"
    assert result["checks"]["resistance_positive"] is False
    assert result["checks"]["resistance_non_decreasing"] is False
    assert result["checks"]["inductance_non_increasing"] is False
    assert result["checks"]["twin_conductor_symmetry"] is False
