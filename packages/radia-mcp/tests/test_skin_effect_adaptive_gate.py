import json, math
from radia_mcp.radia_ngsolve.server import skin_effect_adaptive_energy_loss_gate as mcp_gate
from radia_mcp.radia_ngsolve.skin_effect_adaptive_gate import skin_effect_adaptive_energy_loss_gate


def _case():
    f=1000.0; i=1+0j; z=-0.02-0.4j; v=z*i; p=v*i.conjugate(); w=abs(z.imag)/(2*2*math.pi*f)
    psi=v/(-1j*2*math.pi*f)
    c=lambda x:{"real":x.real,"imag":x.imag}
    rows=[{"mesh_cells":n,"energy_j":w*(1+d),"loss_w":abs(p.real)*(1-e),"adaptive_error":a}
          for n,d,e,a in [(100,0.02,.25,.2),(200,.01,.15,.12),(400,.002,.08,.08),(800,.001,.04,.04),(1600,0,0,.02)]]
    return dict(frequency_hz=f,current=c(i),voltage=c(v),impedance=c(z),power=c(p),flux_linkage=c(psi),total_energy_j=w,total_loss_w=abs(p.real),adaptive_rows=rows)


def test_skin_effect_gate_accepts_closed_port_identities_and_dispatches():
    case=_case(); result=skin_effect_adaptive_energy_loss_gate(**case)
    assert result["status"]=="ok"
    assert json.loads(mcp_gate(**case))["status"]=="ok"


def test_skin_effect_gate_rejects_stale_loss_and_nonincreasing_mesh():
    case=_case(); case["adaptive_rows"][-1]["loss_w"]*=0.7; case["adaptive_rows"][-1]["mesh_cells"]=700
    result=skin_effect_adaptive_energy_loss_gate(**case)
    assert result["status"]=="needs_attention"
    assert result["checks"]["mesh_cells_strictly_increase"] is False
    assert result["checks"]["final_loss_matches_adaptive_history"] is False
