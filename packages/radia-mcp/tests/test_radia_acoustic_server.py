import asyncio, json
from radia_mcp.radia_acoustic import acoustic_capabilities, acoustic_usage, cq_grid_gate, fsi_preflight_gate

def test_capabilities_are_ngsolve_first_and_not_matlab():
    c=acoustic_capabilities()
    assert c["owner"]=="radia.acoustics" and c["matlab_runtime"] is False
    assert "ngsolve.bem" in c["numerical_backends"] and "fsi" in c["apis"] and "cq" in c["apis"]

def test_usage_covers_bem_fsi_and_cq():
    assert "HelmholtzSL" in acoustic_usage("ngsolve_bem")
    assert "VectorH1" in acoustic_usage("fsi")
    assert "s=delta(zeta)/dt" in acoustic_usage("cq")

def test_fsi_preflight_rejects_bad_geometry_and_material():
    assert fsi_preflight_gate(wavenumber=2.0)["ok"]
    bad=fsi_preflight_gate(wavenumber=2.0,c_longitudinal=1.0,c_transverse=1.0,radius_deviation=.1)
    assert not bad["ok"] and not bad["checks"]["positive_lame_lambda"] and not bad["checks"]["spherical_dtn_geometry"]

def test_cq_grid_pins_convention_and_symmetry():
    g=cq_grid_gate(num_time=16,time_step=.1,method="BDF2")
    assert g["ok"] and g["convention"]=="s=delta(zeta)/dt; kappa=i*s/c"
    assert not cq_grid_gate(num_time=2,time_step=.1)["ok"]

def test_server_registers_four_domain_tools():
    from radia_mcp.radia_acoustic.server import mcp, radia_acoustic_capabilities
    names={t.name for t in asyncio.run(mcp.list_tools())}
    assert {"radia_acoustic_usage","radia_acoustic_capabilities","radia_acoustic_fsi_preflight","radia_acoustic_cq_grid"} <= names
    assert json.loads(radia_acoustic_capabilities())["owner"]=="radia.acoustics"
