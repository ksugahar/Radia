from radia_mcp.matlab import matlab_simulink_library_contract


def test_simulink_library_registration_and_ltspice_scope():
    contract = matlab_simulink_library_contract()
    assert contract["status"] == "ready"
    assert contract["schema"] == "radia-mcp.matlab-simulink-library/v2"
    assert contract["registration_code"] == [
        "radia.simulink.buildLibrary",
        "sl_refresh_customizations",
    ]
    assert "LTspice/LTspice Circuit" in contract["blocks"]
    assert contract["blocks"][:5] == [
        "Applications/Electromagnet",
        "Applications/PCB PEEC",
        "Applications/Motor",
        "Applications/Stream Function",
        "Applications/Induction Heating",
    ]
    assert contract["applications"]["initial_backend"] == "python-headless-cli"
    assert contract["applications"]["per_step_python"] == "forbidden"
    assert "optional" in contract["applications"]["mex_policy"]
    assert "retired for every application" in contract["applications"]["notebook_policy"]
    preflight = contract["applications"]["mesh_preflight"]
    assert preflight["checker"] == "check-vol"
    assert preflight["label_contract_schema"] == "radia.vol-label-contract.v1"
    assert preflight["report_schema"] == "cubit-mesh-export.vol-check.v1"
    assert contract["ltspice"]["executable"] == "LTspice.exe"
    assert contract["ltspice"]["legacy_ltc_versions"] == "not supported"
