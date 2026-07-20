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
    assert "IH temporarily" in contract["applications"]["notebook_policy"]
    assert contract["ltspice"]["executable"] == "LTspice.exe"
    assert contract["ltspice"]["legacy_ltc_versions"] == "not supported"
