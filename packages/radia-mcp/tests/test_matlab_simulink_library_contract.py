from radia_mcp.matlab import matlab_simulink_library_contract


def test_simulink_library_registration_and_ltspice_scope():
    contract = matlab_simulink_library_contract()
    assert contract["status"] == "ready"
    assert contract["registration_code"] == [
        "radia.simulink.buildLibrary",
        "sl_refresh_customizations",
    ]
    assert "LTspice/LTspice Circuit" in contract["blocks"]
    assert contract["ltspice"]["executable"] == "LTspice.exe"
    assert contract["ltspice"]["legacy_ltc_versions"] == "not supported"
