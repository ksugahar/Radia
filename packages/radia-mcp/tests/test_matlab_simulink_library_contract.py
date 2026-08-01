from radia_mcp.matlab import matlab_simulink_library_contract


def test_simulink_library_registration_and_ltspice_scope():
    contract = matlab_simulink_library_contract()
    assert contract["status"] == "ready"
    assert contract["schema"] == "radia-mcp.matlab-simulink-library/v3"
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
    assert contract["ltspice"]["simulink_runtime"].startswith("Level-2 MATLAB")
    assert "Standalone MEX function ABI" in contract["ltspice"]["native_boundary"]
    assert contract["ltspice"]["legacy_ltc_versions"] == "not supported"
    assert contract["kicad_ltspice"]["workflow"] == [
        "radia.kicad.exportSpiceNetlist",
        "radia.kicad.prepareLTspice",
        "radia.kicad.buildLTspiceBlock",
    ]
    assert "not written back" in contract["kicad_ltspice"]["reverse_sync"]
    assert "Material Models/Material Dictionary" in contract["blocks"]
    assert "Coupling/Winding Dictionary" in contract["blocks"]
    assert "Applications/Field Study" in contract["blocks"]
    assert "Coupling/Field Study Configuration" in contract["blocks"]
    materials = contract["material_dictionary"]
    assert materials["mesh_format"] == "Netgen .vol"
    assert materials["runtime_bus"] == "RadiaMaterialBus"
    assert materials["per_step_dictionary_lookup"] is False
    assert materials["per_step_strings"] is False
    coupling = contract["circuit_field_coupling"]
    assert coupling["native_backend"].startswith("exact-ZOH")
    assert coupling["detailed_circuit_backend"] == "LTspice interval coupling"
    assert coupling["winding_bus"] == "RadiaWindingBus"
    assert coupling["command_bus"] == "RadiaMachineCommandBus"
    assert coupling["response_bus"] == "RadiaMachineResponseBus"
    assert coupling["mechanical_owner"] == "Simulink or Simscape"
    study = contract["field_study"]
    assert study["runtime_bus"] == "RadiaStudyBus"
    assert study["mesh_format"] == "Netgen .vol"
    assert study["per_step_python"] is False
    assert "steady_heat" in study["physics"]
    assert "harmonic_eddy" in study["physics"]
    assert study["harmonic_eddy_operator"].startswith("(K + j*omega")
