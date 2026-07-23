from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_native_sfunctions_have_lifecycle_hooks_and_error_status():
    for name in ("radia_ih_eddy_sfun.cpp", "radia_ih_thermal_sfun.cpp"):
        text = (ROOT / "src" / "radia" / "simulink" / name).read_text(encoding="utf-8")
        assert "mdlStart" in text
        assert "mdlTerminate" in text
        assert "ssSetErrorStatus" in text
        assert "SS_OPTION_EXCEPTION_FREE_CODE" in text


def test_thermal_sfunction_updates_once_per_discrete_step_and_transports_rotation():
    text = (
        ROOT / "src" / "radia" / "simulink" / "radia_ih_thermal_sfun.cpp"
    ).read_text(encoding="utf-8")
    outputs = text.split("static void mdlOutputs", 1)[1].split("#define MDL_UPDATE", 1)[0]
    update = text.split("static void mdlUpdate", 1)[1].split("static void mdlTerminate", 1)[0]
    assert "advance_thermal" not in outputs
    assert "advance_thermal" in update
    assert "transport_periodic" in update
    assert "ssSetInputPortDirectFeedThrough(S, port, 0)" in text
    assert "ssSetDWorkUsedAsDState(S, 0, 1)" in text
    assert "mdlInitializeConditions" in text


def test_eddy_sfunction_rotates_heat_and_rejects_implicit_angle_support():
    text = (
        ROOT / "src" / "radia" / "simulink" / "radia_ih_eddy_sfun.cpp"
    ).read_text(encoding="utf-8")
    assert "transport_periodic" in text
    assert "-(angle - c->angle_origin_rad)" in text
    assert "rotation_mode is 'none'" in text
    assert "heat_cell_weights" in text
    assert "does not yet implement nonlinear BH iteration" in text
    assert "must contain finite values" in text


def test_vol_checker_is_a_native_simulink_preflight_dependency():
    text = (ROOT / "matlab" / "+radia" / "+simulink" / "validateVolFiles.m").read_text(encoding="utf-8")
    assert "check-vol" in text
    assert "VolCheckFailed" in text
    assert "--report-json" in text
