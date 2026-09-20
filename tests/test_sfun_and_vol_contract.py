from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_level2_sfunctions_own_native_handle_lifecycle():
    package = ROOT / "matlab" / "+radia" / "+simulink"
    for name, prefix in (("ihEddySFunction.m", "ih.eddy"),
                         ("ihThermalSFunction.m", "ih.thermal")):
        text = (package / name).read_text(encoding="utf-8")
        assert 'RegBlockMethod("Start"' in text
        assert 'RegBlockMethod("Terminate"' in text
        assert f"'{prefix}.create'" in text
        assert f"'{prefix}.destroy'" in text
        assert "native_handle_low" in text and "native_handle_high" in text
    thermal = (package / "ihThermalSFunction.m").read_text(encoding="utf-8")
    assert 'RegBlockMethod("Update"' in thermal
    assert "DirectFeedthrough = false" in thermal


def test_thermal_sfunction_updates_once_per_discrete_step_and_transports_rotation():
    wrapper = (ROOT / "matlab" / "+radia" / "+simulink" / "ihThermalSFunction.m").read_text(encoding="utf-8")
    runtime = (ROOT / "src" / "radia" / "simulink" / "radia_ih_runtime.cpp").read_text(encoding="utf-8")
    outputs = wrapper.split("function outputs", 1)[1].split("function update", 1)[0]
    update = wrapper.split("function update", 1)[1].split("function terminate", 1)[0]
    assert "ih.thermal.update" not in outputs
    assert "ih.thermal.update" in update
    assert "advance_thermal" in runtime
    assert "transport_periodic" in runtime
    assert "initializeConditions" in wrapper


def test_eddy_sfunction_rotates_heat_and_rejects_implicit_angle_support():
    text = (ROOT / "src" / "radia" / "simulink" / "radia_ih_runtime.cpp").read_text(encoding="utf-8")
    assert "transport_periodic" in text
    assert "-(angle - config_.angle_origin_rad)" in text
    assert "changing angle requires periodic rotation" in text


def test_vol_checker_is_a_native_simulink_preflight_dependency():
    text = (ROOT / "matlab" / "+radia" / "+simulink" / "validateVolFiles.m").read_text(encoding="utf-8")
    assert "check-vol" in text
    assert "VolCheckFailed" in text
    assert "--report-json" in text


def test_native_ih_rotation_is_applied_to_both_fields():
    runtime = (ROOT / "src" / "radia" / "simulink" / "radia_ih_runtime.cpp").read_text(encoding="utf-8")
    # Exactly two rotations: temperature into the source frame on the way in,
    # and heat back out of it.  A third one used to rotate the stored thermal
    # state as well, which counted the motion twice because that state and the
    # incoming heat are both already in workpiece coordinates; it was removed
    # on 2026-09-16 while this >= 3 anchor stayed behind.  Pinning the count
    # exactly is what makes a re-introduced double rotation fail here.
    assert runtime.count("transport_periodic(") == 2
    # Both rotations are taken from the configured origin, temperature in and
    # heat back out.  The incremental form is what produced the third rotation
    # and the double count, so its absence is the contract now.
    assert "angle - config_.angle_origin_rad" in runtime
    assert "-(angle - config_.angle_origin_rad)" in runtime
    assert "angle - state_.previous_angle_rad" not in runtime


def test_native_ih_recompute_policy_distinguishes_current_and_material_changes():
    eddy = (ROOT / "src" / "radia" / "simulink" / "radia_ih_runtime.cpp").read_text(encoding="utf-8")
    assert "!config_.matrix_temperature_slope.empty() && temperature_changed" in eddy
    assert "heat *= current * current" in eddy
    assert "current_requires_solve" not in eddy


def test_legacy_ih_lut_and_lumped_interfaces_are_removed():
    legacy = (
        "buildIHControlModel.m",
        "evaluateIHEddyHeatDensityLUT.m",
        "evaluateIHPowerLUT.m",
        "ihPlantSFunction.m",
        "makeIHEddyHeatDensityLUT.m",
        "makeIHPlant.m",
        "makeIHPowerLUT.m",
        "simulateIHDrive.m",
        "simulateIHWaveform.m",
    )
    package = ROOT / "matlab" / "+radia" / "+simulink"
    assert all(not (package / name).exists() for name in legacy)
    assert not (ROOT / "matlab" / "radia_ih_plant_sfunction.m").exists()
    assert not (ROOT / "matlab" / "+radia" / "+rl" / "makeIHEnvironment.m").exists()
