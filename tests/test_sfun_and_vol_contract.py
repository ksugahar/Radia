from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_native_sfunctions_have_lifecycle_hooks_and_error_status():
    for name in ("radia_ih_eddy_sfun.cpp", "radia_ih_thermal_sfun.cpp"):
        text = (ROOT / "src" / "radia" / "simulink" / name).read_text(encoding="utf-8")
        assert "mdlStart" in text
        assert "mdlTerminate" in text
        assert "ssSetErrorStatus" in text
        assert "SS_OPTION_EXCEPTION_FREE_CODE" in text


def test_vol_checker_is_a_native_simulink_preflight_dependency():
    text = (ROOT / "matlab" / "+radia" / "+simulink" / "validateVolFiles.m").read_text(encoding="utf-8")
    assert "check-vol" in text
    assert "VolCheckFailed" in text
    assert "--report-json" in text
