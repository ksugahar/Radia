from pathlib import Path


SOURCE = (Path(__file__).resolve().parents[1] / "validation_test" /
          "surface_rom_team10" / "transient.py")


def test_legacy_surface_rom_fom_uses_shared_true_residual_gate():
    source = SOURCE.read_text(encoding="utf-8")
    assert "from radia._residual_gate import RELATIVE_LIMIT, check_true_residual" in source
    assert "relative = check_true_residual(" in source
    assert "reference_norm=linear_reference_norm" in source
    assert "max_relative_linear_residual=step_linear_residual" in source
    assert "linear_true_residual_limit=RELATIVE_LIMIT" in source
