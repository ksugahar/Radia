from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]


def test_matlab_circuit_field_adapter_reuses_native_state_space_mex():
    maker = ROOT / "matlab" / "+radia" / "+simulink" / "makeCircuitFieldStateSpace.m"
    builder = ROOT / "matlab" / "+radia" / "+simulink" / "buildCircuitFieldStateSpaceModel.m"
    maker_text = maker.read_text(encoding="utf-8")
    builder_text = builder.read_text(encoding="utf-8")

    assert "sourceMatrix.' * (K \\ sourceMatrix)" in maker_text
    assert '"radia.circuit-field.state-space.v1"' in maker_text
    assert '"per_branch"' in maker_text
    assert '"radia_state_space_mex_sfunction"' in maker_text
    assert '"radia_state_space_mex_sfunction"' in builder_text
    assert "python_per_step" in maker_text and "false" in maker_text
    assert "field_factorization_per_step" in builder_text
