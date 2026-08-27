from radia_mcp.streamfunction.streamfunction_knowledge import (
    get_streamfunction_documentation,
)


def test_streamfunction_optuna_documents_production_outer_loop():
    text = get_streamfunction_documentation("optuna")
    normalized = " ".join(text.split())

    assert "radia.stream.OptunaRunner" in text
    assert "radia.optuna.createStudy" in text
    assert "Stream Function Optuna" in text
    assert "one complete Stream Function application" in normalized
    assert "never once per Simulink timestep" in normalized
    assert "Keep `aca_eps` fixed" in text
    assert "C++ ACA+ factorization" in normalized
