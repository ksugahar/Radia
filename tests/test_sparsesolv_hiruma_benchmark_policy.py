"""Static contract for the published Hiruma AMS benchmark formulation."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "src/ext/sparsesolv/examples/hiruma/bench_compact_ams.py"


def test_hiruma_epsilon_shift_is_preconditioner_only():
    source = BENCHMARK.read_text(encoding="utf-8")
    system_body = source.split("a = BilinearForm(fes)", 1)[1].split(
        "f = LinearForm(fes)", 1
    )[0]
    surrogate_body = source.split("a_real = BilinearForm(fes_real)", 1)[1].split(
        "a_real.Assemble()", 1
    )[0]

    assert "preconditioner_eps" not in system_body
    assert "preconditioner_eps * nu_cf * u_r * v_r * dx" in surrogate_body
    assert '"system_epsilon_regularization": 0.0' in source
    assert '"epsilon_placement": "preconditioner_only"' in source
