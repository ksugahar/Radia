"""Static contract for the published Hiruma AMS benchmark formulation."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_DIR = ROOT / "validation_test/sparsesolv/hiruma"
BENCHMARK = BENCHMARK_DIR / "bench_compact_ams.py"


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


def test_hiruma_validation_fixture_and_check_report_are_locked():
    import hashlib
    import json

    baseline = json.loads(
        (BENCHMARK_DIR / "compact_ams_baseline.json").read_text(encoding="utf-8")
    )
    report = json.loads(
        (BENCHMARK_DIR / "mesh1_2.5T.vol-check.json").read_text(encoding="utf-8")
    )
    mesh = BENCHMARK_DIR / baseline["fixture"]["file"]

    assert hashlib.sha256(mesh.read_bytes()).hexdigest() == baseline["fixture"]["sha256"]
    assert report["passed"] is True
    assert report["labels"]["passed"] is True
    assert report["boundary_domain_ownership"]["passed"] is True
    assert report["mesh"]["n_elements"] == baseline["fixture"]["ne"]
    assert report["mesh"]["n_points"] == baseline["fixture"]["nv"]
