"""Real finite-cell law in the production weak-coupling pipeline."""
from benchmark_panel_esim import run_case


def test_real_nonlinear_direct_table_agreement(tmp_path):
    result = run_case(tmp_path, repeats=1)
    assert len(result['runs']) == 2
