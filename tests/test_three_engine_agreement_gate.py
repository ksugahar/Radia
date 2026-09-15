"""One-percent core policy; historical solver acceptance is not rewritten."""
import ast
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUNNERS = (
    'validation_test/esrf_three_engine/run_coil_yoke_three_engine.py',
    'validation_test/c_type_three_engine/run_three_engine.py',
    'validation_test/esrf_three_engine/run_hybrid_undulator_three_engine.py',
)


@pytest.mark.parametrize('path,expected', zip(RUNNERS, (0.01, 0.03, 0.05)))
def test_default_core_tolerance(path, expected):
    tree = ast.parse((ROOT / path).read_text(encoding='utf-8'))
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Attribute) and node.func.attr == 'add_argument'
             and node.args and isinstance(node.args[0], ast.Constant)
             and node.args[0].value == '--relative-rms-tolerance']
    assert len(calls) == 1
    assert next(ast.literal_eval(k.value) for k in calls[0].keywords
                if k.arg == 'default') == expected


@pytest.mark.parametrize('path', ['AGENTS.md', 'CLAUDE.md'])
def test_policy_declares_scope(path):
    text = (ROOT / path).read_text(encoding='utf-8')
    assert 'ESRF coil-yoke HDiv/FEM validation defaults to 1 % relative RMS on its declared core stencil' in text


def test_posthoc_core_agreement_preserves_historical_threshold_and_raw_error():
    path = ROOT / ('validation_test/esrf_three_engine/results/candidate_59b094d8/'
                   'three_engine_case6_bdm1_bonus12.json')
    result = json.loads(path.read_text(encoding='utf-8'))
    assert result['relative_rms_tolerance'] == 0.03
    assert max(p['relative_rms'] for p in result['pairwise_core'].values()) == pytest.approx(0.007818616219259033)
    assert max(p['relative_rms'] for p in result['pairwise_core'].values()) < 0.01
    assert max(p['relative_rms'] for p in result['pairwise_raw'].values()) == pytest.approx(0.01407951932211685)
    assert max(p['relative_rms'] for p in result['pairwise_raw'].values()) > 0.01
