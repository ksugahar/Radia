import importlib.util
import re
from pathlib import Path
SPEC = importlib.util.spec_from_file_location('local_push_checks', Path(__file__).parents[1]/'tools/local_push_checks.py')
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)

def test_positive_conflict_not_negative_fixture():
    p = re.compile('forbidden')
    assert m.positive_literal_conflicts("assert 'forbidden' in source", p) == [1]
    assert m.positive_literal_conflicts("assert 'forbidden' not in source\nfixture = 'forbidden'", p) == []

def test_mcp_cpu_dependency_and_unresolved_dynamic_fail_closed():
    assert m.minimal_import_findings('from radia import api')
    assert m.minimal_import_findings('def test():\n import radia')
    assert m.minimal_import_findings('importlib.import_module(name)')
    assert m.minimal_import_findings("importlib.import_module('radia.api')")
    assert m.minimal_import_findings("fixture = 'import radia'\nimport numpy") == []

def test_compact_no_native_execution():
    assert m.COMPACT_TESTS == ['tests/policy', 'tests/test_ci_execution_policy.py', 'tests/test_sparsesolv_ci_contract.py']

def test_exact_candidate_rejects_mutable_head(tmp_path, monkeypatch):
    monkeypatch.setattr(m, 'git', lambda *args: 'another-sha')
    import pytest
    with pytest.raises(ValueError, match='exact candidate'):
        m.plan(tmp_path, 'base', 'a'*40)

def test_existing_router_preserves_cpu_owner():
    import sys
    root = Path(__file__).parents[1]
    sys.path.insert(0, str(root/'tools'))
    try:
        impact = m.load(root/'tools/sparsesolv_ci_impact.py', 'test_local_impact')
        data = impact.manifest()
        required, _ = impact.lane_required('radia', ['src/radia/kelvin_solver.py'], lambda _: ('', ''), data)
        assert required
        required, _ = impact.lane_required('radia-mcp-matrix', ['src/radia/kelvin_solver.py'], lambda _: ('', ''), data)
        assert not required
    finally:
        sys.path.pop(0)

def test_rename_retains_old_lane_input(tmp_path):
    import subprocess
    import sys
    root = Path(__file__).parents[1]
    sys.path.insert(0, str(root/'tools'))
    try:
        impact = m.load(root/'tools/sparsesolv_ci_impact.py', 'rename_impact')
        subprocess.run(['git', 'init', '-q', str(tmp_path)], check=True)
        def git(*args):
            return subprocess.check_output(['git', '-c', 'user.name=Test', '-c', 'user.email=test@example.invalid', *args], cwd=tmp_path, text=True).strip()
        old = tmp_path/'src/radia/solver.py'
        old.parent.mkdir(parents=True)
        old.write_text('unchanged numerical source\n')
        git('add', '.')
        git('commit', '-qm', 'base')
        base = git('rev-parse', 'HEAD')
        git('mv', 'src/radia/solver.py', 'retired.txt')
        git('commit', '-qm', 'rename')
        paths = impact.changed_paths(git('diff', '--name-status', '-z', '-M', f'{base}...HEAD'))
        assert paths == {'src/radia/solver.py', 'retired.txt'}
        assert impact.lane_required('radia', paths, lambda _: ('', ''), impact.manifest())[0]
    finally:
        sys.path.pop(0)
