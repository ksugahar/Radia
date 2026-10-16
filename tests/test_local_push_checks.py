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
