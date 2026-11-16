"""Guards: either validation-record reference input must run repository record checks.

Assumes: main PR/push events, the actual unfiltered policy workflow, and the
repo-contracts runner contract. These are source gates, not numerical validation.
"""
import importlib.util
import shlex
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
INPUTS = ('tools/validation_bibparse.py', 'validation_records/references.bib')
RECORD_TEST = 'tests/test_validation_records.py'


def tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'tools' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


contract = tool('ci_distribution_contract')
tier = tool('run_test_tier')


def workflow(name):
    return yaml.load((ROOT / '.github/workflows' / name).read_text(encoding='utf-8'), Loader=yaml.BaseLoader)


def require_record_route(document, paths, changed):
    for event in ('push', 'pull_request'):
        projection = contract.workflow_events(document)[event]
        assert 'main' in projection.get('branches', ['main'])
        assert 'branches-ignore' not in projection
        assert contract.selected(changed, projection), (event, changed)
    job = document['jobs']['repository-contracts']
    assert 'if' not in job and job.get('continue-on-error', 'false') == 'false'
    matches = []
    for step in job['steps']:
        for line in step.get('run', '').splitlines():
            command = shlex.split(line)
            if command[:2] == ['python', 'tools/run_test_tier.py']:
                matches.append((step, command))
    assert len(matches) == 1
    step, command = matches[0]
    assert 'if' not in step and step.get('continue-on-error', 'false') == 'false'
    assert command[2:] == ['--profile', 'repo-contracts']
    assert RECORD_TEST in paths


@pytest.mark.parametrize('changed', INPUTS)
def test_reference_list_changes_reach_always_on_records(changed):
    assert (ROOT / changed).is_file(), changed
    document = workflow('policy-lint.yml')
    paths, _ = tier.load_profile('repo-contracts')
    require_record_route(document, paths, changed)
