"""Guards: CLI plans cannot omit real native/CLI consumers or hide unknown inputs.

Assumes: immutable, committed miniature repositories with independently declared
consumer chains; no solver/GPU execution. This checks the public plan contract,
not numerical accuracy or the completeness of the real repository's closures.
"""
import json
from pathlib import Path
import subprocess
import sys

import pytest

CLI = Path(__file__).resolve().parents[1] / 'tools/ci_numerical_plan.py'


def git(root, *args):
    return subprocess.check_output(['git', *args], cwd=root, stderr=subprocess.PIPE).decode().strip()


@pytest.fixture
def repository(tmp_path):
    root = tmp_path / 'source'; root.mkdir()
    git(root, 'init', '-q'); git(root, 'config', 'user.email', 'contract@example.invalid')
    git(root, 'config', 'user.name', 'Public source contract')
    files = {
        'tests/ams.py': 'from radia import sparsesolv_ngsolve\n',
        'tests/viewer.py': 'from radia.panels import open_gmsh\n',
        'tests/output.py': 'from radia.gmsh_post_export import write\n',
        'tests/native_contract.py': 'x = 1\n',
        'tests/conftest.py': '',
        'validation_test/motor.py': "import subprocess\ncmd = ['python', 'calc_motor.py']\nsubprocess.run(cmd)\n",
        'validation_test/ih.py': 'from radia.ih import solve\n',
        'src/radia/__init__.py': '', 'src/radia/panels/__init__.py': '',
        'src/radia/panels/open_gmsh.py': 'def show(): return 1\n',
        'src/radia/panels/calc_motor.py': 'from radia.gmsh_post_export import write\n',
        'src/radia/gmsh_post_export.py': 'import gmsh\ndef write(): return 1\n',
        'src/radia/ih.py': 'def solve(): return 1\n',
        'pyproject.toml': '[project]\nname="contract"\nversion="1"\ndependencies=["gmsh>=5.0.0.dev0", "ngsolve==6.2.2608"]\n',
    }
    for name, content in files.items():
        path = root / name; path.parent.mkdir(parents=True, exist_ok=True); path.write_text(content)
    manifest = {'schema': 'radia.test-tier-manifest.v1', 'profiles': {
        'sparsesolv': {'paths': ['tests/ams.py']},
        'solver-numerics': {'paths': ['tests/viewer.py', 'tests/output.py']},
        'solver-heavy': {'paths': ['validation_test/motor.py', 'validation_test/ih.py']},
        'native-build-contracts': {'paths': ['tests/native_contract.py']},
    }, 'impact_rules': {}}
    (root / 'tests/test_tier_manifest.json').write_text(json.dumps(manifest))
    git(root, 'add', '.'); git(root, 'commit', '-qm', 'Independent consumer graph')
    return root, git(root, 'rev-parse', 'HEAD')


def change(root, path, text):
    file = root / path; file.parent.mkdir(parents=True, exist_ok=True); file.write_text(text)
    git(root, 'add', '.'); git(root, 'commit', '-qm', 'Contract input change')


def invoke(root, base, tmp_path, event='pull_request', expected=0):
    report = tmp_path / 'plan.json'
    process = subprocess.run([sys.executable, str(CLI), '--repo', str(root), '--base', base,
                              '--event', event, '--out', str(report)], capture_output=True, text=True)
    assert process.returncode == expected, (process.stdout, process.stderr, report.read_text() if report.exists() else '')
    data = json.loads(report.read_text())
    assert data['qualification'] is False and data['execution_wired'] is False
    return data


def test_viewer_cli_plan_does_not_select_ams_or_heavy(repository, tmp_path):
    root, base = repository
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    plan = invoke(root, base, tmp_path)
    assert plan['ams_paths'] == [] and plan['heavy_case_paths'] == []
    assert plan['python_test_paths'] == ['tests/viewer.py']
    assert plan['unclassified_inputs'] == []


def test_plan_uses_committed_blobs_despite_worktree_replacement(repository, tmp_path):
    root, base = repository
    change(root, 'src/radia/panels/open_gmsh.py', '# UTF-8: 球\ndef show(): return "\\n"\n')
    committed = invoke(root, base, tmp_path)
    (root / 'src/radia/panels/open_gmsh.py').write_text('invalid Python (')
    observed = invoke(root, base, tmp_path)
    assert observed == committed
    assert observed['ams_paths'] == [] and observed['heavy_case_paths'] == []


def test_writer_cli_plan_keeps_motor_subprocess_consumer(repository, tmp_path):
    root, base = repository
    change(root, 'src/radia/gmsh_post_export.py', 'import gmsh\ndef write(): return 2\n')
    plan = invoke(root, base, tmp_path)
    assert plan['ams_paths'] == []
    assert plan['heavy_case_paths'] == ['validation_test/motor.py']
    assert 'tests/output.py' in plan['python_test_paths']
    motor = next(row for row in plan['cases'] if row['path'] == 'validation_test/motor.py')
    assert 'src/radia/panels/calc_motor.py' in motor['dependencies']


def test_native_change_cli_plan_selects_ams_and_build(repository, tmp_path):
    root, base = repository
    change(root, 'src/ext/sparsesolv/ams.h', '// native input changed\n')
    plan = invoke(root, base, tmp_path)
    assert plan['ams_paths'] == ['tests/ams.py']
    assert 'sparsesolv_ngsolve' in plan['native_build_targets']
    assert set(plan['heavy_case_paths']) == {'validation_test/motor.py', 'validation_test/ih.py'}


def test_actual_ams_consumer_change_selects_it(repository, tmp_path):
    root, base = repository
    change(root, 'tests/ams.py', 'from radia import sparsesolv_ngsolve\n# new assertion\n')
    assert invoke(root, base, tmp_path)['ams_paths'] == ['tests/ams.py']


def test_ih_change_does_not_discard_its_case(repository, tmp_path):
    root, base = repository
    change(root, 'src/radia/ih.py', 'def solve(): return 2\n')
    assert invoke(root, base, tmp_path)['heavy_case_paths'] == ['validation_test/ih.py']


@pytest.mark.parametrize('event', ['push', 'workflow_dispatch', 'release', 'schedule'])
def test_non_pr_events_preserve_full_execution_plan(repository, tmp_path, event):
    root, base = repository
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    plan = invoke(root, base, tmp_path, event)
    assert plan['ams_paths'] and len(plan['heavy_case_paths']) == 2


def test_unknown_path_cannot_be_a_false_exemption(repository, tmp_path):
    root, base = repository
    change(root, 'src/radia/unclassified.py', 'x=1\n')
    plan = invoke(root, base, tmp_path)
    assert plan['unclassified_inputs'] == ['src/radia/unclassified.py']
    assert plan['ams_paths'] and len(plan['heavy_case_paths']) == 2


def test_dynamic_unknown_case_is_conservatively_selected(repository, tmp_path):
    root, base = repository
    change(root, 'validation_test/ih.py', "import importlib\nname='variable'\nimportlib.import_module(name)\n")
    dynamic_base = git(root, 'rev-parse', 'HEAD')
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    plan = invoke(root, dynamic_base, tmp_path)
    assert 'validation_test/ih.py' in plan['heavy_case_paths']
    assert next(row['unknowns'] for row in plan['cases'] if row['path'] == 'validation_test/ih.py')


def test_unbounded_subprocess_cannot_hide_a_consumer(repository, tmp_path):
    root, _ = repository
    change(root, 'validation_test/ih.py', 'import subprocess\nsubprocess.run(get_runtime_command())\n')
    base = git(root, 'rev-parse', 'HEAD')
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    plan = invoke(root, base, tmp_path)
    row = next(r for r in plan['cases'] if r['path'] == 'validation_test/ih.py')
    assert row['selected'] and any('subprocess' in reason for reason in row['unknowns'])


def test_parent_initializer_change_is_not_missed(repository, tmp_path):
    root, base = repository
    change(root, 'src/radia/__init__.py', 'import numpy\nnumpy.seterr(all="raise")\n')
    plan = invoke(root, base, tmp_path)
    assert plan['ams_paths'] and len(plan['heavy_case_paths']) == 2


@pytest.mark.parametrize('dependency,ams', [('gmsh', False), ('ngsolve', True)])
def test_dependency_components_do_not_treat_gmsh_as_solver(repository, tmp_path, dependency, ams):
    root, base = repository
    text = (root / 'pyproject.toml').read_text()
    text = text.replace('5.0.0.dev0', '5.0.0.dev20261024') if dependency == 'gmsh' else text.replace('6.2.2608', '6.2.9999')
    change(root, 'pyproject.toml', text)
    plan = invoke(root, base, tmp_path)
    assert bool(plan['ams_paths']) == ams
    assert 'validation_test/motor.py' in plan['heavy_case_paths']


def test_mixed_viewer_native_change_takes_union(repository, tmp_path):
    root, base = repository
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    change(root, 'src/core/native.h', '// native change\n')
    assert invoke(root, base, tmp_path)['ams_paths']


def test_rename_preserves_old_and_new_dependency_sides(repository, tmp_path):
    root, base = repository
    git(root, 'mv', 'src/radia/ih.py', 'src/radia/ih_new.py')
    git(root, 'commit', '-qm', 'Rename dependency')
    plan = invoke(root, base, tmp_path)
    assert {'src/radia/ih.py', 'src/radia/ih_new.py'} <= set(plan['changed_paths'])
    assert 'validation_test/ih.py' in plan['heavy_case_paths']


def test_invalid_base_is_blocked_not_empty_success(repository, tmp_path):
    root, _ = repository
    assert invoke(root, 'no-such-commit', tmp_path, expected=2)['status'] == 'blocked'


def test_empty_profile_is_blocked(repository, tmp_path):
    root, base = repository
    path = root / 'tests/test_tier_manifest.json'; data = json.loads(path.read_text())
    data['profiles']['sparsesolv']['paths'] = []
    change(root, 'tests/test_tier_manifest.json', json.dumps(data))
    assert invoke(root, base, tmp_path, expected=2)['status'] == 'blocked'
