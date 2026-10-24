"""Guards: CLI plans cannot omit real native/CLI consumers or hide unknown inputs.

Assumes: immutable, committed miniature repositories with independently declared
consumer chains; no solver/GPU execution. This checks the public plan contract,
not numerical accuracy or the completeness of the real repository's closures.
"""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

CLI = Path(__file__).resolve().parents[1] / 'tools/ci_numerical_plan.py'


def miniature_environment(**runner_context):
    env = os.environ.copy()
    for key in ('GITHUB_EVENT_NAME', 'GITHUB_EVENT_PATH'):
        env.pop(key, None)
    env.update(runner_context)
    return env


@pytest.fixture
def foreign_runner_context(monkeypatch, tmp_path):
    """A real outer CI event must not define a miniature repository's identity."""
    event = tmp_path / 'outer-runner-event.json'
    event.write_text(json.dumps({'before': 'f' * 40}), encoding='utf-8')
    monkeypatch.setenv('GITHUB_EVENT_NAME', 'push')
    monkeypatch.setenv('GITHUB_EVENT_PATH', str(event))


def git(root, *args):
    return subprocess.check_output(['git', *args], cwd=root, stderr=subprocess.PIPE).decode('utf-8').strip()


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
        path = root / name; path.parent.mkdir(parents=True, exist_ok=True); path.write_text(content, encoding='utf-8')
    manifest = {'schema': 'radia.test-tier-manifest.v1', 'profiles': {
        'sparsesolv': {'paths': ['tests/ams.py']},
        'solver-numerics': {'paths': ['tests/viewer.py', 'tests/output.py']},
        'solver-heavy': {'paths': ['validation_test/motor.py', 'validation_test/ih.py']},
        'native-build-contracts': {'paths': ['tests/native_contract.py']},
    }, 'impact_rules': {}}
    (root / 'tests/test_tier_manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    git(root, 'add', '.'); git(root, 'commit', '-qm', 'Independent consumer graph')
    return root, git(root, 'rev-parse', 'HEAD')


def change(root, path, text):
    file = root / path; file.parent.mkdir(parents=True, exist_ok=True); file.write_text(text, encoding='utf-8')
    git(root, 'add', '.'); git(root, 'commit', '-qm', 'Contract input change')


def invoke(root, base, tmp_path, event='pull_request', expected=0):
    report = tmp_path / 'plan.json'
    process = subprocess.run([sys.executable, str(CLI), '--repo', str(root), '--base', base,
                              '--event', event, '--out', str(report)], capture_output=True, text=True,
                             encoding='utf-8', env=miniature_environment())
    assert process.returncode == expected, (process.stdout, process.stderr, report.read_text(encoding='utf-8') if report.exists() else '')
    data = json.loads(report.read_text(encoding='utf-8'))
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
    (root / 'src/radia/panels/open_gmsh.py').write_text('invalid Python (', encoding='utf-8')
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
    text = (root / 'pyproject.toml').read_text(encoding='utf-8')
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
    path = root / 'tests/test_tier_manifest.json'; data = json.loads(path.read_text(encoding='utf-8'))
    data['profiles']['sparsesolv']['paths'] = []
    change(root, 'tests/test_tier_manifest.json', json.dumps(data))
    assert invoke(root, base, tmp_path, expected=2)['status'] == 'blocked'


def sealed_declaration(root, source, dependencies=(), entries=()):
    import hashlib
    def sha(path):
        return hashlib.sha256(subprocess.check_output(['git', 'show', 'HEAD:' + path], cwd=root)).hexdigest()
    binding = dict(sha256=[sha(source)], resolves=['dynamic import needs declaration: ' + source],
                   dependencies=list(dependencies), reason='Sealed finite generated-code fixture contract',
                   entries={path: [sha(path)] for path in entries})
    data = dict(schema='radia.numerical-inputs.v1', bindings={source: binding}, consumer_inputs={})
    change(root, 'tools/ci_numerical_inputs.json', json.dumps(data))


def test_sealed_finite_loader_can_resolve_but_source_change_cannot(repository, tmp_path):
    root, _ = repository
    change(root, 'tests/conftest.py', '__import__(selected_upstream_name)\n')
    sealed_declaration(root, 'tests/conftest.py')
    base = git(root, 'rev-parse', 'HEAD')
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    assert invoke(root, base, tmp_path)['ams_paths'] == []
    change(root, 'tests/conftest.py', '__import__(different_runtime_name)\n')
    assert invoke(root, base, tmp_path)['ams_paths']


def test_declared_cli_dependency_is_selected_and_case_scope_is_sealed(repository, tmp_path):
    root, _ = repository
    change(root, 'validation_test/ih.py', '__import__(runtime_name)\n')
    sealed_declaration(root, 'validation_test/ih.py', ['src/radia/ih.py'], ['validation_test/ih.py'])
    base = git(root, 'rev-parse', 'HEAD')
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    assert 'validation_test/ih.py' not in invoke(root, base, tmp_path)['heavy_case_paths']
    change(root, 'src/radia/ih.py', 'def solve(): return 42\n')
    assert 'validation_test/ih.py' in invoke(root, base, tmp_path)['heavy_case_paths']


def test_docstring_is_data_but_executed_text_stays_unresolved(repository, tmp_path):
    root, _ = repository
    change(root, 'src/radia/ih.py', '\"\"\"Example:\n    import omitted_example\nprose is not Python\"\"\"\ndef solve(): return 1\n')
    base = git(root, 'rev-parse', 'HEAD')
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    assert invoke(root, base, tmp_path)['heavy_case_paths'] == []
    change(root, 'src/radia/ih.py', 'exec(runtime_program)\ndef solve(): return 1\n')
    assert 'validation_test/ih.py' in invoke(root, base, tmp_path)['heavy_case_paths']


@pytest.mark.parametrize('damage', ['selection', 'checkout'])
def test_runtime_cli_rejects_tampered_plan_or_source(repository, tmp_path, damage, foreign_runner_context):
    root, base = repository
    change(root, 'src/radia/ih.py', 'def solve(): return 2\n')
    value = invoke(root, base, tmp_path)
    runtime = CLI.with_name('ci_numerical_runtime.py')
    command = [sys.executable, str(runtime), '--repo', str(root), '--plan', str(tmp_path / 'plan.json'),
               '--profile', 'solver-heavy', '--out', str(tmp_path / 'verified.json')]
    valid = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=miniature_environment())
    assert valid.returncode == 0, (valid.stdout, valid.stderr,
                                 (tmp_path / 'verified.json').read_text(encoding='utf-8'))
    assert json.loads((tmp_path / 'verified.json').read_text(encoding='utf-8'))['paths'] == ['validation_test/ih.py']
    if damage == 'selection':
        for row in value['cases']: row['selected'] = False
        (tmp_path / 'plan.json').write_text(json.dumps(value), encoding='utf-8')
    else:
        (root / 'src/radia/ih.py').write_text('def solve(): return 999\n', encoding='utf-8')
    rejected = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=miniature_environment())
    assert rejected.returncode == 2, (rejected.stdout, rejected.stderr)


def test_workflow_uses_distinct_plan_outputs_and_verified_execution():
    import yaml
    root = CLI.parents[1]
    workflow = yaml.safe_load((root / '.github/workflows/sparsesolv.yml').read_text(encoding='utf-8'))
    steps = workflow['jobs']['ams-regression']['steps']
    expected = {'Verify source provenance and run required AMS tests': 'ams_required',
                'Run required NGSolve numerical contracts': 'python_required',
                'Validate solver-heavy source or verify trusted prior evidence': 'heavy_required'}
    for name, output in expected.items():
        step = next(step for step in steps if step.get('name') == name)
        assert output in step['if'] and '--numerical-plan' in step['run']
    heavy = next(step for step in steps if step.get('id') == 'heavy-validation')
    assert "-eq 'pull_request'" not in heavy['run']
    assert 'force_validation' in heavy['run']



def test_cli_execution_targets_current_registry_not_retired_files(repository, tmp_path, foreign_runner_context):
    root, base = repository
    path = root / 'tests/test_tier_manifest.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    data['profiles']['solver-heavy']['paths'].remove('validation_test/ih.py')
    (root / 'validation_test/ih.py').unlink()
    change(root, 'tests/test_tier_manifest.json', json.dumps(data))
    value = invoke(root, base, tmp_path)
    assert value['retired_case_paths']['solver-heavy'] == ['validation_test/ih.py']
    assert value['heavy_case_paths'] == ['validation_test/motor.py']
    assert all(row['path'] != 'validation_test/ih.py' for row in value['cases'])
    runtime = CLI.with_name('ci_numerical_runtime.py')
    process = subprocess.run([sys.executable, str(runtime), '--repo', str(root),
                              '--plan', str(tmp_path / 'plan.json'), '--profile', 'solver-heavy',
                              '--out', str(tmp_path / 'verified.json')],
                             capture_output=True, text=True, encoding='utf-8', env=miniature_environment())
    assert process.returncode == 0, (process.stdout, process.stderr,
                                   (tmp_path / 'verified.json').read_text(encoding='utf-8'))
    assert json.loads((tmp_path / 'verified.json').read_text(encoding='utf-8'))['paths'] == ['validation_test/motor.py']


def test_runtime_cli_still_enforces_explicit_runner_event_and_base(repository, tmp_path):
    root, base = repository
    change(root, 'src/radia/ih.py', 'def solve(): return 2\n')
    invoke(root, base, tmp_path, event='push')
    event = tmp_path / 'matching-runner-event.json'
    event.write_text(json.dumps({'before': base}), encoding='utf-8')
    command = [sys.executable, str(CLI.with_name('ci_numerical_runtime.py')), '--repo', str(root),
               '--plan', str(tmp_path / 'plan.json'), '--profile', 'solver-heavy',
               '--out', str(tmp_path / 'verified.json')]
    env = miniature_environment(GITHUB_EVENT_NAME='push', GITHUB_EVENT_PATH=str(event))
    process = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=env)
    assert process.returncode == 0, (tmp_path / 'verified.json').read_text(encoding='utf-8')
    event.write_text(json.dumps({'before': 'f' * 40}), encoding='utf-8')
    process = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=env)
    report = json.loads((tmp_path / 'verified.json').read_text(encoding='utf-8'))
    assert process.returncode == 2 and 'base differs' in report['error'], report
    assert report['error'] in process.stderr, process.stderr
    env['GITHUB_EVENT_NAME'] = 'pull_request'
    process = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=env)
    report = json.loads((tmp_path / 'verified.json').read_text(encoding='utf-8'))
    assert process.returncode == 2 and 'event differs' in report['error'], report
    assert report['error'] in process.stderr, process.stderr


def test_cli_case_reads_joined_documentation_input(repository, tmp_path):
    root, _ = repository
    document = 'docs/induction_heating/FIN_PEEC_PRODUCTION_BOUNDARY.md'
    change(root, document, 'approved boundary\n')
    change(root, 'validation_test/ih.py', 'from pathlib import Path\n'
           'ROOT = Path(__file__).resolve().parents[1]\n'
           'doc = ROOT / "docs" / "induction_heating" / "FIN_PEEC_PRODUCTION_BOUNDARY.md"\n'
           'text = doc.read_text(encoding="utf-8")\n')
    base = git(root, 'rev-parse', 'HEAD')
    change(root, document, 'changed boundary\n')
    value = invoke(root, base, tmp_path)
    row = next(row for row in value['cases'] if row['path'] == 'validation_test/ih.py')
    assert row['selected'] and document in row['dependencies']


@pytest.mark.parametrize('program', [
    'from subprocess import run\nrun(["python", "-c", dynamic_code])\n',
    'from subprocess import Popen as launch\nlaunch(["python", "-c", dynamic_code])\n',
    'import subprocess as sp\nsp.check_output(dynamic_command)\n',
    'import os\nos.system(dynamic_command)\n',
    'import os as operating\noperating.popen(dynamic_command)\n',
])
def test_cli_aliased_process_or_shell_entry_is_not_a_proved_closure(repository, tmp_path, program):
    root, _ = repository
    change(root, 'validation_test/ih.py', program)
    base = git(root, 'rev-parse', 'HEAD')
    change(root, 'src/radia/panels/open_gmsh.py', 'def show(): return 2\n')
    value = invoke(root, base, tmp_path)
    row = next(row for row in value['cases'] if row['path'] == 'validation_test/ih.py')
    assert row['selected'] and row['unknowns']


def test_cli_records_static_third_party_dependency_for_heavy_environment(repository, tmp_path):
    root, base = repository
    change(root, 'validation_test/ih.py', 'from matplotlib import pyplot\n')
    value = invoke(root, base, tmp_path)
    row = next(row for row in value['cases'] if row['path'] == 'validation_test/ih.py')
    assert 'package:matplotlib' in row['components']


@pytest.mark.parametrize('path', [
    'packages/radia-simulink/src/radia_simulink/simulink/contracts/ih_workpiece_v1.json',
    '.github/workflows/sparsesolv.yml',
])
def test_cli_fallback_input_invalidates_trusted_heavy_receipt(repository, tmp_path, path):
    import importlib.util
    root, _ = repository
    change(root, path, 'original\n')
    base = git(root, 'rev-parse', 'HEAD')
    value = invoke(root, base, tmp_path)
    spec = importlib.util.spec_from_file_location('actual_source_receipt_gate', CLI.with_name('validation_source_gate.py'))
    gate = importlib.util.module_from_spec(spec); spec.loader.exec_module(gate)
    original = gate.source_fingerprint(root, {'packages': [('numpy', 'test')]}, value)[0]
    change(root, path, 'changed\n')
    value = invoke(root, base, tmp_path)
    assert len(value['heavy_case_paths']) == 2
    assert value['unclassified_inputs'] or value['global_selection']
    changed = gate.source_fingerprint(root, {'packages': [('numpy', 'test')]}, value)[0]
    assert changed != original, 'fallback-causing source must invalidate prior main evidence'
