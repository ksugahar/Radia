from pathlib import Path
from fnmatch import fnmatchcase

import yaml


def test_ltspice_data_ci_is_scoped_and_runs_the_real_matlab_test():
    root = Path(__file__).resolve().parents[1]
    workflow = yaml.load(
        (root / '.github/workflows/ltspice-data-safety.yml').read_text(),
        Loader=yaml.BaseLoader,
    )
    events = workflow['on']
    assert events['push']['paths'] == events['pull_request']['paths']
    assert 'tests/matlab/test_ltspice_data_safety.m' in events['push']['paths']
    assert 'tests/matlab/fixtures/**' not in events['push']['paths']
    job = workflow['jobs']['matlab-data-contracts']
    assert job['runs-on'] == ['self-hosted', 'Windows', 'X64', '100']
    assert job['timeout-minutes'] == '5'
    assert job['env']['PYTHONIOENCODING'] == 'utf-8'
    code = job['steps'][1]['run']
    assert "runtests('tests/matlab/test_ltspice_data_safety.m')" in code
    assert 'all([r.Passed])' in code and '~any([r.Incomplete])' in code
    assert 'start_matlab' in code and 'eng.quit()' in code
    assert 'child.wait(timeout=180)' in code
    assert "'/PID', str(child.pid), '/T', '/F'" in code
    assert 'Build.ps1' not in code and 'pip install' not in code


def test_optuna_ci_does_not_build_for_unrelated_raw_fixtures():
    root = Path(__file__).resolve().parents[1]
    workflow = yaml.load(
        (root / '.github/workflows/radia-optuna.yml').read_text(),
        Loader=yaml.BaseLoader,
    )
    for event in ('push', 'pull_request'):
        patterns = workflow['on'][event]['paths']
        def selected(path):
            included = False
            for pattern in patterns:
                if fnmatchcase(path, pattern.lstrip('!')):
                    included = not pattern.startswith('!')
            return included

        assert not selected('tests/matlab/fixtures/ltspice_numdgt7_binary.raw')
        assert not selected('tests/matlab/fixtures/kicad_ltspice_bridge.kicad_sch')
        # The Optuna/LTspice objective really consumes this shared fixture.
        assert selected('tests/matlab/fixtures/ltspice_rc.cir')
        assert selected('tests/matlab/fixtures/generate_sobol_direction_numbers.py')
        for fixture in (root / 'tests/matlab/fixtures').glob('*optuna*'):
            assert selected(fixture.relative_to(root).as_posix())


def test_markdown_inputs_are_explicitly_scoped():
    root = Path(__file__).resolve().parents[1]
    for name, included, excluded in (
        ('radia-optuna.yml', ['packages/radia-optuna/README.md',
          'packages/radia-optuna/MATLAB_README.md', 'packages/radia-optuna/THIRD_PARTY_NOTICES.md',
          'packages/radia-optuna/OPTUNA_SIMULINK_LAB.md'], ['packages/radia-optuna/HANDOVER.md']),
        ('ltspice-data-safety.yml', [], ['tests/matlab/fixtures/ltspice_numdgt7_README.md']),
        ('radia-gpu.yml', ['docs/ltspice/README.md', 'packages/radia-gpu/README.md'],
          ['packages/radia-gpu/HANDOVER.md']),
    ):
        workflow = yaml.load((root / '.github/workflows' / name).read_text(), Loader=yaml.BaseLoader)
        for event in ('push', 'pull_request'):
            patterns = workflow['on'][event]['paths']
            assert '!**/*.md' in patterns
            def selected(path):
                result = False
                for pattern in patterns:
                    if fnmatchcase(path, pattern.lstrip('!')):
                        result = not pattern.startswith('!')
                return result
            assert all(selected(path) for path in included)
            assert not any(selected(path) for path in excluded)
