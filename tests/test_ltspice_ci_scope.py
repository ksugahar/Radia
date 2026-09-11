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
    assert job['runs-on'] == ['self-hosted', 'Windows', 'X64', 'mdx']
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
            return any(fnmatchcase(path, pattern) for pattern in patterns)

        assert not selected('tests/matlab/fixtures/ltspice_numdgt7_binary.raw')
        assert not selected('tests/matlab/fixtures/kicad_ltspice_bridge.kicad_sch')
        # The Optuna/LTspice objective really consumes this shared fixture.
        assert selected('tests/matlab/fixtures/ltspice_rc.cir')
        assert selected('tests/matlab/fixtures/generate_sobol_direction_numbers.py')
        for fixture in (root / 'tests/matlab/fixtures').glob('*optuna*'):
            assert selected(fixture.relative_to(root).as_posix())
