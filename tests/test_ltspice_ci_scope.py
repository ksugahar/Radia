from pathlib import Path

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
    code = job['steps'][1]['run']
    assert "runtests('tests/matlab/test_ltspice_data_safety.m')" in code
    assert 'all([r.Passed])' in code and '~any([r.Incomplete])' in code
    assert 'start_matlab' in code and 'eng.quit()' in code
    assert 'child.wait(timeout=180)' in code
    assert "'/PID', str(child.pid), '/T', '/F'" in code
    assert 'Build.ps1' not in code and 'pip install' not in code
