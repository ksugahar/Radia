"""Guards: CI cannot accumulate output in fixed host paths or projection scratch.

Assumes: literal path declarations in workflow YAML/run blocks, runtime paths
rooted by runner/workspace variables, and the CI projection allocation entry.
This static contract is not an arbitrary shell/Python data-flow sandbox.
"""
import importlib.util
from pathlib import Path
import re

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def fixed_paths(value):
    # Ignore URLs, runtime expressions/variables and full-line comments, not
    # arbitrary code after '#'. Relative workspace outputs and cmd /d /c remain.
    value = '\n'.join(line for line in value.splitlines() if not line.lstrip().startswith('#'))
    value = re.sub(r'https?://[^\s\"\']+', 'URL', value)
    value = re.sub(r'(?:\$\{\{.*?\}\}|\$\{[^}]+\}|\$env:[A-Za-z_][A-Za-z_0-9]*|\$[A-Za-z_][A-Za-z_0-9]*)(?:[/\\][^\s\"\';]*)?', 'RUNTIME', value)
    windows = re.findall(r'\b[A-Za-z]:[/\\][^\s\"\']*|\\\\[A-Za-z0-9_.-]+\\[^\s\"\']+', value)
    unix = re.findall(r'(?<![A-Za-z0-9_.$:/])/(?:[A-Za-z0-9_.-]+/)+[^\s\"\';]*', value)
    return windows + [path for path in unix if path != '/dev/null']


def strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from strings(child)


def test_all_workflow_literals_use_runner_or_workspace_scratch():
    violations = []
    for path in sorted((ROOT / '.github/workflows').glob('*.yml')):
        data = yaml.load(path.read_text(encoding='utf-8'), Loader=yaml.BaseLoader)
        violations.extend((path.name, found) for value in strings(data) for found in fixed_paths(value))
    assert not violations, f'Fixed host scratch paths; use RUNNER_TEMP/workspace: {violations}'


@pytest.mark.parametrize('value', [
    r'RADIA_CI_OUTPUT_DIR: Z:\scratch\radia-123',
    'python build.py --out D:/scratch/result',
    r'$out = "E:\builds\scratch"',
    'mkdir -p /tmp/result',
    'python build.py --out /scratch/fixed/result',
    'New-Item -Path "' + chr(92) * 2 + 'example.invalid' + chr(92) + 'share' + chr(92) + 'outputs"',
])
def test_fixed_host_scratch_cannot_hide_in_an_environment_or_command(value):
    assert fixed_paths(value), value


@pytest.mark.parametrize('value', [
    'python build.py --out "$env:RUNNER_TEMP/result/log.json"',
    '${{ runner.temp }}/result/**/*.json',
    '${{ github.workspace }}/build/report.json',
    'python build.py --out build/report.json',
    'cmd.exe /d /c build.cmd',
    'https://example.org/docs/reference',
    'command > /dev/null',
])
def test_runner_paths_workspace_paths_and_command_switches_are_not_fixed_scratch(value):
    assert not fixed_paths(value), value


def projection_builder():
    spec = importlib.util.spec_from_file_location('projection_scratch', ROOT / 'packages/radia-simulink/tools/build_projection_cuda.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_projection_ci_build_is_created_under_actual_runner_temp(monkeypatch, tmp_path):
    module = projection_builder()
    monkeypatch.setenv('GITHUB_ACTIONS', 'true')
    monkeypatch.setenv('RUNNER_TEMP', str(tmp_path))
    output = module.create_build_directory()
    assert output.parent == tmp_path and output.is_dir()
    assert output.name.startswith('projection-cuda-')


@pytest.mark.parametrize('value', [None, 'relative-runner-temp'])
def test_projection_ci_missing_or_relative_runner_temp_fails_loudly(monkeypatch, value):
    module = projection_builder()
    monkeypatch.setenv('GITHUB_ACTIONS', 'true')
    if value is None:
        monkeypatch.delenv('RUNNER_TEMP', raising=False)
    else:
        monkeypatch.setenv('RUNNER_TEMP', value)
    with pytest.raises(RuntimeError, match='RUNNER_TEMP'):
        module.create_build_directory()
