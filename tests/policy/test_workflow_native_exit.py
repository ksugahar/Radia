"""Guards: early Windows native failures cannot be overwritten by later success.

Assumes: line-level PowerShell audit recognizers, explicit shell/default routing;
this is source policy, not arbitrary PowerShell control-flow certification.
"""
import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('native_exit_audit',ROOT/'tools/audit_workflow_native_exit.py')
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

@pytest.mark.parametrize('body',[
 'python first.py\npython last.py',
 '& $py first.py\n& $py last.py',
 'git fetch\nWrite-Host $LASTEXITCODE\npython last.py',
 'python first.py\n$rc = $LASTEXITCODE\npython last.py',
 'python first.py | Tee-Object log\npython last.py',
 'python first.py `\n --option value\npython last.py',
])
def test_discarded_status_is_rejected(body):
    assert audit.unchecked(body)

@pytest.mark.parametrize('guard',[
 'if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }',
 'if ($LASTEXITCODE -ne 0) { throw "failed" }',
 'if ($LASTEXITCODE -ne 0) {\n Write-Host "failure"\n exit $LASTEXITCODE\n}',
])
def test_propagation_before_next_native_is_accepted(guard):
    assert not audit.unchecked('python first.py\n'+guard+'\npython last.py')

def test_here_string_payload_is_not_a_native_command():
    assert len(audit.native_lines("@'\npython payload.py\n'@ | python -\n")) == 1

def test_single_native_step_preserves_runner_exit_propagation():
    assert not audit.unchecked('python only.py')

def test_shell_and_defaults_are_respected():
    assert audit.powershell({}, {'runs-on':'windows-latest'}, {})
    assert audit.powershell({}, {'runs-on':'${{ matrix.runner }}'}, {})
    assert not audit.powershell({}, {'runs-on':'ubuntu-latest'}, {})
    assert not audit.powershell({}, {'runs-on':'windows-latest'}, {'shell':'bash'})
    assert audit.powershell({'defaults':{'run':{'shell':'pwsh'}}}, {'runs-on':'ubuntu-latest'}, {})


def test_all_windows_run_blocks_preserve_recognized_native_status():
    assert audit.audit(ROOT) == []


def test_inline_commands_cannot_bypass_status_checks():
    assert audit.unchecked('python first.py; python last.py')
    assert not audit.unchecked('python -c "print(1); print(2)"')
