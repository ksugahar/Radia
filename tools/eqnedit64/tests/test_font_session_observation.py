"""Exercise the PowerShell acceptance gate without starting any native process."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "build/test_font_session.ps1"


@pytest.mark.parametrize("scenario,code,status,launches", [
    ("healthy", 0, "PASS", 2),
    ("empty-query", 0, "PASS", 2),
    ("missing", 2, "INCONCLUSIVE", 0),
    ("log-denied", 2, "INCONCLUSIVE", 0),
    ("active-log-denied", 2, "INCONCLUSIVE", 2),
    ("control-crash", 2, "INCONCLUSIVE", 0),
    ("pid-changed", 2, "INCONCLUSIVE", 1),
    ("pid-disappeared", 2, "INCONCLUSIVE", 1),
    ("active-crash", 1, "FAIL", 2),
    ("app-failed", 1, "FAIL", 1),
])
def test_observation_contract(tmp_path, scenario, code, status, launches):
    pwsh = shutil.which("pwsh")
    if pwsh is None:
        pytest.skip("PowerShell is required for mocked gate tests")
    app = tmp_path / "not-an-executable.txt"
    app.write_text("never execute this file", encoding="utf-8")
    output = tmp_path / "observation.json"
    harness = tmp_path / "mock.ps1"
    harness.write_text(r'''
param($Gate, $App, $Output, $Scenario)
$env:EQNEDIT64_ISOLATED_TEST_SESSION='1'
$global:launches=0
$global:queries=0
function Get-Process {
    param($Id, $Name, $ErrorAction)
    if ($Id) { return [pscustomobject]@{SessionId=99} }
    if ($Scenario -eq 'missing') { return }
    if ($Scenario -eq 'pid-disappeared' -and $global:launches) { return }
    $number=100
    if ($Scenario -eq 'pid-changed' -and $global:launches) { $number=101 }
    return [pscustomobject]@{SessionId=99; Id=$number}
}
function Get-WinEvent {
    [CmdletBinding()] param($FilterHashtable)
    $global:queries++
    if ($Scenario -eq 'log-denied') { throw 'Access denied to Application log' }
    if ($Scenario -eq 'active-log-denied' -and $global:queries -eq 2) {
        throw 'Application log became unavailable'
    }
    if ($Scenario -eq 'empty-query') {
        $record=[Management.Automation.ErrorRecord]::new(
            [Exception]::new('No events'), 'NoMatchingEventsFound',
            [Management.Automation.ErrorCategory]::ObjectNotFound, $null)
        $PSCmdlet.ThrowTerminatingError($record)
    }
    if (($Scenario -eq 'control-crash' -and $global:queries -eq 1) -or
        ($Scenario -eq 'active-crash' -and $global:queries -eq 2)) {
        return [pscustomobject]@{TimeCreated=$FilterHashtable.StartTime;
            Message='fontdrvhost.exe access violation'; RecordId=42}
    }
}
function Start-Process {
    param($FilePath, $ArgumentList, $WindowStyle, [switch]$Wait, [switch]$PassThru)
    $global:launches++
    if ($WindowStyle -ne 'Hidden') { throw 'Native window requested' }
    $code=0
    if ($Scenario -eq 'app-failed') { $code=17 }
    return [pscustomobject]@{ExitCode=$code}
}
function Start-Sleep { param($Seconds) }
& $Gate -AppPath $App -Iterations 2 -OutputPath $Output
exit $LASTEXITCODE
''', encoding="utf-8")
    run = subprocess.run([pwsh, "-NoProfile", "-File", str(harness), str(SCRIPT),
                          str(app), str(output), scenario], capture_output=True,
                         text=True, timeout=20)
    assert run.returncode == code, run.stdout + run.stderr
    result = json.loads(output.read_text(encoding="utf-8-sig"))
    assert result["status"] == status
    assert len(result["iterations"]) == launches
    assert result["app_sha256"].lower() == hashlib.sha256(app.read_bytes()).hexdigest()
    assert result["session_id"] == 99
    assert result["reason"]
