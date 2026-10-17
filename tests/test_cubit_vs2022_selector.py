"""Exercise Cubit's compiler selector without an installed compiler or build."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('version,exit_code,toolset_exists,accepted', [
    ('17.14.0', 0, True, True), ('18.11.0', 0, True, False),
    ('19.0.0', 0, True, False), (None, 0, True, False), ('18.11.0', 1, True, False), ('17.14.0', 0, False, False),
])
def test_cubit_selector_requires_vs2022(version, exit_code, toolset_exists, accepted):
    pwsh = shutil.which('pwsh')
    if not pwsh:
        pytest.skip('PowerShell selector requires pwsh')
    instances = [] if version is None else [dict(
        installationVersion=version, installationPath='synthetic-vs-install')]
    # Extract only the production function with the PowerShell parser, so no
    # dependency discovery, environment mutation, cleanup or native build runs.
    source = str(ROOT / 'src/cubit_plugin/cubit_build.ps1').replace("'", "''")
    payload = json.dumps(instances).replace("'", "''")
    script = f"""
$ErrorActionPreference='Stop'
$tokens=$null; $errors=$null
$ast=[System.Management.Automation.Language.Parser]::ParseFile('{source}',[ref]$tokens,[ref]$errors)
if ($errors.Count) {{ throw 'PowerShell parse errors' }}
$fn=$ast.Find({{param($n) $n -is [System.Management.Automation.Language.FunctionDefinitionAst] -and $n.Name -eq 'Find-CubitVisualStudio2022'}},$true)
. ([scriptblock]::Create($fn.Extent.Text))
function FakeVswhere {{
    $global:observedArgs=@($args)
    $global:LASTEXITCODE={exit_code}
    '{payload}'
}}
try {{
    function Test-Path {{ return ${str(toolset_exists).lower()} }}
    $value=Find-CubitVisualStudio2022 'FakeVswhere'
    $ok=$true
}} catch {{ $value=$_.Exception.Message; $ok=$false }}
@{{accepted=$ok;value=$value;args=$global:observedArgs}} | ConvertTo-Json -Compress
"""
    result = subprocess.run([pwsh, '-NoProfile', '-NonInteractive', '-Command', script],
                            capture_output=True, text=True, check=True, timeout=30)
    data = json.loads(result.stdout)
    assert data['accepted'] is accepted, data
    args = data['args']
    assert args[args.index('-version') + 1] == '[17,18)'
    assert args[args.index('-requires') + 1] == 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64'
    if accepted:
        assert data['value'] == 'synthetic-vs-install'
    else:
        assert 'Visual Studio Build Tools 2022' in data['value'] or 'MSVC toolset' in data['value']

@pytest.mark.parametrize('matching', [True, False])
def test_cubit_rejects_stale_compiler_cache(tmp_path, matching):
    pwsh = shutil.which('pwsh')
    if not pwsh:
        pytest.skip('PowerShell selector requires pwsh')
    compiler = tmp_path / 'selected' / 'cl.exe'
    cached = compiler if matching else tmp_path / 'other' / 'cl.exe'
    (tmp_path / 'CMakeCache.txt').write_text(
        f'CMAKE_CXX_COMPILER:FILEPATH={cached}\nCMAKE_C_COMPILER:FILEPATH={cached}\n')
    source = str(ROOT / 'src/cubit_plugin/cubit_build.ps1').replace("'", "''")
    script = f"""
$tokens=$null;$errors=$null
$ast=[System.Management.Automation.Language.Parser]::ParseFile('{source}',[ref]$tokens,[ref]$errors)
$fn=$ast.Find({{param($n) $n -is [System.Management.Automation.Language.FunctionDefinitionAst] -and $n.Name -eq 'Assert-CubitCompilerCache'}},$true)
. ([scriptblock]::Create($fn.Extent.Text))
try {{ Assert-CubitCompilerCache '{tmp_path}' '{compiler}'; 'accepted' }} catch {{ 'rejected' }}
"""
    p = subprocess.run([pwsh, '-NoProfile', '-NonInteractive', '-Command', script],
                       capture_output=True, text=True, check=True, timeout=30)
    assert p.stdout.strip() == ('accepted' if matching else 'rejected')
