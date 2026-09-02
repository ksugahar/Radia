param(
    [Parameter(Mandatory = $true)]
    [string]$AppPath
)

$ErrorActionPreference = 'Stop'
$app = (Resolve-Path -LiteralPath $AppPath).Path

function Invoke-EqneditCaptured {
    param([string[]]$Arguments)
    $info = [Diagnostics.ProcessStartInfo]::new()
    $info.FileName = $app
    $info.UseShellExecute = $false
    $info.RedirectStandardOutput = $true
    $info.RedirectStandardError = $true
    foreach ($argument in $Arguments) {
        [void]$info.ArgumentList.Add($argument)
    }
    $process = [Diagnostics.Process]::Start($info)
    $stdout = $process.StandardOutput.ReadToEnd()
    $stderr = $process.StandardError.ReadToEnd()
    if (-not $process.WaitForExit(10000)) {
        $process.Kill($true)
        throw "Eqnedit64 CLI timed out: $($Arguments -join ' ')"
    }
    [pscustomobject]@{
        ExitCode = $process.ExitCode
        Stdout = $stdout
        Stderr = $stderr
    }
}

$help = Invoke-EqneditCaptured @('--help')
if ($help.ExitCode -ne 0 -or
    -not $help.Stdout.Contains('Eqnedit64.exe <入力> <出力>') -or
    -not $help.Stdout.Contains('clipboard-png') -or
    $help.Stderr) {
    throw "Redirected --help contract failed: $($help | ConvertTo-Json -Compress)"
}

$version = Invoke-EqneditCaptured @('--version')
if ($version.ExitCode -ne 0 -or
    $version.Stdout -notmatch '数式エディタ64 [0-9]+\.[0-9]+\.[0-9]+' -or
    -not $version.Stdout.Contains('ビルド:') -or
    $version.Stderr) {
    throw "Redirected --version contract failed: $($version | ConvertTo-Json -Compress)"
}

$invalid = Invoke-EqneditCaptured @('input.tex', 'output.svg')
if ($invalid.ExitCode -ne 94 -or $invalid.Stdout -or
    -not $invalid.Stderr.Contains('clipboard-png')) {
    throw "Invalid-output diagnostic contract failed: $($invalid | ConvertTo-Json -Compress)"
}

Write-Host '[OK] Eqnedit64 CLI help, version, and errors are redirectable.'
