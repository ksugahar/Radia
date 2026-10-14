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
foreach ($removed in @('--texclip', '--clipboard-tex-to-png',
    '--copy-tex', '--copy-tex-file', '--copy-google-slides', '--copy-google-slides-file',
    '--copy-png', '--copy-png-file', '--render-png', '--render-emf',
    '--render-png-file', '--render-emf-file', '--copy-unknown', '--render-unknown')) {
    $result = Invoke-EqneditCaptured @($removed, 'input.tex', 'output.png')
    if ($result.ExitCode -ne 94 -or $result.Stdout -or
        -not $result.Stderr.Contains('Removed compatibility option') -or
        -not $result.Stderr.Contains('clipboard-png')) {
        throw "Removed option was not rejected with a migration diagnostic: $removed"
    }
}
$removedPng = Invoke-EqneditCaptured @('input.tex', 'png')
if ($removedPng.ExitCode -ne 94 -or -not $removedPng.Stderr.Contains('Removed output alias png')) {
    throw 'Ambiguous png clipboard alias was not rejected.'
}
foreach ($destination in @('office', 'powerpoint', 'word', 'slides', 'google-slides',
    'clipboard-png', 'output.png', 'output.emf', 'output.pptxmath')) {
    $missing = Invoke-EqneditCaptured @('missing-cli-contract-input.tex', $destination)
    if ($missing.ExitCode -ne 82) {
        throw "Current destination no longer reaches the checked input reader: $destination"
    }
}
Write-Host '[OK] Removed compatibility inputs reject; current destinations remain supported.'
New-Item -ItemType Directory -Path 'C:\temp' -Force | Out-Null
$invalidColorPath = Join-Path 'C:\temp' ('eqnedit64-invalid-color-' + [guid]::NewGuid().ToString('N') + '.tex')
try {
    [IO.File]::WriteAllText($invalidColorPath, '\color{unsupported}{x}')
    $invalidColor = Invoke-EqneditCaptured @($invalidColorPath, 'office')
    if ($invalidColor.ExitCode -ne 97 -or -not $invalidColor.Stderr.Contains('Unsupported colour name')) {
        throw 'Unsupported colour CLI input failed silently.'
    }
} finally {
    Remove-Item -LiteralPath $invalidColorPath -ErrorAction SilentlyContinue
}
