# Load the freshly built plugin into a disposable headless Cubit 2026.8 and
# check that both commands run.  Never copies into the Cubit install, so it
# is safe while another user has the deployed plugins open.
#   pwsh -File smoke.ps1 [-Plugin build\my_cubit_plugin.dll]
# Passes only on exit code 0, the exact expected log lines and the exact
# file contents.  On failure the work directory (journal, log, output) is kept.
param(
    [string]$CubitInstallDir = $(if ($env:CUBIT_INSTALL_DIR) { $env:CUBIT_INSTALL_DIR } else { 'C:\Program Files\Coreform Cubit 2026.8' }),
    [string]$Plugin = (Join-Path $PSScriptRoot 'build\my_cubit_plugin.dll'),
    [int]$TimeoutSeconds = 300,
    # brick x 1, size 0.25 -> 5^3 nodes, 4^3 hexes, one user block
    [int]$ExpectedNodes = 125,
    [int]$ExpectedBlocks = 1,
    [int]$ExpectedElements = 64
)
$ErrorActionPreference = 'Stop'
if (-not (Test-Path -LiteralPath $Plugin -PathType Leaf)) { throw "Plugin not found: $Plugin" }

# 2026.8 only: the console launcher is bin\cubitx.exe and the loader scans
# *.dll in CUBIT_PLUGIN_DIR.  Older releases used coreform_cubit.com with
# -commandplugindir; this smoke does not cover them.
$bin = Join-Path $CubitInstallDir 'bin'
$console = Join-Path $bin 'cubitx.exe'
if (-not (Test-Path -LiteralPath $console)) { throw "cubitx.exe not found in $bin (this smoke supports Cubit 2026.8 only)" }

$work = Join-Path ([IO.Path]::GetTempPath()) ("cubit_plugin_smoke_" + [guid]::NewGuid().ToString('N'))
$plugins = Join-Path $work 'plugins'
New-Item -ItemType Directory -Force -Path $plugins | Out-Null
$p = $null
try {
    Copy-Item -LiteralPath $Plugin -Destination (Join-Path $plugins ([IO.Path]::GetFileNameWithoutExtension($Plugin) + '.dll'))
    # A non-ASCII directory proves the UTF-8 -> UTF-16 path conversion.
    $outDir = Join-Path $work ([string][char]0x65E5 + [char]0x672C + [char]0x8A9E)   # 日本語
    New-Item -ItemType Directory -Force -Path $outDir | Out-Null
    $out = Join-Path $outDir 'counts.txt'
    $jou = Join-Path $work 'driver.jou'
    @(
        'reset'
        'brick x 1'
        'volume 1 size 0.25'
        'mesh volume 1'
        'block 1 volume 1'
        'myplugin hello "smoke" count 2'
        ('export nodecount "{0}" overwrite' -f $out)
        'exit 0'
    ) | Set-Content -LiteralPath $jou -Encoding utf8

    $psi = [Diagnostics.ProcessStartInfo]::new($console)
    foreach ($a in '-nographics', '-batch', '-nojournal', '-noinitfile', $jou) { $psi.ArgumentList.Add($a) }
    $psi.WorkingDirectory = $work
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.RedirectStandardInput = $true   # closed below: Cubit must not wait on stdin
    $psi.Environment['CUBIT_PLUGIN_DIR'] = $plugins
    $acis = Join-Path $CubitInstallDir 'acis\code\bin'
    if (Test-Path $acis) { $psi.Environment['PATH'] = "$acis;$($psi.Environment['PATH'])" }

    $deadline = [DateTime]::UtcNow.AddSeconds($TimeoutSeconds)
    $p = [Diagnostics.Process]::Start($psi)
    $p.StandardInput.Close()
    $stdout = $p.StandardOutput.ReadToEndAsync()
    $stderr = $p.StandardError.ReadToEndAsync()
    $remaining = { [Math]::Max(0, [int]($deadline - [DateTime]::UtcNow).TotalMilliseconds) }
    $timedOut = -not $p.WaitForExit((& $remaining))
    if ($timedOut) {
        $p.Kill($true)
        [void]$p.WaitForExit(30000)
    }
    # The pipe drains are bounded too: a descendant holding a pipe open must
    # not hang the smoke after the launcher exits.
    $drained = [Threading.Tasks.Task]::WaitAll(@($stdout, $stderr), [Math]::Max(5000, (& $remaining)))
    $log = $(if ($stdout.IsCompleted) { $stdout.Result } else { '<stdout not drained>' }) + "`n" +
           $(if ($stderr.IsCompleted) { $stderr.Result } else { '<stderr not drained>' })
    Set-Content -LiteralPath (Join-Path $work 'cubit.log') -Value $log -Encoding utf8
    if ($timedOut) { throw "Cubit timed out after $TimeoutSeconds s" }
    if (-not $drained) { throw 'Cubit output pipes did not close' }
    $exit = $p.ExitCode

    $problems = @()
    if ($exit -ne 0) { $problems += "exit code $exit" }
    # Cubit may prefix a plugin's message with a run of '*' (seen on the
    # export verb); allow only that before the exact text.
    $hello = ([regex]::Matches($log, '(?m)^\**MYPLUGIN_HELLO smoke\r?$')).Count
    if ($hello -ne 2) { $problems += "hello lines $hello (expected 2)" }
    $expectedLine = "MYPLUGIN_EXPORT nodes=$ExpectedNodes blocks=$ExpectedBlocks elements=$ExpectedElements"
    if ($log -notmatch "(?m)^\**$([regex]::Escape($expectedLine))\r?$") { $problems += "missing '$expectedLine'" }
    $expectedFile = "nodes $ExpectedNodes`nblocks $ExpectedBlocks`nblock_elements $ExpectedElements`n"
    if (-not (Test-Path -LiteralPath $out -PathType Leaf)) {
        $problems += 'output file missing'
    } elseif (([IO.File]::ReadAllText($out) -replace "`r`n", "`n") -cne $expectedFile) {
        $problems += 'output file contents differ'
    }
    if ($problems) {
        Write-Host $log
        throw ("Smoke failed: " + ($problems -join '; '))
    }
    Write-Host "SMOKE OK exit=0 hello=2 $expectedLine"
    Remove-Item -LiteralPath $work -Recurse -Force
} catch {
    Write-Host "Work directory kept: $work"
    throw
} finally {
    if ($p) { $p.Dispose() }
}
