[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$')]
    [string]$Repository,
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[A-Za-z0-9_.-]+$')]
    [string]$ReleaseId
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
Add-Type -AssemblyName System.IO.Compression.FileSystem

if ($PSVersionTable.PSVersion.Major -lt 7) {
    throw 'start_jit_runner.ps1 requires PowerShell 7 or newer.'
}
if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    throw 'GitHub CLI is required.'
}

$runnerRelease = gh api repos/actions/runner/releases/latest |
    ConvertFrom-Json
if ($LASTEXITCODE -ne 0 -or -not $runnerRelease.tag_name) {
    throw 'Could not query the current GitHub Actions runner release.'
}
$runnerVersion = ([string]$runnerRelease.tag_name).TrimStart('v')
$assetName = "actions-runner-win-x64-$runnerVersion.zip"
$asset = $runnerRelease.assets |
    Where-Object { $_.name -ceq $assetName } |
    Select-Object -First 1
if (-not $asset) {
    throw "The runner release has no Windows x64 asset: $assetName"
}

$cacheRoot = Join-Path 'C:\temp\eqnedit64-actions-runner-cache' $runnerVersion
$listener = Join-Path $cacheRoot 'bin\Runner.Listener.exe'
if (-not (Test-Path -LiteralPath $listener -PathType Leaf)) {
    $archive = Join-Path 'C:\temp' $assetName
    $stage = Join-Path 'C:\temp' (
        'eqnedit64-actions-runner-stage-' + [Guid]::NewGuid().ToString('N'))
    try {
        Invoke-WebRequest -Uri $asset.browser_download_url -OutFile $archive
        if ($asset.digest -and ([string]$asset.digest).StartsWith('sha256:')) {
            $expected = ([string]$asset.digest).Substring(7).ToUpperInvariant()
            $actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $archive).Hash
            if ($actual -cne $expected) {
                throw "Actions runner archive checksum mismatch: $actual != $expected"
            }
        }
        New-Item -ItemType Directory -Path $stage | Out-Null
        [IO.Compression.ZipFile]::ExtractToDirectory($archive, $stage)
        if (-not (Test-Path -LiteralPath (
                    Join-Path $stage 'bin\Runner.Listener.exe') -PathType Leaf)) {
            throw 'The downloaded Actions runner archive is incomplete.'
        }
        $cacheParent = Split-Path -Parent $cacheRoot
        New-Item -ItemType Directory -Path $cacheParent -Force | Out-Null
        if (Test-Path -LiteralPath $cacheRoot) {
            throw "Incomplete runner cache already exists: $cacheRoot"
        }
        [IO.Directory]::Move($stage, $cacheRoot)
        $stage = $null
    } finally {
        if (Test-Path -LiteralPath $archive -PathType Leaf) {
            Remove-Item -LiteralPath $archive -Force
        }
        if ($stage -and (Test-Path -LiteralPath $stage -PathType Container)) {
            Remove-Item -LiteralPath $stage -Recurse -Force
        }
    }
}

$runnerName = "LAB-eqnedit64-jit-$ReleaseId"
$requestPath = Join-Path 'C:\temp' (
    'eqnedit64-jit-request-' + [Guid]::NewGuid().ToString('N') + '.json')
try {
    $request = [ordered]@{
        name = $runnerName
        runner_group_id = 1
        labels = @('self-hosted', 'Windows', 'X64', 'windows-radia')
        work_folder = '_work'
    } | ConvertTo-Json
    [IO.File]::WriteAllText(
        $requestPath, $request, [Text.UTF8Encoding]::new($false))
    $response = gh api --method POST `
        "repos/$Repository/actions/runners/generate-jitconfig" `
        --input $requestPath | ConvertFrom-Json
    if ($LASTEXITCODE -ne 0 -or -not $response.encoded_jit_config -or
        -not $response.runner.id) {
        throw 'GitHub did not return a JIT runner configuration.'
    }
} finally {
    if (Test-Path -LiteralPath $requestPath -PathType Leaf) {
        Remove-Item -LiteralPath $requestPath -Force
    }
}

$stdout = Join-Path 'C:\temp' "$runnerName.stdout.log"
$stderr = Join-Path 'C:\temp' "$runnerName.stderr.log"
$process = Start-Process -FilePath $listener `
    -ArgumentList @('run', '--jitconfig', $response.encoded_jit_config) `
    -WorkingDirectory $cacheRoot -WindowStyle Hidden -PassThru `
    -RedirectStandardOutput $stdout -RedirectStandardError $stderr

[pscustomobject]@{
    ProcessId = $process.Id
    RunnerId = [int64]$response.runner.id
    RunnerName = $runnerName
    RunnerVersion = $runnerVersion
    Stdout = $stdout
    Stderr = $stderr
}
