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

$cacheRoot = 'C:\temp\eqnedit64-actions-runner-cache'
New-Item -ItemType Directory -Path $cacheRoot -Force | Out-Null
$cachedArchive = Join-Path $cacheRoot $assetName
if (-not (Test-Path -LiteralPath $cachedArchive -PathType Leaf)) {
    $stageArchive = Join-Path $cacheRoot (
        $assetName + '.' + [Guid]::NewGuid().ToString('N') + '.new')
    try {
        Invoke-WebRequest -Uri $asset.browser_download_url -OutFile $stageArchive
        if ($asset.digest -and ([string]$asset.digest).StartsWith('sha256:')) {
            $expected = ([string]$asset.digest).Substring(7).ToUpperInvariant()
            $actual = (Get-FileHash -Algorithm SHA256 `
                -LiteralPath $stageArchive).Hash
            if ($actual -cne $expected) {
                throw "Actions runner archive checksum mismatch: $actual != $expected"
            }
        }
        [IO.File]::Move($stageArchive, $cachedArchive, $false)
        $stageArchive = $null
    } finally {
        if ($stageArchive -and
            (Test-Path -LiteralPath $stageArchive -PathType Leaf)) {
            Remove-Item -LiteralPath $stageArchive -Force
        }
    }
}
if ($asset.digest -and ([string]$asset.digest).StartsWith('sha256:')) {
    $expected = ([string]$asset.digest).Substring(7).ToUpperInvariant()
    $actual = (Get-FileHash -Algorithm SHA256 `
        -LiteralPath $cachedArchive).Hash
    if ($actual -cne $expected) {
        throw "Cached Actions runner archive checksum mismatch: $actual != $expected"
    }
}

$runnerRoot = Join-Path 'C:\temp' (
    "eqnedit64-actions-runner-$ReleaseId-" + [Guid]::NewGuid().ToString('N'))
if (-not ([IO.Path]::GetFullPath($runnerRoot)).StartsWith(
        'C:\temp\eqnedit64-actions-runner-',
        [StringComparison]::OrdinalIgnoreCase)) {
    throw "Unsafe JIT runner directory: $runnerRoot"
}
New-Item -ItemType Directory -Path $runnerRoot | Out-Null
try {
    [IO.Compression.ZipFile]::ExtractToDirectory($cachedArchive, $runnerRoot)
} catch {
    Remove-Item -LiteralPath $runnerRoot -Recurse -Force
    throw
}
$listener = Join-Path $runnerRoot 'bin\Runner.Listener.exe'
if (-not (Test-Path -LiteralPath $listener -PathType Leaf)) {
    Remove-Item -LiteralPath $runnerRoot -Recurse -Force
    throw 'The downloaded Actions runner archive is incomplete.'
}

$runnerName = "LAB-eqnedit64-jit-$ReleaseId"
$requestPath = Join-Path 'C:\temp' (
    'eqnedit64-jit-request-' + [Guid]::NewGuid().ToString('N') + '.json')
$response = $null
$encodedConfig = $null
$registeredRunnerId = $null
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
    if ($response -and
        $response.PSObject.Properties['encoded_jit_config']) {
        $encodedConfig = [string]$response.encoded_jit_config
    }
    if ($response -and $response.PSObject.Properties['runner'] -and
        $response.runner.PSObject.Properties['id']) {
        $registeredRunnerId = [int64]$response.runner.id
    }
    if ($LASTEXITCODE -ne 0 -or -not $encodedConfig -or
        -not $registeredRunnerId) {
        throw 'GitHub did not return a JIT runner configuration.'
    }
} catch {
    if ($registeredRunnerId) {
        gh api --method DELETE `
            "repos/$Repository/actions/runners/$registeredRunnerId" `
            1>$null 2>$null
    }
    Remove-Item -LiteralPath $runnerRoot -Recurse -Force
    throw
} finally {
    if (Test-Path -LiteralPath $requestPath -PathType Leaf) {
        Remove-Item -LiteralPath $requestPath -Force
    }
}

$stdout = Join-Path 'C:\temp' "$runnerName.stdout.log"
$stderr = Join-Path 'C:\temp' "$runnerName.stderr.log"
$process = $null
try {
    $process = Start-Process -FilePath $listener `
        -ArgumentList @('run', '--jitconfig', $encodedConfig) `
        -WorkingDirectory $runnerRoot -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr

    $onlineDeadline = (Get-Date).ToUniversalTime().AddSeconds(30)
    $runner = $null
    while ((Get-Date).ToUniversalTime() -lt $onlineDeadline) {
        if ($process.HasExited) {
            throw "JIT runner exited before registration became online: $stderr"
        }
        $runnerText = gh api `
            "repos/$Repository/actions/runners/$registeredRunnerId" `
            2>$null
        if ($LASTEXITCODE -eq 0 -and $runnerText) {
            $runner = $runnerText | ConvertFrom-Json
            if ($runner.status -ceq 'online') { break }
        }
        Start-Sleep -Seconds 2
    }
    if (-not $runner -or $runner.status -cne 'online') {
        throw 'JIT runner did not become online within 30 seconds.'
    }
} catch {
    if ($process -and -not $process.HasExited) {
        $process | Stop-Process -Force
        $null = $process.WaitForExit(5000)
    }
    gh api --method DELETE `
        "repos/$Repository/actions/runners/$registeredRunnerId" `
        1>$null 2>$null
    Remove-Item -LiteralPath $runnerRoot -Recurse -Force
    throw
}

[pscustomobject]@{
    Process = $process
    ProcessId = $process.Id
    RunnerId = $registeredRunnerId
    RunnerName = $runnerName
    RunnerVersion = $runnerVersion
    RunnerRoot = $runnerRoot
    Stdout = $stdout
    Stderr = $stderr
}
