[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$SourceExe,
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[0-9a-fA-F]{40}$')]
    [string]$SourceSha,
    [string]$Destination = 'O:\Eqnedit64.exe',
    [string]$ManifestPath = 'O:\Eqnedit64.handtest.json',
    [string]$ReleaseManifestPath = 'O:\Eqnedit64.release.json',
    [string]$LastReleaseManifestPath = 'O:\Eqnedit64.last-release.json',
    [switch]$WhatIf
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if ($PSVersionTable.PSVersion.Major -lt 7) {
    throw 'sync_handtest_to_o.ps1 requires PowerShell 7 or newer.'
}
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw 'Git is required.'
}

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path
$headSha = (git -C $repoRoot rev-parse HEAD).Trim().ToLowerInvariant()
$normalizedSourceSha = $SourceSha.ToLowerInvariant()
if ($headSha -cne $normalizedSourceSha) {
    throw "Hand-test source is not HEAD: HEAD=$headSha requested=$normalizedSourceSha"
}
$trackedStatus = git -C $repoRoot status --porcelain --untracked-files=no
if ($LASTEXITCODE -ne 0 -or $trackedStatus) {
    throw 'Tracked hand-test source is dirty; commit it before updating O:.'
}

$resolvedSource = (Resolve-Path -LiteralPath $SourceExe).Path
$actualVersion = (Get-Item -LiteralPath $resolvedSource).VersionInfo.ProductVersion
if ([string]::IsNullOrWhiteSpace($actualVersion)) {
    throw 'Source EXE has no ProductVersion.'
}
$shortSha = (git -C $repoRoot rev-parse --short $normalizedSourceSha).Trim()
$stampPath = Join-Path $repoRoot 'tools\eqnedit64\src\build_stamp.h'
$stamp = Get-Content -LiteralPath $stampPath -Raw
if ($stamp -notmatch ('EQNEDIT64_BUILD_COMMIT\s+"' +
        [regex]::Escape($shortSha) + '"')) {
    throw "EXE build stamp header is not the hand-test commit $shortSha"
}
$exeText = [Text.Encoding]::ASCII.GetString(
    [IO.File]::ReadAllBytes($resolvedSource))
if (-not $exeText.Contains($shortSha, [StringComparison]::Ordinal)) {
    throw "EXE binary does not contain the hand-test build stamp $shortSha"
}

$signature = Get-AuthenticodeSignature -LiteralPath $resolvedSource
$signer = if ($signature.SignerCertificate) {
    $signature.SignerCertificate.Subject
} else {
    ''
}
if ($signature.Status -ne 'Valid' -or $signer -cne 'CN=ksugahar') {
    throw "Source EXE signature is invalid: $($signature.Status) / $signer"
}
$sourceHash = (Get-FileHash -Algorithm SHA256 `
    -LiteralPath $resolvedSource).Hash
$branch = (git -C $repoRoot branch --show-current).Trim()

$fullDestination = [IO.Path]::GetFullPath($Destination)
$fullManifest = [IO.Path]::GetFullPath($ManifestPath)
$fullReleaseManifest = [IO.Path]::GetFullPath($ReleaseManifestPath)
$fullLastReleaseManifest = [IO.Path]::GetFullPath($LastReleaseManifestPath)
$expectedNames = @{
    $fullDestination = 'Eqnedit64.exe'
    $fullManifest = 'Eqnedit64.handtest.json'
    $fullReleaseManifest = 'Eqnedit64.release.json'
    $fullLastReleaseManifest = 'Eqnedit64.last-release.json'
}
foreach ($entry in $expectedNames.GetEnumerator()) {
    if ([IO.Path]::GetFileName($entry.Key) -cne $entry.Value) {
        throw "Unexpected O: hand-test path: $($entry.Key)"
    }
}
$destinationDirectory = Split-Path -Parent $fullDestination
if (-not (Test-Path -LiteralPath $destinationDirectory -PathType Container)) {
    throw "Destination directory is unavailable: $destinationDirectory"
}
foreach ($path in @(
        $fullManifest, $fullReleaseManifest, $fullLastReleaseManifest)) {
    if ((Split-Path -Parent $path) -cne $destinationDirectory) {
        throw 'All O: hand-test files must be in the destination directory.'
    }
}

$previousRelease = $null
if (Test-Path -LiteralPath $fullReleaseManifest -PathType Leaf) {
    $previousRelease = Get-Content -LiteralPath $fullReleaseManifest -Raw |
        ConvertFrom-Json
} elseif (Test-Path -LiteralPath $fullLastReleaseManifest -PathType Leaf) {
    $previousRelease = Get-Content -LiteralPath $fullLastReleaseManifest -Raw |
        ConvertFrom-Json
}
$previousTag = if ($previousRelease -and $previousRelease.tag) {
    [string]$previousRelease.tag
} else {
    ''
}

$transactionId = '{0}.{1}' -f $PID, [Guid]::NewGuid().ToString('N')
$stageExe = Join-Path $destinationDirectory `
    ('.Eqnedit64.exe.{0}.new' -f $transactionId)
$stageManifest = Join-Path $destinationDirectory `
    ('.Eqnedit64.handtest.json.{0}.new' -f $transactionId)
$backupExe = Join-Path $destinationDirectory `
    ('.Eqnedit64.exe.{0}.backup' -f $transactionId)
$backupManifest = Join-Path $destinationDirectory `
    ('.Eqnedit64.handtest.json.{0}.backup' -f $transactionId)
$backupRelease = Join-Path $destinationDirectory `
    ('.Eqnedit64.release.json.{0}.backup' -f $transactionId)
$backupLastRelease = Join-Path $destinationDirectory `
    ('.Eqnedit64.last-release.json.{0}.backup' -f $transactionId)
$hadDestination = Test-Path -LiteralPath $fullDestination -PathType Leaf
$hadManifest = Test-Path -LiteralPath $fullManifest -PathType Leaf
$hadRelease = Test-Path -LiteralPath $fullReleaseManifest -PathType Leaf
$hadLastRelease = Test-Path -LiteralPath $fullLastReleaseManifest -PathType Leaf

try {
    Copy-Item -LiteralPath $resolvedSource -Destination $stageExe
    $stageHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $stageExe).Hash
    if ($stageHash -cne $sourceHash) {
        throw "Staged O: copy changed bytes: $stageHash != $sourceHash"
    }

    $manifest = [ordered]@{
        schema = 'eqnedit64.o-handtest.v1'
        state = 'handtest'
        version = $actualVersion
        source_sha = $normalizedSourceSha
        source_branch = $branch
        exe_sha256 = $sourceHash
        signer = $signer
        previous_release_tag = $previousTag
        prepared_at = (Get-Date).ToUniversalTime().ToString('o')
    }
    $manifestJson = $manifest | ConvertTo-Json
    [IO.File]::WriteAllText(
        $stageManifest, $manifestJson + "`n", [Text.UTF8Encoding]::new($false))

    if (-not $WhatIf) {
        if ($hadDestination) {
            [IO.File]::Copy($fullDestination, $backupExe, $false)
        }
        if ($hadManifest) {
            [IO.File]::Copy($fullManifest, $backupManifest, $false)
        }
        if ($hadRelease) {
            [IO.File]::Copy($fullReleaseManifest, $backupRelease, $false)
        }
        if ($hadLastRelease) {
            [IO.File]::Copy(
                $fullLastReleaseManifest, $backupLastRelease, $false)
        }

        try {
            if ($hadRelease) {
                [IO.File]::Copy(
                    $fullReleaseManifest, $fullLastReleaseManifest, $true)
                Remove-Item -LiteralPath $fullReleaseManifest -Force
            }
            [IO.File]::Move($stageExe, $fullDestination, $true)
            $stageExe = $null
            [IO.File]::Move($stageManifest, $fullManifest, $true)
            $stageManifest = $null

            $destinationHash = (Get-FileHash -Algorithm SHA256 `
                -LiteralPath $fullDestination).Hash
            $destinationSignature = Get-AuthenticodeSignature `
                -LiteralPath $fullDestination
            $destinationVersion = (Get-Item `
                -LiteralPath $fullDestination).VersionInfo.ProductVersion
            $recorded = Get-Content -LiteralPath $fullManifest -Raw |
                ConvertFrom-Json
            if ($destinationHash -cne $sourceHash -or
                $destinationVersion -cne $actualVersion -or
                $destinationSignature.Status -ne 'Valid' -or
                -not $destinationSignature.SignerCertificate -or
                $destinationSignature.SignerCertificate.Subject -cne
                    'CN=ksugahar' -or
                $recorded.schema -cne 'eqnedit64.o-handtest.v1' -or
                $recorded.state -cne 'handtest' -or
                $recorded.source_sha -cne $normalizedSourceSha -or
                $recorded.exe_sha256 -cne $sourceHash -or
                (Test-Path -LiteralPath $fullReleaseManifest -PathType Leaf)) {
                throw 'O: hand-test post-copy verification failed.'
            }
        } catch {
            $updateFailure = $_
            try {
                $restoreItems = @(
                    [pscustomobject]@{ Path = $fullDestination
                        Backup = $backupExe; Had = $hadDestination },
                    [pscustomobject]@{ Path = $fullManifest
                        Backup = $backupManifest; Had = $hadManifest },
                    [pscustomobject]@{ Path = $fullReleaseManifest
                        Backup = $backupRelease; Had = $hadRelease },
                    [pscustomobject]@{ Path = $fullLastReleaseManifest
                        Backup = $backupLastRelease; Had = $hadLastRelease }
                )
                foreach ($restore in $restoreItems) {
                    if ($restore.Had) {
                        [IO.File]::Copy(
                            $restore.Backup, $restore.Path, $true)
                    } elseif (Test-Path -LiteralPath $restore.Path `
                            -PathType Leaf) {
                        Remove-Item -LiteralPath $restore.Path -Force
                    }
                }
            } catch {
                throw ("O: hand-test update failed and rollback also failed. " +
                    "Update: $updateFailure Rollback: $_")
            }
            throw $updateFailure
        }
    }

    [pscustomobject]@{
        Updated = -not $WhatIf
        State = 'handtest'
        SourceSha = $normalizedSourceSha
        Destination = $fullDestination
        Manifest = $fullManifest
        Version = $actualVersion
        Sha256 = $sourceHash
        Signature = [string]$signature.Status
        Signer = $signer
        PreviousReleaseTag = $previousTag
    }
} finally {
    foreach ($temporary in @(
            $stageExe, $stageManifest, $backupExe, $backupManifest,
            $backupRelease, $backupLastRelease)) {
        if ($temporary -and (Test-Path -LiteralPath $temporary -PathType Leaf)) {
            Remove-Item -LiteralPath $temporary -Force
        }
    }
}
