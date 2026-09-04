[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^eqnedit64-v\d+\.\d+\.\d+$')]
    [string]$Tag,
    [Parameter(Mandatory = $true)]
    [string]$SourceExe,
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[0-9a-fA-F]{40}$')]
    [string]$SourceSha,
    [ValidatePattern('^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$')]
    [string]$Repository = 'ksugahar/Radia',
    [string]$Destination = 'O:\Eqnedit64.exe',
    [string]$ManifestPath = 'O:\Eqnedit64.release.json',
    [string]$HandTestManifestPath = 'O:\Eqnedit64.handtest.json',
    [switch]$WhatIf
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$dryRun = [bool]$WhatIf
# CmdletBinding's caller may set the automatic preference even though this
# script owns an explicit dry-run switch. Keep staging I/O active so the
# preflight can hash and verify its temporary copy.
$WhatIfPreference = $false

if ($PSVersionTable.PSVersion.Major -lt 7) {
    throw 'sync_to_o.ps1 requires PowerShell 7 or newer.'
}
foreach ($command in @('git', 'gh')) {
    if (-not (Get-Command $command -ErrorAction SilentlyContinue)) {
        throw "$command is required."
    }
}

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path
$headSha = (git -C $repoRoot rev-parse HEAD).Trim()
$originMainText = gh api "repos/$Repository/commits/main" --jq .sha
if ($LASTEXITCODE -ne 0 -or
    [string]::IsNullOrWhiteSpace([string]$originMainText)) {
    throw 'Could not query the pushed main commit through GitHub.'
}
$originMainSha = ([string]$originMainText).Trim().ToLowerInvariant()
$normalizedSourceSha = $SourceSha.ToLowerInvariant()
if ($headSha -cne $normalizedSourceSha -or
    $originMainSha -cne $normalizedSourceSha) {
    throw ("Release source is not the exact pushed main commit: " +
        "HEAD=$headSha origin/main=$originMainSha requested=$normalizedSourceSha")
}
$trackedStatus = git -C $repoRoot status --porcelain --untracked-files=no
if ($LASTEXITCODE -ne 0 -or $trackedStatus) {
    throw 'Tracked release source is dirty.'
}

$remoteTags = gh api "repos/$Repository/git/matching-refs/tags/$Tag" |
    ConvertFrom-Json
if ($LASTEXITCODE -ne 0) { throw "Could not query remote tag $Tag" }
if (@($remoteTags).Count -gt 0) {
    throw "Release tag already exists; O: must be prepared before tag push: $Tag"
}

$resolvedSource = (Resolve-Path -LiteralPath $SourceExe).Path
$expectedVersion = $Tag.Substring('eqnedit64-v'.Length)
$actualVersion = (Get-Item -LiteralPath $resolvedSource).VersionInfo.ProductVersion
if ($actualVersion -cne $expectedVersion) {
    throw "Tag/EXE version mismatch: $Tag != $actualVersion"
}

$shortSha = (git -C $repoRoot rev-parse --short $normalizedSourceSha).Trim()
$stampPath = Join-Path $repoRoot 'tools\eqnedit64\src\build_stamp.h'
$stamp = Get-Content -LiteralPath $stampPath -Raw
if ($stamp -notmatch ('EQNEDIT64_BUILD_COMMIT\s+"' +
        [regex]::Escape($shortSha) + '"')) {
    throw "EXE build stamp is not the pushed main commit $shortSha"
}
$exeText = [Text.Encoding]::ASCII.GetString(
    [IO.File]::ReadAllBytes($resolvedSource))
if (-not $exeText.Contains($shortSha, [StringComparison]::Ordinal)) {
    throw "EXE binary does not contain the pushed main build stamp $shortSha"
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

$fullDestination = [IO.Path]::GetFullPath($Destination)
$fullManifest = [IO.Path]::GetFullPath($ManifestPath)
$fullHandTestManifest = [IO.Path]::GetFullPath($HandTestManifestPath)
if ([IO.Path]::GetFileName($fullDestination) -cne 'Eqnedit64.exe') {
    throw "Destination must name Eqnedit64.exe: $fullDestination"
}
if ([IO.Path]::GetFileName($fullManifest) -cne 'Eqnedit64.release.json') {
    throw "Manifest must name Eqnedit64.release.json: $fullManifest"
}
if ([IO.Path]::GetFileName($fullHandTestManifest) -cne
        'Eqnedit64.handtest.json') {
    throw "Hand-test manifest must name Eqnedit64.handtest.json: $fullHandTestManifest"
}
$destinationDirectory = Split-Path -Parent $fullDestination
if (-not (Test-Path -LiteralPath $destinationDirectory -PathType Container)) {
    throw "Destination directory is unavailable: $destinationDirectory"
}
if ((Split-Path -Parent $fullManifest) -cne $destinationDirectory) {
    throw 'EXE and release manifest must be in the same directory.'
}
if ((Split-Path -Parent $fullHandTestManifest) -cne $destinationDirectory) {
    throw 'EXE and hand-test manifest must be in the same directory.'
}

$transactionId = '{0}.{1}' -f $PID, [Guid]::NewGuid().ToString('N')
$stageExe = Join-Path $destinationDirectory `
    ('.Eqnedit64.exe.{0}.new' -f $transactionId)
$stageManifest = Join-Path $destinationDirectory `
    ('.Eqnedit64.release.json.{0}.new' -f $transactionId)
$backupExe = Join-Path $destinationDirectory `
    ('.Eqnedit64.exe.{0}.backup' -f $transactionId)
$backupManifest = Join-Path $destinationDirectory `
    ('.Eqnedit64.release.json.{0}.backup' -f $transactionId)
$backupHandTestManifest = Join-Path $destinationDirectory `
    ('.Eqnedit64.handtest.json.{0}.backup' -f $transactionId)
$hadDestination = Test-Path -LiteralPath $fullDestination -PathType Leaf
$hadManifest = Test-Path -LiteralPath $fullManifest -PathType Leaf
$hadHandTestManifest = Test-Path -LiteralPath $fullHandTestManifest -PathType Leaf

try {
    Copy-Item -LiteralPath $resolvedSource -Destination $stageExe
    $stageHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $stageExe).Hash
    if ($stageHash -cne $sourceHash) {
        throw "Staged O: copy changed bytes: $stageHash != $sourceHash"
    }

    $manifest = [ordered]@{
        schema = 'eqnedit64.o-release.v1'
        tag = $Tag
        version = $expectedVersion
        source_sha = $normalizedSourceSha
        exe_sha256 = $sourceHash
        signer = $signer
        prepared_at = (Get-Date).ToUniversalTime().ToString('o')
    }
    $manifestJson = $manifest | ConvertTo-Json
    [IO.File]::WriteAllText(
        $stageManifest, $manifestJson + "`n", [Text.UTF8Encoding]::new($false))

    if (-not $dryRun) {
        if ($hadDestination) {
            [IO.File]::Copy($fullDestination, $backupExe, $false)
        }
        if ($hadManifest) {
            [IO.File]::Copy($fullManifest, $backupManifest, $false)
        }
        if ($hadHandTestManifest) {
            [IO.File]::Copy(
                $fullHandTestManifest, $backupHandTestManifest, $false)
        }

        try {
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
                $destinationVersion -cne $expectedVersion -or
                $destinationSignature.Status -ne 'Valid' -or
                -not $destinationSignature.SignerCertificate -or
                $destinationSignature.SignerCertificate.Subject -cne 'CN=ksugahar' -or
                $recorded.schema -cne 'eqnedit64.o-release.v1' -or
                $recorded.source_sha -cne $normalizedSourceSha -or
                $recorded.exe_sha256 -cne $sourceHash) {
                throw 'O: post-copy verification failed.'
            }
            if ($hadHandTestManifest) {
                Remove-Item -LiteralPath $fullHandTestManifest -Force
            }
        } catch {
            $updateFailure = $_
            try {
                if ($hadDestination) {
                    [IO.File]::Copy($backupExe, $fullDestination, $true)
                } elseif (Test-Path -LiteralPath $fullDestination -PathType Leaf) {
                    Remove-Item -LiteralPath $fullDestination -Force
                }
                if ($hadManifest) {
                    [IO.File]::Copy($backupManifest, $fullManifest, $true)
                } elseif (Test-Path -LiteralPath $fullManifest -PathType Leaf) {
                    Remove-Item -LiteralPath $fullManifest -Force
                }
                if ($hadHandTestManifest) {
                    [IO.File]::Copy(
                        $backupHandTestManifest, $fullHandTestManifest, $true)
                } elseif (Test-Path -LiteralPath $fullHandTestManifest `
                        -PathType Leaf) {
                    Remove-Item -LiteralPath $fullHandTestManifest -Force
                }
            } catch {
                throw ("O: update failed and rollback also failed. " +
                    "Update: $updateFailure Rollback: $_")
            }
            throw $updateFailure
        }
    }

    [pscustomobject]@{
        Updated = -not $dryRun
        Tag = $Tag
        SourceSha = $normalizedSourceSha
        Destination = $fullDestination
        Manifest = $fullManifest
        Version = $expectedVersion
        Sha256 = $sourceHash
        Signature = [string]$signature.Status
        Signer = $signer
    }
} finally {
    foreach ($temporary in @(
            $stageExe, $stageManifest, $backupExe, $backupManifest,
            $backupHandTestManifest)) {
        if ($temporary -and (Test-Path -LiteralPath $temporary -PathType Leaf)) {
            Remove-Item -LiteralPath $temporary -Force
        }
    }
}
