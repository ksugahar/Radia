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
    [string]$Destination = 'O:\Eqnedit64.exe',
    [string]$ManifestPath = 'O:\Eqnedit64.release.json',
    [string]$HandTestManifestPath = 'O:\Eqnedit64.handtest.json',
    [switch]$WhatIf
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if ($PSVersionTable.PSVersion.Major -lt 7) {
    throw 'sync_to_o.ps1 requires PowerShell 7 or newer.'
}
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw 'Git is required.'
}

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path
git -C $repoRoot fetch origin main --quiet
if ($LASTEXITCODE -ne 0) { throw 'git fetch origin main failed.' }
$headSha = (git -C $repoRoot rev-parse HEAD).Trim()
$originMainSha = (git -C $repoRoot rev-parse origin/main).Trim()
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

$remoteTag = git -C $repoRoot ls-remote --tags origin "refs/tags/$Tag"
if ($LASTEXITCODE -ne 0) { throw "Could not query remote tag $Tag" }
if ($remoteTag) {
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

    if (-not $WhatIf) {
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

    # O: is the human hand-test entry point and cannot be read by CI: the
    # runner service is NETWORK SERVICE, which sees neither a per-user mapped
    # drive nor - these machines being a workgroup, not a domain - the
    # laboratory SMB share.  Stage the same two files as release assets so the
    # tag CI can fetch exactly what was just written to O:.
    $staged = $null
    if (-not $WhatIf) {
        $stagingDir = Join-Path ([IO.Path]::GetTempPath()) (
            'eqnedit64-staging-{0}' -f $transactionId)
        New-Item -ItemType Directory -Path $stagingDir -Force | Out-Null
        try {
            $stagedExe = Join-Path $stagingDir ("Eqnedit64-$expectedVersion.exe")
            $stagedManifest = Join-Path $stagingDir (
                "Eqnedit64-$expectedVersion.release.json")
            [IO.File]::Copy($fullDestination, $stagedExe, $true)
            [IO.File]::Copy($fullManifest, $stagedManifest, $true)
            $uploader = Join-Path $repoRoot 'tools\upload_release_asset.py'
            & python $uploader --repo ksugahar/Radia --tag eqnedit64-staging `
                --create-if-missing $stagedExe $stagedManifest
            if ($LASTEXITCODE -ne 0) {
                throw @"
O: was updated but the release assets could not be staged for CI.
Do not push the tag yet: the tag CI reads the signed executable from the
eqnedit64-staging release, not from O:.  Re-run this script, or upload
$stagedExe and $stagedManifest with tools/upload_release_asset.py.
"@
            }
            $staged = @("Eqnedit64-$expectedVersion.exe",
                        "Eqnedit64-$expectedVersion.release.json")
        } finally {
            if (Test-Path -LiteralPath $stagingDir) {
                Remove-Item -LiteralPath $stagingDir -Recurse -Force
            }
        }
    }

    [pscustomobject]@{
        Updated = -not $WhatIf
        Tag = $Tag
        SourceSha = $normalizedSourceSha
        Destination = $fullDestination
        Manifest = $fullManifest
        StagedAssets = $staged
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
