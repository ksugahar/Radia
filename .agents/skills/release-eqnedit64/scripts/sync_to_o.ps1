[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^eqnedit64-v\d+\.\d+\.\d+$')]
    [string]$Tag,
    [string]$Repository = 'ksugahar/Radia',
    [string]$Destination = 'O:\Eqnedit64.exe',
    [switch]$WhatIf
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if ($PSVersionTable.PSVersion.Major -lt 7) {
    throw 'sync_to_o.ps1 requires PowerShell 7 or newer.'
}
if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    throw 'GitHub CLI (gh) is required.'
}

$expectedVersion = $Tag.Substring('eqnedit64-v'.Length)
$releaseText = gh release view $Tag --repo $Repository `
    --json isDraft,isPrerelease,tagName,url,assets
if ($LASTEXITCODE -ne 0) {
    throw "GitHub Release was not found: $Repository $Tag"
}
$release = $releaseText | ConvertFrom-Json
if ($release.isDraft -or $release.isPrerelease -or $release.tagName -cne $Tag) {
    throw "GitHub Release is not a final exact-tag release: $Tag"
}
$assetNames = @()
foreach ($asset in $release.assets) {
    $assetNames += [string]$asset.name
}
foreach ($required in @('Eqnedit64.exe', 'SHA256SUMS.txt')) {
    if ($required -cnotin $assetNames) {
        throw "GitHub Release is missing $required"
    }
}

$downloadRoot = Join-Path 'C:\temp' (
    'eqnedit64-o-sync-{0}-{1}' -f $PID,
    [Guid]::NewGuid().ToString('N'))
$stage = $null
New-Item -ItemType Directory -Path $downloadRoot | Out-Null

try {
    gh release download $Tag --repo $Repository `
        --pattern Eqnedit64.exe --pattern SHA256SUMS.txt `
        --dir $downloadRoot
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to download release assets for $Tag"
    }

    $downloadedExe = Join-Path $downloadRoot 'Eqnedit64.exe'
    $checksums = Join-Path $downloadRoot 'SHA256SUMS.txt'
    $line = Get-Content -LiteralPath $checksums |
        Where-Object { $_ -match 'Eqnedit64\.exe\s*$' } |
        Select-Object -First 1
    if (-not $line -or
        $line -notmatch '^(?<hash>[0-9a-fA-F]{64})\s+\*?Eqnedit64\.exe\s*$') {
        throw 'Published SHA256SUMS.txt has no valid Eqnedit64.exe entry.'
    }
    $expectedHash = $Matches.hash.ToUpperInvariant()
    $downloadedHash = (Get-FileHash -Algorithm SHA256 `
        -LiteralPath $downloadedExe).Hash
    if ($downloadedHash -cne $expectedHash) {
        throw "Published checksum mismatch: $downloadedHash != $expectedHash"
    }

    $downloadedVersion = (Get-Item -LiteralPath $downloadedExe).VersionInfo.ProductVersion
    if ($downloadedVersion -cne $expectedVersion) {
        throw "Tag/EXE version mismatch: $Tag != $downloadedVersion"
    }
    $downloadedSignature = Get-AuthenticodeSignature -LiteralPath $downloadedExe
    $downloadedSigner = if ($downloadedSignature.SignerCertificate) {
        $downloadedSignature.SignerCertificate.Subject
    } else {
        ''
    }
    if ($downloadedSignature.Status -ne 'Valid' -or
        $downloadedSigner -cne 'CN=ksugahar') {
        throw ("Published EXE signature is invalid: {0} / {1}" -f
            $downloadedSignature.Status, $downloadedSigner)
    }

    $fullDestination = [IO.Path]::GetFullPath($Destination)
    if ([IO.Path]::GetFileName($fullDestination) -cne 'Eqnedit64.exe') {
        throw "Destination must name Eqnedit64.exe: $fullDestination"
    }
    $destinationDirectory = Split-Path -Parent $fullDestination
    if (-not (Test-Path -LiteralPath $destinationDirectory -PathType Container)) {
        throw "Destination directory is unavailable: $destinationDirectory"
    }

    $stage = Join-Path $destinationDirectory `
        ('.Eqnedit64.exe.{0}.{1}.new' -f $PID,
         [Guid]::NewGuid().ToString('N'))
    Copy-Item -LiteralPath $downloadedExe -Destination $stage
    $stageHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $stage).Hash
    if ($stageHash -cne $expectedHash) {
        throw "Staged O: copy changed bytes: $stageHash != $expectedHash"
    }

    if (-not $WhatIf) {
        [IO.File]::Move($stage, $fullDestination, $true)
        $stage = $null

        $destinationHash = (Get-FileHash -Algorithm SHA256 `
            -LiteralPath $fullDestination).Hash
        $destinationSignature = Get-AuthenticodeSignature `
            -LiteralPath $fullDestination
        $destinationVersion = (Get-Item `
            -LiteralPath $fullDestination).VersionInfo.ProductVersion
        if ($destinationHash -cne $expectedHash -or
            $destinationVersion -cne $expectedVersion -or
            $destinationSignature.Status -ne 'Valid' -or
            -not $destinationSignature.SignerCertificate -or
            $destinationSignature.SignerCertificate.Subject -cne 'CN=ksugahar') {
            throw 'O: post-copy verification failed.'
        }

        [pscustomobject]@{
            Updated = $true
            Tag = $Tag
            ReleaseUrl = $release.url
            Destination = $fullDestination
            Version = $destinationVersion
            Sha256 = $destinationHash
            Signature = [string]$destinationSignature.Status
            Signer = $destinationSignature.SignerCertificate.Subject
        }
    } else {
        [pscustomobject]@{
            Updated = $false
            Tag = $Tag
            ReleaseUrl = $release.url
            Destination = $fullDestination
            Version = $downloadedVersion
            Sha256 = $downloadedHash
            Signature = [string]$downloadedSignature.Status
            Signer = $downloadedSigner
        }
    }
} finally {
    if ($stage -and (Test-Path -LiteralPath $stage -PathType Leaf)) {
        Remove-Item -LiteralPath $stage -Force
    }
    $safePrefix = [IO.Path]::GetFullPath('C:\temp\eqnedit64-o-sync-')
    $resolvedDownloadRoot = [IO.Path]::GetFullPath($downloadRoot)
    if ($resolvedDownloadRoot.StartsWith(
            $safePrefix, [StringComparison]::OrdinalIgnoreCase) -and
        (Test-Path -LiteralPath $resolvedDownloadRoot -PathType Container)) {
        Remove-Item -LiteralPath $resolvedDownloadRoot -Recurse -Force
    }
}
