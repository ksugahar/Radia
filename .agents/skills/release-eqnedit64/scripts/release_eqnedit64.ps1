# One command for the mechanical half of an Eqnedit64 release: the steps from
# "main is green" to "the tag is pushed".  It does not replace the judgement
# steps before it - the hand test, the Fable review, the merge to main - which
# stay in SKILL.md.
#
# What it runs on LAB is COMPILATION AND SIGNING ONLY.  The signing key is
# non-exportable, so the executable that ships can only be produced here; the
# GitHub-hosted CI build is an unsigned test artifact.  It starts no GUI,
# renders no equation, and registers no font, so it stays clear of the failure
# that breaks the interactive session's fonts.  Every private-font suite
# belongs to the isolated CI session, never to this script.
#
#   pwsh -File .agents/skills/release-eqnedit64/scripts/release_eqnedit64.ps1 `
#     -Tag eqnedit64-v3.0.16
#
# -WhatIf performs every check and the build, reports what would happen, and
# changes neither O: nor the staging release nor the remote tag.

[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^eqnedit64-v\d+\.\d+\.\d+$')]
    [string]$Tag,
    # Reuse an existing dist\Eqnedit64.exe instead of rebuilding.  The identity
    # checks below still run, so a stale binary is rejected, not trusted.
    [switch]$SkipBuild,
    [switch]$WhatIf
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if ($PSVersionTable.PSVersion.Major -lt 7) {
    throw 'release_eqnedit64.ps1 requires PowerShell 7 or newer.'
}
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw 'Git is required.'
}

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path
$eqnedit = Join-Path $repoRoot 'tools\eqnedit64'
$version = $Tag.Substring('eqnedit64-v'.Length)

function Step([string]$Message) { Write-Host "[release] $Message" }

# --- 1. the source must be exactly the pushed main commit ------------------
Step 'checking the source against origin/main'
git -C $repoRoot fetch origin main --quiet
if ($LASTEXITCODE -ne 0) { throw 'git fetch origin main failed.' }
$headSha = (git -C $repoRoot rev-parse HEAD).Trim()
$originMain = (git -C $repoRoot rev-parse origin/main).Trim()
if ($headSha -cne $originMain) {
    throw "HEAD is not the pushed main commit: HEAD=$headSha origin/main=$originMain"
}
if (git -C $repoRoot status --porcelain --untracked-files=no) {
    throw 'Tracked release source is dirty.'
}

# --- 2. the tag must not exist yet ----------------------------------------
git -C $repoRoot fetch origin --tags --quiet
$remoteTag = git -C $repoRoot ls-remote --tags origin "refs/tags/$Tag"
if ($LASTEXITCODE -ne 0) { throw "Could not query remote tag $Tag" }
if ($remoteTag) {
    throw @"
Release tag $Tag already exists on the remote. The staging assets and O: must
be prepared BEFORE the tag is pushed, so a tag that already exists means the
order was broken. Delete the remote tag and re-run, or release a new version.
"@
}

# --- 3. the tree must already declare this version ------------------------
$versionHeader = Join-Path $eqnedit 'src\eqnedit64_version.h'
$declared = Select-String -LiteralPath $versionHeader `
    -Pattern '#define EQNEDIT64_VERSION_TEXT "([0-9]+\.[0-9]+\.[0-9]+)"'
if (-not $declared) { throw "Could not read the version from $versionHeader" }
$declaredVersion = $declared.Matches[0].Groups[1].Value
if ($declaredVersion -cne $version) {
    throw ("The tree declares $declaredVersion but the tag says $version. " +
        'Bump every version site and push that commit to main first.')
}
Step "source $headSha declares $declaredVersion"

# --- 4. compile and sign (no GUI, no render, no font registration) --------
$dist = Join-Path $eqnedit 'dist\Eqnedit64.exe'
if ($SkipBuild) {
    Step 'skipping the build on request; verifying the existing dist binary'
} else {
    Step 'compiling and signing on LAB'
    & cmd.exe /d /c (Join-Path $eqnedit 'build\build_eqnedt64.bat')
    if ($LASTEXITCODE -ne 0) {
        throw "build_eqnedt64.bat failed with exit code $LASTEXITCODE."
    }
}
if (-not (Test-Path -LiteralPath $dist -PathType Leaf)) {
    throw "The portable release executable is missing: $dist"
}

# --- 5. the binary must be this version, signed by the release identity ----
$item = Get-Item -LiteralPath $dist
$signature = Get-AuthenticodeSignature -LiteralPath $dist
$sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $dist).Hash
$productVersion = $item.VersionInfo.ProductVersion
$signer = if ($signature.SignerCertificate) {
    $signature.SignerCertificate.Subject
} else { '<unsigned>' }
if ($productVersion -cne $version) {
    throw "The built executable reports $productVersion, not $version."
}
if ($signature.Status -ne 'Valid' -or $signer -cne 'CN=ksugahar') {
    throw "The built executable is not validly signed: status=$($signature.Status) signer=$signer"
}
Step "verified $version  sha256=$sha256  signer=$signer"

# --- 6. keep whatever O: held, so a hand-test candidate is not lost -------
$backup = $null
$oneDrive = 'C:\Users\Administrator\OneDrive'
if (Test-Path -LiteralPath $oneDrive) {
    $existing = Get-ChildItem -LiteralPath $oneDrive -Filter 'Eqnedit64*' `
        -File -ErrorAction SilentlyContinue
    if ($existing) {
        $backup = Join-Path 'C:\temp' (
            'eqnedit64-o-backup-{0}' -f (Get-Date -Format 'yyyyMMdd-HHmmss'))
        if ($WhatIf) {
            Step "would back up $($existing.Count) O: file(s) to $backup"
        } else {
            New-Item -ItemType Directory -Path $backup -Force | Out-Null
            $existing | Copy-Item -Destination $backup -Force
            Step "backed up $($existing.Count) O: file(s) to $backup"
        }
    }
}

# --- 7. stage to O: and to the release CI reads ---------------------------
Step 'staging to O: and the eqnedit64-staging release'
$syncArgs = @(
    '-Tag', $Tag,
    '-SourceExe', $dist,
    '-SourceSha', $headSha
)
if ($WhatIf) { $syncArgs += '-WhatIf' }
# sync_to_o.ps1 throws on every failure it detects, and $ErrorActionPreference
# is Stop, so a returned object means it succeeded.  Do not read $LASTEXITCODE
# here: it belongs to whatever native command that script happened to run last.
$sync = & (Join-Path $PSScriptRoot 'sync_to_o.ps1') @syncArgs
$sync | Format-List | Out-String | Write-Host

# --- 8. the tag goes last -------------------------------------------------
if ($WhatIf) {
    Step "would create and push the annotated tag $Tag at $headSha"
} else {
    Step "creating and pushing $Tag"
    git -C $repoRoot tag -a $Tag $headSha -m "Eqnedit64 $version"
    if ($LASTEXITCODE -ne 0) { throw "Could not create the tag $Tag" }
    git -C $repoRoot push origin $Tag
    if ($LASTEXITCODE -ne 0) {
        throw ("The tag was created locally but not pushed. " +
            "Push it with: git push origin $Tag")
    }
}

[pscustomobject]@{
    Released = -not $WhatIf
    Tag = $Tag
    Version = $version
    SourceSha = $headSha
    Executable = $dist
    Sha256 = $sha256
    Signer = $signer
    OBackup = $backup
    NextStep = 'Watch the tag CI, then the GitHub Release and the PyPI wheels.'
}
