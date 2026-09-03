[CmdletBinding(SupportsShouldProcess = $true, ConfirmImpact = 'Medium')]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^\d+\.\d+\.\d+$')]
    [string]$Version,
    [ValidatePattern('^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$')]
    [string]$Repository = 'ksugahar/Radia',
    [ValidateRange(10, 90)]
    [int]$TimeoutMinutes = 45
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if ($PSVersionTable.PSVersion.Major -lt 7) {
    throw 'publish.ps1 requires PowerShell 7 or newer.'
}
foreach ($command in @('git', 'gh', 'python', 'cmd.exe')) {
    if (-not (Get-Command $command -ErrorAction SilentlyContinue)) {
        throw "$command is required."
    }
}
gh auth status --hostname github.com 1>$null 2>$null
if ($LASTEXITCODE -ne 0) { throw 'GitHub CLI is not authenticated.' }

$env:GIT_CONFIG_COUNT = '1'
$env:GIT_CONFIG_KEY_0 = 'safe.directory'
$env:GIT_CONFIG_VALUE_0 = '*'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path
$tag = "eqnedit64-v$Version"
$deadline = (Get-Date).ToUniversalTime().AddMinutes($TimeoutMinutes)

function Invoke-Git {
    param([Parameter(ValueFromRemainingArguments = $true)][string[]]$Arguments)
    $output = git -C $repoRoot @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "git failed: $($Arguments -join ' ')"
    }
    return $output
}

function Get-GhJson {
    param([Parameter(Mandatory = $true)][string[]]$Arguments)
    $text = gh @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "gh failed: $($Arguments -join ' ')"
    }
    if (-not $text) { return $null }
    return $text | ConvertFrom-Json
}

function Wait-Until {
    param(
        [Parameter(Mandatory = $true)][scriptblock]$Probe,
        [Parameter(Mandatory = $true)][string]$Description
    )
    while ((Get-Date).ToUniversalTime() -lt $deadline) {
        $result = & $Probe
        if ($null -ne $result) { return $result }
        Start-Sleep -Seconds 15
    }
    throw "Timed out waiting for $Description."
}

function Assert-MainChecks {
    param([Parameter(Mandatory = $true)][string]$Sha)
    $runs = Get-GhJson @('run', 'list', '--repo', $Repository,
        '--commit', $Sha, '--event', 'push', '--limit', '30', '--json',
        'databaseId,name,status,conclusion,headSha,createdAt')
    foreach ($name in @('Eqnedit64', 'Policy Lint')) {
        $run = $runs | Where-Object {
            $_.name -ceq $name -and $_.headSha -ceq $Sha
        } | Sort-Object createdAt -Descending | Select-Object -First 1
        if (-not $run -or $run.status -cne 'completed' -or
            $run.conclusion -cne 'success') {
            throw "Required main check is not green for ${Sha}: $name"
        }
    }
}

function Assert-VersionFiles {
    $expectations = [ordered]@{
        'tools\eqnedit64\CMakeLists.txt' =
            "project(Eqnedit64 VERSION $Version LANGUAGES CXX)"
        'tools\eqnedit64\src\eqnedit64_version.h' =
            "#define EQNEDIT64_VERSION_TEXT `"$Version`""
        'packages\eqnedit64\pyproject.toml' = "version = `"$Version`""
        'packages\eqnedit64\src\eqnedit64\__init__.py' =
            "__version__ = `"$Version`""
        'packages\eqnedit64\tests\verify_installed.py' =
            "assert eqnedit64.__version__ == `"$Version`""
    }
    foreach ($entry in $expectations.GetEnumerator()) {
        $path = Join-Path $repoRoot $entry.Key
        if (-not (Test-Path -LiteralPath $path -PathType Leaf) -or
            -not ([IO.File]::ReadAllText($path).Contains(
                    $entry.Value, [StringComparison]::Ordinal))) {
            throw "Release version is not synchronized in $($entry.Key)"
        }
    }
    $parts = $Version.Split('.')
    $header = [IO.File]::ReadAllText((Join-Path $repoRoot `
        'tools\eqnedit64\src\eqnedit64_version.h'))
    $nativeMacros = @(
        "#define EQNEDIT64_VERSION_MAJOR $($parts[0])",
        "#define EQNEDIT64_VERSION_MINOR $($parts[1])",
        "#define EQNEDIT64_VERSION_PATCH $($parts[2])",
        "#define EQNEDIT64_VERSION_TUPLE $($parts[0]),$($parts[1]),$($parts[2]),0",
        "#define EQNEDIT64_VERSION_TEXT_W L`"$Version`""
    )
    foreach ($macro in $nativeMacros) {
        if (-not $header.Contains($macro, [StringComparison]::Ordinal)) {
            throw "The native version macro is not synchronized: $macro"
        }
    }
}

function Get-RemoteTagSha {
    $refs = Get-GhJson @('api',
        "repos/$Repository/git/matching-refs/tags/$tag")
    if (@($refs).Count -eq 0) { return $null }
    if (@($refs).Count -ne 1 -or $refs[0].ref -cne "refs/tags/$tag") {
        throw "Remote tag lookup was ambiguous: $tag"
    }
    $commit = Get-GhJson @('api', "repos/$Repository/commits/$tag")
    return ([string]$commit.sha).ToLowerInvariant()
}

function Publish-AnnotatedTag {
    param([Parameter(Mandatory = $true)][string]$Sha)
    $name = ([string](Invoke-Git config user.name)).Trim()
    $email = ([string](Invoke-Git config user.email)).Trim()
    if (-not $name -or -not $email) {
        throw 'Git user.name and user.email are required for the release tag.'
    }
    $tagRequestPath = Join-Path 'C:\temp' (
        'eqnedit64-tag-' + [Guid]::NewGuid().ToString('N') + '.json')
    $refRequestPath = Join-Path 'C:\temp' (
        'eqnedit64-ref-' + [Guid]::NewGuid().ToString('N') + '.json')
    try {
        $tagRequest = [ordered]@{
            tag = $tag
            message = "Eqnedit64 $Version"
            object = $Sha
            type = 'commit'
            tagger = [ordered]@{
                name = $name
                email = $email
                date = (Get-Date).ToUniversalTime().ToString('o')
            }
        } | ConvertTo-Json -Depth 3
        [IO.File]::WriteAllText(
            $tagRequestPath, $tagRequest, [Text.UTF8Encoding]::new($false))
        $tagObject = Get-GhJson @('api', '--method', 'POST',
            "repos/$Repository/git/tags", '--input', $tagRequestPath)
        if (-not $tagObject.sha) {
            throw 'GitHub did not create the annotated tag object.'
        }
        $refRequest = [ordered]@{
            ref = "refs/tags/$tag"
            sha = [string]$tagObject.sha
        } | ConvertTo-Json
        [IO.File]::WriteAllText(
            $refRequestPath, $refRequest, [Text.UTF8Encoding]::new($false))
        $null = Get-GhJson @('api', '--method', 'POST',
            "repos/$Repository/git/refs", '--input', $refRequestPath)
    } finally {
        foreach ($path in @($tagRequestPath, $refRequestPath)) {
            if (Test-Path -LiteralPath $path -PathType Leaf) {
                Remove-Item -LiteralPath $path -Force
            }
        }
    }
}

function Assert-ORelease {
    param([Parameter(Mandatory = $true)][string]$Sha)
    $exe = 'O:\Eqnedit64.exe'
    $manifestPath = 'O:\Eqnedit64.release.json'
    if (-not (Test-Path -LiteralPath $exe -PathType Leaf) -or
        -not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) {
        throw 'O: release EXE or manifest is missing.'
    }
    $manifest = Get-Content -LiteralPath $manifestPath -Raw |
        ConvertFrom-Json
    $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $exe).Hash
    $signature = Get-AuthenticodeSignature -LiteralPath $exe
    $productVersion = (Get-Item -LiteralPath $exe).VersionInfo.ProductVersion
    if ($manifest.schema -cne 'eqnedit64.o-release.v1' -or
        $manifest.tag -cne $tag -or $manifest.version -cne $Version -or
        $manifest.source_sha -cne $Sha -or
        $manifest.exe_sha256 -cne $hash -or
        $productVersion -cne $Version -or
        $signature.Status -ne 'Valid' -or
        -not $signature.SignerCertificate -or
        $signature.SignerCertificate.Subject -cne 'CN=ksugahar') {
        throw 'O: release gate is not the exact signed release.'
    }
    return [pscustomobject]@{
        Exe = $exe
        Manifest = $manifestPath
        Hash = $hash
        Signature = [string]$signature.Status
        Signer = $signature.SignerCertificate.Subject
    }
}

function Wait-ForRun {
    param(
        [Parameter(Mandatory = $true)][string]$Workflow,
        [Parameter(Mandatory = $true)][string]$Sha,
        [Parameter(Mandatory = $true)][datetime]$NotBefore,
        [string]$Branch
    )
    return Wait-Until -Description "$Workflow run" -Probe {
        $arguments = @('run', 'list', '--repo', $Repository, '--workflow',
            $Workflow, '--limit', '20', '--json',
            'databaseId,status,conclusion,headSha,createdAt,url')
        if ($Branch) { $arguments += @('--branch', $Branch) }
        $runs = Get-GhJson $arguments
        $run = $runs | Where-Object {
            $_.headSha -ceq $Sha -and
            ([datetime]$_.createdAt).ToUniversalTime() -ge $NotBefore
        } | Sort-Object createdAt | Select-Object -First 1
        if ($run) { return $run }
        return $null
    }
}

function Watch-Run {
    param([Parameter(Mandatory = $true)]$Run)
    while ((Get-Date).ToUniversalTime() -lt $deadline) {
        $state = Get-GhJson @('run', 'view', ([string]$Run.databaseId),
            '--repo', $Repository, '--json', 'status,conclusion,url')
        if ($state.status -ceq 'completed') {
            if ($state.conclusion -cne 'success') {
                gh run view ([string]$Run.databaseId) --repo $Repository `
                    --log-failed
                throw "GitHub Actions run failed: $($state.url)"
            }
            return $state
        }
        Start-Sleep -Seconds 15
    }
    throw "Timed out waiting for GitHub Actions run: $($Run.url)"
}

function Test-PublicReleaseExists {
    try {
        gh release view $tag --repo $Repository 1>$null 2>$null
        if ($LASTEXITCODE -ne 0) { return $false }
        $null = Invoke-RestMethod -Uri "https://pypi.org/pypi/eqnedit64/$Version/json"
        return $true
    } catch {
        return $false
    }
}

function Verify-PublicRelease {
    param([Parameter(Mandatory = $true)]$ORelease)
    $verifyRoot = Join-Path 'C:\temp' (
        "eqnedit64-$Version-public-verify-" +
        (Get-Date).ToUniversalTime().ToString('yyyyMMddTHHmmssZ'))
    New-Item -ItemType Directory -Path $verifyRoot | Out-Null
    gh release download $tag --repo $Repository --pattern 'Eqnedit64.exe' `
        --pattern 'SHA256SUMS.txt' --dir $verifyRoot
    if ($LASTEXITCODE -ne 0) { throw 'GitHub Release download failed.' }
    $releaseExe = Join-Path $verifyRoot 'Eqnedit64.exe'
    $releaseHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $releaseExe).Hash
    if ($releaseHash -cne $ORelease.Hash) {
        throw 'GitHub Release EXE differs from O:.'
    }
    $checksumLine = Get-Content -LiteralPath (
        Join-Path $verifyRoot 'SHA256SUMS.txt') |
        Where-Object { $_ -match 'Eqnedit64\.exe$' } |
        Select-Object -First 1
    if (-not $checksumLine -or
        (($checksumLine -split '\s+')[0]).ToUpperInvariant() -cne
            $ORelease.Hash) {
        throw 'GitHub Release checksum does not describe the O: EXE.'
    }

    $pypi = Wait-Until -Description "PyPI eqnedit64 $Version" -Probe {
        try {
            return Invoke-RestMethod -Uri `
                "https://pypi.org/pypi/eqnedit64/$Version/json"
        } catch {
            return $null
        }
    }
    $wheels = @($pypi.urls | Where-Object { $_.packagetype -ceq 'bdist_wheel' })
    if ($wheels.Count -ne 4) {
        throw "Expected four PyPI wheels, found $($wheels.Count)."
    }
    foreach ($pythonTag in @('cp310', 'cp311', 'cp312', 'cp313')) {
        if (-not ($wheels.filename -match $pythonTag)) {
            throw "PyPI wheel is missing: $pythonTag"
        }
    }

    Add-Type -AssemblyName System.IO.Compression
    foreach ($wheel in $wheels) {
        $wheelPath = Join-Path $verifyRoot $wheel.filename
        Invoke-WebRequest -Uri $wheel.url -OutFile $wheelPath
        $wheelHash = (Get-FileHash -Algorithm SHA256 `
            -LiteralPath $wheelPath).Hash
        if ($wheelHash -cne ([string]$wheel.digests.sha256).ToUpperInvariant()) {
            throw "PyPI wheel checksum mismatch: $($wheel.filename)"
        }
        $archive = [IO.Compression.ZipFile]::OpenRead($wheelPath)
        try {
            $entry = $archive.Entries | Where-Object {
                $_.FullName -ceq 'eqnedit64/Eqnedit64.exe'
            } | Select-Object -First 1
            if (-not $entry) {
                throw "Wheel has no bundled EXE: $($wheel.filename)"
            }
            $embedded = Join-Path $verifyRoot `
                ($wheel.filename + '.Eqnedit64.exe')
            $input = $entry.Open()
            $output = [IO.File]::Create($embedded)
            try { $input.CopyTo($output) } finally {
                $output.Dispose()
                $input.Dispose()
            }
            $embeddedHash = (Get-FileHash -Algorithm SHA256 `
                -LiteralPath $embedded).Hash
            if ($embeddedHash -cne $ORelease.Hash) {
                throw "Wheel EXE differs from O:: $($wheel.filename)"
            }
        } finally {
            $archive.Dispose()
        }
    }
    return [pscustomobject]@{
        GitHubRelease = "https://github.com/$Repository/releases/tag/$tag"
        PyPI = "https://pypi.org/project/eqnedit64/$Version/"
        VerificationDirectory = $verifyRoot
        WheelCount = $wheels.Count
    }
}

Assert-VersionFiles
$sha = ([string](Invoke-Git rev-parse HEAD)).Trim().ToLowerInvariant()
$originMain = ([string](gh api "repos/$Repository/commits/main" `
    --jq .sha)).Trim().ToLowerInvariant()
if ($LASTEXITCODE -ne 0 -or -not $originMain) {
    throw 'Could not query the pushed main commit through GitHub.'
}
if ($sha -cne $originMain) {
    throw "Release source must equal origin/main: HEAD=$sha origin/main=$originMain"
}
$status = Invoke-Git status --porcelain
if ($status) { throw 'Release source is dirty.' }
Assert-MainChecks -Sha $sha

$remoteTagSha = Get-RemoteTagSha
if ($remoteTagSha -and $remoteTagSha -cne $sha) {
    throw "Remote tag points elsewhere: $tag -> $remoteTagSha, expected $sha"
}

if (-not $remoteTagSha) {
    Push-Location $repoRoot
    try {
        cmd.exe /d /c tools\eqnedit64\build\build_eqnedt64.bat
        if ($LASTEXITCODE -ne 0) { throw 'Eqnedit64 build failed.' }
        pwsh -NoProfile -File tools\eqnedit64\tests\test_cli_output.ps1 `
            -AppPath tools\eqnedit64\dist\Eqnedit64.exe
        if ($LASTEXITCODE -ne 0) { throw 'Eqnedit64 CLI verification failed.' }
    } finally {
        Pop-Location
    }
    if (-not $PSCmdlet.ShouldProcess(
            "$Repository $tag", 'Update O:, tag, and publish')) {
        & (Join-Path $PSScriptRoot 'sync_to_o.ps1') -Tag $tag `
            -SourceExe (Join-Path $repoRoot 'tools\eqnedit64\dist\Eqnedit64.exe') `
            -SourceSha $sha -WhatIf
        return
    }
    $sync = & (Join-Path $PSScriptRoot 'sync_to_o.ps1') -Tag $tag `
        -SourceExe (Join-Path $repoRoot 'tools\eqnedit64\dist\Eqnedit64.exe') `
        -SourceSha $sha -Repository $Repository
    $oRelease = Assert-ORelease -Sha $sha

    $tagPushTime = (Get-Date).ToUniversalTime().AddSeconds(-5)
    Publish-AnnotatedTag -Sha $sha
    $publishedTagSha = Get-RemoteTagSha
    if ($publishedTagSha -cne $sha) {
        throw "Published tag did not resolve to exact main: $publishedTagSha"
    }
} else {
    if (-not $PSCmdlet.ShouldProcess(
            "$Repository $tag", 'Resume and verify publication')) {
        return
    }
    $oRelease = Assert-ORelease -Sha $sha
    $tagPushTime = [datetime]'2000-01-01T00:00:00Z'
}

$jit = $null
try {
    if (-not (Test-PublicReleaseExists)) {
        $jit = & (Join-Path $PSScriptRoot 'start_jit_runner.ps1') `
            -Repository $Repository -ReleaseId ($tag -replace '[^A-Za-z0-9_.-]', '-')
        $tagRun = Wait-ForRun -Workflow 'Eqnedit64' -Sha $sha `
            -NotBefore $tagPushTime -Branch $tag
        $null = Watch-Run -Run $tagRun
        $releaseRun = Wait-ForRun -Workflow 'release-eqnedit64-pypi.yml' `
            -Sha $sha -NotBefore ([datetime]$tagRun.createdAt)
        $null = Watch-Run -Run $releaseRun
    }
    $public = Verify-PublicRelease -ORelease $oRelease
} finally {
    if ($jit) {
        $process = Get-Process -Id $jit.ProcessId -ErrorAction SilentlyContinue
        if ($process) {
            $process | Stop-Process -Force
            $null = $process.WaitForExit(5000)
        }
        gh api --method DELETE `
            "repos/$Repository/actions/runners/$($jit.RunnerId)" 1>$null 2>$null
        $runnerRoot = [IO.Path]::GetFullPath([string]$jit.RunnerRoot)
        if (-not $runnerRoot.StartsWith(
                'C:\temp\eqnedit64-actions-runner-',
                [StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to remove unsafe JIT runner directory: $runnerRoot"
        }
        if (Test-Path -LiteralPath $runnerRoot -PathType Container) {
            Remove-Item -LiteralPath $runnerRoot -Recurse -Force
        }
    }
}

[pscustomobject]@{
    Version = $Version
    Tag = $tag
    SourceSha = $sha
    Destination = $oRelease.Exe
    Manifest = $oRelease.Manifest
    Sha256 = $oRelease.Hash
    Signature = $oRelease.Signature
    Signer = $oRelease.Signer
    GitHubRelease = $public.GitHubRelease
    PyPI = $public.PyPI
    Wheels = $public.WheelCount
    VerificationDirectory = $public.VerificationDirectory
}
