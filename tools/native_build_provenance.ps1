function Invoke-NativeProvenanceGit {
    param(
        [Parameter(Mandatory)][string]$RepoRoot,
        [Parameter(Mandatory)][string[]]$Arguments
    )

    Set-StrictMode -Version Latest

    $output = @(& git -C $RepoRoot @Arguments 2>&1)
    if ($LASTEXITCODE -ne 0) {
        throw "git $($Arguments -join ' ') failed in ${RepoRoot}: $($output -join [Environment]::NewLine)"
    }
    return ($output -join "`n")
}

function Get-NativeBuildSourceIdentity {
    param([Parameter(Mandatory)][string]$RepoRoot)

    Set-StrictMode -Version Latest
    # The full build may refresh this tracked distribution artifact. It is an
    # output, not an input to _radia_pybind or radia_mex, so exclude it from
    # the source-state race check.
    $sourcePathspec = @(
        "--", ".",
        ":(exclude)packages/cubit-mesh-export/src/cubit_mesh_export/cubit_mesh_export.ccm"
    )
    $commit = (Invoke-NativeProvenanceGit -RepoRoot $RepoRoot `
        -Arguments @("rev-parse", "HEAD")).Trim().ToLowerInvariant()
    if ($commit -notmatch '^[0-9a-f]{40}$') {
        throw "git rev-parse HEAD returned an invalid commit: $commit"
    }
    $status = Invoke-NativeProvenanceGit -RepoRoot $RepoRoot `
        -Arguments (@("status", "--porcelain=v1", "--untracked-files=all") + $sourcePathspec)
    $diff = Invoke-NativeProvenanceGit -RepoRoot $RepoRoot `
        -Arguments (@("diff", "--binary", "HEAD") + $sourcePathspec)
    $stateBytes = [Text.Encoding]::UTF8.GetBytes($status + "`0" + $diff)
    # This fingerprints Git's status plus tracked binary diff. For untracked
    # paths the status records names, not contents; such a build remains dirty
    # and is never eligible for an accepted manifest.
    $changeFingerprint = [Convert]::ToHexString(
        [Security.Cryptography.SHA256]::HashData($stateBytes)
    ).ToLowerInvariant()
    return [pscustomobject]@{
        source_commit = $commit
        source_dirty = [bool]($status.Length -gt 0)
        source_change_fingerprint_sha256 = $changeFingerprint
    }
}

function Clear-NativeBuildProvenance {
    param([Parameter(Mandatory)][string]$BinaryPath)

    Set-StrictMode -Version Latest
    $manifestPath = "$BinaryPath.build.json"
    if (Test-Path -LiteralPath $manifestPath) {
        Remove-Item -LiteralPath $manifestPath -Force
    }
}

function Write-NativeBuildProvenance {
    param(
        [Parameter(Mandatory)][string]$BinaryPath,
        [Parameter(Mandatory)][string]$RepoRoot,
        [Parameter(Mandatory)]$StartIdentity
    )

    Set-StrictMode -Version Latest
    $manifestPath = "$BinaryPath.build.json"
    $temporaryPath = "$manifestPath.tmp-$PID-$([Guid]::NewGuid().ToString('N'))"
    Clear-NativeBuildProvenance -BinaryPath $BinaryPath
    try {
        if (-not (Test-Path -LiteralPath $BinaryPath -PathType Leaf)) {
            throw "Cannot record native build provenance; binary is missing: $BinaryPath"
        }
        $endIdentity = Get-NativeBuildSourceIdentity -RepoRoot $RepoRoot
        foreach ($field in @(
            "source_commit", "source_dirty", "source_change_fingerprint_sha256"
        )) {
            if ($endIdentity.$field -ne $StartIdentity.$field) {
                throw "Source identity changed during native build: $field"
            }
        }
        if ($StartIdentity.source_dirty) {
            Write-Warning (
                "Native build completed, but provenance was not issued because " +
                "the source checkout was dirty: $BinaryPath"
            )
            return
        }
        $binary = Get-Item -LiteralPath $BinaryPath
        $manifest = [ordered]@{
            schema = "radia.native-build-provenance.v1"
            source_commit = $StartIdentity.source_commit
            source_dirty = $StartIdentity.source_dirty
            source_change_fingerprint_sha256 = `
                $StartIdentity.source_change_fingerprint_sha256
            binary_name = $binary.Name
            binary_bytes = $binary.Length
            binary_sha256 = (Get-FileHash -LiteralPath $binary.FullName `
                -Algorithm SHA256).Hash.ToLowerInvariant()
            generated_at_utc = [DateTimeOffset]::UtcNow.ToString("o")
        }
        $json = $manifest | ConvertTo-Json -Depth 3
        [IO.File]::WriteAllText($temporaryPath, $json + [Environment]::NewLine,
            [Text.UTF8Encoding]::new($false))
        Move-Item -LiteralPath $temporaryPath -Destination $manifestPath -Force
        Write-Host "  Native provenance: $manifestPath" -ForegroundColor Cyan
    } catch {
        if (Test-Path -LiteralPath $temporaryPath) {
            Remove-Item -LiteralPath $temporaryPath -Force
        }
        Clear-NativeBuildProvenance -BinaryPath $BinaryPath
        throw
    }
}
