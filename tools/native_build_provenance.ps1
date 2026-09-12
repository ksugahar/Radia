Set-StrictMode -Version Latest

function Invoke-NativeProvenanceGit {
    param(
        [Parameter(Mandatory)][string]$RepoRoot,
        [Parameter(Mandatory)][string[]]$Arguments
    )

    $output = @(& git -C $RepoRoot @Arguments 2>&1)
    if ($LASTEXITCODE -ne 0) {
        throw "git $($Arguments -join ' ') failed in ${RepoRoot}: $($output -join [Environment]::NewLine)"
    }
    return ($output -join "`n")
}

function Get-NativeBuildSourceIdentity {
    param([Parameter(Mandatory)][string]$RepoRoot)

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
    $stateHash = [Convert]::ToHexString(
        [Security.Cryptography.SHA256]::HashData($stateBytes)
    ).ToLowerInvariant()
    return [pscustomobject]@{
        source_commit = $commit
        source_dirty = [bool]($status.Length -gt 0)
        source_state_sha256 = $stateHash
    }
}

function Clear-NativeBuildProvenance {
    param([Parameter(Mandatory)][string]$BinaryPath)

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

    $manifestPath = "$BinaryPath.build.json"
    $temporaryPath = "$manifestPath.tmp-$PID-$([Guid]::NewGuid().ToString('N'))"
    Clear-NativeBuildProvenance -BinaryPath $BinaryPath
    try {
        if (-not (Test-Path -LiteralPath $BinaryPath -PathType Leaf)) {
            throw "Cannot record native build provenance; binary is missing: $BinaryPath"
        }
        $endIdentity = Get-NativeBuildSourceIdentity -RepoRoot $RepoRoot
        foreach ($field in @("source_commit", "source_dirty", "source_state_sha256")) {
            if ($endIdentity.$field -ne $StartIdentity.$field) {
                throw "Source identity changed during native build: $field"
            }
        }
        $binary = Get-Item -LiteralPath $BinaryPath
        $manifest = [ordered]@{
            schema = "radia.native-build-provenance.v1"
            source_commit = $StartIdentity.source_commit
            source_dirty = $StartIdentity.source_dirty
            source_state_sha256 = $StartIdentity.source_state_sha256
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
