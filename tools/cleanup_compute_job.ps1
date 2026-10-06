param(
    [Parameter(Mandatory)][string]$JobRoot,
    [Parameter(Mandatory)][string]$RecoveryArchive,
    [Parameter(Mandatory)][string]$RecoveredArchiveSHA256,
    [Parameter(Mandatory)][string]$EvidenceCommit,
    [string[]]$StagingFiles = @()
)
$ErrorActionPreference = 'Stop'
if ($EvidenceCommit -notmatch '^[0-9a-f]{40}$') { throw 'Verified evidence commit SHA required' }
$tempRoot = [IO.Path]::GetFullPath('C:\temp')
$resolvedJob = [IO.Path]::GetFullPath($JobRoot).TrimEnd('\')
if ($resolvedJob -eq $tempRoot -or -not $resolvedJob.StartsWith($tempRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Job escapes C:\temp' }
$targets = @($resolvedJob, $RecoveryArchive) + $StagingFiles
foreach ($target in $targets) {
    $fullTarget = [IO.Path]::GetFullPath($target)
    if (-not $fullTarget.StartsWith($tempRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw "Escaping target: $target" }
    $item = Get-Item -LiteralPath $fullTarget
    if ($item.Attributes -band [IO.FileAttributes]::ReparsePoint) { throw "Linked target: $target" }
    if ($item.PSIsContainer) {
        $links = Get-ChildItem -LiteralPath $fullTarget -Force -Recurse | Where-Object { $_.Attributes -band [IO.FileAttributes]::ReparsePoint }
        if ($links) { throw "Links/junctions inside target: $target" }
    }
}
$hash = (Get-FileHash -LiteralPath $RecoveryArchive -Algorithm SHA256).Hash.ToLowerInvariant()
if ($hash -ne $RecoveredArchiveSHA256) { throw 'Recovery archive changed since verified recovery' }
$active = Get-CimInstance Win32_Process -Filter "Name='python.exe'" | Where-Object { $_.CommandLine -and $_.CommandLine.Contains($resolvedJob, [StringComparison]::OrdinalIgnoreCase) }
if ($active) { throw 'Owned compute process still running' }
foreach ($target in $targets) { Remove-Item -LiteralPath $target -Recurse -Force }
$remaining = @($targets | Where-Object { Test-Path -LiteralPath $_ })
if ($remaining.Count) { throw 'Cleanup did not remove all owned targets' }
[pscustomobject]@{ status='removed'; host=$env:COMPUTERNAME; targets=$targets; evidence_commit=$EvidenceCommit; archive_sha256=$hash } | ConvertTo-Json -Depth 4
