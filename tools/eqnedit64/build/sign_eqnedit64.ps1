param(
    [Parameter(Mandatory = $true)]
    [string]$Path,
    [string]$Subject = 'CN=ksugahar'
)

$ErrorActionPreference = 'Stop'
$resolved = (Resolve-Path -LiteralPath $Path).Path
$certificate = Get-ChildItem -LiteralPath 'Cert:\CurrentUser\My' | Where-Object {
    $_.Subject -eq $Subject -and $_.HasPrivateKey -and
    $_.NotAfter -gt (Get-Date) -and
    $_.EnhancedKeyUsageList.ObjectId -contains '1.3.6.1.5.5.7.3.3'
} | Sort-Object NotAfter -Descending | Select-Object -First 1
if (-not $certificate) {
    throw 'Eqnedit64 developer certificate is missing. Run build\setup_developer_signing.ps1 once.'
}

$candidates = @()
$kitsBin = 'C:\Program Files (x86)\Windows Kits\10\bin'
if (Test-Path -LiteralPath $kitsBin) {
    $candidates += Get-ChildItem -LiteralPath $kitsBin -Directory |
        Sort-Object Name -Descending |
        ForEach-Object { Join-Path $_.FullName 'x64\signtool.exe' }
}
$candidates +=
    'C:\Program Files (x86)\Windows Kits\10\App Certification Kit\signtool.exe'
$signTool = $candidates | Where-Object { Test-Path -LiteralPath $_ } |
    Select-Object -First 1
if (-not $signTool) { throw 'signtool.exe was not found in the Windows SDK.' }

& $signTool sign /fd SHA256 /sha1 $certificate.Thumbprint /s My `
    /d 'Eqnedit64 - TeX Equation Editor' $resolved | Out-Host
if ($LASTEXITCODE -ne 0) { throw "signtool failed with exit code $LASTEXITCODE." }

$signature = Get-AuthenticodeSignature -LiteralPath $resolved
if ($signature.Status -ne 'Valid' -or
    $signature.SignerCertificate.Thumbprint -ne $certificate.Thumbprint) {
    throw "Authenticode verification failed: $($signature.Status) $($signature.StatusMessage)"
}
Write-Host ('[OK] developer-signed {0} ({1})' -f $resolved,
    $certificate.Thumbprint)
