$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot '../../../.agents/skills/release-eqnedit64/scripts/assert_release_signature.ps1')
function Get-AuthenticodeSignature { return $script:fixture }
$approved = '7E80D4365CFDBA7B22B33C35EB7CD6A41050E065'
foreach ($sample in @(
    @('Valid', 'CN=ksugahar', $approved, $true),
    @('Valid', 'CN=ksugahar', 'BAEBD216FF82C03F57B2796C42DC8F4D083E07E1', $true),
    @('Valid', 'CN=ksugahar', ('0' * 40), $false),
    @('Valid', 'CN=other', $approved, $false),
    @('UnknownError', 'CN=ksugahar', $approved, $false),
    @('NotTrusted', 'CN=ksugahar', $approved, $false),
    @('HashMismatch', 'CN=ksugahar', $approved, $false),
    @('NotSigned', 'CN=ksugahar', $approved, $false))) {
    $script:fixture = [pscustomobject]@{
        Status = $sample[0]
        SignerCertificate = [pscustomobject]@{ Subject=$sample[1]; Thumbprint=$sample[2] }
    }
    $accepted = $false
    try { $null = Assert-EqneditReleaseSignature 'fixture.exe'; $accepted = $true } catch {}
    if ($accepted -ne $sample[3]) { throw "Signer policy failed: $($sample[0..2] -join ', ')" }
}
Write-Host 'PASS: pinned signer policy rejects same-name keys and invalid signatures'
