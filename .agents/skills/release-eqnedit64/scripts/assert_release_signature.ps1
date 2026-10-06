# Approved release certificates, including the historical LAB signer.
function Assert-EqneditReleaseSignature([string]$Path, [switch]$AllowUntrustedRoot) {
    $approved = @(
        '7E80D4365CFDBA7B22B33C35EB7CD6A41050E065',
        'BAEBD216FF82C03F57B2796C42DC8F4D083E07E1'
    )
    $signature = Get-AuthenticodeSignature -LiteralPath $Path
    # Pin the certificate before adding its public root on a disposable
    # hosted runner. Re-evaluate Authenticode afterwards and require Valid.
    if ($AllowUntrustedRoot -and $signature.SignerCertificate -and
        $signature.SignerCertificate.Thumbprint -in $approved -and
        $signature.SignerCertificate.Subject -ceq 'CN=ksugahar') {
        $store = [Security.Cryptography.X509Certificates.X509Store]::new('Root', 'CurrentUser')
        try {
            $store.Open('ReadWrite')
            $store.Add($signature.SignerCertificate)
        } finally { $store.Close() }
        $signature = Get-AuthenticodeSignature -LiteralPath $Path
    }
    if (-not $signature.SignerCertificate -or
        $signature.SignerCertificate.Subject -cne 'CN=ksugahar' -or
        $signature.SignerCertificate.Thumbprint -notin $approved -or
        $signature.Status -ne 'Valid') {
        throw "Unapproved Eqnedit64 signature: $($signature.Status)"
    }
    return $signature
}
