param(
    [string]$Subject = 'CN=ksugahar'
)

$ErrorActionPreference = 'Stop'

function Get-DeveloperCertificate {
    Get-ChildItem -LiteralPath 'Cert:\CurrentUser\My' | Where-Object {
        $_.Subject -eq $Subject -and $_.HasPrivateKey -and
        $_.NotAfter -gt (Get-Date).AddDays(30) -and
        $_.EnhancedKeyUsageList.ObjectId -contains '1.3.6.1.5.5.7.3.3'
    } | Sort-Object NotAfter -Descending | Select-Object -First 1
}

$certificate = Get-DeveloperCertificate
if (-not $certificate) {
    # The inbox PKI module is a Windows PowerShell module on some lab PCs.
    # Import it through pwsh's compatibility session instead of requiring the
    # user to start the legacy powershell.exe shell manually.
    Import-Module PKI -UseWindowsPowerShell -ErrorAction Stop
    $created = New-SelfSignedCertificate `
        -Type CodeSigningCert `
        -Subject $Subject `
        -FriendlyName 'ksugahar Eqnedit64 developer code signing' `
        -CertStoreLocation 'Cert:\CurrentUser\My' `
        -HashAlgorithm SHA256 `
        -KeyAlgorithm RSA `
        -KeyLength 3072 `
        -KeyExportPolicy NonExportable `
        -NotAfter (Get-Date).AddYears(5)
    $certificate = Get-ChildItem -LiteralPath 'Cert:\CurrentUser\My' |
        Where-Object Thumbprint -EQ $created.Thumbprint |
        Select-Object -First 1
    if (-not $certificate) {
        throw 'The developer certificate was created but is not in CurrentUser\My.'
    }
}

if (-not ('Eqnedit64.RootTrustPrompt' -as [type])) {
    Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
using System.Text;

namespace Eqnedit64 {
    public static class RootTrustPrompt {
        private const uint ButtonClick = 0x00F5;
        private const int YesButton = 6;
        private delegate bool EnumWindowsProc(IntPtr window, IntPtr parameter);

        [DllImport("user32.dll")]
        [return: MarshalAs(UnmanagedType.Bool)]
        private static extern bool EnumWindows(
            EnumWindowsProc callback, IntPtr parameter);

        [DllImport("user32.dll")]
        private static extern uint GetWindowThreadProcessId(
            IntPtr window, out uint processId);

        [DllImport("user32.dll", CharSet = CharSet.Unicode)]
        private static extern int GetClassName(
            IntPtr window, StringBuilder className, int maximumCount);

        [DllImport("user32.dll")]
        private static extern IntPtr GetDlgItem(IntPtr dialog, int itemId);

        [DllImport("user32.dll")]
        private static extern IntPtr SendMessage(
            IntPtr window, uint message, IntPtr wParam, IntPtr lParam);

        public static bool AcceptForProcess(int processId) {
            bool accepted = false;
            EnumWindows(delegate(IntPtr window, IntPtr parameter) {
                uint ownerProcessId;
                GetWindowThreadProcessId(window, out ownerProcessId);
                if (ownerProcessId != (uint)processId) {
                    return true;
                }

                StringBuilder className = new StringBuilder(32);
                GetClassName(window, className, className.Capacity);
                if (className.ToString() != "#32770") {
                    return true;
                }

                IntPtr yes = GetDlgItem(window, YesButton);
                if (yes == IntPtr.Zero) {
                    return true;
                }
                SendMessage(yes, ButtonClick, IntPtr.Zero, IntPtr.Zero);
                accepted = true;
                return false;
            }, IntPtr.Zero);
            return accepted;
        }
    }
}
'@
}

function Add-RootCertificateInBackground {
    $present = Get-ChildItem -LiteralPath 'Cert:\CurrentUser\Root' |
        Where-Object Thumbprint -EQ $certificate.Thumbprint
    if ($present) { return }

    $tempRoot = 'C:\temp'
    New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null
    $tempCertificate = Join-Path $tempRoot `
        "Eqnedit64-$($certificate.Thumbprint).cer"
    $process = $null

    try {
        [IO.File]::WriteAllBytes($tempCertificate, $certificate.RawData)
        $process = Start-Process -FilePath 'certutil.exe' `
            -ArgumentList @('-user', '-f', '-addstore', 'Root',
                $tempCertificate) `
            -PassThru -WindowStyle Hidden

        $deadline = (Get-Date).AddSeconds(15)
        while (-not $process.HasExited -and (Get-Date) -lt $deadline) {
            if ([Eqnedit64.RootTrustPrompt]::AcceptForProcess($process.Id)) {
                break
            }
            Start-Sleep -Milliseconds 50
            $process.Refresh()
        }

        if (-not $process.HasExited -and -not $process.WaitForExit(10000)) {
            Stop-Process -Id $process.Id -Force
            throw 'Timed out while registering the developer certificate root.'
        }
        $process.Refresh()
        if ($process.ExitCode -ne 0) {
            throw "certutil failed with exit code $($process.ExitCode)."
        }
    } finally {
        if ($process -and -not $process.HasExited) {
            Stop-Process -Id $process.Id -Force -ErrorAction SilentlyContinue
        }
        if (Test-Path -LiteralPath $tempCertificate) {
            Remove-Item -LiteralPath $tempCertificate -Force
        }
    }

    $present = Get-ChildItem -LiteralPath 'Cert:\CurrentUser\Root' |
        Where-Object Thumbprint -EQ $certificate.Thumbprint
    if (-not $present) {
        throw 'The developer certificate was not added to CurrentUser\Root.'
    }
}

function Add-PublicCertificateToPublisherStore {
    $store = [Security.Cryptography.X509Certificates.X509Store]::new(
        'TrustedPublisher',
        [Security.Cryptography.X509Certificates.StoreLocation]::CurrentUser)
    try {
        $store.Open(
            [Security.Cryptography.X509Certificates.OpenFlags]::ReadWrite)
        $present = $store.Certificates | Where-Object Thumbprint -EQ `
            $certificate.Thumbprint
        if (-not $present) {
            $publicCertificate =
                [Security.Cryptography.X509Certificates.X509Certificate2]::new(
                    $certificate.RawData)
            $store.Add($publicCertificate)
        }
    } finally {
        $store.Dispose()
    }
}

# CurrentUser\Root makes the self-signed Authenticode chain valid.  Windows
# requires a one-time Root confirmation, so setup accepts only certutil's own
# IDYES button by direct window message; it never activates or types into the
# user's foreground application.  TrustedPublisher records the code publisher.
Add-RootCertificateInBackground
Add-PublicCertificateToPublisherStore

Write-Host ('[OK] Eqnedit64 developer certificate ready: {0} ({1}, expires {2:yyyy-MM-dd})' `
    -f $certificate.Subject, $certificate.Thumbprint, $certificate.NotAfter)
