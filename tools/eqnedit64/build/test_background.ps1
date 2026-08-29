$ErrorActionPreference = 'Stop'

$projectRoot = Split-Path -Parent $PSScriptRoot
$app = Join-Path $PSScriptRoot 'Eqnedit64.exe'
$portableApp = Join-Path $projectRoot 'dist\Eqnedit64.exe'
$unitTest = Join-Path $PSScriptRoot 'test_tex_document.exe'
$operationOutput = 'C:\temp\Eqnedit64-operation-test.tex'
$operationLogOutput = 'C:\temp\Eqnedit64-operation-log-test.log'
$usabilityBundleOutput = 'C:\temp\Eqnedit64-usability-bundle-test.json'

if (-not (Test-Path -LiteralPath $app)) {
    throw "Eqnedit64.exe is missing. Run build_eqnedt64.bat first."
}
if (-not (Test-Path -LiteralPath $portableApp)) {
    throw "Portable Eqnedit64.exe is missing. Run build_eqnedt64.bat first."
}
$releaseFiles = @(Get-ChildItem -LiteralPath (Split-Path -Parent $portableApp) -File)
if ($releaseFiles.Count -ne 1 -or $releaseFiles[0].Name -ne 'Eqnedit64.exe') {
    throw "dist must contain only Eqnedit64.exe: $($releaseFiles.Name -join ', ')"
}
if ((Get-FileHash -Algorithm SHA256 -LiteralPath $app).Hash -ne
    (Get-FileHash -Algorithm SHA256 -LiteralPath $portableApp).Hash) {
    throw 'Portable executable differs from the tested build executable.'
}
$signature = Get-AuthenticodeSignature -LiteralPath $portableApp
if ($signature.Status -ne 'Valid' -or
    $signature.SignerCertificate.Subject -ne
        'CN=ksugahar' -or
    -not ($signature.SignerCertificate.EnhancedKeyUsageList.ObjectId -contains
        '1.3.6.1.5.5.7.3.3')) {
    throw "Portable executable has no valid Eqnedit64 developer signature: $($signature.Status)"
}

& $unitTest
if ($LASTEXITCODE -ne 0) { throw "TeX document tests failed: $LASTEXITCODE" }

$self = Start-Process -FilePath $app -ArgumentList '--self-test' `
    -WindowStyle Hidden -Wait -PassThru
if ($self.ExitCode -ne 0) { throw "Eqnedit64 self-test failed: $($self.ExitCode)" }

$statusLayout = Start-Process -FilePath $app -ArgumentList '--status-layout-test' `
    -WindowStyle Hidden -Wait -PassThru
if ($statusLayout.ExitCode -ne 0) {
    throw "Eqnedit64 hidden status-layout test failed: $($statusLayout.ExitCode)"
}

$visualScale = Start-Process -FilePath $app -ArgumentList '--visual-scale-test' `
    -WindowStyle Hidden -Wait -PassThru
if ($visualScale.ExitCode -ne 0) {
    throw "Eqnedit64 hidden visual-scale test failed: $($visualScale.ExitCode)"
}

$uiInteraction = Start-Process -FilePath $app -ArgumentList '--ui-interaction-test' `
    -WindowStyle Hidden -Wait -PassThru
if ($uiInteraction.ExitCode -ne 0) {
    throw "Eqnedit64 hidden UI-interaction test failed: $($uiInteraction.ExitCode)"
}

$operation = Start-Process -FilePath $app `
    -ArgumentList @('--operation-test', $operationOutput) `
    -WindowStyle Hidden -Wait -PassThru
if ($operation.ExitCode -ne 0) {
    throw "Eqnedit64 operation test failed: $($operation.ExitCode)"
}

$operationLog = Start-Process -FilePath $app `
    -ArgumentList @('--operation-log-test', $operationLogOutput) `
    -WindowStyle Hidden -Wait -PassThru
if ($operationLog.ExitCode -ne 0) {
    throw "Eqnedit64 operation-log test failed: $($operationLog.ExitCode)"
}

$logText = [IO.File]::ReadAllText($operationLogOutput, [Text.UTF8Encoding]::new($true))
$logRequirements = @(
    'Eqnedit64 operation log v2',
    'debug.start',
    'test.keyboard',
    'user.marker',
    'caret=',
    'selection=',
    'dirty=',
    'elapsed_ms=',
    'delta_ms=',
    'focus=',
    'input_style=',
    'alignment=',
    'zoom_percent=',
    'equation_mode=',
    'latex=x^{2}',
    'debug.stop'
)
foreach ($required in $logRequirements) {
    if (-not $logText.Contains($required)) {
        throw "Operation log is missing [$required]:`n$logText"
    }
}

& (Join-Path $PSScriptRoot 'analyze_last_operation_log.ps1') `
    -LogPath $operationLogOutput -OutputPath $usabilityBundleOutput `
    -Privacy structure
$usabilityBundleText = [IO.File]::ReadAllText(
    $usabilityBundleOutput, [Text.UTF8Encoding]::new($false))
$usabilityBundle = $usabilityBundleText | ConvertFrom-Json
if ($usabilityBundle.schema -ne 'eqnedit64.usability-review-bundle.v1' -or
    $usabilityBundle.operation_log_version -ne 2 -or
    $usabilityBundle.privacy -ne 'structure' -or
    $usabilityBundle.summary.event_count -ne 4 -or
    $usabilityBundle.summary.explicit_marker_count -ne 1 -or
    $usabilityBundle.candidates[0].detector -ne 'explicit_user_marker') {
    throw "Usability review bundle contract is invalid:`n$usabilityBundleText"
}
if ($usabilityBundleText.Contains('x^{2}')) {
    throw 'Structure-privacy usability bundle leaked raw equation text.'
}

$saved = [IO.File]::ReadAllText($operationOutput, [Text.UTF8Encoding]::new($false))
$requirements = @(
    '\begin{equation}',
    '\begin{aligned}',
    'F',
    '= ma',
    '\\',
    'E',
    '= mc^{2}',
    '\end{aligned}',
    '\end{equation}'
)
foreach ($required in $requirements) {
    if (-not $saved.Contains($required)) {
        throw "Operation output is missing [$required]:`n$saved"
    }
}

$bytes = [IO.File]::ReadAllBytes($app)
$peOffset = [BitConverter]::ToInt32($bytes, 0x3c)
$machine = [BitConverter]::ToUInt16($bytes, $peOffset + 4)
if ($machine -ne 0x8664) { throw ('Executable is not x64: 0x{0:X4}' -f $machine) }
$optionalHeader = $peOffset + 24
$optionalMagic = [BitConverter]::ToUInt16($bytes, $optionalHeader)
if ($optionalMagic -ne 0x20b) { throw 'Executable is not PE32+.' }
$dllCharacteristics = [BitConverter]::ToUInt16($bytes, $optionalHeader + 0x46)
if (($dllCharacteristics -band 0x4000) -eq 0) {
    throw 'Executable does not enable Control Flow Guard.'
}

# Resolve the PE debug directory and require an RSDS CodeView record. /Fd by
# itself produces no usable symbols; this catches accidental removal of /Zi or
# /DEBUG before a crash leaves another unresolvable minidump.
$debugRva = [BitConverter]::ToUInt32($bytes, $optionalHeader + 112 + 6 * 8)
$debugSize = [BitConverter]::ToUInt32($bytes, $optionalHeader + 112 + 6 * 8 + 4)
$sectionCount = [BitConverter]::ToUInt16($bytes, $peOffset + 6)
$optionalSize = [BitConverter]::ToUInt16($bytes, $peOffset + 20)
$sectionTable = $optionalHeader + $optionalSize
$debugOffset = $null
for ($i = 0; $i -lt $sectionCount; $i++) {
    $section = $sectionTable + 40 * $i
    $virtualSize = [BitConverter]::ToUInt32($bytes, $section + 8)
    $virtualAddress = [BitConverter]::ToUInt32($bytes, $section + 12)
    $rawSize = [BitConverter]::ToUInt32($bytes, $section + 16)
    $rawPointer = [BitConverter]::ToUInt32($bytes, $section + 20)
    $span = [Math]::Max($virtualSize, $rawSize)
    if ($debugRva -ge $virtualAddress -and
        $debugRva -lt $virtualAddress + $span) {
        $debugOffset = $rawPointer + ($debugRva - $virtualAddress)
        break
    }
}
if ($null -eq $debugOffset -or $debugSize -lt 28) {
    throw 'Executable has no readable PE debug directory.'
}
$hasCodeView = $false
for ($at = [int]$debugOffset; $at + 28 -le $debugOffset + $debugSize; $at += 28) {
    $type = [BitConverter]::ToUInt32($bytes, $at + 12)
    $dataSize = [BitConverter]::ToUInt32($bytes, $at + 16)
    $dataPointer = [BitConverter]::ToUInt32($bytes, $at + 24)
    if ($type -eq 2 -and $dataSize -ge 24 -and
        $dataPointer + 4 -le $bytes.Length -and
        [Text.Encoding]::ASCII.GetString($bytes, $dataPointer, 4) -eq 'RSDS') {
        $hasCodeView = $true
        break
    }
}
if (-not $hasCodeView) { throw 'Executable has no RSDS CodeView debug record.' }
$commit = (& git -C $projectRoot rev-parse --short HEAD).Trim()
foreach ($symbol in @('Eqnedit64.pdb', 'Eqnedit64.map')) {
    $symbolPath = Join-Path $projectRoot ("symbols\$commit\$symbol")
    if (-not (Test-Path -LiteralPath $symbolPath -PathType Leaf)) {
        throw "Archived release symbol is missing: $symbolPath"
    }
}

$version = (Get-Item -LiteralPath $app).VersionInfo
if ($version.ProductName -ne 'Eqnedit64' -or $version.ProductVersion -notlike '3.0.3*') {
    throw "Version resource is missing or invalid."
}

$asciiImage = [Text.Encoding]::ASCII.GetString($bytes).ToUpperInvariant()
foreach ($forbiddenDependency in @('VCRUNTIME', 'MSVCP', 'PYTHON3', 'EQNEDIT_CORE.PYD')) {
    if ($asciiImage.Contains($forbiddenDependency)) {
        throw "Portable executable contains forbidden runtime dependency: $forbiddenDependency"
    }
}
$unicodeImage = [Text.Encoding]::Unicode.GetString($bytes)
if (-not $unicodeImage.Contains('Copyright (c) 2019-2025 The Bootstrap Authors') -or
    -not $unicodeImage.Contains('Permission is hereby granted, free of charge')) {
    throw 'Portable executable does not contain the required embedded icon license notice.'
}

foreach ($portableMode in @('--self-test', '--status-layout-test',
        '--visual-scale-test', '--ui-interaction-test')) {
    $portable = Start-Process -FilePath $portableApp -ArgumentList $portableMode `
        -WorkingDirectory (Split-Path -Parent $portableApp) `
        -WindowStyle Hidden -Wait -PassThru
    if ($portable.ExitCode -ne 0) {
        throw "Portable executable failed $portableMode`: $($portable.ExitCode)"
    }
}

Remove-Item -LiteralPath $operationOutput -Force
Remove-Item -LiteralPath $operationLogOutput -Force
Remove-Item -LiteralPath $usabilityBundleOutput -Force

Write-Host 'PASS: TeX paste/envelope unit tests'
Write-Host 'PASS: hidden Eqnedit64 self-test'
Write-Host 'PASS: hidden native status-bar parts, font-height, and simple-mode test'
Write-Host 'PASS: offscreen 96/120/144/192 dpi equation, caret, selection, and font visibility'
Write-Host 'PASS: hidden WM_CHAR, WM_KEYDOWN, WM_COMMAND, source, and status interaction test'
Write-Host 'PASS: raw TeX, PowerPoint-native Office Math, EMF/DIBV5, and Google Slides 300 dpi/24 pt payload generation'
Write-Host 'PASS: background command-sequence save test'
Write-Host 'PASS: background operation-log and F12 marker test'
Write-Host 'PASS: privacy-aware LLM usability bundle from the real operation log'
Write-Host 'PASS: x64 PE and version resources'
Write-Host 'PASS: Control Flow Guard, RSDS CodeView record, and commit-archived PDB/map'
Write-Host 'PASS: valid ksugahar developer signature'
Write-Host 'PASS: single-file portable release, embedded license, static runtime, and isolated hidden startup'
