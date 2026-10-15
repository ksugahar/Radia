param(
    [int]$ProcessId,
    [ValidateSet('radia_cuda_ic.mexw64','radia_cuda_ams.mexw64','cusparse64_12.dll','cublas64_12.dll')][string]$ModuleName,
    [string]$OutputPath
)
$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$modules=@(Get-Process -Id $ProcessId -Module | Where-Object ModuleName -eq $ModuleName)
if($modules.Count -ne 1){throw 'Expected one loaded numerical MEX module'}
$module=$modules[0]
$hash=[Security.Cryptography.SHA256]::Create()
try {$digest=([BitConverter]::ToString($hash.ComputeHash([IO.File]::ReadAllBytes($module.FileName)))).Replace('-','').ToLower()}
finally {$hash.Dispose()}
@{path=$module.FileName;sha256=$digest}|ConvertTo-Json -Depth 4|Set-Content -LiteralPath $OutputPath -Encoding UTF8
