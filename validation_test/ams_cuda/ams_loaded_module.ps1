param([int]$ProcessId,[string]$OutputPath)
$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$modules=@(Get-Process -Id $ProcessId -Module | Where-Object ModuleName -eq 'radia_cuda_ams.mexw64')
if($modules.Count -ne 1){throw 'Expected exactly one loaded AMS MEX module'}
$module=$modules[0]
@{path=$module.FileName;sha256=([BitConverter]::ToString([Security.Cryptography.SHA256]::Create().ComputeHash([IO.File]::ReadAllBytes($module.FileName)))).Replace('-','').ToLower()} | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $OutputPath -Encoding UTF8
