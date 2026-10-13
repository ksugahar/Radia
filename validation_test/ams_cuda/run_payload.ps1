param([Parameter(Mandatory=$true)][string]$Job,
      [string]$MatlabRoot='C:/Program Files/MATLAB/R2026a')
$ErrorActionPreference='Stop'
$jobRoot=(Resolve-Path -LiteralPath $Job).Path
$manifest=Get-Content -LiteralPath "$jobRoot/manifest.json" -Raw | ConvertFrom-Json
foreach($item in $manifest.PSObject.Properties){
    if((Get-FileHash -LiteralPath (Join-Path $jobRoot $item.Name) -Algorithm SHA256).Hash.ToLower() -ne $item.Value){throw "Payload hash mismatch: $($item.Name)"}
}
# Child-only environment; neither installed Python nor NGSolve package paths.
$env:PATH="$jobRoot/runtime;$MatlabRoot/bin/win64;$env:SystemRoot/System32;$env:SystemRoot"
$env:OMP_NUM_THREADS='1';$env:MKL_NUM_THREADS='1';$env:MKL_THREADING_LAYER='SEQUENTIAL'
$expression="addpath('$($jobRoot.Replace('\','/').Replace("'","''"))');accept_payload('$($jobRoot.Replace('\','/').Replace("'","''"))')"
$process=Start-Process -FilePath "$MatlabRoot/bin/matlab.exe" -ArgumentList @('-nodesktop','-nosplash','-singleCompThread','-batch',('"'+$expression+'"'),'-logfile',('"'+$jobRoot+'/matlab.log"')) -WindowStyle Hidden -PassThru
$process.PriorityClass='BelowNormal'
$process.WaitForExit()
if($process.ExitCode -ne 0){throw "MATLAB payload failed: $($process.ExitCode)"}
