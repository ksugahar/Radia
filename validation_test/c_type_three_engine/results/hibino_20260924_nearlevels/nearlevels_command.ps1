$env:OMP_NUM_THREADS='8'
$env:MKL_NUM_THREADS='8'
$env:OPENBLAS_NUM_THREADS='8'
Set-Location C:/temp/hdiv-nearlevels-20260924
$env:HDIV_LEVEL='coarse'
./venv/Scripts/python.exe -u compare_near_levels_hdiv.py 2>&1 | Tee-Object -FilePath coarse.log
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
$env:HDIV_LEVEL='medium'
./venv/Scripts/python.exe -u compare_near_levels_hdiv.py 2>&1 | Tee-Object -FilePath medium.log
exit $LASTEXITCODE
