$env:OMP_NUM_THREADS='8'
$env:MKL_NUM_THREADS='8'
$env:OPENBLAS_NUM_THREADS='8'
Set-Location C:/temp/hdiv-bdm2levels-20260924
$env:HDIV_LEVEL='medium'
./venv/Scripts/python.exe -u bdm2_levels.py 2>&1 | Tee-Object -FilePath medium.log
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
$env:HDIV_LEVEL='fine'
./venv/Scripts/python.exe -u bdm2_levels.py 2>&1 | Tee-Object -FilePath fine.log
exit $LASTEXITCODE
