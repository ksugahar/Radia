$env:OMP_NUM_THREADS='8'
$env:MKL_NUM_THREADS='8'
$env:OPENBLAS_NUM_THREADS='8'
Set-Location C:/temp/hdiv-fullfiner-20260924
./venv/Scripts/python.exe -u run_finer_three.py 2>&1 | Tee-Object -FilePath run.log
exit $LASTEXITCODE
