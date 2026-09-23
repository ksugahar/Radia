$env:OMP_NUM_THREADS='8'
$env:MKL_NUM_THREADS='8'
$env:OPENBLAS_NUM_THREADS='8'
$env:HDIV_LEVEL='finer'
Set-Location C:/temp/hdiv-bdm2finer-20260924
./venv/Scripts/python.exe -u bdm2_finer.py 2>&1 | Tee-Object -FilePath run.log
exit $LASTEXITCODE
