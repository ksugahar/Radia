$env:OMP_NUM_THREADS='8'
$env:MKL_NUM_THREADS='8'
$env:OPENBLAS_NUM_THREADS='8'
Set-Location C:/temp/hdiv-sourcequad-20260924
./venv/Scripts/python.exe -u compare_source_quadrature_hdiv.py 2>&1 | Tee-Object -FilePath run.log
exit $LASTEXITCODE
