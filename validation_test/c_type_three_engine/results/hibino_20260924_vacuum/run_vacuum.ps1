$ErrorActionPreference='Stop'
Set-Location C:/temp/hdiv-vacuum-20260924
if (Test-Path vacuum.json) { throw 'Prior result exists' }
$env:OMP_NUM_THREADS='8'
$env:MKL_NUM_THREADS='8'
$env:OPENBLAS_NUM_THREADS='8'
python -m venv --system-site-packages venv
if ($LASTEXITCODE -ne 0) { throw 'venv failed' }
./venv/Scripts/python.exe -m pip install --no-deps ./radia-5.0.0-cp312-cp312-win_amd64.whl
if ($LASTEXITCODE -ne 0) { throw 'wheel failed' }
./venv/Scripts/python.exe -m py_compile vacuum_omega.py
if ($LASTEXITCODE -ne 0) { throw 'compile failed' }
./venv/Scripts/python.exe -u vacuum_omega.py 2>&1 | Tee-Object -FilePath run.log
exit $LASTEXITCODE
