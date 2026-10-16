"""Test artifacts must be classified separately before any MATLAB process starts."""
import json
from pathlib import Path
import subprocess
import sys
import pytest

@pytest.mark.parametrize('test',[False,True])
def test_ih_explicit_build_classification(test):
    script=Path(__file__).parents[1]/'tools/build_ih_cuda.py'
    command=[sys.executable,str(script),'--describe-build-mode']
    if test:command.append('--transaction-test')
    report=json.loads(subprocess.check_output(command,text=True))
    assert report['transaction_test_build'] is test
    assert report['distribution_eligible'] is (not test)
    assert report['output_directory_name']==('build-ih-transaction-test' if test else 'build-ih-cuda')
