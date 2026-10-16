"""Test-only lifecycle injection must never select production artifact output."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT=Path(__file__).resolve().parents[1]/'tools/build_ams_cuda.py'

@pytest.mark.parametrize('flag',[None,'0','1'])
def test_lifecycle_build_is_explicit_and_separate(flag):
    env=os.environ.copy();env.pop('RADIA_AMS_GEOMETRY_LIFECYCLE_TEST',None)
    if flag is not None:env['RADIA_AMS_GEOMETRY_LIFECYCLE_TEST']=flag
    run=subprocess.run([sys.executable,str(SCRIPT),'--describe-build-mode'],env=env,capture_output=True,text=True,check=True)
    report=json.loads(run.stdout)
    test=flag=='1'
    assert report['geometry_lifecycle_test_build']==test
    assert report['distribution_eligible'] is (not test)
    assert report['output_directory_name']==('build-ams-geometry-test' if test else 'build-ams-cuda')
    assert bool(report['cuda_test_definitions'])==test

def test_unknown_lifecycle_flag_does_not_silently_build_production():
    env=os.environ.copy();env['RADIA_AMS_GEOMETRY_LIFECYCLE_TEST']='true'
    run=subprocess.run([sys.executable,str(SCRIPT),'--describe-build-mode'],env=env,capture_output=True,text=True)
    assert run.returncode!=0
    assert 'must be 0 or 1' in run.stderr
"""Packaging rejects lifecycle-test and unclassified AMS receipts."""
import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location('ams_distribution_verify', Path(__file__).parents[1] / 'verify_wheel.py')
verify = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)


@pytest.mark.parametrize('receipt,accepted', [
    ({'geometry_lifecycle_test_build': False, 'distribution_eligible': True}, True),
    ({'geometry_lifecycle_test_build': True, 'distribution_eligible': True}, False),
    ({'geometry_lifecycle_test_build': True, 'distribution_eligible': False}, False),
    ({'geometry_lifecycle_test_build': False, 'distribution_eligible': False}, False),
    ({'geometry_lifecycle_test_build': False}, False),
    ({}, False),
    ({'distribution_eligible': 1}, False),
    ({'distribution_eligible': 'true'}, False),
])
def test_ams_distribution_requires_explicit_production_eligibility(receipt, accepted):
    if accepted:
        verify.verify_ams_distribution_eligibility(receipt)
    else:
        with pytest.raises(ValueError):
            verify.verify_ams_distribution_eligibility(receipt)
