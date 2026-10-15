"""Reject incomplete/nonfinite IH CPU evidence before release qualification."""
import importlib.util
from pathlib import Path
import pytest

spec=importlib.util.spec_from_file_location('ih_installed_acceptance',Path(__file__).resolve().parents[3]/'validation_test/ih_cuda/run_installed.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)

@pytest.mark.parametrize('fault',[None,'nan','negative','boolean','missing','residual','rejected'])
def test_ih_cpu_qualification(fault):
    cpu=dict(accepted=True,relative_errors=dict(temperature=0.,heat=1e-12,impedance=1e-12,electromagnetic=1e-12),true_relative_residual=1e-12)
    if fault=='nan':cpu['relative_errors']['heat']=float('nan')
    if fault=='negative':cpu['relative_errors']['heat']=-1.
    if fault=='boolean':cpu['true_relative_residual']=True
    if fault=='missing':cpu['relative_errors'].pop('impedance')
    if fault=='residual':cpu['true_relative_residual']=1e-5
    if fault=='rejected':cpu['accepted']=False
    if fault:
        with pytest.raises(ValueError):runner.check_cpu_release_gate(cpu)
    else:runner.check_cpu_release_gate(cpu)
